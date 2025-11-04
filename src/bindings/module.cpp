// module.cpp — nanobind vector-only API (fixed ndarray ctor)
// Build notes:
//  - Requires nanobind and Python Development.Module
//  - Optionally compile with OpenMP by defining _OPENMP and linking OpenMP::OpenMP_CXX
//  - Expects qss::chemical_potential_B(...) and related headers in include/qss/

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/vector.h>

#include <array>
#include <cstdint>    // int64_t for strides
#include <cstring>
#include <stdexcept>

#ifdef _OPENMP
  #include <omp.h>
#endif

// Project headers
#include "qss/types.hpp"
#include "qss/distributions.hpp"
#include "qss/solver.hpp"
#include "qss/solver_params.hpp"
#include "qss/chemical_potential.hpp"

namespace nb = nanobind;
using nb::arg;

// ----------------- small utilities -----------------
static qss::SolverParams make_params(double tol, int max_iters,
                                     double expand_factor, int max_expands,
                                     double a0, double b0) {
    qss::SolverParams p;
    p.tol = tol; p.max_iters = max_iters;
    p.expand_factor = expand_factor; p.max_expands = max_expands;
    p.a0 = a0; p.b0 = b0;
    return p;
}

// Allocate a C-contiguous 4D ndarray owned by a capsule
static nb::ndarray<double> make_c_contig_4d(nb::ssize_t s0, nb::ssize_t s1,
                                            nb::ssize_t s2, nb::ssize_t s3) {
    const size_t N = (size_t)s0 * (size_t)s1 * (size_t)s2 * (size_t)s3;
    double* data = new double[N];
    nb::capsule owner(data, [](void* p) noexcept { delete[] (double*)p; });

    // nanobind wants shape as size_t[] and strides (bytes) as int64_t[]
    std::array<size_t, 4> shape{ (size_t)s0, (size_t)s1, (size_t)s2, (size_t)s3 };
    std::array<int64_t, 4> strides;
    strides[3] = (int64_t)sizeof(double);
    strides[2] = (int64_t)s3 * strides[3];
    strides[1] = (int64_t)s2 * strides[2];
    strides[0] = (int64_t)s1 * strides[1];

    return nb::ndarray<double>((void*)data, 4, shape.data(), owner, strides.data());
}

// Strict 1-D double view (C- or non-C-contig is fine: we use stride)
struct Vec1D {
    const double* base {nullptr};
    nb::ssize_t stride_elems {0};
    size_t n {0};
    inline double at(size_t i) const { return *(base + (ptrdiff_t)i * stride_elems); }
};

static Vec1D as_vec1d(nb::ndarray<nb::ro, double> a, const char* name) {
    if (a.ndim() != 1) throw std::runtime_error(std::string(name) + " must be 1-D");
    Vec1D v;
    v.base = a.data();
    v.stride_elems = a.stride(0) / (nb::ssize_t)sizeof(double);
    v.n = (size_t)a.shape(0);
    if (v.n == 0) throw std::runtime_error(std::string(name) + " cannot be empty");
    return v;
}

// DistType vector: accept list/tuple/np array of enums/ints, 1-D only
static std::vector<qss::DistType> as_dist_vec(nb::handle h, const char* name) {
    if (!nb::isinstance<nb::iterable>(h))
        throw std::runtime_error(std::string(name) + " must be a 1-D sequence of DistType");
    nb::list seq = nb::cast<nb::list>(nb::module_::import_("builtins").attr("list")(h));
    if (seq.size() == 0) throw std::runtime_error(std::string(name) + " cannot be empty");
    std::vector<qss::DistType> out;
    out.reserve((size_t)seq.size());
    for (nb::handle it : seq) out.push_back(nb::cast<qss::DistType>(it));
    return out;
}

// ----------------- module -----------------
NB_MODULE(qss, m) {
    m.doc() = "Quantum statistics solver bindings (simple vector-only API)";

    nb::enum_<qss::DistType>(m, "DistType")
        .value("FD",     qss::DistType::FD)
        .value("FLB",    qss::DistType::FLB)
        .value("NFLB",   qss::DistType::NFLB)
        .value("NFLB_i", qss::DistType::NFLB_i);

    // distros, n_assumed, beta, H_a are all 1-D vectors (no scalars, no broadcasting)
    m.def("chemical_potential_B_grid",
        [](nb::handle distros_h,
           nb::ndarray<nb::ro, double> n_assumed_1d,
           nb::ndarray<nb::ro, double> beta_1d,
           nb::ndarray<nb::ro, double> H_a_T_1d,
           nb::ndarray<const double, nb::ndim<1>, nb::c_contig> kx, // len nx
           nb::ndarray<const double, nb::ndim<1>, nb::c_contig> ky, // len ny
           nb::ndarray<const double, nb::ndim<2>, nb::c_contig> E,  // (ny, nx)
           double g,
           double muB_eV_per_T,
           double tol, int max_iters,
           double expand_factor, int max_expands,
           double a0, double b0,
           int threads // 0 => default
        ) {
            if (E.ndim() != 2) throw std::runtime_error("E must be 2-D");
            const size_t ny = (size_t)E.shape(0), nx = (size_t)E.shape(1);
            if ((size_t)kx.size() != nx || (size_t)ky.size() != ny)
                throw std::runtime_error("shape(E) must be (len(ky), len(kx))");

            // Copy kx, ky, E into qss containers
            qss::Vec vx(nx), vy(ny);
            std::memcpy(vx.data(), kx.data(), nx * sizeof(double));
            std::memcpy(vy.data(), ky.data(), ny * sizeof(double));
            qss::Mat Em(ny, qss::Vec(nx));
            for (size_t j = 0; j < ny; ++j)
                std::memcpy(Em[j].data(), E.data() + j * nx, nx * sizeof(double));

            // Strict 1-D inputs
            std::vector<qss::DistType> distros = as_dist_vec(distros_h, "distros");
            Vec1D n  = as_vec1d(n_assumed_1d, "n_assumed");
            Vec1D bt = as_vec1d(beta_1d,      "beta");
            Vec1D hh = as_vec1d(H_a_T_1d,     "H_a_T");

            const size_t ND = distros.size();
            const size_t NN = n.n, NB = bt.n, NH = hh.n;

            auto out = make_c_contig_4d((nb::ssize_t)ND, (nb::ssize_t)NN,
                                        (nb::ssize_t)NB, (nb::ssize_t)NH);
            double* __restrict pout = out.data();

            qss::SolverParams params = make_params(
                tol, max_iters, expand_factor, max_expands, a0, b0);

            nb::gil_scoped_release release;

            #ifdef _OPENMP
            if (threads > 0) omp_set_num_threads(threads);
            #pragma omp parallel for collapse(4) schedule(static)
            #endif
            for (ptrdiff_t id = 0; id < (ptrdiff_t)ND; ++id)
            for (ptrdiff_t inx= 0; inx < (ptrdiff_t)NN; ++inx)
            for (ptrdiff_t ib = 0; ib < (ptrdiff_t)NB; ++ib)
            for (ptrdiff_t ih = 0; ih < (ptrdiff_t)NH; ++ih) {
                const qss::DistType dist = distros[(size_t)id];
                const double n_assumed   = n.at((size_t)inx);
                const double beta        = bt.at((size_t)ib);
                const double H_a_T       = hh.at((size_t)ih);

                const double mu = qss::chemical_potential_B(
                    dist, beta, n_assumed, vx, vy, Em, H_a_T, g, muB_eV_per_T, params);

                const size_t flat = (size_t)id * (NN*NB*NH)
                                  + (size_t)inx * (NB*NH)
                                  + (size_t)ib *  NH
                                  + (size_t)ih;
                pout[flat] = mu;
            }

            return out; // shape: (ND, NN, NB, NH)
        },
        arg("distros"),
        arg("n_assumed"),
        arg("beta"),
        arg("H_a_T"),
        arg("kx"),
        arg("ky"),
        arg("E"),
        arg("g") = 2.0,
        arg("muB_eV_per_T") = qss::QSS_MU_B_EV_PER_T,
        arg("tol") = 1e-6,
        arg("max_iters") = 300,
        arg("expand_factor") = 1.5,
        arg("max_expands") = 200,
        arg("a0") = -10.0,
        arg("b0") = 10.0,
        arg("threads") = 0
    );
}
