#include "bindings.hpp"
#include "qss/distributions.hpp"
#include "qss/types.hpp"
#include <cstring>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

namespace nb = nanobind;
using Out1D = nb::ndarray<nb::numpy, double, nb::ndim<1>>;
using Vec1 = nb::ndarray<nb::numpy, const double, nb::ndim<1>, nb::c_contig>;

static Out1D make_out_1d(size_t n)
{
    double *data = new double[n];
    nb::capsule owner(data, [](void *p) noexcept { delete[] (double *)p; });
    return Out1D(data, {n}, owner, {1});
}

// ---- 2-D matrix variant: returns ndarray with SAME SHAPE (rows, cols) ----
using Mat2In = nb::ndarray<nb::numpy, const double, nb::ndim<2>, nb::c_contig>;
using Mat2Out = nb::ndarray<nb::numpy, double, nb::ndim<2>, nb::c_contig>;

// Helper: allocate a C-contiguous 2D float64 NumPy array of shape (rows, cols)
// Strides are in *elements* for this constructor.
static Mat2Out make_out_2d(size_t rows, size_t cols)
{
    double *data = new double[rows * cols];
    nb::capsule owner(data, [](void *p) noexcept { delete[] (double *)p; });
    // C-order: stride along row = cols, stride along col = 1
    return Mat2Out(data, {rows, cols}, owner, {(int64_t)cols, 1});
}

namespace qss_bindings
{

void bind_distributions(nb::module_ &m)
{
    nb::enum_<qss::DistType>(m, "DistType")
        .value("FD", qss::DistType::FD)
        .value("FLB", qss::DistType::FLB)
        .value("NFL_REWORKED", qss::DistType::NFL_REWORKED)
        .value("NFLB", qss::DistType::NFLB)
        .value("SSL", qss::DistType::SSL)
        .value("NFL", qss::DistType::NFL)
        .value("HUBBARD", qss::DistType::HUBBARD);
    // Scalar variant: returns float | None
    m.def(
        "distribution",
        [](qss::DistType dist, double mu, double kBT, double delta, int spin_shift, double U, double e_k, bool n_kd = false)
            -> std::optional<double> {
            return qss::distribution(e_k, mu, kBT, delta, spin_shift, U, dist, n_kd);
        },
        "dist",
        "mu",
        "kBT",
        "delta",
        "spin_shift",
        "U",
        "e_k",
        "n_kd",
        "Scalar distribution value at energy e_k (float or None).");

    // Vector variant: returns 1-D numpy.ndarray[float64]
    m.def(
        "distribution",
        [](qss::DistType dist, double mu, double kBT, double delta, int spin_shift, double U, Vec1 E, bool n_kd = false) -> Out1D {
            const size_t n = (size_t)E.shape(0);
            if (n == 0)
                return make_out_1d(0);

            qss::Vec Ev(n);
            std::memcpy(Ev.data(), E.data(), n * sizeof(double));

            auto res = qss::distribution(Ev, mu, kBT, delta, spin_shift, U, dist, n_kd);
            if (!res)
                throw std::runtime_error("distribution: invalid DistType");

            Out1D out = make_out_1d(n);
            std::memcpy(out.data(), res->data(), n * sizeof(double));
            return out;
        },
        "dist",
        "mu",
        "kBT",
        "delta",
        "spin_shift",
        "U",
        "E",
        "n_kd",
        "Vectorized distribution over a 1-D energy array E.");

    m.def(
        "distribution",
        [](qss::DistType dist, double mu, double kBT, double delta, int spin_shift, double U, Mat2In E, bool n_kd = false) -> Mat2Out {
            const size_t rows = E.shape(0);
            const size_t cols = E.shape(1);

            // Copy input
            qss::Mat mat_in(rows, qss::Vec(cols));
            const double *in_ptr = E.data();
            for (size_t i = 0; i < rows; ++i)
            {
                std::memcpy(mat_in[i].data(), in_ptr + i * cols,
                            cols * sizeof(double));
            }

            // Call C++ implementation
            auto res = qss::distribution(mat_in, mu, kBT, delta, spin_shift, U, dist, n_kd);
            if (!res)
                throw std::runtime_error("distribution: invalid DistType");

            // Copy output
            Mat2Out out = make_out_2d(rows, cols);
            double *out_ptr = out.data();
            for (size_t i = 0; i < rows; ++i)
            {
                std::memcpy(out_ptr + i * cols, (*res)[i].data(),
                            cols * sizeof(double));
            }
            return out;
        },
        nb::arg("dist"),
        nb::arg("mu"),
        nb::arg("kBT"),
        nb::arg("delta"),
        nb::arg("spin_shift"),
        nb::arg("U"),
        nb::arg("E"),
        nb::arg("n_kd") = false,
        "Elementwise distribution over a 2-D float64 matrix E (C-contiguous).\n"
        "Returns an array with the same shape (rows, cols).");
}

} // namespace qss_bindings
