// src/bindings/module.cpp (excerpt)
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/vector.h>
#include <stdexcept>
#include <cstring>
#include "qss/chemical_potential.hpp"
#include "qss/solver_params.hpp"
#include "qss/types.hpp"

namespace nb = nanobind;
using namespace nb::literals;

static qss::SolverParams make_params(
    double tol, int max_iters, double expand_factor, int max_expands,
    double a0, double b0) {
  qss::SolverParams p;
  p.tol = tol; p.max_iters = max_iters;
  p.expand_factor = expand_factor; p.max_expands = max_expands;
  p.a0 = a0; p.b0 = b0;
  return p;
}

NB_MODULE(qss, m) {
  nb::enum_<qss::DistType>(m, "DistType")
      .value("FD", qss::DistType::FD)
      .value("FLB", qss::DistType::FLB)
      .value("NFLB", qss::DistType::NFLB)
      .value("NFLB_i", qss::DistType::NFLB_i);

  // Scalar beta -> mu (B-based)
  m.def("chemical_potential_B",
    [](qss::DistType dist,
       double beta, double n_assumed,
       nb::ndarray<const double, nb::c_contig> kx,
       nb::ndarray<const double, nb::c_contig> ky,
       nb::ndarray<const double, nb::c_contig> E,
       double H_a_T, double g, double muB_eV_per_T,
       double tol, int max_iters,
       double expand_factor, int max_expands,
       double a0, double b0) {
      if (E.ndim()!=2) throw std::runtime_error("E must be 2D");
      const size_t ny = E.shape(0), nx = E.shape(1);
      if (kx.size()!=nx || ky.size()!=ny)
        throw std::runtime_error("shape(E) must be (len(ky), len(kx))");

      qss::Vec vx(nx), vy(ny);
      std::memcpy(vx.data(), kx.data(), nx*sizeof(double));
      std::memcpy(vy.data(), ky.data(), ny*sizeof(double));
      qss::Mat Em(ny, qss::Vec(nx));
      for (size_t j=0;j<ny;++j)
        std::memcpy(Em[j].data(), E.data()+j*nx, nx*sizeof(double));

      qss::SolverParams p = make_params(tol, max_iters, expand_factor, max_expands, a0, b0);

      double mu;
      {
        nb::gil_scoped_release _;
        mu = qss::chemical_potential_B(dist, beta, n_assumed, vx, vy, Em, H_a_T, g, muB_eV_per_T, p);
      }
      return mu;
    },
    "dist"_a, "beta"_a, "n_assumed"_a, "kx"_a, "ky"_a, "E"_a,
    "H_a_T"_a, "g"_a = 2.0, "muB_eV_per_T"_a = qss::QSS_MU_B_EV_PER_T,
    "tol"_a = 1e-6, "max_iters"_a = 300,
    "expand_factor"_a = 1.5, "max_expands"_a = 200,
    "a0"_a = -10.0, "b0"_a = 10.0
  );

  // Vector betas -> mus (B-based)
  m.def("chemical_potential_vec_B",
    [](qss::DistType dist,
       nb::ndarray<const double, nb::c_contig> betas, 
       nb::ndarray<const double, nb::c_contig> n_assumed,
       nb::ndarray<const double, nb::c_contig> kx,
       nb::ndarray<const double, nb::c_contig> ky,
       nb::ndarray<const double, nb::c_contig> E,
       nb::ndarray<const double, nb::c_contig> H_a_T, 
       double g, double muB_eV_per_T,
       double tol, int max_iters,
       double expand_factor, int max_expands,
       double a0, double b0) {
      if (betas.ndim()!=1) throw std::runtime_error("betas must be 1D");
      if (E.ndim()!=2) throw std::runtime_error("E must be 2D");
      const size_t ny = E.shape(0), nx = E.shape(1);
      if (kx.size()!=nx || ky.size()!=ny)
        throw std::runtime_error("shape(E) must be (len(ky), len(kx))");

      qss::Vec vx(nx), vy(ny);
      std::memcpy(vx.data(), kx.data(), nx*sizeof(double));
      std::memcpy(vy.data(), ky.data(), ny*sizeof(double));
      qss::Mat Em(ny, qss::Vec(nx));
      for (size_t j=0;j<ny;++j)
        std::memcpy(Em[j].data(), E.data()+j*nx, nx*sizeof(double));

      qss::Vec vbetas(betas.size());
      std::memcpy(vbetas.data(), betas.data(), betas.size()*sizeof(double));

      qss::SolverParams p = make_params(tol, max_iters, expand_factor, max_expands, a0, b0);

      qss::Vec mus;
      {
        nb::gil_scoped_release _;
        mus = qss::chemical_potential_B(dist, vbetas, n_assumed, vx, vy, Em, H_a_T, g, muB_eV_per_T, p);
      }

  // nanobind can convert std::vector<double> directly to a numpy array
  return nb::cast(mus);
    },
    "dist"_a, "betas"_a, "n_assumed"_a, "kx"_a, "ky"_a, "E"_a,
    "H_a_T"_a, "g"_a = 2.0, "muB_eV_per_T"_a = qss::QSS_MU_B_EV_PER_T,
    "tol"_a = 1e-6, "max_iters"_a = 300,
    "expand_factor"_a = 1.5, "max_expands"_a = 200,
    "a0"_a = -5000.0, "b0"_a = 1000.0
  );
}
