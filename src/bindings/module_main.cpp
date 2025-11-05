#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include "bindings.hpp"

namespace nb = nanobind;

NB_MODULE(qss, m) {
  // Ensure NumPy is available for ndarray construction
  nb::module_::import_("numpy");
  m.doc() = "Quantum statistics solver bindings";

  // Call binder functions implemented in other .cpp files
  qss_bindings::bind_distributions(m);
  qss_bindings::bind_chemical_potential(m);
}
