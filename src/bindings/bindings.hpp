#pragma once
#include <nanobind/nanobind.h>

namespace qss_bindings
{
void bind_distributions(nanobind::module_ &m);
void bind_chemical_potential(nanobind::module_ &m);
} // namespace qss_bindings
