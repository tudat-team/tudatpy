#include <pybind11/pybind11.h>
#include "astro/expose_astro.h"

namespace py = pybind11;

void add_astro_to_kernel(py::module_& m)
{
    tudatpy::astro::expose_astro(m);
}