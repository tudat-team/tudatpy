#include <pybind11/pybind11.h>
#include "constants/expose_constants.h"

namespace py = pybind11;

void add_constants_to_kernel(py::module_& m)
{
    tudatpy::constants::expose_constants(m);
}