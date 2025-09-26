#include <pybind11/pybind11.h>
#include "math/expose_math.h"

namespace py = pybind11;

void add_math_to_kernel(py::module_& m)
{
    tudatpy::math::expose_math(m);
}