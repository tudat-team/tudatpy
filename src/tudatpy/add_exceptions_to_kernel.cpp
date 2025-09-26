#include <pybind11/pybind11.h>
#include "exceptions/expose_exceptions.h"

namespace py = pybind11;

void add_exceptions_to_kernel(py::module_& m)
{
    tudatpy::exceptions::expose_exceptions(m);
}