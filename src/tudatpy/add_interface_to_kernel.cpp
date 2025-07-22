#include <pybind11/pybind11.h>
#include "interface/expose_interface.h"

namespace py = pybind11;

void add_interface_to_kernel(py::module_& m)
{
    tudatpy::interface::expose_interface(m);
}