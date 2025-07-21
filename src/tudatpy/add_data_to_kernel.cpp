#include <pybind11/pybind11.h>
#include "data/expose_data.h"

namespace py = pybind11;

void add_data_to_kernel(py::module_& m)
{
    tudatpy::data::expose_data(m);
}