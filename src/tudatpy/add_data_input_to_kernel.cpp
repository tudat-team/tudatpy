#include <pybind11/pybind11.h>

#include "data_input/expose_data_input.h"

namespace py = pybind11;

void add_data_input_to_kernel( py::module_& m )
{
    tudatpy::data_input::expose_data_input( m );
}
