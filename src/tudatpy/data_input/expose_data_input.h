#ifndef TUDATPY_EXPOSE_DATA_INPUT_H
#define TUDATPY_EXPOSE_DATA_INPUT_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{

namespace data_input
{

void expose_data_input( py::module& m );

}  // namespace data_input

}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_DATA_INPUT_H
