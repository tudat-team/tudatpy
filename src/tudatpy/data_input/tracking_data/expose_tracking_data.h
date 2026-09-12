#ifndef TUDATPY_EXPOSE_DATA_INPUT_TRACKING_DATA_H
#define TUDATPY_EXPOSE_DATA_INPUT_TRACKING_DATA_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{

namespace data_input
{

namespace tracking_data
{

void expose_tracking_data( py::module& m );

}  // namespace tracking_data

}  // namespace data_input

}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_DATA_INPUT_TRACKING_DATA_H
