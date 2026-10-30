#ifndef TUDATPY_EXPOSE_DATA_INPUT_TRACKING_DATA_LANDMARKS_H
#define TUDATPY_EXPOSE_DATA_INPUT_TRACKING_DATA_LANDMARKS_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{

namespace data_input
{

namespace tracking_data
{

namespace landmarks
{

void expose_landmarks( py::module& m );

}  // namespace landmarks

}  // namespace tracking_data

}  // namespace data_input

}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_DATA_INPUT_TRACKING_DATA_LANDMARKS_H
