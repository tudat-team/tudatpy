#ifndef TUDATPY_EXPOSE_DATA_INPUT_TRACKING_DATA_FDETS_H
#define TUDATPY_EXPOSE_DATA_INPUT_TRACKING_DATA_FDETS_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{

namespace data_input
{

namespace tracking_data
{

namespace fdets
{

void expose_fdets( py::module& m );

}  // namespace fdets

}  // namespace tracking_data

}  // namespace data_input

}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_DATA_INPUT_TRACKING_DATA_FDETS_H
