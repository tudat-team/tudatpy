#ifndef TUDATPY_EXPOSE_DATA_ACCESS_TRACKING_H
#define TUDATPY_EXPOSE_DATA_ACCESS_TRACKING_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{

namespace data_access
{

namespace tracking
{

void expose_tracking( py::module& m );

}  // namespace tracking

}  // namespace data_access

}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_DATA_ACCESS_TRACKING_H
