#ifndef TUDATPY_EXPOSE_DATA_ACCESS_TRACKING_IFMS_H
#define TUDATPY_EXPOSE_DATA_ACCESS_TRACKING_IFMS_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{

namespace data_access
{

namespace tracking
{

namespace ifms
{

void expose_ifms( py::module& m );

}  // namespace ifms

}  // namespace tracking

}  // namespace data_access

}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_DATA_ACCESS_TRACKING_IFMS_H
