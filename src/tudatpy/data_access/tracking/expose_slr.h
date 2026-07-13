#ifndef TUDATPY_EXPOSE_DATA_ACCESS_TRACKING_SLR_H
#define TUDATPY_EXPOSE_DATA_ACCESS_TRACKING_SLR_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{

namespace data_access
{

namespace tracking
{

namespace slr
{

void expose_slr( py::module& m );

}  // namespace slr

}  // namespace tracking

}  // namespace data_access

}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_DATA_ACCESS_TRACKING_SLR_H
