#ifndef TUDATPY_EXPOSE_DATA_ACCESS_TRACKING_FDETS_H
#define TUDATPY_EXPOSE_DATA_ACCESS_TRACKING_FDETS_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{

namespace data_access
{

namespace tracking
{

namespace fdets
{

void expose_fdets( py::module& m );

}  // namespace fdets

}  // namespace tracking

}  // namespace data_access

}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_DATA_ACCESS_TRACKING_FDETS_H
