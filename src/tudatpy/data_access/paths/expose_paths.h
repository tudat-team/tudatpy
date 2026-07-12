#ifndef TUDATPY_EXPOSE_DATA_ACCESS_PATHS_H
#define TUDATPY_EXPOSE_DATA_ACCESS_PATHS_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{

namespace data_access
{

namespace paths
{

void expose_paths( py::module& m );

}  // namespace paths

}  // namespace data_access

}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_DATA_ACCESS_PATHS_H
