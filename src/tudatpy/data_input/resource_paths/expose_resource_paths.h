#ifndef TUDATPY_EXPOSE_DATA_INPUT_RESOURCE_PATHS_H
#define TUDATPY_EXPOSE_DATA_INPUT_RESOURCE_PATHS_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{

namespace data_input
{

namespace resource_paths
{

void expose_resource_paths( py::module& m );

}  // namespace resource_paths

}  // namespace data_input

}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_DATA_INPUT_RESOURCE_PATHS_H
