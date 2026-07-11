#ifndef TUDATPY_EXPOSE_DATA_ACCESS_H
#define TUDATPY_EXPOSE_DATA_ACCESS_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{

namespace data_access
{

void expose_data_access( py::module& m );

}  // namespace data_access

}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_DATA_ACCESS_H
