#ifndef TUDATPY_EXPOSE_DATA_ACCESS_ENVIRONMENT_ILRS_H
#define TUDATPY_EXPOSE_DATA_ACCESS_ENVIRONMENT_ILRS_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{
namespace data_access
{
namespace environment
{
namespace ilrs
{

void expose_ilrs( py::module& m );

}  // namespace ilrs
}  // namespace environment
}  // namespace data_access
}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_DATA_ACCESS_ENVIRONMENT_ILRS_H
