#ifndef TUDATPY_EXPOSE_DATA_INPUT_ENVIRONMENT_DATA_ILRS_H
#define TUDATPY_EXPOSE_DATA_INPUT_ENVIRONMENT_DATA_ILRS_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace ilrs
{

void expose_ilrs( py::module& m );

}  // namespace ilrs
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_DATA_INPUT_ENVIRONMENT_DATA_ILRS_H
