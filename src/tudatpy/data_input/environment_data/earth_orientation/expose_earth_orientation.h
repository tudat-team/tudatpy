#ifndef TUDATPY_EXPOSE_DATA_INPUT_EARTH_ORIENTATION_H
#define TUDATPY_EXPOSE_DATA_INPUT_EARTH_ORIENTATION_H
#include <pybind11/pybind11.h>
namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace earth_orientation
{
void expose_earth_orientation( pybind11::module& module );
}
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
#endif
