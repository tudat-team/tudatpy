#ifndef TUDATPY_EXPOSE_STATION_LOADING_H
#define TUDATPY_EXPOSE_STATION_LOADING_H
#include <pybind11/pybind11.h>
namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace station_loading
{
void expose_station_loading( pybind11::module& module );
}
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
#endif
