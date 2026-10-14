#ifndef TUDATPY_EXPOSE_STATION_HISTORY_H
#define TUDATPY_EXPOSE_STATION_HISTORY_H
#include <pybind11/pybind11.h>
namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace station_history
{
void expose_station_history( pybind11::module& m );
}
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
#endif
