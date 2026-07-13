#ifndef TUDATPY_EXPOSE_DATA_ACCESS_ENVIRONMENT_SPACE_WEATHER_H
#define TUDATPY_EXPOSE_DATA_ACCESS_ENVIRONMENT_SPACE_WEATHER_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{
namespace data_access
{
namespace environment
{
namespace space_weather
{

void expose_space_weather( py::module& m );

}  // namespace space_weather
}  // namespace environment
}  // namespace data_access
}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_DATA_ACCESS_ENVIRONMENT_SPACE_WEATHER_H
