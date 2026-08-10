#ifndef TUDATPY_EXPOSE_DATA_INPUT_ENVIRONMENT_DATA_SPACE_WEATHER_H
#define TUDATPY_EXPOSE_DATA_INPUT_ENVIRONMENT_DATA_SPACE_WEATHER_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace space_weather
{

void expose_space_weather( py::module& m );

}  // namespace space_weather
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_DATA_INPUT_ENVIRONMENT_DATA_SPACE_WEATHER_H
