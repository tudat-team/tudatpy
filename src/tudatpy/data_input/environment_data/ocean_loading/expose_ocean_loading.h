#ifndef TUDATPY_EXPOSE_OCEAN_LOADING_H
#define TUDATPY_EXPOSE_OCEAN_LOADING_H

#include <pybind11/pybind11.h>

namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace ocean_loading
{
void expose_ocean_loading( pybind11::module& module );
}
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy

#endif
