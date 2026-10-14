#ifndef TUDATPY_EXPOSE_ATMOSPHERIC_LOADING_H
#define TUDATPY_EXPOSE_ATMOSPHERIC_LOADING_H
#include <pybind11/pybind11.h>
namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace atmospheric_loading
{
void expose_atmospheric_loading( pybind11::module& module );
}
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
#endif
