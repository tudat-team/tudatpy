#ifndef TUDATPY_EXPOSE_GRAVITY_DATA_H
#define TUDATPY_EXPOSE_GRAVITY_DATA_H

#include <pybind11/pybind11.h>

namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace gravity_field
{

//! Register gravity provider records, conventions and filename readers.
void expose_gravity_field( pybind11::module& module );

}  // namespace gravity_field
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy

#endif
