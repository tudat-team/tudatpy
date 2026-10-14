#ifndef TUDATPY_EXPOSE_TERRESTRIAL_FRAME_TRANSFORMATIONS_H
#define TUDATPY_EXPOSE_TERRESTRIAL_FRAME_TRANSFORMATIONS_H
#include <pybind11/pybind11.h>
namespace tudatpy
{
namespace dynamics
{
namespace environment_setup
{
namespace terrestrial_frame_transformations
{
void expose_terrestrial_frame_transformations( pybind11::module& module );
}
}  // namespace environment_setup
}  // namespace dynamics
}  // namespace tudatpy
#endif
