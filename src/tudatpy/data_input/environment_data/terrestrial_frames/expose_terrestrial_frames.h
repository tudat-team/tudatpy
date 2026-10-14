#ifndef TUDATPY_EXPOSE_TERRESTRIAL_FRAMES_H
#define TUDATPY_EXPOSE_TERRESTRIAL_FRAMES_H
#include <pybind11/pybind11.h>
namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace terrestrial_frames
{
void expose_terrestrial_frames( pybind11::module& module );
}
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
#endif
