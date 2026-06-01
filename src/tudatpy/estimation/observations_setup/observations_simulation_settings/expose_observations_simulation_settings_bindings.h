#ifndef TUDATPY_ESTIMATION_OBSERVATIONS_SETUP_OBSERVATION_SIMULATION_SETTINGS_BINDINGS_H
#define TUDATPY_ESTIMATION_OBSERVATIONS_SETUP_OBSERVATION_SIMULATION_SETTINGS_BINDINGS_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{
namespace estimation
{
namespace observations_setup
{
namespace observations_simulation_settings
{

void expose_observation_simulation_settings_core_bindings( py::module& m );
void expose_observation_simulation_settings_factory_bindings( py::module& m );

}  // namespace observations_simulation_settings
}  // namespace observations_setup
}  // namespace estimation
}  // namespace tudatpy

#endif
