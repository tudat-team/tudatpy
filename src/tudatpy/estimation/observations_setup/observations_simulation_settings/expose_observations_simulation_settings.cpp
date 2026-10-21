/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_observations_simulation_settings.h"

#include "expose_observations_simulation_settings_bindings.h"

namespace py = pybind11;

namespace tudatpy
{
namespace estimation
{
namespace observations_setup
{

namespace observations_simulation_settings
{

void expose_observation_simulation_settings_types( py::module& m )
{
    expose_observation_simulation_settings_core_bindings( m );
}

void expose_observations_simulation_settings( py::module& m )
{
    expose_observation_simulation_settings_factory_bindings( m );
}

}  // namespace observations_simulation_settings
}  // namespace observations_setup
}  // namespace estimation
}  // namespace tudatpy
