/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_observations_setup.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"

namespace py = pybind11;

namespace tudatpy
{
namespace estimation
{
namespace observations_setup
{

void expose_observations_setup( py::module& m )
{

    auto ancillary_settings = m.def_submodule( "ancillary_settings" );
    ancillary_settings::expose_ancillary_settings( ancillary_settings );

    auto observations_dependent_variables = m.def_submodule( "observations_dependent_variables" );
    observations_dependent_variables::expose_observations_dependent_variables( observations_dependent_variables );

    auto observations_simulation_settings = m.def_submodule( "observations_simulation_settings" );
    observations_simulation_settings::expose_observations_simulation_settings( observations_simulation_settings );

    auto observations_wrapper = m.def_submodule( "observations_wrapper" );
    observations_wrapper::expose_observations_wrapper( observations_wrapper );

    auto random_noise = m.def_submodule( "random_noise" );
    random_noise::expose_random_noise( random_noise );

    auto viability = m.def_submodule( "viability" );
    viability::expose_viability( viability );

}

}  // namespace observations_setup
}  // namespace estimation
}  // namespace tudatpy
