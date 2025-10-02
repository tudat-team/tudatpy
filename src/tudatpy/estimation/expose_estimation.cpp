/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_estimation.h"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/basics/timeType.h"

namespace py = pybind11;
namespace tudatpy
{

namespace estimation
{

void expose_estimation( py::module &m )
{
    auto observable_models_setup_submodule = m.def_submodule( "observable_models_setup" );
    observable_models_setup::expose_observable_models_setup( observable_models_setup_submodule );

    auto observable_models_submodule = m.def_submodule( "observable_models" );
    observable_models::expose_observable_models( observable_models_submodule );

    auto observations_submodule = m.def_submodule( "observations" );
    observations::expose_observations( observations_submodule );

    auto observations_setup_submodule = m.def_submodule( "observations_setup" );
    observations_setup::expose_observations_setup( observations_setup_submodule );

    auto estimation_analysis_submodule = m.def_submodule( "estimation_analysis" );
    estimation_analysis::expose_estimation_analysis_estimator( estimation_analysis_submodule );
    estimation_analysis::expose_estimation_analysis_ephemeris_fit( estimation_analysis_submodule );
    estimation_analysis::expose_estimation_analysis( estimation_analysis_submodule );


};

}  // namespace estimation
}  // namespace tudatpy
