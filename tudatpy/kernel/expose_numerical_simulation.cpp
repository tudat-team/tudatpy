/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "docstrings.h"
#include "scalarTypes.h"

#include "expose_numerical_simulation.h"

#include "expose_numerical_simulation/expose_environment_setup.h"
#include "expose_numerical_simulation/expose_estimation_setup.h"
#include "expose_numerical_simulation/expose_propagation_setup.h"

#include "expose_numerical_simulation/expose_environment.h"
#include "expose_numerical_simulation/expose_estimation.h"
#include "expose_numerical_simulation/expose_propagation.h"

#include "tudat/basics/timeType.h"
#include "tudat/astro/basic_astro/dateTime.h"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace tp = tudat::propagators;
namespace tss = tudat::simulation_setup;
namespace tni = tudat::numerical_integrators;
namespace tep = tudat::estimatable_parameters;
namespace tom = tudat::observation_models;
namespace tba = tudat::basic_astrodynamics;

namespace tudatpy {


namespace numerical_simulation {

void expose_numerical_simulation(py::module &m) {


    auto environment_submodule = m.def_submodule("environment");
    environment::expose_environment(environment_submodule);

    auto propagation_submodule = m.def_submodule("propagation");
    propagation::expose_propagation(propagation_submodule);

    auto estimation_submodule = m.def_submodule("estimation");
    estimation::expose_estimation2(estimation_submodule);
    estimation::expose_estimation(estimation_submodule);
    estimation::expose_single_observation_set(estimation_submodule);
    estimation::expose_observation_collection(estimation_submodule);
    estimation::expose_propagated_covariance(estimation_submodule);

    auto environment_setup_submodule = m.def_submodule("environment_setup");
    environment_setup::expose_environment_setup(environment_setup_submodule);

    auto propagation_setup_submodule = m.def_submodule("propagation_setup");
    propagation_setup::expose_propagation_setup(propagation_setup_submodule);

    auto estimation_setup_submodule = m.def_submodule("estimation_setup");
    estimation_setup::expose_estimation_setup(estimation_setup_submodule);

};

}// namespace numerical_simulation
}// namespace tudatpy
