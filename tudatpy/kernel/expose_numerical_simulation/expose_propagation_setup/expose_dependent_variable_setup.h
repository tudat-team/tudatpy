/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATBUNDLE_EXPOSE_DEPENDENT_VARIABLE_SETUP_H
#define TUDATBUNDLE_EXPOSE_DEPENDENT_VARIABLE_SETUP_H

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "tudat/simulation/propagation_setup/accelerationSettings.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"
#include "tudat/simulation/propagation_setup/createEnvironmentUpdater.h"
#include "tudat/simulation/propagation_setup/createMassRateModels.h"
#include "tudat/simulation/propagation_setup/createStateDerivativeModel.h"
#include "tudat/simulation/propagation_setup/createTorqueModel.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/simulation/propagation_setup/environmentUpdater.h"
#include "tudat/simulation/propagation_setup/propagationCR3BPFullProblem.h"
#include "tudat/simulation/propagation_setup/propagationLambertTargeterFullProblem.h"
#include "tudat/simulation/propagation_setup/propagationOutput.h"
#include "tudat/simulation/propagation_setup/propagationOutputSettings.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/simulation/propagation_setup/propagationTermination.h"
#include "tudat/simulation/propagation_setup/propagationTerminationSettings.h"
#include "tudat/simulation/propagation_setup/setNumericallyIntegratedStates.h"
#include "tudat/simulation/propagation_setup/torqueSettings.h"

namespace py = pybind11;

namespace tudatpy {
namespace numerical_simulation {
namespace propagation_setup {
namespace dependent_variable {

    void expose_dependent_variable_setup(py::module &m);

}// namespace dependent_variable
}// namespace propagation_setup
}// namespace numerical_simulation
}// namespace tudatpy

#endif //TUDATBUNDLE_EXPOSE_DEPENDENT_VARIABLE_SETUP_H
