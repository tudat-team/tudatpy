//
// Created by Filippo Oggionni on 20/07/21.
//

#ifndef TUDATBUNDLE_EXPOSE_MASS_RATE_SETUP_H
#define TUDATBUNDLE_EXPOSE_MASS_RATE_SETUP_H

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
#include "tudat/simulation/propagation_setup/createThrustModelGuidance.h"
#include "tudat/simulation/propagation_setup/createTorqueModel.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/simulation/propagation_setup/environmentUpdater.h"
#include "tudat/simulation/propagation_setup/propagationCR3BPFullProblem.h"
#include "tudat/simulation/propagation_setup/propagationLambertTargeterFullProblem.h"
#include "tudat/simulation/propagation_setup/propagationOutput.h"
#include "tudat/simulation/propagation_setup/propagationOutputSettings.h"
#include "tudat/simulation/propagation_setup/propagationPatchedConicFullProblem.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/simulation/propagation_setup/propagationTermination.h"
#include "tudat/simulation/propagation_setup/propagationTerminationSettings.h"
#include "tudat/simulation/propagation_setup/setNumericallyIntegratedStates.h"
#include "tudat/simulation/propagation_setup/thrustSettings.h"
#include "tudat/simulation/propagation_setup/torqueSettings.h"

namespace tudatpy {
namespace simulation {
namespace propagation_setup {

    void expose_mass_rate_setup(py::module &m);

}// namespace propagation_setup
}// namespace simulation
}// namespace tudatpy

#endif //TUDATBUNDLE_EXPOSE_MASS_RATE_SETUP_H
