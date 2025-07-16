/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_observation.h"

#include <pybind11/functional.h>

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"
#include "tudat/simulation/estimation_setup/processOdfFile.h"
#include "tudat/simulation/estimation_setup/processTrackingTxtFile.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"

namespace tss = tudat::simulation_setup;
namespace tom = tudat::observation_models;
namespace tuc = tudat::unit_conversions;
namespace ti = tudat::interpolators;



namespace tudatpy
{
namespace numerical_simulation
{
namespace estimation_setup
{
namespace observation
{

void expose_observation_setup( py::module& m )
{


}


}  // namespace observation
}  // namespace estimation_setup
}  // namespace numerical_simulation
}  // namespace tudatpy