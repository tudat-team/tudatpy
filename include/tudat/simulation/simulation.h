/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SIMULATIONHEADER_H
#define TUDAT_SIMULATIONHEADER_H

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/basic_astro/keplerPropagator.h"
#include "tudat/astro/basic_astro/sphericalStateConversions.h"
#include "tudat/astro/basic_astro/stateRepresentationConversions.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"

#include "tudat/io/basicInputOutput.h"
#include "tudat/io/mapTextFileReader.h"

#include "tudat/basics/basicTypedefs.h"

#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/filters/createFilter.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"
#include "tudat/math/interpolators/createInterpolator.h"

#ifdef TUDAT_BUILD_WITH_SPICE_INTERFACE

#include "tudat/interface/spice/spiceInterface.h"

#endif

#include "tudat/simulation/environment/body.h"
#include "tudat/simulation/environment/createBodies.h"
#include "tudat/simulation/environment/defaultBodies.h"
#include "tudat/simulation/estimation/createEstimatableParameters.h"
#include "tudat/simulation/estimation/estimatableParameterSettings.h"
#include "tudat/simulation/propagation/accelerationSettings.h"
#include "tudat/simulation/propagation/propagationSettings.h"
#include "tudat/simulation/propagation/propagationOutputSettings.h"
#include "tudat/simulation/propagation/propagationTerminationSettings.h"
#include "tudat/simulation/estimation/createNumericalSimulator.h"
#include "tudat/simulation/propagation/thrustSettings.h"
#include "tudat/simulation/propagation/createMassRateModels.h"

#endif // TUDAT_SIMULATIONHEADER_H
