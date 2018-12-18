/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_FULLPROPAGATIONLAMBERTTARGETER_H
#define TUDAT_FULLPROPAGATIONLAMBERTTARGETER_H

#include <string>
#include <boost/bind.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationSettings.h"
#include "Tudat/Astrodynamics/ReferenceFrames/aerodynamicAngleCalculator.h"
#include "Tudat/Astrodynamics/Propagators/nBodyCowellStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/nBodyEnckeStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/nBodyGaussKeplerStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/nBodyGaussModifiedEquinoctialStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/nBodyUnifiedStateModelQuaternionsStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/nBodyUnifiedStateModelModifiedRodriguesParametersStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/nBodyUnifiedStateModelExponentialMapStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/rotationalMotionStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/rotationalMotionQuaternionsStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/rotationalMotionModifiedRodriguesParametersStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/rotationalMotionExponentialMapStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/bodyMassStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/customStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/stateDerivativeCircularRestrictedThreeBodyProblem.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include "Tudat/Astrodynamics/MissionSegments/lambertTargeter.h"
#include "Tudat/Astrodynamics/MissionSegments/lambertTargeterIzzo.h"
#include "Tudat/Astrodynamics/MissionSegments/lambertRoutines.h"

namespace tudat
{

namespace propagators
{

//! Function to directly setup CR3BP bodyMap
simulation_setup::NamedBodyMap setupBodyMapLambertTargeter(
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate);


//! Function to directly setup CR3BP acceleration map
basic_astrodynamics::AccelerationMap setupAccelerationMapLambertTargeter(
        const std::string& nameCentralBody,
        const std::string& nameBodyToPropagate,
        const std::vector< std::string >& bodiesToPropagate,
        const std::vector< std::string >& centralBodies,
        const simulation_setup::NamedBodyMap& bodyMap );


//! Function to propagate the full dynamics problem and compare it with Lambert targeter
void propagateLambertTargeterAndFullProblem(
        const Eigen::Vector3d& cartesianPositionAtDeparture,
        const Eigen::Vector3d& cartesianPositionAtArrival,
        const double timeOfFlight,
        simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::vector< std::string >& bodiesToPropagate,
        const std::vector< std::string >& centralBodies,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::map< double, Eigen::Vector6d >& lambertTargeterResult,
        const std::map< double, Eigen::Vector6d >& fullProblemResult);

std::pair< Eigen::Vector6d, Eigen::Vector6d > getDifferenceFullPropagationWrtLambertTargeterAtDepartureAndArrival(
        const Eigen::Vector3d& cartesianPositionAtDeparture,
        const Eigen::Vector3d& cartesianPositionAtArrival,
        const double timeOfFlight,
        simulation_setup::NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::vector< std::string >& bodiesToPropagate,
        const std::vector< std::string >& centralBody,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings);


} // namespace propagators

} // namespace tudat

#endif // TUDAT_FULLPROPAGATIONLAMBERTTARGETER_H
