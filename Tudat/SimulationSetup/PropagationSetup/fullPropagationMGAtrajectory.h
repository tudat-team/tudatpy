/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_FULLPROPAGATIONMGATRAJECTORY_H
#define TUDAT_FULLPROPAGATIONMGATRAJECTORY_H

#include <string>
#include <boost/bind.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationSettings.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include "Tudat/Astrodynamics/MissionSegments/lambertTargeter.h"
#include "Tudat/Astrodynamics/MissionSegments/lambertTargeterIzzo.h"
#include "Tudat/Astrodynamics/MissionSegments/lambertRoutines.h"
#include "Tudat/Astrodynamics/TrajectoryDesign/trajectory.h"

namespace tudat
{

namespace propagators
{

void fullPropagationMGA(
        simulation_setup::NamedBodyMap& bodyMap,
        const int numberOfLegs,
        const std::vector< std::string >& transferBodyOrder,
        const std::vector< std::string >& bodiesAndManoeuvresOrder,
        const std::vector< std::string >& centralBody,
        const std::vector< std::string >& bodyToPropagate,
        const std::vector< transfer_trajectories::TransferLegType>& legTypeVector,
        const std::vector< ephemerides::EphemerisPointer >& ephemerisVector,
        const Eigen::VectorXd& gravitationalParameterVector,
        const Eigen::VectorXd& trajectoryVariableVector,
        const double centralBodyGravitationalParameter,
        const Eigen::VectorXd& minimumPericenterRadiiVector,
        const Eigen::VectorXd& semiMajorAxesVector,
        const Eigen::VectorXd& eccentricitiesVector,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        std::map< int, std::map< double, Eigen::Vector6d > >& lambertTargeterResultForEachLeg,
        std::map< int, std::map< double, Eigen::Vector6d > >& fullProblemResultForEachLeg);


std::map< int, std::pair< Eigen::Vector6d, Eigen::Vector6d > > getDifferenceFullPropagationWrtLambertTargeterMGA(
        simulation_setup::NamedBodyMap& bodyMap,
        const int numberOfLegs,
        const std::vector< std::string >& transferBodyOrder,
        const std::vector< std::string >& bodiesAndManoeuvresOrder,
        const std::vector< std::string >& centralBody,
        const std::vector< std::string >& bodyToPropagate,
        const std::vector< transfer_trajectories::TransferLegType >& legTypeVector,
        const std::vector< ephemerides::EphemerisPointer >& ephemerisVector,
        const Eigen::VectorXd& gravitationalParameterVector,
        const Eigen::VectorXd& trajectoryVariableVector,
        const double centralBodyGravitationalParameter,
        const Eigen::VectorXd& minimumPericenterRadiiVector,
        const Eigen::VectorXd& semiMajorAxesVector,
        const Eigen::VectorXd& eccentricitiesVector,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings);


} // namespace propagators

} // namespace tudat

#endif // TUDAT_FULLPROPAGATIONMGATRAJECTORY_H
