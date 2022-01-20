/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROPAGATIONLOWTHRUSTPROBLEM_H
#define TUDAT_PROPAGATIONLOWTHRUSTPROBLEM_H

#include "tudat/astro/low_thrust/lowThrustLeg.h"
#include "tudat/simulation/propagation_setup/accelerationSettings.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"
#include "tudat/simulation/environment_setup/body.h"

namespace tudat
{

namespace simulation_setup
{

void computeLowThrustLegSemiAnalyticalAndFullPropagation(
        const std::shared_ptr< low_thrust_trajectories::LowThrustLeg > lowThrustLeg,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >,
        std::shared_ptr< propagators::PropagatorSettings< double > > >& propagatorSettings,
        std::map< double, Eigen::VectorXd >& fullPropagationResults,
        std::map< double, Eigen::Vector6d >& semiAnalyticalResults,
        std::map< double, Eigen::VectorXd >& dependentVariablesHistory );


basic_astrodynamics::AccelerationMap retrieveLowThrustAccelerationMap(
        const std::shared_ptr< low_thrust_trajectories::LowThrustLeg > lowThrustLeg,
        const simulation_setup::SystemOfBodies& bodies,
        const std::string& bodyToPropagate,
        const std::string& centralBody,
        const std::function< double( const double ) > specificImpulseFunction,
        const double lowThrustLegInitialTime );


//! Define appropriate translational state propagator settings for the full propagation.
std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > createLowThrustTranslationalStatePropagatorSettings(
        const std::shared_ptr< low_thrust_trajectories::LowThrustLeg > lowThrustLeg,
        const std::string& bodyToPropagate,
        const std::string& centralBody,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave );


//! Define appropriate propagator settings for the full propagation.
std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >,
std::shared_ptr< propagators::PropagatorSettings< double > > > createLowThrustPropagatorSettings(
        const std::shared_ptr< low_thrust_trajectories::LowThrustLeg > lowThrustLeg,
        const double bodyMassAtMidPoint,
        const simulation_setup::SystemOfBodies& bodies,
        const std::string& bodyToPropagate,
        const std::string& centralBody,
        const std::function< double( const double ) > specificImpulseFunction,
        const basic_astrodynamics::AccelerationMap perturbingAccelerationsMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings,
        const std::shared_ptr< propagators::DependentVariableSaveSettings >& dependentVariablesToSave );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_PROPAGATIONLOWTHRUSTPROBLEM_H
