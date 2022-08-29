/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/acceleration_partials/aerodynamicAccelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

//! Function for updating partial w.r.t. the bodies' positions
void AerodynamicAccelerationPartial::update( const double currentTime )
{
    Eigen::Vector6d nominalState = vehicleStateGetFunction_( );
    Eigen::Vector6d perturbedState;

    // Compute state partial by numerical difference
    Eigen::Vector3d upperturbedAcceleration, downperturbedAcceleration;
    for( unsigned int i = 0; i < 6; i++ )
    {
        // Perturb state upwards
        perturbedState = nominalState;
        perturbedState( i ) += bodyStatePerturbations_( i );

        // Update environment/acceleration to perturbed state.
        flightConditions_->resetCurrentTime( );
        aerodynamicAcceleration_->resetCurrentTime( );
        vehicleStateSetFunction_( perturbedState );
        flightConditions_->updateConditions( currentTime );
        aerodynamicAcceleration_->updateMembers( currentTime );

        // Retrieve perturbed acceleration.
        upperturbedAcceleration = aerodynamicAcceleration_->getAcceleration( );

        // Perturb state downwards
        perturbedState = nominalState;
        perturbedState( i ) -= bodyStatePerturbations_( i );

        // Update environment/acceleration to perturbed state.
        flightConditions_->resetCurrentTime( );
        aerodynamicAcceleration_->resetCurrentTime( );
        vehicleStateSetFunction_( perturbedState );
        flightConditions_->updateConditions( currentTime );
        aerodynamicAcceleration_->updateMembers( currentTime );

        // Retrieve perturbed acceleration.
        downperturbedAcceleration = aerodynamicAcceleration_->getAcceleration( );

        // Compute partial
        currentAccelerationStatePartials_.block( 0, i, 3, 1 ) =
                ( upperturbedAcceleration - downperturbedAcceleration ) / ( 2.0 * bodyStatePerturbations_( i ) );
    }

    // Reset environment/acceleration mode to nominal conditions
    flightConditions_->resetCurrentTime( );
    aerodynamicAcceleration_->resetCurrentTime( );

    vehicleStateSetFunction_( nominalState );
    flightConditions_->updateConditions( currentTime );
    aerodynamicAcceleration_->updateMembers( currentTime );

    currentTime_ = currentTime;
}

} // namespace acceleration_partials

} // namespace tudat
