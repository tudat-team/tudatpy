/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/aerodynamicAccelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{


//! Function for updating partial w.r.t. the bodies' positions
void AerodynamicAccelerationPartial::update( const double currentTime )
{
    Eigen::Vector6d nominalState = vehicleStateGetFunction_( );
    Eigen::Vector6d perturbedState;

    Eigen::Vector3d upperturbedAcceleration, downperturbedAcceleration;
    for( unsigned int i = 0; i < 6; i++ )
    {
        perturbedState = nominalState;
        perturbedState( i ) += bodyStatePerturbations_( i );

        flightConditions_->resetCurrentTime( TUDAT_NAN );
        aerodynamicAcceleration_->resetTime( TUDAT_NAN );
        vehicleStateSetFunction_( perturbedState );

        flightConditions_->updateConditions( currentTime );
        aerodynamicAcceleration_->updateMembers( currentTime );
        upperturbedAcceleration = aerodynamicAcceleration_->getAcceleration( );

        perturbedState = nominalState;
        perturbedState( i ) -= bodyStatePerturbations_( i );

        flightConditions_->resetCurrentTime( TUDAT_NAN );
        aerodynamicAcceleration_->resetTime( TUDAT_NAN );
        vehicleStateSetFunction_( perturbedState );

        flightConditions_->updateConditions( currentTime );
        aerodynamicAcceleration_->updateMembers( currentTime );
        downperturbedAcceleration = aerodynamicAcceleration_->getAcceleration( );

        currentAccelerationStatePartials_.block( 0, i, 3, 1 ) =
                ( upperturbedAcceleration - downperturbedAcceleration ) / ( 2.0 * bodyStatePerturbations_( i ) );
    }

    flightConditions_->resetCurrentTime( TUDAT_NAN );
    aerodynamicAcceleration_->resetTime( TUDAT_NAN );
    vehicleStateSetFunction_( nominalState );

    flightConditions_->updateConditions( currentTime );
    aerodynamicAcceleration_->updateMembers( currentTime );
}


} // namespace acceleration_partials

} // namespace tudat
