/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <iostream>
#include "tudat/astro/LowThrustTrajectories/hybridMethodModel.h"
#include "tudat/math/quadrature/createNumericalQuadrature.h"
#include "tudat/astro/basic_astro/modifiedEquinoctialElementConversions.h"

namespace tudat
{
namespace low_thrust_trajectories
{

using namespace orbital_element_conversions;

//! Retrieve hybrid method acceleration model (including thrust and central gravity acceleration)
basic_astrodynamics::AccelerationMap HybridMethodModel::getLowThrustTrajectoryAccelerationMap( )
{
    throw std::runtime_error(
            "HybridMethodModel::getLowThrustTrajectoryAccelerationMap used the removed deprecated thrust direction/magnitude "
            "settings interface. Configure thrust with an engine model and body-fixed-direction rotation model instead." );
}

//! Propagate the spacecraft trajectory to time-of-flight.
Eigen::Vector6d HybridMethodModel::propagateTrajectory( )
{
    Eigen::Vector6d propagatedState = propagateTrajectory( 0.0, timeOfFlight_, stateAtDeparture_, initialSpacecraftMass_ );
    return propagatedState;
}

//! Propagate the spacecraft trajectory to a given time.
Eigen::Vector6d HybridMethodModel::propagateTrajectory( double initialTime,
                                                        double finalTime,
                                                        Eigen::Vector6d initialState,
                                                        double initialMass )
{
    throw std::runtime_error(
            "HybridMethodModel propagation used the removed deprecated thrust direction/magnitude settings interface. Configure thrust "
            "with an engine model and body-fixed-direction rotation model instead." );
}

//! Propagate the trajectory to set of epochs.
std::map< double, Eigen::Vector6d > HybridMethodModel::propagateTrajectory( std::vector< double > epochs,
                                                                            std::map< double, Eigen::Vector6d >& propagatedTrajectory )
{
    // Initialise propagated state.
    Eigen::Vector6d propagatedState = stateAtDeparture_;

    // Initialise mass of the spacecraft at departure.
    bodies_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );
    double currentMass = initialSpacecraftMass_;

    for( int epochIndex = 0; epochIndex < epochs.size( ); epochIndex++ )
    {
        double currentTime = epochs[ epochIndex ];
        if( epochIndex > 0 )
        {
            if( currentTime < epochs[ epochIndex - 1 ] )
            {
                throw std::runtime_error(
                        "Error when propagating trajectory with hybrid method, epochs at which the trajectory should be "
                        "computed are not in increasing order." );
            }
        }
        if( ( currentTime < 0.0 ) || ( currentTime > timeOfFlight_ ) )
        {
            throw std::runtime_error(
                    "Error when propagating trajectory with hybrid method, epochs at which the trajectory should be "
                    "computed are not constrained between 0.0 and timeOfFlight." );
        }

        if( epochIndex == 0 )
        {
            if( currentTime > 0.0 )
            {
                propagatedState = propagateTrajectory( 0.0, currentTime, propagatedState, currentMass );
                currentMass = bodies_[ bodyToPropagate_ ]->getBodyMass( );
            }
            propagatedTrajectory[ currentTime ] = propagatedState;
        }
        else
        {
            propagatedState = propagateTrajectory( epochs[ epochIndex - 1 ], currentTime, propagatedState, currentMass );
            currentMass = bodies_[ bodyToPropagate_ ]->getBodyMass( );
            propagatedTrajectory[ currentTime ] = propagatedState;
        }
    }

    bodies_[ centralBody_ ]->setConstantBodyMass( initialSpacecraftMass_ );

    return propagatedTrajectory;
}

//! Return the deltaV associated with the thrust profile of the trajectory.
double HybridMethodModel::computeDeltaV( )
{
    // Compute (constant) mass rate.
    double massRate = -maximumThrust_ / ( specificImpulse_ * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION );

    // Compute time during which the engine was switched on.
    double engineSwitchedOnDuration = ( massAtTimeOfFlight_ - initialSpacecraftMass_ ) / massRate;

    std::shared_ptr< numerical_quadrature::QuadratureSettings< double > > quadratureSettings_ =
            std::make_shared< numerical_quadrature::GaussianQuadratureSettings< double > >( 0.0, 16 );

    // Thrust acceleration function to use quadrature.
    // Define thrust acceleration as a function of time (to be integrated to compute the associated deltaV).
    std::function< double( const double ) > thrustAcceleration = [ = ]( const double currentTime ) {
        // Compute current mass.
        double currentMass = initialSpacecraftMass_ + massRate * currentTime;

        // Compute and return current thrust acceleration.
        double currentThrustAcceleration = maximumThrust_ / currentMass;

        return currentThrustAcceleration;
    };

    // Create numerical quadrature from quadrature settings.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( thrustAcceleration, quadratureSettings_, engineSwitchedOnDuration );

    // Compute deltaV analytically.
    double deltaV = -specificImpulse_ * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION *
            std::log( 1.0 + massRate / initialSpacecraftMass_ * engineSwitchedOnDuration );

    return deltaV;
}

}  // namespace low_thrust_trajectories
}  // namespace tudat
