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
#include "Tudat/Astrodynamics/LowThrustDirectMethods/simsFlanaganLeg.h"
#include "Tudat/Mathematics/NumericalQuadrature/createNumericalQuadrature.h"

namespace tudat
{
namespace low_thrust_direct_methods
{

int SimsFlanaganLeg::convertTimeToLegSegment( double currentTime )
{
    int indexSegment;

    double segmentDurationForwardPropagation = timeAtMatchPoint_ / numberSegmentsForwardPropagation_;
    double segmentDurationBackwardPropagation = timeAtMatchPoint_ / numberSegmentsBackwardPropagation_;

    if ( currentTime <= timeAtMatchPoint_ )
    {
        indexSegment = currentTime / segmentDurationForwardPropagation;
    }
    else if ( currentTime == timeOfFlight_ )
    {
        indexSegment = numberSegments_ - 1;
    }
    else
    {
        indexSegment = numberSegmentsForwardPropagation_ + ( currentTime - timeAtMatchPoint_ ) / segmentDurationBackwardPropagation;
    }

    return indexSegment;
}

std::shared_ptr< simulation_setup::ThrustAccelerationSettings > SimsFlanaganLeg::getConstantThrustAccelerationSettingsPerSegment(
        unsigned int indexSegment )
{
    // Define (constant) thrust magnitude function.
    std::function< double( const double ) > thrustMagnitudeFunction = [ = ]( const double currentTime )
    {
        return maximumThrust_ * throttles_[ indexSegment ].norm();
    };

    // Define thrust magnitude settings from thrust magnitude function.
    std::shared_ptr< simulation_setup::FromFunctionThrustMagnitudeSettings > thrustMagnitudeSettings
            = std::make_shared< simulation_setup::FromFunctionThrustMagnitudeSettings >(
                thrustMagnitudeFunction, specificImpulseFunction_ );

    // Define thrust direction function (constant over one leg segment).
    std::function< Eigen::Vector3d( const double ) > thrustDirectionFunction = [ = ]( const double currentTime )
    {
        return throttles_[ indexSegment ].normalized();
    };


    // Define thrust direction settings from the direction of thrust acceleration retrieved from the shaping method.
    std::shared_ptr< simulation_setup::CustomThrustDirectionSettings > thrustDirectionSettings =
            std::make_shared< simulation_setup::CustomThrustDirectionSettings >( thrustDirectionFunction );

    // Define thrust acceleration settings.
    std::shared_ptr< simulation_setup::ThrustAccelerationSettings > thrustAccelerationSettings =
            std::make_shared< simulation_setup::ThrustAccelerationSettings >(
                thrustDirectionSettings, thrustMagnitudeSettings );

    return thrustAccelerationSettings;
}

basic_astrodynamics::AccelerationMap SimsFlanaganLeg::getAccelerationModelPerSegment( unsigned int indexSegment )
{

    // Define acceleration settings.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > accelerationsSettings;

    // Add point-mass gravitational acceleration from central body.
    accelerationsSettings[ centralBody_ ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );

    // Retrieve thrust acceleration settings.
    accelerationsSettings[ bodyToPropagate_ ].push_back( getConstantThrustAccelerationSettingsPerSegment( indexSegment ) );

    // Create acceleration map.
    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ bodyToPropagate_ ] = accelerationsSettings;

    // Create the acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap_, accelerationMap, std::vector< std::string >{ bodyToPropagate_ }, std::vector< std::string >{ centralBody_ } );

    return accelerationModelMap;

}


std::shared_ptr< simulation_setup::ThrustAccelerationSettings > SimsFlanaganLeg::getThrustAccelerationSettingsFullLeg( )
{

    // Define thrust magnitude function.
    std::function< double( const double ) > thrustMagnitudeFunction = [ = ]( const double currentTime )
    {
        int indexSegment = convertTimeToLegSegment( currentTime );

        return maximumThrust_ * throttles_[ indexSegment ].norm();
    };

    // Define thrust magnitude settings from thrust magnitude function.
    std::shared_ptr< simulation_setup::FromFunctionThrustMagnitudeSettings > thrustMagnitudeSettings
            = std::make_shared< simulation_setup::FromFunctionThrustMagnitudeSettings >(
                thrustMagnitudeFunction, specificImpulseFunction_ );

    // Define thrust direction function.
    std::function< Eigen::Vector3d( const double ) > thrustDirectionFunction = [ = ]( const double currentTime )
    {
        int indexSegment = convertTimeToLegSegment( currentTime );

        return throttles_[ indexSegment ].normalized();
    };


    // Define thrust direction settings from the direction of thrust acceleration retrieved from the shaping method.
    std::shared_ptr< simulation_setup::CustomThrustDirectionSettings > thrustDirectionSettings =
            std::make_shared< simulation_setup::CustomThrustDirectionSettings >( thrustDirectionFunction );

    // Define thrust acceleration settings.
    std::shared_ptr< simulation_setup::ThrustAccelerationSettings > thrustAccelerationSettings =
            std::make_shared< simulation_setup::ThrustAccelerationSettings >(
                thrustDirectionSettings, thrustMagnitudeSettings );

    return thrustAccelerationSettings;
}

basic_astrodynamics::AccelerationMap SimsFlanaganLeg::getAccelerationModelFullLeg( )
{
    // Define acceleration settings.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > accelerationsSettings;

    // Add point-mass gravitational acceleration from central body.
    accelerationsSettings[ centralBody_ ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );

    // Retrieve thrust acceleration settings.
    accelerationsSettings[ bodyToPropagate_ ].push_back( getThrustAccelerationSettingsFullLeg( ) );

    // Create acceleration map.
    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ bodyToPropagate_ ] = accelerationsSettings;

    // Create the acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap_, accelerationMap, std::vector< std::string >{ bodyToPropagate_ }, std::vector< std::string >{ centralBody_ } );

    return accelerationModelMap;
}

//! Propagate the mass from departure to a given segment.
double SimsFlanaganLeg::propagateMassToSegment( int indexSegment )
{
    double currentMass = initialSpacecraftMass_;

    for ( unsigned int currentSegment = 0 ; currentSegment < indexSegment ; currentSegment++ )
    {

        // Compute time at half of the current segment.
        double currentTime = timesAtNodes_[ currentSegment ] +
                ( timesAtNodes_[ currentSegment + 1 ] - timesAtNodes_[ currentSegment ] ) / 2.0;

        // Compute segment duration.
        double segmentDuration;
        if ( currentSegment < numberSegmentsForwardPropagation_ )
        {
            segmentDuration = timeAtMatchPoint_ / numberSegmentsForwardPropagation_;
        }
        else
        {
            segmentDuration = timeAtMatchPoint_ / numberSegmentsBackwardPropagation_;
        }

        // Compute current deltaV vector.
        Eigen::Vector3d currentDeltaV;
        for ( int i = 0 ; i < 3 ; i++ )
        {
            currentDeltaV[ i ] = maximumThrust_ / currentMass * segmentDuration * throttles_[ currentSegment ][ i ];
        }

        // Update mass of the spacecraft.
        currentMass *= std::exp( - currentDeltaV.norm() /
                          ( specificImpulseFunction_( currentTime ) * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION ) );
    }

    // Return mass of the spacecraft at half of the required segment.
    return currentMass;
}


//! Propagate the spacecraft trajectory from departure to match point (forward propagation).
void SimsFlanaganLeg::propagateForwardFromDepartureToMatchPoint( )
{
    // Initialise current state at the leg departure.
    Eigen::Vector6d currentState = stateAtDeparture_;

    currentState = propagateTrajectoryForward( 0.0, timeAtMatchPoint_, stateAtDeparture_, timeAtMatchPoint_ / numberSegmentsForwardPropagation_ );

    // Set state vector at match point once the forward propagation is over.
    stateAtMatchPointFromForwardPropagation_ = currentState;

    // Set mass of the spacecraft at match point after the forward propagation.
    massAtMatchPointFromForwardPropagation_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass( );


}

//! Propagate the spacecraft trajectory from arrival to match point (backward propagation).
void SimsFlanaganLeg::propagateBackwardFromArrivalToMatchPoint( )
{
    // Initialise current state at the leg arrival.
    Eigen::Vector6d currentState = stateAtArrival_;

    // Initialise mass of the spacecraft at the end of the leg.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( propagateMassToSegment( numberSegments_ ) );

    currentState = propagateTrajectoryBackward( timeOfFlight_, timeAtMatchPoint_, stateAtArrival_, timeAtMatchPoint_ / numberSegmentsBackwardPropagation_ );

    // Set state vector at match point once the backward propagation is over.
    stateAtMatchPointFromBackwardPropagation_ = currentState;

    // Set mass of the spacecraft at match point after the backward propagation.
    massAtMatchPointFromBackwardPropagation_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass( );

}


//! Propagate the trajectory inside one segment.
Eigen::Vector6d SimsFlanaganLeg::propagateInsideForwardSegment( double initialTime, double finalTime, double segmentDuration, Eigen::Vector6d initialState )
{

    // Compute time elapsed since start of the current leg segment.
    double timeElapsedCurrentSegment = finalTime - initialTime;

    // Compute current segment.
    int currentSegment = convertTimeToLegSegment( initialTime );

    // Initialise propagated state.
    Eigen::Vector6d propagatedState = initialState;

    // Update mass of the vehicle.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( propagateMassToSegment( currentSegment ) );


    // Propagate from start of the current leg segment to the required propagation final time.
    if ( ( finalTime <= ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ) )
         || ( initialTime > ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ) ) )
    {
        // Directly propagate the current state over the time elapsed since the start of the current leg segment.
        propagatedState = orbital_element_conversions::convertKeplerianToCartesianElements(
                    orbital_element_conversions::propagateKeplerOrbit(
                        orbital_element_conversions::convertCartesianToKeplerianElements(
                            propagatedState, centralBodyGravitationalParameter_ ),
                        timeElapsedCurrentSegment, centralBodyGravitationalParameter_ ), centralBodyGravitationalParameter_ );
    }
    else
    {
        // First propagate the current state to half of the current segment.
        if ( ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ) - initialTime > 0.0 )
        {
            propagatedState = orbital_element_conversions::convertKeplerianToCartesianElements(
                        orbital_element_conversions::propagateKeplerOrbit(
                            orbital_element_conversions::convertCartesianToKeplerianElements(
                                propagatedState, centralBodyGravitationalParameter_ ),
                            ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ) - initialTime, centralBodyGravitationalParameter_ ),
                        centralBodyGravitationalParameter_ );
        }

        // Compute the deltaV that needs to be applied at half of the current segment.
        Eigen::Vector3d deltaVvector;
        for ( unsigned int i = 0 ; i < 3 ; i++ )
        {
            deltaVvector[ i ] = maximumThrust_ / bodyMap_[ bodyToPropagate_ ]->getBodyMass()
                    * segmentDuration * throttles_[ currentSegment ][ i ];
        }

        // Add the deltaV vector to the state propagated until half of the current segment.
        propagatedState.segment( 3, 3 ) += deltaVvector;

        // Propagate the updated state from half of the current segment to required propagation final time.
        if ( finalTime - ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ) > 0.0  )
        {
            propagatedState = orbital_element_conversions::convertKeplerianToCartesianElements(
                        orbital_element_conversions::propagateKeplerOrbit(
                            orbital_element_conversions::convertCartesianToKeplerianElements(
                                propagatedState, centralBodyGravitationalParameter_ ),
                            finalTime - ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ), centralBodyGravitationalParameter_ ),
                        centralBodyGravitationalParameter_ );
        }


        // Update mass of the vehicle.
        bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( propagateMassToSegment( currentSegment + 1 ) );

    }

    return propagatedState;
}


//! Propagate the trajectory inside one segment.
Eigen::Vector6d SimsFlanaganLeg::propagateInsideBackwardSegment( double initialTime, double finalTime, double segmentDuration, Eigen::Vector6d initialState )
{
    // Compute time elapsed since start of the current leg segment.
    double timeElapsedCurrentSegment = finalTime - initialTime;

    // Compute current segment.
    int currentSegment = convertTimeToLegSegment( finalTime );

    // Initialise propagated state.
    Eigen::Vector6d propagatedState = initialState;


    // Propagate from start of the current leg segment to the required propagation final time.
    if ( ( finalTime >= ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ) )
         || ( initialTime < ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ) ) )
    {
        // Directly propagate the current state over the time elapsed since the start of the current leg segment.
        propagatedState = orbital_element_conversions::convertKeplerianToCartesianElements(
                    orbital_element_conversions::propagateKeplerOrbit(
                        orbital_element_conversions::convertCartesianToKeplerianElements(
                            propagatedState, centralBodyGravitationalParameter_ ),
                        timeElapsedCurrentSegment, centralBodyGravitationalParameter_ ), centralBodyGravitationalParameter_ );
    }
    else
    {
        // First propagate the current state to half of the current segment.
        if ( initialTime - ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ) > 0.0 )
        {
            propagatedState = orbital_element_conversions::convertKeplerianToCartesianElements(
                        orbital_element_conversions::propagateKeplerOrbit(
                            orbital_element_conversions::convertCartesianToKeplerianElements(
                                propagatedState, centralBodyGravitationalParameter_ ),
                            ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ) - initialTime,
                            centralBodyGravitationalParameter_ ), centralBodyGravitationalParameter_ );
        }

        // Update mass of the vehicle.
        bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( propagateMassToSegment( currentSegment ) );

        // Compute the deltaV that needs to be applied at half of the current segment.
        Eigen::Vector3d deltaVvector;
        for ( unsigned int i = 0 ; i < 3 ; i++ )
        {
            deltaVvector[ i ] = maximumThrust_ / bodyMap_[ bodyToPropagate_ ]->getBodyMass()
                    * segmentDuration * throttles_[ currentSegment ][ i ];
        }

        // Add the deltaV vector to the state propagated until half of the current segment.
        propagatedState.segment( 3, 3 ) -= deltaVvector;

        // Propagate the updated state from half of the current segment to required propagation final time.
        if ( ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ) - finalTime > 0.0  )
        {
            propagatedState = orbital_element_conversions::convertKeplerianToCartesianElements(
                        orbital_element_conversions::propagateKeplerOrbit(
                            orbital_element_conversions::convertCartesianToKeplerianElements(
                                propagatedState, centralBodyGravitationalParameter_ ),
                            finalTime - ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ), centralBodyGravitationalParameter_ ),
                        centralBodyGravitationalParameter_ );
        }



    }

    return propagatedState;
}


//! Propagate the trajectory forward to given time.
Eigen::Vector6d SimsFlanaganLeg::propagateTrajectoryForward( double initialTime, double finalTime, Eigen::Vector6d initialState, double segmentDuration )
{
    // Compute index of leg segment which corresponds to the propagation initial time.
    int initialSegment = convertTimeToLegSegment( initialTime );

    // Compute index of leg segment which corresponds to the propagation final time.
    int finalSegment = convertTimeToLegSegment( finalTime );

    // Initialise propagated state.
    Eigen::Vector6d propagatedState = initialState;

    // Initialise current segment.
    int currentSegment = initialSegment;


    if ( initialSegment == finalSegment )
    {
        // Directly propagate to final time inside the current segment.
        if ( finalTime - initialTime > 0 )
        {
            propagatedState = propagateInsideForwardSegment( initialTime, finalTime, segmentDuration, propagatedState );
        }
    }
    else
    {
        // First propagate to beginning of next segment.
        propagatedState = propagateInsideForwardSegment( initialTime, timesAtNodes_[ initialSegment + 1 ], segmentDuration, propagatedState );

        // Update current segment.
        currentSegment++;

        // Compute number of segments the trajectory should be propagated through before reaching the segment corresponding
        // to the final propagation time.
        int propagatedSegments = ( finalTime - timesAtNodes_[ currentSegment ] ) / segmentDuration;

        // Propagate through these segments.
        for ( int i = 0 ; i < propagatedSegments ; i++ )
        {
            propagatedState = propagateInsideForwardSegment( timesAtNodes_[ currentSegment ], timesAtNodes_[ currentSegment + 1 ],
                    segmentDuration, propagatedState );
            currentSegment++;
        }

        // Propagate inside last segment.
        if ( ( finalTime - timesAtNodes_[ currentSegment ] ) > 0.0 )
        {
            propagatedState = propagateInsideForwardSegment( timesAtNodes_[ currentSegment ], finalTime, segmentDuration, propagatedState );
            currentSegment++;
        }

    }

    return propagatedState;

}


//! Propagate the trajectory backward to given time.
Eigen::Vector6d SimsFlanaganLeg::propagateTrajectoryBackward( double initialTime, double finalTime, Eigen::Vector6d initialState, double segmentDuration )
{
    // Compute index of leg segment which corresponds to the propagation initial time.
    int initialSegment = convertTimeToLegSegment( initialTime );

    // Compute index of leg segment which corresponds to the propagation final time.
    int finalSegment = convertTimeToLegSegment( finalTime );

    // Initialise propagated state.
    Eigen::Vector6d propagatedState = initialState;

    // Initialise current segment.
    int currentSegment = initialSegment;


    if ( initialSegment == finalSegment )
    {
        // Directly propagate to final time inside the current segment.
        if ( finalTime - initialTime < 0 )
        {
            propagatedState = propagateInsideBackwardSegment( initialTime, finalTime, segmentDuration, propagatedState );
        }
    }
    else
    {

        // First propagate to beginning of next segment.
        propagatedState = propagateInsideBackwardSegment( initialTime, timesAtNodes_[ initialSegment ], segmentDuration, propagatedState );

        // Update current segment.
        currentSegment--;

        // Compute number of segments the trajectory should be propagated through before reaching the segment corresponding
        // to the final propagation time.
        int propagatedSegments = ( timesAtNodes_[ currentSegment + 1 ] - finalTime ) / segmentDuration;

        // Propagate through these segments.
        for ( int i = 0 ; i < propagatedSegments ; i++ )
        {
            propagatedState = propagateInsideBackwardSegment( timesAtNodes_[ currentSegment + 1 ], timesAtNodes_[ currentSegment ],
                    segmentDuration, propagatedState );
            currentSegment--;
        }

        // Propagate inside last segment.
        if ( ( timesAtNodes_[ currentSegment + 1 ] - finalTime ) > 0.0 )
        {
            propagatedState = propagateInsideBackwardSegment( timesAtNodes_[ currentSegment + 1 ], finalTime, segmentDuration, propagatedState );
            currentSegment--;
        }

    }

    return propagatedState;
}


//! Propagate the trajectory to set of epochs.
std::map< double, Eigen::Vector6d > SimsFlanaganLeg::propagateTrajectoryForward(
        std::vector< double > epochs, std::map< double, Eigen::Vector6d >& propagatedTrajectory, Eigen::Vector6d initialState,
        double initialMass, double initialTime, double segmentDuration )
{
        // Initialise propagated state.
        Eigen::Vector6d propagatedState = initialState;

        // Initialise mass of the spacecraft at departure.
        bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialMass );


        for ( int epochIndex = 0 ; epochIndex < epochs.size( ) ; epochIndex++ )
        {
            double currentTime = epochs[ epochIndex ];
            if ( epochIndex > 0 )
            {
                if ( currentTime < epochs[ epochIndex - 1 ] )
                {
                    throw std::runtime_error( "Error when propagating trajectory in Sims-Flanagan, epochs at which the trajectory should be "
                                              "computed are not in increasing order." );
                }
            }
            if ( ( currentTime < 0.0 ) || ( currentTime > timeOfFlight_ ) )
            {
                throw std::runtime_error( "Error when propagating trajectory in Sims-Flanagan, epochs at which the trajectory should be "
                                          "computed are not constrained between 0.0 and timeOfFlight." );
            }


            if ( epochIndex == 0 )
            {
                if ( currentTime > initialTime )
                {
                    propagatedState = propagateTrajectoryForward( initialTime, currentTime, propagatedState,
                                                                  timeAtMatchPoint_ / numberSegmentsForwardPropagation_ );
                }
                propagatedTrajectory[ currentTime ] = propagatedState;
            }
            else
            {
                propagatedState = propagateTrajectoryForward( epochs[ epochIndex - 1 ], currentTime,
                        propagatedState, segmentDuration );
                propagatedTrajectory[ currentTime ] = propagatedState;
            }

        }

        bodyMap_[ centralBody_ ]->setConstantBodyMass( initialMass );

        return propagatedTrajectory;
}



//! Propagate the trajectory to set of epochs.
std::map< double, Eigen::Vector6d > SimsFlanaganLeg::propagateTrajectoryBackward(
        std::vector< double > epochs, std::map< double, Eigen::Vector6d >& propagatedTrajectory,
        Eigen::Vector6d initialState, double initialMass, double initialTime, double segmentDuration )
{

        // Initialise propagated state.
        Eigen::Vector6d propagatedState = initialState;

        // Initialise mass of the spacecraft at departure.
        bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialMass );

        for ( int epochIndex = 0 ; epochIndex < epochs.size() ; epochIndex++ )
        {
            double currentTime = epochs[ epochIndex ];
            if ( epochIndex > 0 )
            {
                if ( currentTime > epochs[ epochIndex - 1 ] )
                {
                    throw std::runtime_error( "Error when propagating trajectory backward in Sims-Flanagan, epochs at which the trajectory should be "
                                              "computed are not in decreasing order." );
                }
            }
            if ( ( currentTime < 0.0 ) || ( currentTime > timeOfFlight_ ) )
            {
                throw std::runtime_error( "Error when propagating trajectory backward in Sims-Flanagan, epochs at which the trajectory should be "
                                          "computed are not constrained between 0.0 and timeOfFlight." );
            }


            if ( epochIndex == 0 )
            {
                if ( currentTime < initialTime )
                {
                    propagatedState = propagateTrajectoryBackward( initialTime, currentTime, propagatedState,
                                                                   timeAtMatchPoint_ / numberSegmentsForwardPropagation_ );
                }
                propagatedTrajectory[ currentTime ] = propagatedState;
            }
            else
            {
                propagatedState = propagateTrajectoryBackward( epochs[ epochIndex - 1 ], currentTime, propagatedState, segmentDuration );
                propagatedTrajectory[ currentTime ] = propagatedState;
            }

        }

        return propagatedTrajectory;
}


} // namespace low_thrust_direct_methods
} // namespace tudat
