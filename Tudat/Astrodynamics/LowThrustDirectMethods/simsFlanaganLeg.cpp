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
        return maximumThrust_ /*/ bodyMap_[ bodyToPropagate_ ]->getBodyMass()*/ * throttles_[ indexSegment ].norm();
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

//    double segmentDurationForwardPropagation = timeAtMatchPoint_ / numberSegmentsForwardPropagation_;
//    double segmentDurationBackwardPropagation = timeAtMatchPoint_ / numberSegmentsBackwardPropagation_;

    // Define thrust magnitude function.
    std::function< double( const double ) > thrustMagnitudeFunction = [ = ]( const double currentTime )
    {
        int indexSegment = convertTimeToLegSegment( currentTime );
//        if ( currentTime < timeAtMatchPoint_ )
//        {
//            indexSegment = currentTime / segmentDurationForwardPropagation;
//        }
//        else if ( currentTime == timeOfFlight_ )
//        {
//            indexSegment = numberSegments_ - 1;
//        }
//        else
//        {
//            indexSegment = numberSegmentsForwardPropagation_ + ( currentTime - timeAtMatchPoint_ ) / segmentDurationBackwardPropagation;
//        }

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
//        if ( currentTime < timeAtMatchPoint_ )
//        {
//            indexSegment = currentTime / segmentDurationForwardPropagation;
//        }
//        else if ( currentTime == timeOfFlight_ )
//        {
//            indexSegment = numberSegments_ - 1;
//        }
//        else
//        {
//            indexSegment = numberSegmentsForwardPropagation_ + ( currentTime - timeAtMatchPoint_ ) / segmentDurationBackwardPropagation;
//        }

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

//! Propagate the mass from departure to arrival for low order solution (return mass of the spacecraft at end of the leg).
double SimsFlanaganLeg::propagateMassToLegArrival( )
{
    double currentMass = initialSpacecraftMass_;

    for ( unsigned int currentSegment = 0 ; currentSegment < numberSegments_ ; currentSegment++ )
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

    // Return mass of the spacecraft at the end of the leg.
    return currentMass;
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

//    Eigen::Vector6d testWithoutThrust = orbital_element_conversions::convertKeplerianToCartesianElements(
//                orbital_element_conversions::propagateKeplerOrbit(
//                orbital_element_conversions::convertCartesianToKeplerianElements( currentState, centralBodyGravitationalParameter_ ),
//                timeOfFlight_ / 2.0, centralBodyGravitationalParameter_ ), centralBodyGravitationalParameter_ );


//    // Loop over the leg segments of the forward propagation.
//    for ( unsigned int currentSegment = 0 ; currentSegment < numberSegmentsForwardPropagation_ ; currentSegment++ )
//    {
//        // Forward propagation over one segment of the leg.
//        currentState = propagateForwardSegment( currentSegment, currentState );
//    }
    currentState = propagateTrajectoryForward( 0.0, timeAtMatchPoint_, stateAtDeparture_, timeAtMatchPoint_ / numberSegmentsForwardPropagation_ );

    // Set state vector at match point once the forward propagation is over.
    stateAtMatchPointFromForwardPropagation_ = currentState;

    // Set mass of the spacecraft at match point after the forward propagation.
    massAtMatchPointFromForwardPropagation_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass( );

    std::cout << "END OF FORWARD PROPAGATION: " << "\n\n";
    std::cout << "state at match point: " << stateAtMatchPointFromForwardPropagation_ << "\n\n";
    std::cout << "mass at match point: " << massAtMatchPointFromForwardPropagation_ << "\n\n";

//    std::cout << "test without thrust: " << testWithoutThrust << "\n\n";

}

//! Propagate the spacecraft trajectory from arrival to match point (backward propagation).
void SimsFlanaganLeg::propagateBackwardFromArrivalToMatchPoint( )
{
    // Initialise current state at the leg arrival.
    Eigen::Vector6d currentState = stateAtArrival_;

    // Initialise mass of the spacecraft at the end of the leg.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( propagateMassToLegArrival( ) );

//    Eigen::Vector6d testWithoutThrust = orbital_element_conversions::convertKeplerianToCartesianElements(
//                orbital_element_conversions::propagateKeplerOrbit(
//                orbital_element_conversions::convertCartesianToKeplerianElements( currentState, centralBodyGravitationalParameter_ ),
//                - timeOfFlight_ / 2.0, centralBodyGravitationalParameter_ ), centralBodyGravitationalParameter_ );

//    // Loop over the leg segments of the backward propagation.
//    for ( unsigned int currentSegment = numberSegments_ - 1 ;
//          currentSegment >= numberSegmentsForwardPropagation_ ; currentSegment-- )
//    {
//        // Backward propagation over one segment of the leg.
//        currentState = propagateBackwardSegment( currentSegment, currentState );
//    }
    currentState = propagateTrajectoryBackward( timeOfFlight_, timeAtMatchPoint_, stateAtArrival_, timeAtMatchPoint_ / numberSegmentsBackwardPropagation_ );

    // Set state vector at match point once the backward propagation is over.
    stateAtMatchPointFromBackwardPropagation_ = currentState;

    // Set mass of the spacecraft at match point after the backward propagation.
    massAtMatchPointFromBackwardPropagation_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass( );


    std::cout << "END OF BACKWARD PROPAGATION: " << "\n\n";
    std::cout << "state at match point: " << stateAtMatchPointFromBackwardPropagation_ << "\n\n";
    std::cout << "mass at match point: " << massAtMatchPointFromBackwardPropagation_ << "\n\n";

//    std::cout << "test without thrust: " << testWithoutThrust << "\n\n";

//    std::cout << "mass at the end of the leg (after propagation): " << propagateMassToLegArrival( ) << "n\n";

}


//! Propagate the trajectory inside one segment (low order solution).
Eigen::Vector6d SimsFlanaganLeg::propagateInsideSegment( double initialTime, double finalTime, double segmentDuration, Eigen::Vector6d initialState )
{

    // Compute time elapsed since start of the current leg segment.
    double timeElapsedCurrentSegment = finalTime - initialTime;

    // Compute current segment.
    int currentSegment = convertTimeToLegSegment( initialTime );
    std::cout << "mid-segment: " << ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ) << "\n\n";
    std::cout << "current segment: " << currentSegment << "\n\n";

    // Initialise propagated state.
    Eigen::Vector6d propagatedState = initialState;

    // Update mass of the vehicle.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( propagateMassToSegment( currentSegment ) );


    // Propagate from start of the current leg segment to the required propagation final time.
    if ( /*timeElapsedCurrentSegment <= segmentDuration */
         ( finalTime <= ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ) )
         || ( initialTime > ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ) ) )
    {
        std::cout << "no DELTAV: " << "\n\n";
        // Directly propagate the current state over the time elapsed since the start of the current leg segment.
        propagatedState = orbital_element_conversions::convertKeplerianToCartesianElements(
                    orbital_element_conversions::propagateKeplerOrbit(
                        orbital_element_conversions::convertCartesianToKeplerianElements(
                            propagatedState, centralBodyGravitationalParameter_ ),
                        timeElapsedCurrentSegment, centralBodyGravitationalParameter_ ), centralBodyGravitationalParameter_ );
    }
    else
    {
        std::cout << "detection DELTAV" << "\n\n";

        // First propagate the current state to half of the current segment.
        if ( ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ) - initialTime > 0.0 )
        {
            propagatedState = orbital_element_conversions::convertKeplerianToCartesianElements(
                        orbital_element_conversions::propagateKeplerOrbit(
                            orbital_element_conversions::convertCartesianToKeplerianElements(
                                propagatedState, centralBodyGravitationalParameter_ ),
                            ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ) - initialTime /*segmentDuration / 2.0*/, centralBodyGravitationalParameter_ ), centralBodyGravitationalParameter_ );
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
                            finalTime - ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ) /*timeElapsedCurrentSegment - segmentDuration / 2.0*/, centralBodyGravitationalParameter_ ),
                        centralBodyGravitationalParameter_ );
        }

        // Update mass of the vehicle.
//        bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( bodyMap_[ bodyToPropagate_ ]->getBodyMass() *
//                    std::exp( - deltaVvector.norm() / ( specificImpulseFunction_( initialTime + segmentDuration / 2.0 ) *
//                                                        physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION ) ) );

        // Update mass of the vehicle.
        bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( propagateMassToSegment( currentSegment + 1 ) );

    }

    std::cout << "propagated state: " << propagatedState << "\n\n";

    return propagatedState;
}


//! Propagate the trajectory inside one segment (low order solution).
Eigen::Vector6d SimsFlanaganLeg::propagateInsideBackwardSegment( double initialTime, double finalTime, double segmentDuration, Eigen::Vector6d initialState )
{
    // Compute time elapsed since start of the current leg segment.
    double timeElapsedCurrentSegment = finalTime - initialTime;

    // Compute current segment.
    int currentSegment = convertTimeToLegSegment( finalTime );
    std::cout << "mid-segment: " << ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ) << "\n\n";
    std::cout << "current segment: " << currentSegment << "\n\n";

    // Initialise propagated state.
    Eigen::Vector6d propagatedState = initialState;


    // Propagate from start of the current leg segment to the required propagation final time.
    if ( ( finalTime >= ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ) )
         || ( initialTime < ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ) ) )
    {
        std::cout << "no DELTAV: " << "\n\n";
        // Directly propagate the current state over the time elapsed since the start of the current leg segment.
        propagatedState = orbital_element_conversions::convertKeplerianToCartesianElements(
                    orbital_element_conversions::propagateKeplerOrbit(
                        orbital_element_conversions::convertCartesianToKeplerianElements(
                            propagatedState, centralBodyGravitationalParameter_ ),
                        timeElapsedCurrentSegment, centralBodyGravitationalParameter_ ), centralBodyGravitationalParameter_ );
    }
    else
    {
        std::cout << "detection DELTAV" << "\n\n";

        // First propagate the current state to half of the current segment.
        if ( initialTime - ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ) > 0.0 )
        {
            propagatedState = orbital_element_conversions::convertKeplerianToCartesianElements(
                        orbital_element_conversions::propagateKeplerOrbit(
                            orbital_element_conversions::convertCartesianToKeplerianElements(
                                propagatedState, centralBodyGravitationalParameter_ ),
                            ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ) - initialTime /*segmentDuration / 2.0*/, centralBodyGravitationalParameter_ ), centralBodyGravitationalParameter_ );
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
                            finalTime - ( timesAtNodes_[ currentSegment ] + segmentDuration / 2.0 ) /*timeElapsedCurrentSegment - segmentDuration / 2.0*/, centralBodyGravitationalParameter_ ),
                        centralBodyGravitationalParameter_ );
        }



    }

    std::cout << "propagated state: " << propagatedState << "\n\n";

    return propagatedState;
}


//! Propagate the trajectory to given time (low order solution).
Eigen::Vector6d SimsFlanaganLeg::propagateTrajectory( double initialTime, double finalTime, Eigen::Vector6d initialState )
{
    // Compute index of leg segment which corresponds to the propagation initial time.
    int initialSegment = convertTimeToLegSegment( initialTime );

    // Compute index of leg segment which corresponds to the propagation final time.
    int finalSegment = convertTimeToLegSegment( finalTime );

//    // Re-initialise mass of the spacecraft in the body map.
//    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );

    // Initialise propagated state.
    Eigen::Vector6d propagatedState = initialState; //stateAtDeparture_;

    // Initialise current segment.
    int currentSegment = initialSegment;

    // Compute segment duration.
    double segmentDuration;
    bool areInitialAndFinalSegmentsInSamePropagationPart;
    bool isPropagationForward;
    if ( ( initialTime <= timeAtMatchPoint_ ) && ( finalTime <= timeAtMatchPoint_ ) )
    {
//        std::cout << "FORWARD PROPAGATION: " << "\n\n";
        segmentDuration = timeAtMatchPoint_ / numberSegmentsForwardPropagation_;
        areInitialAndFinalSegmentsInSamePropagationPart = true;
        isPropagationForward = true;
    }
    else if ( ( initialTime >= timeAtMatchPoint_ ) && ( finalTime > timeAtMatchPoint_ ) )
    {
//        std::cout << "BACKWARD PROPAGATION: " << "\n\n";
        segmentDuration = timeAtMatchPoint_ / numberSegmentsBackwardPropagation_;
        areInitialAndFinalSegmentsInSamePropagationPart = true;
        isPropagationForward = false;
    }
    else
    {
//        std::cout << "BOTH FORWARD AND BACKWARD PROPAGATION: " << "\n\n";
        areInitialAndFinalSegmentsInSamePropagationPart = false;
    }


    if ( areInitialAndFinalSegmentsInSamePropagationPart )
    {

        if ( initialSegment == finalSegment /*( finalTime - initialTime ) < segmentDuration*/ )
        {
            std::cout << "propagation over less than one segment: " << "\n\n";
            // Directly propagate to final time inside the current segment.
            propagatedState = propagateInsideSegment( initialTime, finalTime, segmentDuration, propagatedState );
            std::cout << "initial time: " << initialTime << "\n\n";
            std::cout << "final time: " << finalTime << "\n\n";
        }
        else
        {
            std::cout << "propagate through several segments" << "\n\n";
//            std::cout << "initial segment: " << initialSegment << "\n\n";
            // First propagate to beginning of next segment.
            propagatedState = propagateInsideSegment( initialTime, timesAtNodes_[ initialSegment + 1 ], segmentDuration, propagatedState );

            // Update current segment.
            currentSegment++;

            // Compute number of segments the trajectory should be propagated through before reaching the segment corresponding
            // to the final propagation time.
            int propagatedSegments = ( finalTime - timesAtNodes_[ currentSegment ] ) / segmentDuration;
//            std::cout << "propagated segments: " << propagatedSegments << "\n\n";

            // Propagate through these segments.
            for ( int i = 0 ; i < propagatedSegments ; i++ )
            {
                propagatedState = propagateInsideSegment( timesAtNodes_[ currentSegment ], timesAtNodes_[ currentSegment + 1 ],
                        segmentDuration, propagatedState );
                currentSegment++;
            }

            // Propagate inside last segment.
            if ( ( finalTime - timesAtNodes_[ currentSegment ] ) > 0.0 )
            {
                propagatedState = propagateInsideSegment( timesAtNodes_[ currentSegment ], finalTime, segmentDuration, propagatedState );
                currentSegment++;
            }

//            std::cout << "number of segments to match point: " << currentSegment << "\n\n";

        }
    }


    // If initial time corresponds to forward propagation part and final time corresponds to backward propagation part.
    else
    {
        std::cout << "final time in backward propagation: " << "\n\n";
        // Compute number of segments between initial time and time at match point.
        int propagatedSegmentsToMatchPoint = ( timeAtMatchPoint_ - initialTime ) / ( timeAtMatchPoint_ / numberSegmentsForwardPropagation_ );

        // Propagate from initial time to beginning next segment.
        propagatedState = propagateInsideSegment( initialTime, timesAtNodes_[ initialSegment + 1 ],
                timeAtMatchPoint_ / numberSegmentsForwardPropagation_, propagatedState );
        currentSegment++;

        // Propagate from beginning next segment to match point.
        for ( int i = 0 ; i < propagatedSegmentsToMatchPoint ; i++ )
        {
            propagateInsideSegment( timesAtNodes_[ currentSegment ], timesAtNodes_[ currentSegment + 1 ],
                                    timeAtMatchPoint_ / numberSegmentsForwardPropagation_, propagatedState );
            currentSegment++;
        }

        // Compute number of segments between match point and final time.
        int propagatedSegmentsFromMatchPoint = ( finalTime - timeAtMatchPoint_ ) / ( timeAtMatchPoint_ / numberSegmentsBackwardPropagation_ );

        // Propagate from match point to beginning of the segment corresponding to the final time.
        for ( int i = 0 ; i < propagatedSegmentsFromMatchPoint ; i++ )
        {
            propagateInsideSegment( timesAtNodes_[ currentSegment ], timesAtNodes_[ currentSegment + 1 ],
                    timeAtMatchPoint_ / numberSegmentsBackwardPropagation_, propagatedState );
            currentSegment++;
        }

        // Propagate from beginning of the last segment to final time.
        propagatedState = propagateInsideSegment( timesAtNodes_[ currentSegment ], finalTime,
                                                  timeAtMatchPoint_ / numberSegmentsBackwardPropagation_, propagatedState );

    }

    return propagatedState;

}




//! Propagate the trajectory forward to given time (low order solution).
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
        std::cout << "propagation over less than one segment: " << "\n\n";
        // Directly propagate to final time inside the current segment.
        if ( finalTime - initialTime > 0 )
        {
            propagatedState = propagateInsideSegment( initialTime, finalTime, segmentDuration, propagatedState );
        }
    }
    else
    {
        std::cout << "propagate through several segments" << "\n\n";

        // First propagate to beginning of next segment.
        propagatedState = propagateInsideSegment( initialTime, timesAtNodes_[ initialSegment + 1 ], segmentDuration, propagatedState );

        // Update current segment.
        currentSegment++;

        // Compute number of segments the trajectory should be propagated through before reaching the segment corresponding
        // to the final propagation time.
        int propagatedSegments = ( finalTime - timesAtNodes_[ currentSegment ] ) / segmentDuration;
//            std::cout << "propagated segments: " << propagatedSegments << "\n\n";

        // Propagate through these segments.
        for ( int i = 0 ; i < propagatedSegments ; i++ )
        {
            propagatedState = propagateInsideSegment( timesAtNodes_[ currentSegment ], timesAtNodes_[ currentSegment + 1 ],
                    segmentDuration, propagatedState );
            currentSegment++;
        }

        // Propagate inside last segment.
        if ( ( finalTime - timesAtNodes_[ currentSegment ] ) > 0.0 )
        {
            propagatedState = propagateInsideSegment( timesAtNodes_[ currentSegment ], finalTime, segmentDuration, propagatedState );
            currentSegment++;
        }

//            std::cout << "number of segments to match point: " << currentSegment << "\n\n";

    }

    return propagatedState;

}


//! Propagate the trajectory backward to given time (low order solution).
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
        std::cout << "propagation over less than one segment: " << "\n\n";
        // Directly propagate to final time inside the current segment.
        if ( finalTime - initialTime < 0 )
        {
            propagatedState = propagateInsideBackwardSegment( initialTime, finalTime, segmentDuration, propagatedState );
        }
        std::cout << "initial time: " << initialTime << "\n\n";
        std::cout << "final time: " << finalTime << "\n\n";
    }
    else
    {
        std::cout << "propagate through several segments" << "\n\n";

        // First propagate to beginning of next segment.
        propagatedState = propagateInsideBackwardSegment( initialTime, timesAtNodes_[ initialSegment ], segmentDuration, propagatedState );

        // Update current segment.
        currentSegment--;

        // Compute number of segments the trajectory should be propagated through before reaching the segment corresponding
        // to the final propagation time.
        int propagatedSegments = ( timesAtNodes_[ currentSegment + 1 ] - finalTime ) / segmentDuration;
//            std::cout << "propagated segments: " << propagatedSegments << "\n\n";

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

//            std::cout << "number of segments to match point: " << currentSegment << "\n\n";

    }

    return propagatedState;
}



//! Propagate the trajectory to set of epochs (low order solution).
std::map< double, Eigen::Vector6d > SimsFlanaganLeg::propagateTrajectory(
        std::vector< double > epochs,
        std::map< double, Eigen::Vector6d >& propagatedTrajectory )
{

//    // Clear output map.
//    propagatedTrajectory.clear( );

    // Initialise propagated state.
    Eigen::Vector6d propagatedState = stateAtDeparture_;

    // Initialise mass of the spacecraft at departure.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );

//    // Time up to which the trajectory has to be propagated.
//        double finalPropagationTime = epochs[ 0 ];
//        for ( int i = 0 ; i < epochs.size() ; i++ )
//        {
//            if ( epochs[ i ] > finalPropagationTime )
//            {
//                finalPropagationTime = epochs[ i ];
//            }
//        }
    //    if ( ( finalPropagationTime < 0.0 ) || ( finalPropagationTime > timeOfFlight_ ) )
    //    {
    //        throw std::runtime_error( "Error when propagating trajectory in Sims-Flanagan, epochs at which the trajectory should be "
    //                                  "computed are not constrained between 0.0 and timeOfFlight" );
    //    }


    for ( int epochIndex = 0 ; epochIndex < epochs.size() ; epochIndex++ )
    {
        double currentTime = epochs[ epochIndex ];
//        std::cout << "current time: " << currentTime << "\n\n";
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

//        // Compute index of leg segment which corresponds to the propagation final time.
//        int currentSegment = convertTimeToLegSegment( currentTime );

        if ( epochIndex == 0 )
        {
            if ( currentTime > 0.0 )
            {
                propagatedState = propagateTrajectory( 0.0, currentTime, propagatedState );
            }
            propagatedTrajectory[ currentTime ] = propagatedState;
        }
        else
        {
            propagatedState = propagateTrajectory( epochs[ epochIndex - 1 ], currentTime, propagatedState );
            propagatedTrajectory[ currentTime ] = propagatedState;
        }

    }

    return propagatedTrajectory;

}



//! Propagate the trajectory to set of epochs (low order solution).
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
                if ( ( currentTime > initialTime ) /*&& ( currentTime <= timeAtMatchPoint_ )*/ )
                {
                    propagatedState = propagateTrajectoryForward( initialTime, currentTime, propagatedState,
                                                                  timeAtMatchPoint_ / numberSegmentsForwardPropagation_ );
                }
                propagatedTrajectory[ currentTime ] = propagatedState;
            }
            else
            {
                if ( currentTime < timeAtMatchPoint_ )
                {
                    propagatedState = propagateTrajectoryForward( epochs[ epochIndex - 1 ], currentTime,
                            propagatedState, segmentDuration ); //timeAtMatchPoint_ / numberSegmentsForwardPropagation_ );
                }
                if ( currentTime >= timeAtMatchPoint_ )
                {
                    propagatedState = propagateTrajectoryForward( epochs[ epochIndex - 1 ], currentTime,
                            propagatedState, segmentDuration ); // timeAtMatchPoint_ / numberSegmentsBackwardPropagation_ /*segmentDuration*/ );
                }
                propagatedTrajectory[ currentTime ] = propagatedState;
            }

        }

        bodyMap_[ centralBody_ ]->setConstantBodyMass( initialMass );

        return propagatedTrajectory;
}



//! Propagate the trajectory to set of epochs (low order solution).
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
    //        std::cout << "current time: " << currentTime << "\n\n";
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
                if ( ( currentTime < initialTime ) /*&& ( currentTime <= timeAtMatchPoint_ )*/ )
                {
                    propagatedState = propagateTrajectoryBackward( initialTime, currentTime, propagatedState, timeAtMatchPoint_ / numberSegmentsForwardPropagation_ );
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
