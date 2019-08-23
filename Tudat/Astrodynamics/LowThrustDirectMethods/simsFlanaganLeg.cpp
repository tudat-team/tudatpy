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


Eigen::Vector6d SimsFlanaganLeg::propagateForwardSegment( unsigned int indexSegment, Eigen::Vector6d initialState )
{

    // Compute time at half of the current segment.
    double currentTime = timesAtNodes_[ indexSegment ] +
            ( timesAtNodes_[ indexSegment + 1 ] - timesAtNodes_[ indexSegment ] ) / 2.0;

    // Compute segment duration for the forward propagation (from departure to match point).
    double segmentDuration = timeAtMatchPoint_ / numberSegmentsForwardPropagation_;

    // Total deltaV.
    Eigen::Vector3d deltaVvector;
    for ( unsigned int i = 0 ; i < 3 ; i++ )
    {
        deltaVvector[ i ] = maximumThrust_ / bodyMap_[ bodyToPropagate_ ]->getBodyMass() * segmentDuration * throttles_[ indexSegment ][ i ];
    }

    // Propagate to mid-segment.
    Eigen::Vector6d stateAtMidsegment =  orbital_element_conversions::convertKeplerianToCartesianElements(
                orbital_element_conversions::propagateKeplerOrbit(
                orbital_element_conversions::convertCartesianToKeplerianElements( initialState, centralBodyGravitationalParameter_ ),
                segmentDuration / 2.0, centralBodyGravitationalParameter_ ), centralBodyGravitationalParameter_ );

//    std::cout << "FP: current segment: " << indexSegment << "\n\n";
//    std::cout << "deltaV vector: " << deltaVvector.transpose() << "\n\n";
//    std::cout << "FP: state before deltaV: " << stateAtMidsegment.transpose() << "\n\n";
//    std::cout << "FP: current mass: " << bodyMap_[ bodyToPropagate_ ]->getBodyMass() << "\n\n";

    // Add dV at the middle of the segment.
    stateAtMidsegment.segment( 3, 3 ) += deltaVvector;

//    std::cout << "FP: state after deltaV: " << stateAtMidsegment.transpose() << "\n\n";

    // Propagate from mid-segment to end of the segment.
    Eigen::Vector6d stateAtSegmentEnd = orbital_element_conversions::convertKeplerianToCartesianElements(
                orbital_element_conversions::propagateKeplerOrbit(
                orbital_element_conversions::convertCartesianToKeplerianElements( stateAtMidsegment, centralBodyGravitationalParameter_ ),
                segmentDuration / 2.0, centralBodyGravitationalParameter_ ), centralBodyGravitationalParameter_ );

//    std::cout << "FP: state at segment end: " << stateAtSegmentEnd.transpose() << "\n\n";

//    std::cout << "FP: updated mass: " << bodyMap_[ bodyToPropagate_ ]->getBodyMass() *
//               std::exp( - deltaVvector.norm() / ( specificImpulse_ * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION ) ) << "\n\n";

    // Update mass of the vehicle.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( bodyMap_[ bodyToPropagate_ ]->getBodyMass() *
                std::exp( - deltaVvector.norm() / ( specificImpulseFunction_( currentTime ) * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION ) ) );

    // Update total deltaV value.
    totalDeltaV_ += deltaVvector.norm();

//    std::cout << "FP: current total deltaV: " << totalDeltaV_ << "\n\n";

//    /*Eigen::Vector6d */stateAtSegmentEnd = initialState;

    // Return state at end of the leg segment.
    return stateAtSegmentEnd;

}


Eigen::Vector6d SimsFlanaganLeg::propagateBackwardSegment( unsigned int indexSegment, Eigen::Vector6d initialState )
{
    // Compute segment duration for the backward propagation (from arrival to match point).
    double segmentDuration = timeAtMatchPoint_ / numberSegmentsBackwardPropagation_;

    // Total deltaV.
    Eigen::Vector3d deltaVvector;
    for ( unsigned int i = 0 ; i < 3 ; i++ )
    {
        deltaVvector[ i ] = maximumThrust_ / bodyMap_[ bodyToPropagate_ ]->getBodyMass() * segmentDuration
                * throttles_[ indexSegment ][ i ];
    }

    // Propagate backward to mid-segment.
    Eigen::Vector6d stateAtMidsegment = orbital_element_conversions::convertKeplerianToCartesianElements(
                orbital_element_conversions::propagateKeplerOrbit(
                orbital_element_conversions::convertCartesianToKeplerianElements( initialState, centralBodyGravitationalParameter_ ),
                - segmentDuration / 2.0, centralBodyGravitationalParameter_ ), centralBodyGravitationalParameter_ );

//    std::cout << "BP: current segment: " << indexSegment << "\n\n";
//    std::cout << "deltaV vector: " << deltaVvector.transpose() << "\n\n";
//    std::cout << "BP: state before deltaV: " << stateAtMidsegment.transpose() << "\n\n";
//    std::cout << "BP: current mass: " << bodyMap_[ bodyToPropagate_ ]->getBodyMass() << "\n\n";

    // Add dV at the middle of the segment.
    stateAtMidsegment.segment( 3, 3 ) -= deltaVvector;

//    std::cout << "BP: state after deltaV: " << stateAtMidsegment.transpose() << "\n\n";

    // Propagate backward from mid-segment to end of the segment.
    Eigen::Vector6d stateAtSegmentEnd = orbital_element_conversions::convertKeplerianToCartesianElements(
                orbital_element_conversions::propagateKeplerOrbit(
                orbital_element_conversions::convertCartesianToKeplerianElements( stateAtMidsegment, centralBodyGravitationalParameter_ ),
                - segmentDuration / 2.0, centralBodyGravitationalParameter_ ), centralBodyGravitationalParameter_ );

//    std::cout << "BP: state at segment end: " << stateAtSegmentEnd.transpose() << "\n\n";

//    std::cout << "BP: updated mass: " << bodyMap_[ bodyToPropagate_ ]->getBodyMass() /
//               std::exp( - deltaVvector.norm() / ( specificImpulse_ * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION ) ) << "\n\n";

//    std::cout << "BP: updated mass: " << propagateMassToSegment( indexSegment ) << "\n\n";

//    // Update mass of the vehicle.
//    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( bodyMap_[ bodyToPropagate_ ]->getBodyMass() /
//                                   std::exp( - deltaVvector.norm() / ( specificImpulse_ * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION ) ) );

    // Other method to update the mass of the vehicle.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( propagateMassToSegment( indexSegment ) );

    // Update total deltaV value.
    totalDeltaV_ += deltaVvector.norm();

//    std::cout << "BP: current total deltaV: " << totalDeltaV_ << "\n\n";

    // Return state at the end of the leg segment.
    return stateAtSegmentEnd;

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


    // Loop over the leg segments of the forward propagation.
    for ( unsigned int currentSegment = 0 ; currentSegment < numberSegmentsForwardPropagation_ ; currentSegment++ )
    {
        // Forward propagation over one segment of the leg.
        currentState = propagateForwardSegment( currentSegment, currentState );
    }

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

    // Loop over the leg segments of the backward propagation.
    for ( unsigned int currentSegment = numberSegments_ - 1 ;
          currentSegment >= numberSegmentsForwardPropagation_ ; currentSegment-- )
    {
        // Backward propagation over one segment of the leg.
        currentState = propagateBackwardSegment( currentSegment, currentState );
    }

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


//! Propagate the mass from departure to arrival for high order solution (return mass of the spacecraft at end of the leg).
std::map< double, Eigen::VectorXd > SimsFlanaganLeg::propagateMassToLegArrivalHighOrderSolution(
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    // Retrieve acceleration map.
    basic_astrodynamics::AccelerationMap accelerationMap = getAccelerationModelFullLeg( );

    // Create mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModel;
    massRateModel[ bodyToPropagate_ ] = createMassRateModel( bodyToPropagate_, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
                                                       bodyMap_, accelerationMap );

    // Redefine integrator settings.
    integratorSettings->initialTime_ = 0.0;
    integratorSettings->initialTimeStep_ = std::fabs( integratorSettings->initialTimeStep_ );

    // Redefine termination settings.
    std::shared_ptr< propagators::PropagationTerminationSettings > terminationSettings
            = std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight_, true );

    // Create settings for propagating the mass of the vehicle.
    std::shared_ptr< propagators::MassPropagatorSettings< double > > massPropagatorSettings
            = std::make_shared< propagators::MassPropagatorSettings< double > >(
                std::vector< std::string >{ bodyToPropagate_ }, massRateModel,
                ( Eigen::Matrix< double, 1, 1 >( ) << initialSpacecraftMass_ ).finished( ), terminationSettings );

    // Propagate the mass forward over the full leg.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap_, integratorSettings, massPropagatorSettings );

    // Return mass history.
    return dynamicsSimulator.getEquationsOfMotionNumericalSolution();
}


double SimsFlanaganLeg::computeDeltaVperSegmentHighOrderSolution(
        unsigned int indexSegment,
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    double initialTime = timesAtNodes_[ indexSegment ];
    double finalTime = timesAtNodes_[ indexSegment + 1 ];

    std::shared_ptr< numerical_quadrature::QuadratureSettings< double > > quadratureSettings_ =
            std::make_shared< numerical_quadrature::GaussianQuadratureSettings < double > >( initialTime, 16 );

//    // Thrust acceleration function to use quadrature.
//    // Define thrust acceleration as a function of time (to be integrated to compute the associated deltaV).
//    std::function< double( const double ) > thrustAcceleration = [ = ] ( const double currentTime ){

//        // Compute mass history.
//        std::map< double, Eigen::VectorXd > massHistory = propagateMassToLegArrivalHighOrderSolution( integratorSettings );

//        // Interpolate mass to current time.
//        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > massInterpolator =
//                std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >( massHistory, 8 );
//        double currentMass = massInterpolator->interpolate( currentTime )[ 0 ];

//        Eigen::Vector3d accelerationVector;
//        for ( int i = 0 ; i < 3 ; i++ )
//        {
//            accelerationVector[ i ] = maximumThrust_ / currentMass * throttles_[ indexSegment ][ i ];
//        }
//        return accelerationVector.norm();
//    };

    // Thrust acceleration function to use integrator.
    // Define thrust acceleration as a function of time (to be integrated to compute the associated deltaV).
    std::function< Eigen::Vector1d( const double, const Eigen::Vector1d& ) > thrustAcceleration = [ = ] ( const double currentTime,
            const Eigen::Vector1d& independentVariable ){

        // Compute mass history.
        std::map< double, Eigen::VectorXd > massHistory = propagateMassToLegArrivalHighOrderSolution( integratorSettings );

        // Interpolate mass to current time.
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > massInterpolator =
                std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >( massHistory, 8 );
        double currentMass = massInterpolator->interpolate( currentTime )[ 0 ];

        Eigen::Vector3d accelerationVector;
        for ( int i = 0 ; i < 3 ; i++ )
        {
            accelerationVector[ i ] = maximumThrust_ / currentMass * throttles_[ indexSegment ][ i ];
        }
        return ( Eigen::Vector1d() << accelerationVector.norm( ) ).finished( );
    };

//    // Create numerical quadrature from quadrature settings.
//    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
//            numerical_quadrature::createQuadrature( thrustAcceleration, quadratureSettings_, finalTime );

    // Redefine integrator settings.
    integratorSettings->initialTime_ = initialTime;
    integratorSettings->initialTimeStep_ = std::fabs( integratorSettings->initialTimeStep_ );

    // Create numerical integrator.
    std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::Vector1d /*, double, double*/ > > integrator =
            numerical_integrators::createIntegrator< double, Eigen::Vector1d >( thrustAcceleration, Eigen::Vector1d::Zero( ), integratorSettings );
    double deltaV = integrator->integrateTo( finalTime, integratorSettings->initialTimeStep_ )[ 0 ];

//    return quadrature->getQuadrature( );
    return deltaV;
}


//! Propagate the spacecraft trajectory forward over a leg segment (high order solution).
Eigen::Vector6d SimsFlanaganLeg::propagateForwardSegmentHighOrderSolution(
        unsigned int indexSegment, Eigen::Vector6d initialState,
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        propagators::TranslationalPropagatorType propagatorType )
{
    // Retrieve acceleration map.
    basic_astrodynamics::AccelerationMap accelerationMap = getAccelerationModelPerSegment( indexSegment );

    // Create mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModel;
    massRateModel[ bodyToPropagate_ ] = createMassRateModel( bodyToPropagate_, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
                                                       bodyMap_, accelerationMap );

    // Redefine integrator settings.
    integratorSettings->initialTime_ = timesAtNodes_[ indexSegment ];
    integratorSettings->initialTimeStep_ = std::fabs( integratorSettings->initialTimeStep_ );

    // Redefine termination settings.
    std::shared_ptr< propagators::PropagationTerminationSettings > terminationSettings
            = std::make_shared< propagators::PropagationTimeTerminationSettings >( timesAtNodes_[ indexSegment + 1 ], true );

    // Define propagator settings.
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
            ( std::vector< std::string >{ centralBody_ }, accelerationMap, std::vector< std::string >{ bodyToPropagate_ },
              initialState, terminationSettings, propagatorType );


    // Create settings for propagating the mass of the vehicle.
    std::shared_ptr< propagators::MassPropagatorSettings< double > > massPropagatorSettings
            = std::make_shared< propagators::MassPropagatorSettings< double > >(
                std::vector< std::string >{ bodyToPropagate_ }, massRateModel,
                ( Eigen::Matrix< double, 1, 1 >( ) << bodyMap_[ bodyToPropagate_ ]->getBodyMass( ) ).finished( ), terminationSettings );

    // Create list of propagation settings.
    std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
    propagatorSettingsVector.push_back( propagatorSettings );
    propagatorSettingsVector.push_back( massPropagatorSettings );

    // Backward hybrid propagation settings.
    std::shared_ptr< propagators::PropagatorSettings< double > > hybridPropagatorSettings =
            std::make_shared< propagators::MultiTypePropagatorSettings< double > >( propagatorSettingsVector, terminationSettings );

//    std::cout << "FP: mass before propagation: " << bodyMap_[ bodyToPropagate_ ]->getBodyMass() << "\n\n";

    // Propagate the trajectory forward over the leg segment.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap_, integratorSettings, hybridPropagatorSettings );
    Eigen::VectorXd stateAtSegmentEnd = dynamicsSimulator.getEquationsOfMotionNumericalSolution().rbegin()->second;

//    std::cout << "FP: current index: " << indexSegment << "\n\n";
//    std::cout << "FP: state at segment end: " << stateAtSegmentEnd << "\n\n";

    // Update total deltaV value.
    totalDeltaV_ += computeDeltaVperSegmentHighOrderSolution( indexSegment, integratorSettings );

//    std::cout << "FP: current deltaV: " << totalDeltaV_ << "\n\n";
//    std::cout << "FP: mass in body map before forced update: " << bodyMap_[ bodyToPropagate_ ]->getBodyMass() << "\n\n";

    // Update the mass of the spacecraft at the end of the propagation
    // (ensure consistency between stateAtSegmentEnd and mass od the spacecraft in the bodyMap -> interpolation errors?)
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( stateAtSegmentEnd[ 6 ] );

//    std::cout << "FP: mass in body map after forced update: " << bodyMap_[ bodyToPropagate_ ]->getBodyMass() << "\n\n";

    // Return the propagated state vector at the end of the leg segment.
    return stateAtSegmentEnd.segment( 0, 6 );
}


//! Propagate the spacecraft trajectory backward over a leg segment (high order solution).
Eigen::Vector6d SimsFlanaganLeg::propagateBackwardSegmentHighOrderSolution(
        unsigned int indexSegment, Eigen::Vector6d initialState,
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        propagators::TranslationalPropagatorType propagatorType )
{
    // Retrieve acceleration map.
    basic_astrodynamics::AccelerationMap accelerationMap = getAccelerationModelPerSegment( indexSegment );

    // Create mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModel;
    massRateModel[ bodyToPropagate_ ] = createMassRateModel( bodyToPropagate_, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
                                                       bodyMap_, accelerationMap );

    // Redefine integrator settings.
    integratorSettings->initialTime_ = timesAtNodes_[ indexSegment + 1 ];
    integratorSettings->initialTimeStep_ = - std::fabs( integratorSettings->initialTimeStep_ );

    // Redefine termination settings.
    std::shared_ptr< propagators::PropagationTerminationSettings > terminationSettings
            = std::make_shared< propagators::PropagationTimeTerminationSettings >( timesAtNodes_[ indexSegment ], true );

    // Define propagator settings.
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
            ( std::vector< std::string >{ centralBody_ }, accelerationMap, std::vector< std::string >{ bodyToPropagate_ },
              initialState, terminationSettings, propagatorType );


    // Create settings for propagating the mass of the vehicle.
    std::shared_ptr< propagators::MassPropagatorSettings< double > > massPropagatorSettings
            = std::make_shared< propagators::MassPropagatorSettings< double > >(
                std::vector< std::string >{ bodyToPropagate_ }, massRateModel,
                ( Eigen::Matrix< double, 1, 1 >( ) << bodyMap_[ bodyToPropagate_ ]->getBodyMass( ) ).finished( ), terminationSettings );

    // Create list of propagation settings.
    std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
    propagatorSettingsVector.push_back( propagatorSettings );
    propagatorSettingsVector.push_back( massPropagatorSettings );

    // Backward hybrid propagation settings.
    std::shared_ptr< propagators::PropagatorSettings< double > > hybridPropagatorSettings =
            std::make_shared< propagators::MultiTypePropagatorSettings< double > >( propagatorSettingsVector, terminationSettings );

    // Propagate the trajectory backward over the leg segment.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap_, integratorSettings, hybridPropagatorSettings );
    Eigen::VectorXd stateAtSegmentEnd = dynamicsSimulator.getEquationsOfMotionNumericalSolution().begin()->second;

//    std::cout << "BP: current index: " << indexSegment << "\n\n";
//    std::cout << "BP: state at segment end: " << stateAtSegmentEnd << "\n\n";

    // Update total deltaV value.
    totalDeltaV_ += computeDeltaVperSegmentHighOrderSolution( indexSegment, integratorSettings );

//    std::cout << "BP: current deltaV: " << totalDeltaV_ << "\n\n";

    // Update the mass of the spacecraft at the end of the propagation (MIGHT BE UNNECESSARY?)
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( stateAtSegmentEnd[ 6 ] );

//    std::cout << "BP: mass in body map after forced update: " << bodyMap_[ bodyToPropagate_ ]->getBodyMass() << "\n\n";

    // Return the propagated state vector at the end of the leg segment.
    return stateAtSegmentEnd.segment( 0, 6 );

}


//! Propagate the spacecraft trajectory from departure to match point (forward propagation) (high order solution).
void SimsFlanaganLeg::propagateForwardFromDepartureToMatchPointHighOrderSolution(
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        propagators::TranslationalPropagatorType propagatorType )
{
    // Initialise current state at the leg departure.
    Eigen::Vector6d currentState = stateAtDeparture_;

    Eigen::Vector6d testWithoutThrust = orbital_element_conversions::convertKeplerianToCartesianElements(
                orbital_element_conversions::propagateKeplerOrbit(
                    orbital_element_conversions::convertCartesianToKeplerianElements(
                        stateAtDeparture_, centralBodyGravitationalParameter_ ), timeOfFlight_ / 2.0, centralBodyGravitationalParameter_),
                centralBodyGravitationalParameter_ );

    // Loop over the leg segments of the forward propagation.
    for ( unsigned int currentSegment = 0 ; currentSegment < numberSegmentsForwardPropagation_ ; currentSegment++ )
    {
        // Forward propagation over one segment of the leg.
        currentState = propagateForwardSegmentHighOrderSolution( currentSegment, currentState, integratorSettings, propagatorType /*, propagatorSettings*/ );
    }

    // Set state vector at match point once the forward propagation is over.
    stateAtMatchPointFromForwardPropagation_ = currentState;

    // Set mass of the spacecraft at match point after the forward propagation.
    massAtMatchPointFromForwardPropagation_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass( );


    std::cout << "END OF FORWARD PROPAGATION: " << "\n\n";
    std::cout << "state at match point: " << stateAtMatchPointFromForwardPropagation_ << "\n\n";
    std::cout << "mass at match point: " << massAtMatchPointFromForwardPropagation_ << "\n\n";

//    std::cout << "test without thrust: " << testWithoutThrust << "\n\n";

}

//! Propagate the spacecraft trajectory from arrival to match point (backward propagation) (high order solution).
void SimsFlanaganLeg::propagateBackwardFromArrivalToMatchPointHighOrderSolution(
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        propagators::TranslationalPropagatorType propagatorType )
{
    // Initialise the spacecraft mass at leg arrival.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( propagateMassToLegArrivalHighOrderSolution( integratorSettings ).rbegin()->second[ 0 ] );
    std::cout << "mass at the end of the leg (propagation): " << propagateMassToLegArrivalHighOrderSolution( integratorSettings ).rbegin()->second[ 0 ] << "\n\n";

    // Initialise current state at the leg arrival.
    Eigen::Vector6d currentState = stateAtArrival_;

    Eigen::Vector6d testWithoutThrust = orbital_element_conversions::convertKeplerianToCartesianElements(
                orbital_element_conversions::propagateKeplerOrbit(
                    orbital_element_conversions::convertCartesianToKeplerianElements(
                        stateAtArrival_, centralBodyGravitationalParameter_ ), - timeOfFlight_ / 2.0, centralBodyGravitationalParameter_),
                centralBodyGravitationalParameter_ );

    // Loop over the leg segments of the backward propagation.
    for ( unsigned int currentSegment = numberSegments_ - 1 ; currentSegment >= numberSegmentsForwardPropagation_ ; currentSegment-- )
    {
        // Backward propagation over one segment of the leg.
        currentState = propagateBackwardSegmentHighOrderSolution( currentSegment, currentState, integratorSettings, propagatorType /*, propagatorSettings*/ );
    }

    // Set state vector at match point once the backward propagation is over.
    stateAtMatchPointFromBackwardPropagation_ = currentState;

    // Set mass of the spacecraft at match point after the backward propagation.
    massAtMatchPointFromBackwardPropagation_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass( );


    std::cout << "END OF BACKWARD PROPAGATION: " << "\n\n";
    std::cout << "state at match point: " << stateAtMatchPointFromBackwardPropagation_ << "\n\n";
    std::cout << "mass at match point: " << massAtMatchPointFromBackwardPropagation_ << "\n\n";

//    std::cout << "test without thrust: " << testWithoutThrust << "\n\n";

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
        double initialMass, double segmentDuration )
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
//                if ( ( currentTime > 0.0 ) && ( currentTime <= timeAtMatchPoint_ ) )
//                {
//                    propagatedState = propagateTrajectoryForward( 0.0, currentTime, propagatedState, timeAtMatchPoint_ / numberSegmentsForwardPropagation_ );
//                }
//                if ( ( currentTime > 0.0 ) && ( currentTime > timeAtMatchPoint_ ) )
//                {
//                    propagatedState = propagateTrajectoryForward( 0.0, currentTime, propagatedState, timeAtMatchPoint_ / numberSegmentsBackwardPropagation_ );
//                }
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
        Eigen::Vector6d initialState, double initialMass, double segmentDuration )
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
//                if ( ( currentTime > 0.0 ) && ( currentTime <= timeAtMatchPoint_ ) )
//                {
//                    propagatedState = propagateTrajectoryForward( initial, currentTime, propagatedState, timeAtMatchPoint_ / numberSegmentsForwardPropagation_ );
//                }
////                if ( ( currentTime > 0.0 ) && ( currentTime > timeAtMatchPoint_ ) )
////                {
////                    propagatedState = propagateTrajectoryForward( 0.0, currentTime, propagatedState, timeAtMatchPoint_ / numberSegmentsBackwardPropagation_ );
////                }
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






//! Propagate the trajectory to given time (high order solution).
Eigen::Vector6d SimsFlanaganLeg::propagateTrajectoryHighOrderSolution(
        double currentTime,
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        propagators::TranslationalPropagatorType propagatorType )
{
    // Compute index of leg segment which corresponds to the propagation final time.
    int currentSegment = convertTimeToLegSegment( currentTime );

    // Re-initialise mass of the spacecraft in the body map.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );

    // Initialise propagated state.
    Eigen::Vector6d propagatedState = stateAtDeparture_;

    // Propagate over all segments and stop at the beginning of the current one.
    for ( int segment = 0 ; segment < currentSegment ; segment++ )
    {
        propagatedState = propagateForwardSegmentHighOrderSolution( segment, propagatedState, integratorSettings, propagatorType );
    }


    // Retrieve acceleration map.
    basic_astrodynamics::AccelerationMap accelerationMap = getAccelerationModelPerSegment( currentSegment );

    // Create mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModel;
    massRateModel[ bodyToPropagate_ ] = createMassRateModel( bodyToPropagate_,
                                                             std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
                                                       bodyMap_, accelerationMap );

    // Redefine integrator settings.
    integratorSettings->initialTime_ = timesAtNodes_[ currentSegment ];
    integratorSettings->initialTimeStep_ = std::fabs( integratorSettings->initialTimeStep_ );

    // Redefine termination settings.
    std::shared_ptr< propagators::PropagationTerminationSettings > terminationSettings
            = std::make_shared< propagators::PropagationTimeTerminationSettings >( currentTime, true );

    // Define propagator settings.
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
            ( std::vector< std::string >{ centralBody_ }, accelerationMap, std::vector< std::string >{ bodyToPropagate_ },
              propagatedState, terminationSettings, propagatorType );

    // Create settings for propagating the mass of the vehicle.
    std::shared_ptr< propagators::MassPropagatorSettings< double > > massPropagatorSettings
            = std::make_shared< propagators::MassPropagatorSettings< double > >(
                std::vector< std::string >{ bodyToPropagate_ }, massRateModel,
                ( Eigen::Matrix< double, 1, 1 >( ) << bodyMap_[ bodyToPropagate_ ]->getBodyMass( ) ).finished( ), terminationSettings );

    // Create list of propagation settings.
    std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
    propagatorSettingsVector.push_back( propagatorSettings );
    propagatorSettingsVector.push_back( massPropagatorSettings );

    // Backward hybrid propagation settings.
    std::shared_ptr< propagators::PropagatorSettings< double > > hybridPropagatorSettings =
            std::make_shared< propagators::MultiTypePropagatorSettings< double > >( propagatorSettingsVector, terminationSettings );

    // Propagate the trajectory forward over the leg segment.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap_, integratorSettings, hybridPropagatorSettings );
    Eigen::VectorXd stateAtSegmentEnd = dynamicsSimulator.getEquationsOfMotionNumericalSolution().rbegin()->second;

    // Update the mass of the spacecraft at the end of the propagation
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( stateAtSegmentEnd[ 6 ] );

    // Return the propagated state vector at the end of the leg segment.
    return stateAtSegmentEnd.segment( 0, 6 );

}


//! Propagate the trajectory to set of epochs (high order solution).
void SimsFlanaganLeg::propagateTrajectoryHighOrderSolution(
        std::vector< double > epochs,
        std::map< double, Eigen::Vector6d >& propagatedTrajectory,
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        propagators::TranslationalPropagatorType propagatorType )
{
    // Re-initialise mass of the spacecraft in the body map.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );

    // Initialise propagated state.
    Eigen::Vector6d propagatedState = stateAtDeparture_;


    // Retrieve acceleration map.
    basic_astrodynamics::AccelerationMap accelerationMap = getAccelerationModelFullLeg( );

    // Create mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModel;
    massRateModel[ bodyToPropagate_ ] = createMassRateModel( bodyToPropagate_, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
                                                             bodyMap_, accelerationMap );


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

        double initialTime;
        double finalTime;
        if ( epochIndex == 0 )
        {
            initialTime = 0.0;
            finalTime = epochs[ epochIndex ];
        }
        else
        {
            initialTime = epochs[ epochIndex - 1 ];
            finalTime = epochs[ epochIndex ];
        }


        // Redefine integrator settings.
        integratorSettings->initialTime_ = initialTime;
        integratorSettings->initialTimeStep_ = std::fabs( integratorSettings->initialTimeStep_ );

        // Redefine termination settings.
        std::shared_ptr< propagators::PropagationTerminationSettings > terminationSettings
                = std::make_shared< propagators::PropagationTimeTerminationSettings >( finalTime, true );

        // Define propagator settings.
        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                ( std::vector< std::string >{ centralBody_ }, accelerationMap, std::vector< std::string >{ bodyToPropagate_ },
                  propagatedState, terminationSettings, propagatorType );

        // Create settings for propagating the mass of the vehicle.
        std::shared_ptr< propagators::MassPropagatorSettings< double > > massPropagatorSettings
                = std::make_shared< propagators::MassPropagatorSettings< double > >(
                    std::vector< std::string >{ bodyToPropagate_ }, massRateModel,
                    ( Eigen::Matrix< double, 1, 1 >( ) << bodyMap_[ bodyToPropagate_ ]->getBodyMass( ) ).finished( ), terminationSettings );

        // Create list of propagation settings.
        std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
        propagatorSettingsVector.push_back( propagatorSettings );
        propagatorSettingsVector.push_back( massPropagatorSettings );

        // Backward hybrid propagation settings.
        std::shared_ptr< propagators::PropagatorSettings< double > > hybridPropagatorSettings =
                std::make_shared< propagators::MultiTypePropagatorSettings< double > >( propagatorSettingsVector, terminationSettings );

        // Propagate the trajectory forward over the leg segment.
        propagators::SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap_, integratorSettings, hybridPropagatorSettings );
        Eigen::VectorXd stateAtSegmentEnd = dynamicsSimulator.getEquationsOfMotionNumericalSolution().rbegin()->second;

        propagatedState = stateAtSegmentEnd;

        // Update the mass of the spacecraft at the end of the propagation
        bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( stateAtSegmentEnd[ 6 ] );

        propagatedTrajectory[ currentTime ] = stateAtSegmentEnd;

    }

}



} // namespace low_thrust_direct_methods
} // namespace tudat
