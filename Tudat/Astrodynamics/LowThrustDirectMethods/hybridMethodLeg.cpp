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
#include "Tudat/Astrodynamics/LowThrustDirectMethods/hybridMethodLeg.h"
#include "Tudat/Mathematics/NumericalQuadrature/createNumericalQuadrature.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/modifiedEquinoctialElementConversions.h"

namespace tudat
{
namespace low_thrust_direct_methods
{

using namespace orbital_element_conversions;

//Eigen::Vector5d computeOptimalThrustComponents(
//        Eigen::Vector6d& currentState,
//        Eigen::Vector6d& currentCoStates,
//        double mass,
//        const double maximumThrustMagnitude )
//{
//    // Declaring optimal thrust components vector.
//    Eigen::Vector5d optimalThrustComponents;

//    // Retrieve modified equinoctial elements.
//    double p = currentState[ semiParameterIndex ];
//    double f = currentState[ fElementIndex ];
//    double g = currentState[ gElementIndex ];
//    double h = currentState[ hElementIndex ];
//    double k = currentState[ kElementIndex ];
//    double L = currentState[ trueLongitudeIndex ];

//    double w1 = 1.0 + f * std::sin( L ) + g * std::sin( L );
//    double w2 = 1.0 + h * h + k * k;

//    // Compute all required auxiliary variables to compute optimal angle alpha.
//    double lambdap = currentCoStates[ semiParameterIndex ] * ( 2.0 * p ) / w1;
//    double lambdaf1 = currentCoStates[ fElementIndex ] * std::sin( L );
//    double lambdag1 = currentCoStates[ gElementIndex ] * std::cos( L );
//    double lambdaf2 = currentCoStates[ fElementIndex ] * ( ( w1 + 1.0 ) * std::cos( L ) + f ) / w1;
//    double lambdag2 = currentCoStates[ gElementIndex ] * ( ( w1 + 1.0 ) * std::sin( L ) + g ) / w1;

//    // Compute sinus of the optimal value of angle alpha.
//    double sinOptimalAlpha = - ( lambdaf1 - lambdag1 ) /
//            std::sqrt( ( lambdaf1 - lambdag1 ) * ( lambdaf1 - lambdag1 ) +
//                       ( lambdap + lambdaf2 + lambdag2 ) * ( lambdap + lambdaf2 + lambdag2 ) );

//    // Compute cosinus of the optimal value of angle alpha.
//    double cosOptimalAlpha = - ( lambdap + lambdaf2 + lambdag2 ) /
//            std::sqrt( ( lambdaf1 - lambdag1 ) * ( lambdaf1 - lambdag1 ) +
//                       ( lambdap + lambdaf2 + lambdag2 ) * ( lambdap + lambdaf2 + lambdag2 ) );

//    // Compute all required auxiliary variables to compute optial angle beta.
//    lambdap = currentCoStates[ semiParameterIndex ] * ( 2.0 * p ) / w1 * cosOptimalAlpha;
//    lambdaf1 = currentCoStates[ fElementIndex ] * std::sin( L ) * sinOptimalAlpha;
//    lambdag1 = currentCoStates[ gElementIndex ] * std::cos( L ) * sinOptimalAlpha;
//    lambdaf2 = currentCoStates[ fElementIndex ] * ( ( 1.0 + w1 ) * std::cos( L ) + f ) / w1 * cosOptimalAlpha;
//    lambdag2 = currentCoStates[ gElementIndex ] * ( ( 1.0 + w1 ) * std::sin( L ) + g ) / w1 * cosOptimalAlpha;
//    double lambdaf3 = currentCoStates[ fElementIndex ] * ( g / w1 ) * ( h * std::sin( L ) - k * std::cos( L ) );
//    double lambdag3 = currentCoStates[ gElementIndex ] * ( f / w1 ) * ( h * std::sin( L ) - k * std::cos( L ) );
//    double lambdah = currentCoStates[ hElementIndex ] * ( w2 * std::cos( L ) ) / ( 2.0 * w1 );
//    double lambdak = currentCoStates[ kElementIndex ] * ( w2 * std::sin( L ) ) / ( 2.0 * w1 );

//    // Compute cosinus of optimal thrust angle beta.
//    double cosOptimalBeta = - ( - lambdaf3 + lambdag3 + lambdah + lambdak ) /
//            std::sqrt( ( - lambdaf3 + lambdag3 + lambdah + lambdak ) * ( - lambdaf3 + lambdag3 + lambdah + lambdak )
//                       + ( lambdap + lambdaf1 - lambdag1 + lambdaf2 + lambdag2 )
//                       * ( lambdap + lambdaf1 - lambdag1 + lambdaf2 + lambdag2 ) );

//    // Compute sinus of optimal thrust angle beta.
//    double sinOptimalBeta = - ( lambdap + lambdaf1 - lambdag1 + lambdaf2 + lambdag2 ) /
//            std::sqrt( ( - lambdaf3 + lambdag3 + lambdah + lambdak ) * ( - lambdaf3 + lambdag3 + lambdah + lambdak )
//                       + ( lambdap + lambdaf1 - lambdag1 + lambdaf2 + lambdag2 )
//                       * ( lambdap + lambdaf1 - lambdag1 + lambdaf2 + lambdag2 ) );


//    // Switching function for the thrust magnitude.
//    double thrustMagnitudeSwitchingCondition = ( 1.0 / mass ) *
//            ( lambdap * cosOptimalBeta + lambdah * sinOptimalBeta + lambdak * sinOptimalBeta
//            + lambdaf1 * cosOptimalBeta + lambdaf2 * cosOptimalBeta - lambdaf3 * sinOptimalBeta
//            - lambdag1 * cosOptimalBeta + lambdag2 * cosOptimalBeta + lambdag3 * sinOptimalBeta );

//    double thrustMagnitude = 0.0;
//    if ( thrustMagnitudeSwitchingCondition <= 0.0 )
//    {
//        thrustMagnitude = maximumThrustMagnitude;
//    }


//    // Build output vector containing the cosinus and sinus of the optimal thrust angles alpha and beta, respectively, as well as the
//    // current thrust magnitude.
//    optimalThrustComponents[ 0 ] = cosOptimalAlpha;
//    optimalThrustComponents[ 1 ] = sinOptimalAlpha;
//    optimalThrustComponents[ 2 ] = cosOptimalBeta;
//    optimalThrustComponents[ 3 ] = sinOptimalBeta;
//    optimalThrustComponents[ 4 ] = thrustMagnitude;


//    return optimalThrustComponents;

//}

//Eigen::Vector6d computeCurrentCoStates(
//        const Eigen::Vector6d& initialCoStates,
//        const Eigen::Vector6d& finalCoStates,
//        const double currentTime )
//{
//    // Declare current co-states vector.
//    Eigen::Vector6d currentCoStates;

//    for ( int i = 0 ; i < 6 ; i++ )
//    {
//        currentCoStates[ i ] = initialCoStates[ i ]
//                + ( currentTime / timeOfFlight_ ) * ( finalCoStates[ i ] - initialCoStates[ i ] );
//    }

//    return currentCoStates;
//}

//! Retrieve MEE costates-based thrust acceleration.
//std::shared_ptr< propulsion::ThrustAcceleration > HybridMethodLeg::getMEEcostatesBasedThrustAccelerationModel( )
std::shared_ptr< simulation_setup::AccelerationSettings > HybridMethodLeg::getMEEcostatesBasedThrustAccelerationSettings( )
{
    // Define thrust direction settings from the MEE costates.
    std::shared_ptr< simulation_setup::MeeCostateBasedThrustDirectionSettings > thrustDirectionSettings =
            std::make_shared< simulation_setup::MeeCostateBasedThrustDirectionSettings >( bodyToPropagate_, centralBody_,
                                                                                          costatesFunction_ );

    std::function< double ( const double ) > specificImpulseFunction = [ = ]( const double currentTime )
    {
      return specificImpulse_;
    };

    // Define bang-bang thrust magnitude settings based on MEE co-states.
    std::shared_ptr< simulation_setup::FromMeeCostatesBangBangThrustMagnitudeSettings > thrustMagnitudeSettings
            = std::make_shared< simulation_setup::FromMeeCostatesBangBangThrustMagnitudeSettings >(
                maximumThrust_, specificImpulseFunction, costatesFunction_, bodyToPropagate_, centralBody_ );

    // Define thrust acceleration settings.
    std::shared_ptr< simulation_setup::ThrustAccelerationSettings > thrustAccelerationSettings =
            std::make_shared< simulation_setup::ThrustAccelerationSettings >(
                thrustDirectionSettings, thrustMagnitudeSettings );

//    // Create MEE costates-based thrust acceleration model.
//    std::shared_ptr< propulsion::ThrustAcceleration > meeCostatesBasedThrustAccelerationModel =
//            createThrustAcceleratioModel( thrustAccelerationSettings, bodyMap_, bodyToPropagate_ );

    return thrustAccelerationSettings;

}


//! Retrieve hybrid method acceleration model (including thrust and central gravity acceleration)
basic_astrodynamics::AccelerationMap HybridMethodLeg::getAccelerationModel( )
{
    // Acceleration from the central body.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
//    bodyToPropagateAccelerations[ centralBody_ ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
//                                                                basic_astrodynamics::central_gravity ) );
    bodyToPropagateAccelerations[ bodyToPropagate_ ].push_back( getMEEcostatesBasedThrustAccelerationSettings( ) );

    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ bodyToPropagate_ ] = bodyToPropagateAccelerations;

    // Create the acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap_, accelerationMap, std::vector< std::string >{ bodyToPropagate_ }, std::vector< std::string >{ centralBody_ } );

////     Add thurst acceleration.
//    accelerationModelMap[ bodyToPropagate_ ] = getMEEcostatesBasedThrustAccelerationModel( );

    return accelerationModelMap;
}


//! Propagate the spacecraft trajectory.
Eigen::Vector6d HybridMethodLeg::propagateTrajectory(
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings/*,
        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettings*/ )
{
    // Re-initialise integrator settings.
    integratorSettings->initialTime_ = 0.0;

//    propagatorSettings->resetInitialStates( stateAtDeparture_ );
//    propagatorSettings->bodiesToIntegrate_ = std::vector< std::string >{ bodyToPropagate_ };
//    propagatorSettings->centralBodies_ = std::vector< std::string >{ centralBody_ };
//    propagatorSettings->propagator_ = propagators::gauss_modified_equinoctial;


//    // Retrieve bang-bang thrust acceleration model from MEE costates.
//    std::shared_ptr< propulsion::ThrustAcceleration > thrustAcceleration = getMEEcostatesBasedThrustAccelerationModel( );


    // Acceleration from the central body.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
    bodyToPropagateAccelerations[ centralBody_ ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );
    bodyToPropagateAccelerations[ bodyToPropagate_ ].push_back( getMEEcostatesBasedThrustAccelerationSettings( ) );

    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ bodyToPropagate_ ] = bodyToPropagateAccelerations;

    // Create the acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap_, accelerationMap, std::vector< std::string >{ bodyToPropagate_ }, std::vector< std::string >{ centralBody_ } );

//    accelerationModelMap[ bodyToPropagate_ ][ bodyToPropagate_ ].push_back( thrustAcceleration );

    // Create mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModel;
    massRateModel[ bodyToPropagate_ ] = createMassRateModel( bodyToPropagate_, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
                                                       bodyMap_, accelerationModelMap );

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::modified_equinocial_state_dependent_variable,  bodyToPropagate_, centralBody_ ) );
//    dependentVariables.push_back( std::make_shared< propagators::SingleAccelerationDependentVariableSaveSettings >(
//                        basic_astrodynamics::thrust_acceleration, bodyToPropagate_, bodyToPropagate_, 0 ) );
    dependentVariables.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::total_acceleration_dependent_variable, bodyToPropagate_ ) );

    // Create object with list of dependent variables
    std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< propagators::DependentVariableSaveSettings >( dependentVariables );


    bool hasTimeOfFlightBeenReached = false;
    Eigen::Vector6d initialModifiedEquinoctialElements = convertCartesianToModifiedEquinoctialElements( stateAtDeparture_, centralBodyGravitationalParameter_ );
    double initialTrueLongitude = initialModifiedEquinoctialElements[ trueLongitudeIndex ];

    // Compute intermediate true longitude = initialTrueLongitude + pi (modulo 2pi).
    double intermediateTrueLongitude = initialModifiedEquinoctialElements[ trueLongitudeIndex ]
            + mathematical_constants::PI;
    if ( intermediateTrueLongitude > 2.0 * mathematical_constants::PI )
    {
        intermediateTrueLongitude += - 2.0 * mathematical_constants::PI;
    }

    std::map< double, Eigen::VectorXd > results;
    std::map< double, Eigen::VectorXd > dependentVariableSolution;
    Eigen::VectorXd propagationResult; propagationResult.resize( 7 );
    propagationResult.segment( 0, 6 ) = stateAtDeparture_;
    propagationResult[ 6 ] = initialSpacecraftMass_;
    double currentInitialPropagationTime = 0.0;

    // Ensure that the propagation stops when the required time of flight is required.
    std::shared_ptr< propagators::PropagationTimeTerminationSettings > timeOfFlightTerminationSettings
            = std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight_, false );

    std::function< Eigen::Vector6d( ) > vehicleStateFunction = std::bind(
                &simulation_setup::Body::getState, bodyMap_.at( bodyToPropagate_ ) );

    std::function< Eigen::Vector6d( ) > centralBodyStateFunction = std::bind(
                &simulation_setup::Body::getState, bodyMap_.at( centralBody_ ) );

    int numberOfCompletedRevolutions = 0;
    std::map< double, Eigen::VectorXd > stateHistoryAfterOneRevolution;
    std::map< double, Eigen::VectorXd > dependentVariablesHistoryAfterOneRevolution;

    bool isOrbitalAveragingUsed = false;


    while( !hasTimeOfFlightBeenReached )
    {

        if ( ( !isOrbitalAveragingUsed ) || ( numberOfCompletedRevolutions == 0 ) )
        {

        // Create custom termination settings.
        std::function< bool( const double, const std::function< Eigen::Vector6d( ) >,
                             const std::function< Eigen::Vector6d( ) > ) > customTerminationFunction = [ = ]
                ( const double currentTime, const std::function< Eigen::Vector6d( ) > getSpacecraftState,
                const std::function< Eigen::Vector6d( ) > getCentralBodyState )
        {
          bool terminationCondition = false;

//          Eigen::Vector6d initialModifiedEquinoctialElements = convertCartesianToModifiedEquinoctialElements(
//                      stateAtDeparture_, centralBodyGravitationalParameter_ );

          Eigen::Vector6d currentState = getSpacecraftState( ) - getCentralBodyState( );

          Eigen::Vector6d currentModifiedEquinoctialElements = convertCartesianToModifiedEquinoctialElements(
                      currentState, centralBodyGravitationalParameter_ );

//          double initialTrueLongitude = initialModifiedEquinoctialElements[ trueLongitudeIndex ];
          double currentTrueLongitude = currentModifiedEquinoctialElements[ trueLongitudeIndex ];

          if ( ( initialTrueLongitude >= mathematical_constants::PI )
               && ( currentTrueLongitude < initialTrueLongitude ) && ( currentTrueLongitude >= intermediateTrueLongitude ) )
          {
                terminationCondition = true;
          }
          if ( ( initialTrueLongitude < mathematical_constants::PI )
               && ( currentTrueLongitude >= intermediateTrueLongitude ) )
          {
              terminationCondition = true;
          }
          if ( currentTime > timeOfFlight_ )
          {
              terminationCondition = true;
          }
          return terminationCondition;
        };


        std::shared_ptr< propagators::PropagationTerminationSettings > customTerminationSettings =
                std::make_shared< propagators::PropagationCustomTerminationSettings >(
                    std::bind( customTerminationFunction, std::placeholders::_1, vehicleStateFunction, centralBodyStateFunction ) ) ;

//        std::vector< std::shared_ptr< propagators::PropagationTerminationSettings > > terminationSettingsList;

//        // Add time of flight termination condition to the list of termination settings.
//        terminationSettingsList.push_back( timeOfFlightTerminationSettings );




//        std::shared_ptr< propagators::PropagationDependentVariableTerminationSettings > trueLongitudeTerminationCondition
//                = std::make_shared< propagators::PropagationDependentVariableTerminationSettings >(
//                    std::make_shared< propagators::SingleDependentVariableSaveSettings >(
//                    propagators::modified_equinocial_state_dependent_variable, bodyToPropagate_, centralBody_,
//                    trueLongitudeIndex ), intermediateTrueLongitude, false ); //, true, std::make_shared< root_finders::RootFinderSettings >(
////                    root_finders::bisection_root_finder, 1.0E-6, 400 ) );


//        // Add true longitude termination condition to the list of termination settings.
//        terminationSettingsList.push_back( trueLongitudeTerminationCondition );


//        std::shared_ptr< propagators::PropagationTerminationSettings > terminationSettings
//                = std::make_shared< propagators::PropagationHybridTerminationSettings >( terminationSettingsList, 1 );

        // Re-define propagator settings.
        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > translationalStatePropagatorSettings =
                std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >(
                    std::vector< std::string >{ centralBody_ }, accelerationModelMap, std::vector< std::string >{ bodyToPropagate_ },
                    propagationResult.segment( 0, 6 ), customTerminationSettings, propagators::gauss_modified_equinoctial, dependentVariablesToSave/*,
                    propagatorSettings->getDependentVariablesToSave( )*/ );

//        // Create mass rate models
//        std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModel;
//        massRateModel[ bodyToPropagate_ ] = createMassRateModel( bodyToPropagate_, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
//                                                           bodyMap_, accelerationModelMap );

        // Create settings for propagating the mass of the vehicle.
        std::shared_ptr< propagators::MassPropagatorSettings< double > > massPropagatorSettings
                = std::make_shared< propagators::MassPropagatorSettings< double > >(
                    std::vector< std::string >{ bodyToPropagate_ }, massRateModel,
                    ( Eigen::Matrix< double, 1, 1 >( ) << bodyMap_[ bodyToPropagate_ ]->getBodyMass( ) ).finished( ),
                    customTerminationSettings );

        // Create list of propagation settings.
        std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
        propagatorSettingsVector.push_back( translationalStatePropagatorSettings );
        propagatorSettingsVector.push_back( massPropagatorSettings );

        // Hybrid propagation settings.
        std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings =
                std::make_shared< propagators::MultiTypePropagatorSettings< double > >( propagatorSettingsVector, customTerminationSettings,
                                                                                        dependentVariablesToSave );

        integratorSettings->initialTime_ = currentInitialPropagationTime;

        // Propagate the trajectory.
        propagators::SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap_, integratorSettings, propagatorSettings );
        std::map< double, Eigen::VectorXd > intermediateResults = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        propagationResult = intermediateResults.rbegin( )->second;
        currentInitialPropagationTime = intermediateResults.rbegin( )->first;

        // Retrieve dependent variables history.
        std::map< double, Eigen::VectorXd > intermediateDependentVariableSolution = dynamicsSimulator.getDependentVariableHistory( );

        for ( std::map< double, Eigen::VectorXd >::iterator itr = intermediateResults.begin( ) ; itr != intermediateResults.end( ) ; itr++ )
        {
            results[ itr->first ] = itr->second;
            dependentVariableSolution[ itr->first ] = intermediateDependentVariableSolution[ itr->first ];
        }


        if ( currentInitialPropagationTime >= timeOfFlight_ )
        {
            hasTimeOfFlightBeenReached = true;
        }
        else
        {






//        // Modify termination conditions settings.
//        terminationSettingsList.clear( );

//        // Add time of flight termination condition to the list of termination settings.
//        terminationSettingsList.push_back( timeOfFlightTerminationSettings );

//        trueLongitudeTerminationCondition = std::make_shared< propagators::PropagationDependentVariableTerminationSettings >(
//                    std::make_shared< propagators::SingleDependentVariableSaveSettings >(
//                    propagators::modified_equinocial_state_dependent_variable, bodyToPropagate_, centralBody_,
//                    trueLongitudeIndex ), initialModifiedEquinoctialElements[ trueLongitudeIndex ], false ); //, true,
////                    std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0E-6, 400 ) );


            // Create custom termination settings.
/*            std::function< bool( const double, const std::function< Eigen::Vector6d( ) >,
                                 const std::function< Eigen::Vector6d( ) > ) >*/ customTerminationFunction = [ = ]
                    ( const double currentTime, const std::function< Eigen::Vector6d( ) > getSpacecraftState,
                    const std::function< Eigen::Vector6d( ) > getCentralBodyState )
            {
              bool terminationCondition = false;

//              Eigen::Vector6d initialModifiedEquinoctialElements = convertCartesianToModifiedEquinoctialElements(
//                          stateAtDeparture_, centralBodyGravitationalParameter_ );

              Eigen::Vector6d currentState = getSpacecraftState( ) - getCentralBodyState( );

              Eigen::Vector6d currentModifiedEquinoctialElements = convertCartesianToModifiedEquinoctialElements(
                          currentState, centralBodyGravitationalParameter_ );

//              double initialTrueLongitude = initialModifiedEquinoctialElements[ trueLongitudeIndex ];
              double currentTrueLongitude = currentModifiedEquinoctialElements[ trueLongitudeIndex ];

              if ( ( intermediateTrueLongitude >= mathematical_constants::PI )
                   && ( currentTrueLongitude < intermediateTrueLongitude ) && ( currentTrueLongitude >= initialTrueLongitude ) )
              {
                    terminationCondition = true;
              }
              if ( ( intermediateTrueLongitude < mathematical_constants::PI )
                   && ( currentTrueLongitude >= initialTrueLongitude ) )
              {
                  terminationCondition = true;
              }
              if ( currentTime > timeOfFlight_ )
              {
                  terminationCondition = true;
              }
              return terminationCondition;
            };
//            std::shared_ptr< propagators::PropagationTerminationSettings > customTerminationSettings =
//                    std::make_shared< propagators::PropagationCustomTerminationSettings >(
//                        std::bind( &customTerminationFunction, std::placeholders::_1, cartesianStateFunction ) );

//            std::function< Eigen::Vector6d( ) > spacecraftStateFunction =
//                    std::bind( &simulation_setup::Body::getState, bodyMap_.at( bodyToPropagate_ ) );


            /*std::shared_ptr< propagators::PropagationTerminationSettings >*/ customTerminationSettings =
                    std::make_shared< propagators::PropagationCustomTerminationSettings >(
                        std::bind( customTerminationFunction, std::placeholders::_1, vehicleStateFunction, centralBodyStateFunction ) ) ;
//            std::function< bool( const double ) > testFunction =
//                        std::bind( customTerminationFunction, std::placeholders::_1, spacecraftStateFunction /*std::placeholders::_2*/ ) ; // spacecraftStateFunction ) );
//            std::shared_ptr< propagators::PropagationTerminationSettings > customTerminationSettings =
//                    std::make_shared< propagators::PropagationCustomTerminationSettings >( testFunction ) ;

//            std::cout << "test function: " << testFunction( 0.0 ) << "\n\n";


//        // Add true longitude termination condition to the list of termination settings.
//        terminationSettingsList.push_back( trueLongitudeTerminationCondition );


//        terminationSettings = std::make_shared< propagators::PropagationHybridTerminationSettings >( terminationSettingsList, 1 );

        // Re-define propagator settings.
        translationalStatePropagatorSettings = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >(
                    std::vector< std::string >{ centralBody_ }, accelerationModelMap, std::vector< std::string >{ bodyToPropagate_ },
                    propagationResult.segment( 0, 6 ), customTerminationSettings, propagators::gauss_modified_equinoctial, dependentVariablesToSave/*,
                    propagatorSettings->getDependentVariablesToSave( )*/ );

        std::cout << "initial state second part of the propagation: " << propagationResult.segment( 0, 6 ) << "\n\n";


        // Create settings for propagating the mass of the vehicle.
        massPropagatorSettings = std::make_shared< propagators::MassPropagatorSettings< double > >(
                    std::vector< std::string >{ bodyToPropagate_ }, massRateModel,
                    ( Eigen::Matrix< double, 1, 1 >( ) << bodyMap_[ bodyToPropagate_ ]->getBodyMass( ) ).finished( ), customTerminationSettings );

        // Create list of propagation settings.
        propagatorSettingsVector.clear( );
        propagatorSettingsVector.push_back( translationalStatePropagatorSettings );
        propagatorSettingsVector.push_back( massPropagatorSettings );

        // Hybrid propagation settings.
        propagatorSettings = std::make_shared< propagators::MultiTypePropagatorSettings< double > >( propagatorSettingsVector, customTerminationSettings,
                                                                                        dependentVariablesToSave );

        integratorSettings->initialTime_ = currentInitialPropagationTime; //  results.rbegin( )->first;
        std::cout << "initial time second part of the propagation: " << results.rbegin( )->first << "\n\n";

        // Propagate the trajectory.
        propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorSecondPartOrbit( bodyMap_, integratorSettings, propagatorSettings );
        intermediateResults = dynamicsSimulatorSecondPartOrbit.getEquationsOfMotionNumericalSolution( );
        propagationResult = dynamicsSimulatorSecondPartOrbit.getEquationsOfMotionNumericalSolution().rbegin()->second;
//        propagationResult = intermediateResults.rbegin( )->second;
        currentInitialPropagationTime = intermediateResults.rbegin( )->first;

        // Retrieve dependent variables history.
        intermediateDependentVariableSolution = dynamicsSimulatorSecondPartOrbit.getDependentVariableHistory( );

        for ( std::map< double, Eigen::VectorXd >::iterator itr = intermediateResults.begin( ) ; itr != intermediateResults.end( ) ; itr++ )
        {
            results[ itr->first ] = itr->second;
            dependentVariableSolution[ itr->first ] = intermediateDependentVariableSolution[ itr->first ];
        }

        if ( numberOfCompletedRevolutions == 0 )
        {
            stateHistoryAfterOneRevolution = results;
            dependentVariablesHistoryAfterOneRevolution = dependentVariableSolution;
        }

        // Check if time of flight has been reached.
        if ( currentInitialPropagationTime >= timeOfFlight_ )
        {
            hasTimeOfFlightBeenReached = true;
        }
        else
        {
            numberOfCompletedRevolutions++;
        }



//        hasTimeOfFlightBeenReached = true;

        }


    }


    // Use orbital averaging to propagate the trajectory.
    else
    {

        // Compute averaged state derivative.
        Eigen::Vector7d averagedStateDerivative = computeAveragedStateDerivative(
                    stateHistoryAfterOneRevolution, dependentVariablesHistoryAfterOneRevolution );

        // Compute number of steps to propagate until time of flight.
        int numberOfSteps = ( timeOfFlight_ - currentInitialPropagationTime ) / integratorSettings->initialTimeStep_ + 1;

        std::cout << "number of steps until TOF: " << numberOfSteps << "\n\n";
        std::cout << "( timeOfFlight_ - currentInitialPropagationTime ): " << ( timeOfFlight_ - currentInitialPropagationTime ) << "\n\n";

        // Define current state.
        Eigen::Vector7d currentState;
        Eigen::Vector6d currentCartesianState = propagationResult.segment( 0, 6 );
        currentState.segment( 0, 6 ) = convertCartesianToModifiedEquinoctialElements(
                    currentCartesianState, centralBodyGravitationalParameter_ );
        currentState[ 6 ] = propagationResult[ 6 ];

        // Propagate over required number of steps.
        propagationResult = currentState + averagedStateDerivative * integratorSettings->initialTimeStep_ * numberOfSteps;
        std::cout << "propagation results after OA: " << propagationResult << "\n\n";
        std::cout << "averaged state derivative: " << averagedStateDerivative << "\n\n";
        std::cout << "current state: " << currentState << "\n\n";
        std::cout << "delta current state: " << averagedStateDerivative * integratorSettings->initialTimeStep_ * numberOfSteps << "\n\n";

        hasTimeOfFlightBeenReached = true;

    }


    }


        std::cout << "number of completed revolutions: " << numberOfCompletedRevolutions << "\n\n";
        std::cout << "initial time step: " << integratorSettings->initialTimeStep_ << "\n\n";


//    std::vector< std::shared_ptr< propagators::PropagationTerminationSettings > > terminationSettingsList;

//    terminationSettingsList.push_back( std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight_, true ) );


//    double intermediateTrueLongitude = initialModifiedEquinoctialElements[ trueLongitudeIndex ]
//            + mathematical_constants::PI;
//    if ( intermediateTrueLongitude > 2.0 * mathematical_constants::PI )
//    {
//        intermediateTrueLongitude += - 2.0 * mathematical_constants::PI;
//    }



//    std::shared_ptr< propagators::PropagationDependentVariableTerminationSettings > trueLongitudeTerminationCondition
//            = std::make_shared< propagators::PropagationDependentVariableTerminationSettings >(
//                std::make_shared< propagators::SingleDependentVariableSaveSettings >(
//                propagators::modified_equinocial_state_dependent_variable, bodyToPropagate_, centralBody_,
//                trueLongitudeIndex ), intermediateTrueLongitude, false, true, std::make_shared< root_finders::RootFinderSettings >(
//                root_finders::bisection_root_finder, 1.0E-6, 100 ) );

//    terminationSettingsList.push_back( std::make_shared< propagators::PropagationDependentVariableTerminationSettings >(
//                                           std::make_shared< propagators::SingleDependentVariableSaveSettings >(
//                                               propagators::modified_equinocial_state_dependent_variable,
//                                               bodyToPropagate_, centralBody_, trueLongitudeIndex ),
//                                           intermediateTrueLongitude
//                                           /*initialModifiedEquinoctialElements[ trueLongitudeIndex ]*/, false, true,
//                                           std::make_shared< root_finders::RootFinderSettings >(
//                                               root_finders::bisection_root_finder, 1.0E-6, 100 ) ) );



//    std::shared_ptr< propagators::PropagationTerminationSettings > terminationSettings
//            = std::make_shared< propagators::PropagationHybridTerminationSettings >( terminationSettingsList, 1 );

//    // Re-define propagator settings.
//    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > translationalStatePropagatorSettings =
//            std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >(
//                std::vector< std::string >{ centralBody_ }, accelerationModelMap, std::vector< std::string >{ bodyToPropagate_ },
//                stateAtDeparture_, terminationSettings, propagators::gauss_modified_equinoctial, dependentVariablesToSave/*,
//                propagatorSettings->getDependentVariablesToSave( )*/ );

//    // Create mass rate models
//    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModel;
//    massRateModel[ bodyToPropagate_ ] = createMassRateModel( bodyToPropagate_, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
//                                                       bodyMap_, accelerationModelMap );

//    // Create settings for propagating the mass of the vehicle.
//    std::shared_ptr< propagators::MassPropagatorSettings< double > > massPropagatorSettings
//            = std::make_shared< propagators::MassPropagatorSettings< double > >(
//                std::vector< std::string >{ bodyToPropagate_ }, massRateModel,
//                ( Eigen::Matrix< double, 1, 1 >( ) << initialSpacecraftMass_ /*bodyMap_[ bodyToPropagate_ ]->getBodyMass( )*/ ).finished( ),
//                terminationSettings );

//    // Create list of propagation settings.
//    std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
//    propagatorSettingsVector.push_back( translationalStatePropagatorSettings );
//    propagatorSettingsVector.push_back( massPropagatorSettings );

//    // Backward hybrid propagation settings.
//    std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings =
//            std::make_shared< propagators::MultiTypePropagatorSettings< double > >( propagatorSettingsVector, terminationSettings,
//                                                                                    dependentVariablesToSave );

//    // Propagate the trajectory.
//    propagators::SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap_, integratorSettings, propagatorSettings );
//    std::map< double, Eigen::VectorXd > results = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
//    Eigen::VectorXd propagationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution().rbegin()->second;

    // Retrieve state and mass of the spacecraft at the end of the propagation.
    Eigen::Vector6d propagatedState = propagationResult.segment( 0, 6 );
    massAtTimeOfFlight_ = propagationResult[ 6 ];

    // Compute deltaV required by the trajectory.
    totalDeltaV_ = computeTotalDeltaV( );

//    // Retrieve dependent variables history.
//    std::map< double, Eigen::VectorXd > dependentVariableSolution = dynamicsSimulator.getDependentVariableHistory( );

    input_output::writeDataMapToTextFile( dependentVariableSolution,
                                          "dependentVariablesHybridMethod.dat",
                                          "C:/Users/chamb/Documents/Master_2/SOCIS/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( results,
                                          "propagationResultHybridMethod.dat",
                                          "C:/Users/chamb/Documents/Master_2/SOCIS/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    return propagatedState;
}


//! Propagate the spacecraft trajectory to a given time.
Eigen::Vector6d HybridMethodLeg::propagateTrajectory( double initialTime, double finalTime, Eigen::Vector6d initialState, double initialMass,
                                     std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    // Re-initialise integrator settings.
    integratorSettings->initialTime_ = initialTime;

    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialMass );

    std::cout << "initial mass in propagate trajectory to set of epochs: " << bodyMap_[ bodyToPropagate_ ]->getBodyMass() << "\n\n";

//    propagatorSettings->resetInitialStates( stateAtDeparture_ );
//    propagatorSettings->bodiesToIntegrate_ = std::vector< std::string >{ bodyToPropagate_ };
//    propagatorSettings->centralBodies_ = std::vector< std::string >{ centralBody_ };
//    propagatorSettings->propagator_ = propagators::gauss_modified_equinoctial;


//    // Retrieve bang-bang thrust acceleration model from MEE costates.
//    std::shared_ptr< propulsion::ThrustAcceleration > thrustAcceleration = getMEEcostatesBasedThrustAccelerationModel( );


    // Acceleration from the central body.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
    bodyToPropagateAccelerations[ centralBody_ ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );
    bodyToPropagateAccelerations[ bodyToPropagate_ ].push_back( getMEEcostatesBasedThrustAccelerationSettings( ) );

    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ bodyToPropagate_ ] = bodyToPropagateAccelerations;

    // Create the acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap_, accelerationMap, std::vector< std::string >{ bodyToPropagate_ }, std::vector< std::string >{ centralBody_ } );

//    accelerationModelMap[ bodyToPropagate_ ][ bodyToPropagate_ ].push_back( thrustAcceleration );

    // Create mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModel;
    massRateModel[ bodyToPropagate_ ] = createMassRateModel( bodyToPropagate_, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
                                                       bodyMap_, accelerationModelMap );

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::modified_equinocial_state_dependent_variable,  bodyToPropagate_, centralBody_ ) );
//    dependentVariables.push_back( std::make_shared< propagators::SingleAccelerationDependentVariableSaveSettings >(
//                        basic_astrodynamics::thrust_acceleration, bodyToPropagate_, bodyToPropagate_, 0 ) );
    dependentVariables.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::total_acceleration_dependent_variable, bodyToPropagate_ ) );

    // Create object with list of dependent variables
    std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< propagators::DependentVariableSaveSettings >( dependentVariables );


    bool hasTimeOfFlightBeenReached = false;
    Eigen::Vector6d initialModifiedEquinoctialElements = convertCartesianToModifiedEquinoctialElements( stateAtDeparture_, centralBodyGravitationalParameter_ );
    double initialTrueLongitude = initialModifiedEquinoctialElements[ trueLongitudeIndex ];

    // Compute intermediate true longitude = initialTrueLongitude + pi (modulo 2pi).
    double intermediateTrueLongitude = initialModifiedEquinoctialElements[ trueLongitudeIndex ]
            + mathematical_constants::PI;
    if ( intermediateTrueLongitude > 2.0 * mathematical_constants::PI )
    {
        intermediateTrueLongitude += - 2.0 * mathematical_constants::PI;
    }

    std::map< double, Eigen::VectorXd > results;
    std::map< double, Eigen::VectorXd > dependentVariableSolution;
    Eigen::VectorXd propagationResult; propagationResult.resize( 7 );
    propagationResult.segment( 0, 6 ) = initialState;
    propagationResult[ 6 ] = initialMass;
    double currentInitialPropagationTime = initialTime;

    // Ensure that the propagation stops when the required time of flight is required.
    std::shared_ptr< propagators::PropagationTimeTerminationSettings > timeOfFlightTerminationSettings
            = std::make_shared< propagators::PropagationTimeTerminationSettings >( finalTime, false );

    std::function< Eigen::Vector6d( ) > vehicleStateFunction = std::bind(
                &simulation_setup::Body::getState, bodyMap_.at( bodyToPropagate_ ) );

    std::function< Eigen::Vector6d( ) > centralBodyStateFunction = std::bind(
                &simulation_setup::Body::getState, bodyMap_.at( centralBody_ ) );

    int numberOfCompletedRevolutions = 0;
    std::map< double, Eigen::VectorXd > stateHistoryAfterOneRevolution;
    std::map< double, Eigen::VectorXd > dependentVariablesHistoryAfterOneRevolution;

    bool isOrbitalAveragingUsed = false;


    while( !hasTimeOfFlightBeenReached )
    {

        if ( ( !isOrbitalAveragingUsed ) || ( numberOfCompletedRevolutions == 0 ) )
        {

        // Create custom termination settings.
        std::function< bool( const double, const std::function< Eigen::Vector6d( ) >,
                             const std::function< Eigen::Vector6d( ) > ) > customTerminationFunction = [ = ]
                ( const double currentTime, const std::function< Eigen::Vector6d( ) > getSpacecraftState,
                const std::function< Eigen::Vector6d( ) > getCentralBodyState )
        {
          bool terminationCondition = false;

//          Eigen::Vector6d initialModifiedEquinoctialElements = convertCartesianToModifiedEquinoctialElements(
//                      stateAtDeparture_, centralBodyGravitationalParameter_ );

          Eigen::Vector6d currentState = getSpacecraftState( ) - getCentralBodyState( );

          Eigen::Vector6d currentModifiedEquinoctialElements = convertCartesianToModifiedEquinoctialElements(
                      currentState, centralBodyGravitationalParameter_ );

//          double initialTrueLongitude = initialModifiedEquinoctialElements[ trueLongitudeIndex ];
          double currentTrueLongitude = currentModifiedEquinoctialElements[ trueLongitudeIndex ];

          if ( ( initialTrueLongitude >= mathematical_constants::PI )
               && ( currentTrueLongitude < initialTrueLongitude ) && ( currentTrueLongitude >= intermediateTrueLongitude ) )
          {
                terminationCondition = true;
          }
          if ( ( initialTrueLongitude < mathematical_constants::PI )
               && ( currentTrueLongitude >= intermediateTrueLongitude ) )
          {
              terminationCondition = true;
          }
          if ( currentTime >= finalTime )
          {
              terminationCondition = true;
          }
          return terminationCondition;
        };


        std::shared_ptr< propagators::PropagationTerminationSettings > customTerminationSettings =
                std::make_shared< propagators::PropagationCustomTerminationSettings >(
                    std::bind( customTerminationFunction, std::placeholders::_1, vehicleStateFunction, centralBodyStateFunction ) ) ;


        // Re-define propagator settings.
        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > translationalStatePropagatorSettings =
                std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >(
                    std::vector< std::string >{ centralBody_ }, accelerationModelMap, std::vector< std::string >{ bodyToPropagate_ },
                    propagationResult.segment( 0, 6 ), customTerminationSettings, propagators::gauss_modified_equinoctial, dependentVariablesToSave/*,
                    propagatorSettings->getDependentVariablesToSave( )*/ );

//        // Create mass rate models
//        std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModel;
//        massRateModel[ bodyToPropagate_ ] = createMassRateModel( bodyToPropagate_, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
//                                                           bodyMap_, accelerationModelMap );

        // Create settings for propagating the mass of the vehicle.
        std::shared_ptr< propagators::MassPropagatorSettings< double > > massPropagatorSettings
                = std::make_shared< propagators::MassPropagatorSettings< double > >(
                    std::vector< std::string >{ bodyToPropagate_ }, massRateModel,
                    ( Eigen::Matrix< double, 1, 1 >( ) << bodyMap_[ bodyToPropagate_ ]->getBodyMass( ) ).finished( ),
                    customTerminationSettings );

        // Create list of propagation settings.
        std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
        propagatorSettingsVector.push_back( translationalStatePropagatorSettings );
        propagatorSettingsVector.push_back( massPropagatorSettings );

        // Hybrid propagation settings.
        std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings =
                std::make_shared< propagators::MultiTypePropagatorSettings< double > >( propagatorSettingsVector, customTerminationSettings,
                                                                                        dependentVariablesToSave );

        integratorSettings->initialTime_ = currentInitialPropagationTime;

        // Propagate the trajectory.
        propagators::SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap_, integratorSettings, propagatorSettings );
        std::map< double, Eigen::VectorXd > intermediateResults = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        propagationResult = intermediateResults.rbegin( )->second;
        currentInitialPropagationTime = intermediateResults.rbegin( )->first;

        // Retrieve dependent variables history.
        std::map< double, Eigen::VectorXd > intermediateDependentVariableSolution = dynamicsSimulator.getDependentVariableHistory( );

        for ( std::map< double, Eigen::VectorXd >::iterator itr = intermediateResults.begin( ) ; itr != intermediateResults.end( ) ; itr++ )
        {
            results[ itr->first ] = itr->second;
            dependentVariableSolution[ itr->first ] = intermediateDependentVariableSolution[ itr->first ];
        }


        if ( currentInitialPropagationTime >= finalTime )
        {
            hasTimeOfFlightBeenReached = true;
        }
        else
        {


            // Create custom termination settings.
/*            std::function< bool( const double, const std::function< Eigen::Vector6d( ) >,
                                 const std::function< Eigen::Vector6d( ) > ) >*/ customTerminationFunction = [ = ]
                    ( const double currentTime, const std::function< Eigen::Vector6d( ) > getSpacecraftState,
                    const std::function< Eigen::Vector6d( ) > getCentralBodyState )
            {
              bool terminationCondition = false;

//              Eigen::Vector6d initialModifiedEquinoctialElements = convertCartesianToModifiedEquinoctialElements(
//                          stateAtDeparture_, centralBodyGravitationalParameter_ );

              Eigen::Vector6d currentState = getSpacecraftState( ) - getCentralBodyState( );

              Eigen::Vector6d currentModifiedEquinoctialElements = convertCartesianToModifiedEquinoctialElements(
                          currentState, centralBodyGravitationalParameter_ );

//              double initialTrueLongitude = initialModifiedEquinoctialElements[ trueLongitudeIndex ];
              double currentTrueLongitude = currentModifiedEquinoctialElements[ trueLongitudeIndex ];

              if ( ( intermediateTrueLongitude >= mathematical_constants::PI )
                   && ( currentTrueLongitude < intermediateTrueLongitude ) && ( currentTrueLongitude >= initialTrueLongitude ) )
              {
                    terminationCondition = true;
              }
              if ( ( intermediateTrueLongitude < mathematical_constants::PI )
                   && ( currentTrueLongitude >= initialTrueLongitude ) )
              {
                  terminationCondition = true;
              }
              if ( currentTime >= finalTime )
              {
                  terminationCondition = true;
              }
              return terminationCondition;
            };


            /*std::shared_ptr< propagators::PropagationTerminationSettings >*/ customTerminationSettings =
                    std::make_shared< propagators::PropagationCustomTerminationSettings >(
                        std::bind( customTerminationFunction, std::placeholders::_1, vehicleStateFunction, centralBodyStateFunction ) ) ;

        // Re-define propagator settings.
        translationalStatePropagatorSettings = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >(
                    std::vector< std::string >{ centralBody_ }, accelerationModelMap, std::vector< std::string >{ bodyToPropagate_ },
                    propagationResult.segment( 0, 6 ), customTerminationSettings, propagators::gauss_modified_equinoctial, dependentVariablesToSave/*,
                    propagatorSettings->getDependentVariablesToSave( )*/ );

        std::cout << "initial state second part of the propagation: " << propagationResult.segment( 0, 6 ) << "\n\n";


        // Create settings for propagating the mass of the vehicle.
        massPropagatorSettings = std::make_shared< propagators::MassPropagatorSettings< double > >(
                    std::vector< std::string >{ bodyToPropagate_ }, massRateModel,
                    ( Eigen::Matrix< double, 1, 1 >( ) << bodyMap_[ bodyToPropagate_ ]->getBodyMass( ) ).finished( ), customTerminationSettings );

        // Create list of propagation settings.
        propagatorSettingsVector.clear( );
        propagatorSettingsVector.push_back( translationalStatePropagatorSettings );
        propagatorSettingsVector.push_back( massPropagatorSettings );

        // Hybrid propagation settings.
        propagatorSettings = std::make_shared< propagators::MultiTypePropagatorSettings< double > >( propagatorSettingsVector, customTerminationSettings,
                                                                                        dependentVariablesToSave );

        integratorSettings->initialTime_ = currentInitialPropagationTime; //  results.rbegin( )->first;
        std::cout << "initial time second part of the propagation: " << results.rbegin( )->first << "\n\n";

        // Propagate the trajectory.
        propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorSecondPartOrbit( bodyMap_, integratorSettings, propagatorSettings );
        intermediateResults = dynamicsSimulatorSecondPartOrbit.getEquationsOfMotionNumericalSolution( );
        propagationResult = dynamicsSimulatorSecondPartOrbit.getEquationsOfMotionNumericalSolution().rbegin()->second;
//        propagationResult = intermediateResults.rbegin( )->second;
        currentInitialPropagationTime = intermediateResults.rbegin( )->first;

        // Retrieve dependent variables history.
        intermediateDependentVariableSolution = dynamicsSimulatorSecondPartOrbit.getDependentVariableHistory( );

        for ( std::map< double, Eigen::VectorXd >::iterator itr = intermediateResults.begin( ) ; itr != intermediateResults.end( ) ; itr++ )
        {
            results[ itr->first ] = itr->second;
            dependentVariableSolution[ itr->first ] = intermediateDependentVariableSolution[ itr->first ];
        }

        if ( numberOfCompletedRevolutions == 0 )
        {
            stateHistoryAfterOneRevolution = results;
            dependentVariablesHistoryAfterOneRevolution = dependentVariableSolution;
        }

        // Check if time of flight has been reached.
        if ( currentInitialPropagationTime >= finalTime )
        {
            hasTimeOfFlightBeenReached = true;
        }
        else
        {
            numberOfCompletedRevolutions++;
        }



//        hasTimeOfFlightBeenReached = true;

        }


    }


    // Use orbital averaging to propagate the trajectory.
    else
    {

        // Compute averaged state derivative.
        Eigen::Vector7d averagedStateDerivative = computeAveragedStateDerivative(
                    stateHistoryAfterOneRevolution, dependentVariablesHistoryAfterOneRevolution );

        // Compute number of steps to propagate until time of flight.
        int numberOfSteps = ( finalTime - currentInitialPropagationTime ) / integratorSettings->initialTimeStep_ + 1;

        std::cout << "number of steps until TOF: " << numberOfSteps << "\n\n";
        std::cout << "( finalTime - currentInitialPropagationTime ): " << ( finalTime - currentInitialPropagationTime ) << "\n\n";

        // Define current state.
        Eigen::Vector7d currentState;
        Eigen::Vector6d currentCartesianState = propagationResult.segment( 0, 6 );
        currentState.segment( 0, 6 ) = convertCartesianToModifiedEquinoctialElements(
                    currentCartesianState, centralBodyGravitationalParameter_ );
        currentState[ 6 ] = propagationResult[ 6 ];

        // Propagate over required number of steps.
        propagationResult = currentState + averagedStateDerivative * integratorSettings->initialTimeStep_ * numberOfSteps;
        std::cout << "propagation results after OA: " << propagationResult << "\n\n";
        std::cout << "averaged state derivative: " << averagedStateDerivative << "\n\n";
        std::cout << "current state: " << currentState << "\n\n";
        std::cout << "delta current state: " << averagedStateDerivative * integratorSettings->initialTimeStep_ * numberOfSteps << "\n\n";

        hasTimeOfFlightBeenReached = true;

    }


    }


        std::cout << "number of completed revolutions: " << numberOfCompletedRevolutions << "\n\n";




    // Retrieve state and mass of the spacecraft at the end of the propagation.
    Eigen::Vector6d propagatedState = propagationResult.segment( 0, 6 );
//    massAtTimeOfFlight_ = propagationResult[ 6 ];

    // Compute deltaV required by the trajectory.
    totalDeltaV_ = computeTotalDeltaV( );

    std::cout << "current mass in propagate trajectory( initialTime, finalTime, ... ): " <<
                 bodyMap_[ bodyToPropagate_ ]->getBodyMass() << "\n\n";
    std::cout << "current state in propagate trajectory( initialTime, finalTime, ... ): " <<
                 propagatedState << "\n\n";


//    input_output::writeDataMapToTextFile( dependentVariableSolution,
//                                          "dependentVariablesHybridMethod.dat",
//                                          "C:/Users/chamb/Documents/Master_2/SOCIS/",
//                                          "",
//                                          std::numeric_limits< double >::digits10,
//                                          std::numeric_limits< double >::digits10,
//                                          "," );

//    input_output::writeDataMapToTextFile( results,
//                                          "propagationResultHybridMethod.dat",
//                                          "C:/Users/chamb/Documents/Master_2/SOCIS/",
//                                          "",
//                                          std::numeric_limits< double >::digits10,
//                                          std::numeric_limits< double >::digits10,
//                                          "," );

    return propagatedState;
}


//! Propagate the trajectory to set of epochs.
std::map< double, Eigen::Vector6d > HybridMethodLeg::propagateTrajectory(
        std::vector< double > epochs, std::map< double, Eigen::Vector6d >& propagatedTrajectory, Eigen::Vector6d initialState,
        double initialMass, double initialTime, std::shared_ptr<  numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    // Initialise propagated state.
    Eigen::Vector6d propagatedState = initialState;

    // Initialise mass of the spacecraft at departure.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialMass );
    double currentMass = initialMass;

    std::cout << "initial mass in propagate trajectory to set of epochs: " << currentMass << "\n\n";


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
            if ( currentTime >= initialTime )
            {
                propagatedState = propagateTrajectory( initialTime, currentTime, propagatedState, currentMass, integratorSettings );
                currentMass = bodyMap_[ bodyToPropagate_ ]->getBodyMass( );
            }
            propagatedTrajectory[ currentTime ] = propagatedState;
            std::cout << "current time test: " << currentTime << "\n\n";
            std::cout << "initial time test: " << initialTime << "\n\n";
        }
        else
        {
            propagatedState = propagateTrajectory( epochs[ epochIndex - 1 ], currentTime, propagatedState, currentMass, integratorSettings );
            currentMass = bodyMap_[ bodyToPropagate_ ]->getBodyMass( );
            propagatedTrajectory[ currentTime ] = propagatedState;
        }
        std::cout << "current mass in propagateTrajectory( epochs,... ): " << currentMass << "\n\n";
        std::cout << "current state in propagateTrajectory( epochs,... ): " << propagatedState << "\n\n";

    }

    bodyMap_[ centralBody_ ]->setConstantBodyMass( initialMass );

        std::cout << "TEST: " << "\n\n";


    return propagatedTrajectory;
}


//! Return the deltaV associated with the thrust profile of the trajectory.
double HybridMethodLeg::computeTotalDeltaV( )
{

    // Compute (constant) mass rate.
    double massRate = - maximumThrust_ /
            ( specificImpulse_ * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION );

    // Compute time during which the engine was switched on.
    double engineSwitchedOnDuration = ( massAtTimeOfFlight_ - initialSpacecraftMass_ ) / massRate;

    std::shared_ptr< numerical_quadrature::QuadratureSettings< double > > quadratureSettings_ =
            std::make_shared< numerical_quadrature::GaussianQuadratureSettings < double > >( 0.0, 16 );

    // Thrust acceleration function to use quadrature.
    // Define thrust acceleration as a function of time (to be integrated to compute the associated deltaV).
    std::function< double( const double ) > thrustAcceleration = [ = ] ( const double currentTime ){

        // Compute current mass.
        double currentMass = initialSpacecraftMass_ + massRate * currentTime;

        // Compute and return current thrust acceleration.
        double currentThrustAcceleration = maximumThrust_ / currentMass;

        return currentThrustAcceleration;
    };

//    // Thrust acceleration function to use integrator.
//    // Define thrust acceleration as a function of time (to be integrated to compute the associated deltaV).
//    std::function< Eigen::Vector1d( const double, const Eigen::Vector1d& ) > thrustAcceleration = [ = ] ( const double currentTime,
//            const Eigen::Vector1d& independentVariable ){

//        // Compute current mass.
//        double currentMass = initialSpacecraftMass_ + massRate * currentTime;

//        // Compute and return current thrust acceleration.
//        double currentThrustAcceleration = maximumThrust_ / currentMass;

//        return ( Eigen::Vector1d() << currentThrustAcceleration ).finished( );
//    };

    // Create numerical quadrature from quadrature settings.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( thrustAcceleration, quadratureSettings_, engineSwitchedOnDuration );

//    // Redefine integrator settings.
//    integratorSettings->initialTime_ = 0.0;
//    integratorSettings->initialTimeStep_ = std::fabs( integratorSettings->initialTimeStep_ );

//    // Create numerical integrator.
//    std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::Vector1d > > integrator =
//            numerical_integrators::createIntegrator< double, Eigen::Vector1d >( thrustAcceleration, Eigen::Vector1d::Zero( ),
//                                                                                integratorSettings );

//    double deltaV = integrator->integrateTo( engineSwitchedOnDuration, integratorSettings->initialTimeStep_ )[ 0 ];

    // Compute deltaV analytically.
    double deltaV = - specificImpulse_ * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION
            * std::log( 1.0 + massRate / initialSpacecraftMass_ * engineSwitchedOnDuration );
    std::cout << "deltaV computed analytically: " << deltaV << "\n\n";


//    return quadrature->getQuadrature( );
    return deltaV;
}


//! Compute dynamics matrix.
Eigen::Matrix< double, 6, 3 > HybridMethodLeg::computeDynamicsMatrix( Eigen::Vector6d modifiedEquinoctialElements )
{
    // Retrieve each modified equinoctial element independently.
    double p = modifiedEquinoctialElements[ semiParameterIndex ];
    double f = modifiedEquinoctialElements[ fElementIndex ];
    double g = modifiedEquinoctialElements[ gElementIndex ];
    double h = modifiedEquinoctialElements[ hElementIndex ];
    double k = modifiedEquinoctialElements[ kElementIndex ];
    double L = modifiedEquinoctialElements[ trueLongitudeIndex ];

    // Compute auxialiary variables.
    double w1 = 1.0 + f * std::cos( L ) + g * std::sin( L );
    double w2 = 1.0 + h * h + k * k;

    // Compute dynamics matrix.
    Eigen::Matrix< double, 6, 3 > dynamicsMatrix = Eigen::MatrixXd::Zero( 6, 3 );

    dynamicsMatrix( 0, 1 ) = 2.0 * p / w1;
    dynamicsMatrix( 1, 0 ) = std::sin( L );
    dynamicsMatrix( 1, 1 ) = ( 1.0 / w1 ) * ( ( w1 + 1.0 ) * std::cos( L ) + f );
    dynamicsMatrix( 1, 2 ) = ( g / w1 ) * ( h * std::sin( L ) - k * std::cos( L ) );
    dynamicsMatrix( 2, 0 ) = - std::cos( L );
    dynamicsMatrix( 2, 1 ) = ( 1.0 / w1 ) * ( ( w1 + 1.0 ) * std::sin( L ) + g );
    dynamicsMatrix( 2, 2 ) = ( f / w1 ) * ( h * std::sin( L ) - k * std::cos( L ) );
    dynamicsMatrix( 3, 2 ) = ( w2 * std::cos( L ) ) / ( 2.0 * w1 );
    dynamicsMatrix( 4, 2 ) = ( w2 * std::sin( L ) ) / ( 2.0 * w1 );
    dynamicsMatrix( 5, 2 ) = ( 1.0 / w1 ) * ( h * std::sin( L ) - k * std::cos( L ) );

    dynamicsMatrix *= std::sqrt( p / centralBodyGravitationalParameter_ );

    return dynamicsMatrix;


}


//! Compute averaged state derivative.
Eigen::Vector7d HybridMethodLeg::computeAveragedStateDerivative(
        std::map< double, Eigen::VectorXd > stateHistory,
        std::map< double, Eigen::VectorXd > dependentVariableHistory )
{
    std::map< double, Eigen::Vector7d > stateDerivativeMap;

    double timeAfterOneRevolution;

    // Compute state derivative for each epoch.
    for ( std::map< double, Eigen::VectorXd >::iterator itr = dependentVariableHistory.begin( ) ; itr != dependentVariableHistory.end( )-- ; itr++ )
    {
        // Retrieve modified equinoctial elements and inertial perturbing accelerations from set of dependent variables.
        Eigen::Vector6d modifiedEquinoctialElements = itr->second.segment( 0, 6 );
        Eigen::Vector3d perturbingAccelerationsInInertialFrame = itr->second.segment( 6, 3 );

        // Retrieve current mass.
        Eigen::Vector6d cartesianState = stateHistory[ itr->first ].segment( 0, 6 );
        double currentMass = stateHistory[ itr->first ][ 6 ];

        // Retrieve each modified equinoctial element independently.
        double p = modifiedEquinoctialElements[ semiParameterIndex ];
        double f = modifiedEquinoctialElements[ fElementIndex ];
        double g = modifiedEquinoctialElements[ gElementIndex ];
        double h = modifiedEquinoctialElements[ hElementIndex ];
        double k = modifiedEquinoctialElements[ kElementIndex ];
        double L = modifiedEquinoctialElements[ trueLongitudeIndex ];

        // Compute auxialiary variables.
        double w1 = 1.0 + f * std::cos( L ) + g * std::sin( L );
//        double w2 = 1.0 + h * h + k * k;

        // Compute cartesian state.
//        Eigen::Vector6d cartesianState = convertModifiedEquinoctialToCartesianElements( modifiedEquinoctialElements,
//                                                                                        centralBodyGravitationalParameter_, false );

        // Transform perturbing accelerations vector from inertial to RSW frame.
        Eigen::Matrix3d rotationMatrixFromInertialToRsw = reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrix( cartesianState );

        Eigen::Vector3d perturbingAccelerationsInRsw = rotationMatrixFromInertialToRsw * perturbingAccelerationsInInertialFrame;


        // Compute dynamics matrix.
        Eigen::Matrix< double, 6, 3 > dynamicsMatrix = computeDynamicsMatrix( modifiedEquinoctialElements );

        // Compute vector b.
        Eigen::Vector6d b = ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0,
                              std::sqrt( centralBodyGravitationalParameter_ * p ) * ( w1 / p ) * ( w1 / p ) ).finished( );

        // Compute state derivative vector.
        Eigen::Vector6d modifiedEquinoctialElementsDerivatives = dynamicsMatrix * perturbingAccelerationsInRsw + b;
        double massDerivative = - ( perturbingAccelerationsInInertialFrame.norm( ) * currentMass )  /
                ( specificImpulse_ * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION );

        Eigen::Vector7d stateDerivativeVector;
        stateDerivativeVector.segment( 0, 6 ) = modifiedEquinoctialElementsDerivatives;
        stateDerivativeVector[ 6 ] = massDerivative;
//        std::cout << "mass derivative: " << massDerivative << "\n\n";

        stateDerivativeMap[ itr->first ] = stateDerivativeVector;

        timeAfterOneRevolution = itr->first;

    }


    // Compute averaged state derivative.
    Eigen::Vector7d averagedStateDerivative = Eigen::Vector7d::Zero( );

    int counterTest = 0;
    double stepSize;

    for ( std::map< double, Eigen::VectorXd >::iterator nextStepItr = dependentVariableHistory.begin( ) ;
          nextStepItr != dependentVariableHistory.end( )--; nextStepItr++ )
    {
        std::map< double, Eigen::VectorXd >::iterator itr = nextStepItr++;
        if ( nextStepItr != dependentVariableHistory.end( ) )
        {
        averagedStateDerivative += ( nextStepItr->first - itr->first ) * ( stateDerivativeMap[ nextStepItr->first ] + stateDerivativeMap[ itr->first ] );
//        std::cout << "test time stuff: " << nextStepItr->first - itr->first << "\n\n";
//        std::cout << "stateDerivativeMap[ itr->first ]: " << stateDerivativeMap[ itr->first ] << "\n\n";
//        std::cout << "stateDerivativeMap[ nextStepItr->first ]: " << stateDerivativeMap[ nextStepItr->first ] << "\n\n";
        stepSize =  nextStepItr->first - itr->first;
        counterTest++;

        }
        nextStepItr--;
    }

//    std::cout << "counter test: " << counterTest << "\n\n";
//    std::cout << "time after revolution: " << timeAfterOneRevolution / stepSize << "\n\n";
//    std::cout << "time after revolution: " << timeAfterOneRevolution << "\n\n";
//    std::cout << "timeAfterRevolution - initialTime: " << timeAfterOneRevolution - stateHistory.begin( )->first << "\n\n";
//    std::cout << "TEST: " << ( 1.0 / ( 2.0  * ( /*stateHistory.rbegin( )->first*/ - stateHistory.begin( )->first ) ) ) << "\n\n";
    averagedStateDerivative = ( 1.0 / ( 2.0 * ( timeAfterOneRevolution - stateHistory.begin( )->first ) ) ) * averagedStateDerivative;
//            * intermediateSum;

    return averagedStateDerivative;

}


} // namespace low_thrust_direct_methods
} // namespace tudat
