/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <string>
#include <thread>

#include <limits>


#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/desaturationDeltaV.h"


namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_desaturation_deltaV_values_estimation )

//Using declarations.
using namespace tudat::observation_models;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;
using namespace tudat::basic_astrodynamics;

//! Unit test to check if desaturation deltaV values are estimated correctly
BOOST_AUTO_TEST_CASE( test_DesaturationDeltaVsEstimation )
{

    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Mars" );

    // Specify number of observation days.
    int numberOfDaysOfData = 1;

    // Specify initial time
    double initialEphemerisTime = 1.0e7;
    double finalEphemerisTime = initialEphemerisTime + numberOfDaysOfData * 86400.0;

    // Define times and deltaV magnitudes for momentum wheel desaturation maneuvers.
    std::vector< double > thrustMidTimes = { initialEphemerisTime + 1.0 * 3600.0, initialEphemerisTime + 3.0 * 3600.0,
                                             initialEphemerisTime + 5.0 * 3600.0 };
    std::vector< Eigen::Vector3d > deltaVValues =
    { 1.0E-3 * ( Eigen::Vector3d( ) << 0.3, -2.5, 3.4 ).finished( ),
      1.0E-3 * ( Eigen::Vector3d( ) << 2.0, 5.9, -0.5 ).finished( ),
      1.0E-3 * ( Eigen::Vector3d( ) << -1.6, 4.4, -5.8 ).finished( ) };
    double totalManeuverTime = 90.0;
    double maneuverRiseTime = 15.0;


    // Create bodies needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                "ECLIPJ2000", "IAU_Earth", spice_interface::computeRotationQuaternionBetweenFrames(
                    "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
                    initialEphemerisTime, 2.0 * mathematical_constants::PI / ( physical_constants::JULIAN_DAY ) );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( 400.0 );

    // Create aerodynamic coefficient interface settings.
    double referenceArea = 4.0;
    double aerodynamicCoefficient = 1.2;
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >(
                referenceArea, aerodynamicCoefficient * ( Eigen::Vector3d( ) << 1.2, -0.1, -0.4 ).finished( ),
                negative_aerodynamic_frame_coefficients );

    // Create and set aerodynamic coefficients object
    bodies.at( "Vehicle" )->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle", bodies ) );

    // Create radiation pressure settings
    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
    std::shared_ptr< RadiationPressureInterfaceSettings > asterixRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodies.at( "Vehicle" )->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    asterixRadiationPressureSettings, "Vehicle", bodies ) );

    // Create ground stations.
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "Station1" );
    groundStationNames.push_back( "Station2" );
    groundStationNames.push_back( "Station3" );

    createGroundStation( bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ),
                         coordinate_conversions::geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station2", ( Eigen::Vector3d( ) << 0.0, -0.55, 2.0 ).finished( ),
                         coordinate_conversions::geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station3", ( Eigen::Vector3d( ) << 0.0, 0.05, 4.0 ).finished( ),
                         coordinate_conversions::geodetic_position );

    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 8, 8 ) );
    accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                    basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfVehicle[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
                                                    basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::cannon_ball_radiation_pressure ) );
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::aerodynamic ) );
    accelerationsOfVehicle[ "Vehicle" ].push_back( std::make_shared< MomentumWheelDesaturationAccelerationSettings >(
                    thrustMidTimes, deltaVValues, totalManeuverTime, maneuverRiseTime ) );

    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;


    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;
    bodiesToIntegrate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    // Set Keplerian elements for Vehicle.
    Eigen::Vector6d initialStateInKeplerianElements;
    initialStateInKeplerianElements( semiMajorAxisIndex ) = 7200.0E3;
    initialStateInKeplerianElements( eccentricityIndex ) = 0.05;
    initialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    initialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = unit_conversions::convertDegreesToRadians( 235.7 );
    initialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 23.4 );
    initialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

    // Set (perturbed) initial state.
    Eigen::Matrix< double, 6, 1 > systemInitialState = convertKeplerianToCartesianElements(
                initialStateInKeplerianElements, earthGravitationalParameter );

    // Create propagator settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, finalEphemerisTime, cowell );

    // Create integrator settings.
    std::shared_ptr< IntegratorSettings< double > > integratorSettings
            = std::make_shared< RungeKuttaVariableStepSizeSettings< double > >
            ( initialEphemerisTime, 40.0, CoefficientSets::rungeKuttaFehlberg78,
              40.0, 40.0, 1.0, 1.0 );

    // Define link ends.
    std::vector< LinkDefinition > stationReceiverLinkEnds;
    std::vector< LinkDefinition > stationTransmitterLinkEnds;

    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
    {
        LinkDefinition linkEnds;
        linkEnds[ transmitter ] = std::pair< std::string, std::string >( std::make_pair( "Earth", groundStationNames.at( i ) ) );
        linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Vehicle", "" );
        stationTransmitterLinkEnds.push_back( linkEnds );

        linkEnds[ receiver ] = std::pair< std::string, std::string >( std::make_pair( "Earth", groundStationNames.at( i ) ) );
        linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Vehicle", "" );
        stationReceiverLinkEnds.push_back( linkEnds );
    }

    std::map< ObservableType, std::vector< LinkDefinition > > linkEndsPerObservable;
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 1 ] );

    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 1 ] );
    linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ 2 ] );

    linkEndsPerObservable[ angular_position ].push_back( stationReceiverLinkEnds[ 2 ] );
    linkEndsPerObservable[ angular_position ].push_back( stationTransmitterLinkEnds[ 1 ] );


    // Define parameters to be estimated.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialStateParameterSettings< double >( propagatorSettings, bodies );

    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", desaturation_delta_v_values ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", radiation_pressure_coefficient ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", constant_drag_coefficient ) );

    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  2, 0, 2, 2, "Earth", spherical_harmonics_cosine_coefficient_block ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  2, 1, 2, 2, "Earth", spherical_harmonics_sine_coefficient_block ) );

    parameterNames.push_back(  std::make_shared< EstimatableParameterSettings > ( "Earth", rotation_pole_position ) );
    parameterNames.push_back(  std::make_shared< EstimatableParameterSettings > ( "Earth", ground_station_position, "Station1" ) );


    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double >( parameterNames, bodies , propagatorSettings );

    printEstimatableParameterEntries( parametersToEstimate );

    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    for( std::map< ObservableType, std::vector< LinkDefinition > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;


        std::vector< LinkDefinition > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            observationSettingsList.push_back(
                        std::make_shared< ObservationModelSettings >(
                            currentObservable, currentLinkEndsList.at( i ), std::shared_ptr< LightTimeCorrectionSettings >( ) ) );
        }
    }

    // Create orbit determination object.
    OrbitDeterminationManager< double, double > orbitDeterminationManager = OrbitDeterminationManager< double, double >(
                bodies, parametersToEstimate, observationSettingsList, integratorSettings, propagatorSettings );

    // Compute list of observation times.
    std::vector< double > baseTimeList;
    double observationTimeStart = initialEphemerisTime + 1000.0;
    double  observationInterval = 60.0;
    for( int i = 0; i < numberOfDaysOfData; i++ )
    {
        for( unsigned int j = 0; j < 500; j++ )
        {
            baseTimeList.push_back( observationTimeStart + static_cast< double >( i ) * 86400.0 +
                                    static_cast< double >( j ) * observationInterval );
        }
    }


    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
    for( std::map< ObservableType, std::vector< LinkDefinition > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkDefinition > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            measurementSimulationInput.push_back(
                        std::make_shared< TabulatedObservationSimulationSettings< > >(
                            currentObservable, currentLinkEndsList[ i ], baseTimeList, receiver ) );
        }
    }

    // Simulate observations.
    std::shared_ptr< ObservationCollection< > > observationsAndTimes = simulateObservations< double, double >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    // Perturb parameter estimate.
    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< double >( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    Eigen::Matrix< double, Eigen::Dynamic, 1 > parameterPerturbation =
            Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( truthParameters.rows( ) );

    // Perturbe initial state estimate.
    parameterPerturbation.segment( 0, 3 ) = Eigen::Vector3d::Constant( 1.0 );
    parameterPerturbation.segment( 3, 3 ) = Eigen::Vector3d::Constant( 1.E-3 );

    // Perturb deltaVs estimate.
    for ( unsigned int i = 8 ; i < 8 + deltaVValues.size() * 3 ; i++)
    {
        parameterPerturbation[ i ] = 1.0e-3;
    }

    initialParameterEstimate += parameterPerturbation;
    parametersToEstimate->resetParameterValues( initialParameterEstimate );



    // Define estimation input
    std::shared_ptr< EstimationInput< double, double  > > estimationInput = std::make_shared< EstimationInput< double, double > >(
                observationsAndTimes );

    std::map< observation_models::ObservableType, double > weightPerObservable;
    weightPerObservable[ one_way_range ] = 1.0 / ( 1.0 * 1.0 );
    weightPerObservable[ angular_position ] = 1.0 / ( 1.0E-5 * 1.0E-5 );
    weightPerObservable[ one_way_doppler ] = 1.0 / ( 1.0E-11 * 1.0E-11 * physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT  );

    estimationInput->setConstantPerObservableWeightsMatrix( weightPerObservable );
    estimationInput->defineEstimationSettings( true, true, true, true, false );
    estimationInput->setConvergenceChecker(
                std::make_shared< EstimationConvergenceChecker >( 4 ) );

    // Perform estimation
    std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters(
                estimationInput );

    Eigen::VectorXd estimationError = estimationOutput->parameterEstimate_ - truthParameters;
    std::cout <<"estimation error: "<< ( estimationError ).transpose( ) << std::endl;


    // Check if parameters are correctly estimated
    Eigen::VectorXd estimatedParametervalues = estimationOutput->parameterEstimate_;

    // Initial state.
    for( unsigned int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( truthParameters( i ) - estimationOutput->parameterEstimate_( i ) ), 0.1 );
        BOOST_CHECK_SMALL( std::fabs( truthParameters( i + 3 ) - estimationOutput->parameterEstimate_( i + 3 ) ), 1.0E-6 );
    }
    // Radiation pressure and drag coefficients.
    BOOST_CHECK_SMALL( std::fabs( truthParameters( 6 ) - estimationOutput->parameterEstimate_( 6 ) ), 1.0e-4 );
    BOOST_CHECK_SMALL( std::fabs( truthParameters( 7 ) - estimationOutput->parameterEstimate_( 7 ) ), 1.0e-4 );

    // Momentum wheel desaturation deltaV values.
    for ( unsigned int i = 8 ; i < 8 + deltaVValues.size() * 3 ; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( truthParameters( i ) - estimationOutput->parameterEstimate_( i ) ), 1.0E-9 );
    }

    // Gravity field coefficients.
    for ( unsigned int i = 17 ; i < 22 ; i++ ){
        BOOST_CHECK_SMALL( std::fabs( truthParameters( i ) - estimationOutput->parameterEstimate_( i ) ), 1.0E-12 );
    }

    // Earth pole position.
    for ( unsigned int i = 22 ; i < 24 ; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( truthParameters( i ) - estimationOutput->parameterEstimate_( i ) ), 1.0E-12 );
    }

    // Ground station position.
    for ( unsigned int i = 24 ; i < 27 ; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( truthParameters( i ) - estimationOutput->parameterEstimate_( i ) ), 1.0E-6 );
    }

}

}


BOOST_AUTO_TEST_SUITE_END( )

}
