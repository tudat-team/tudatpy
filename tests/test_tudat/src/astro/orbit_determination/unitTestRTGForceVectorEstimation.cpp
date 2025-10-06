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
BOOST_AUTO_TEST_SUITE( test_rtg_force_vector_estimation )

// Using declarations.
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

//! Unit test to check if rtg acceleration model values are estimated correctly
BOOST_AUTO_TEST_CASE( test_RTGForceVectorEstimation )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Specify number of observation days.
    int numberOfDaysOfData = 1;
    // Specify initial time
    double initialEphemerisTime = 1E07;
    double finalEphemerisTime = initialEphemerisTime + numberOfDaysOfData * 86400.0;


    // Create bodies needed in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );

    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
            "ECLIPJ2000",
            "IAU_Earth",
            spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
            initialEphemerisTime,
            2.0 * mathematical_constants::PI / ( physical_constants::JULIAN_DAY ) );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );


    ///////////// DYNAMICAL MODEL SETUP /////////////////////////////////////////////////

    // Define Relevant Epochs
    double referenceEpoch = initialEphemerisTime;

    // Define function describing rotational ephemeris of vehicle
    std::function<Eigen::Matrix3d(double)> timeDependentRotationFunction =
    [](double epoch) {
        double angleRad = 1/70. * epoch * mathematical_constants::PI / 180.0;
        return Eigen::AngleAxisd(angleRad, Eigen::Vector3d::UnitZ()).toRotationMatrix();
    };

    // Assign and update
    bodies.at( "Vehicle" )
            ->setRotationalEphemeris( createRotationModel( std::make_shared< CustomRotationModelSettings >(
                                                                   "ECLIPJ2000",
                                                                   "VehicleFixed",
                                                                   timeDependentRotationFunction,
                                                                   1.0 ),
                                                           "Vehicle",
                                                           bodies ) );
    bodies.at( "Vehicle" )->setCurrentRotationalStateToLocalFrameFromEphemeris( referenceEpoch );

    // Define function describing mass function of vehicle
    double initialVehicleMass = 5000;

    // Define vehicle mass function
    std::function<double(double)> vehicleMassFunction =
        [=](double epoch) {
            //double delta_epoch = epoch - referenceEpoch;
            double delta_test = epoch - initialEphemerisTime - 20000;
            if (delta_test > -1000. && delta_test <= 1000.) {
                return initialVehicleMass - (delta_test+1000.);
            } else if (delta_test <= -1000.) {
                return initialVehicleMass;
            } else {
                return initialVehicleMass - 2000.;
            }
    };

    // Assign and Update
    bodies.at( "Vehicle" )->setBodyMassFunction( vehicleMassFunction );
    bodies.at( "Vehicle" )->updateMass( referenceEpoch );

    Eigen::Vector3d rtgForceVectorValues;
    rtgForceVectorValues << 0.5E-5, 0.5E-5, 0.5E-5;
    double decayScaleFactor = 1.5e-05;       // corresponding to a half-life of approximately half a day

    // Define origin of integration
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );


    /////////////    /////////////////////////////////////////////////

    // Create ground stations.
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "Station1" );
    groundStationNames.push_back( "Station2" );
    groundStationNames.push_back( "Station3" );

    createGroundStation( bodies.at( "Earth" ),
                         "Station1",
                         ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ),
                         coordinate_conversions::geodetic_position );
    createGroundStation( bodies.at( "Earth" ),
                         "Station2",
                         ( Eigen::Vector3d( ) << 0.0, -0.55, 2.0 ).finished( ),
                         coordinate_conversions::geodetic_position );
    createGroundStation( bodies.at( "Earth" ),
                         "Station3",
                         ( Eigen::Vector3d( ) << 0.0, 0.05, 4.0 ).finished( ),
                         coordinate_conversions::geodetic_position );

    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 8, 8 ) );
    accelerationsOfVehicle[ "Vehicle" ].push_back( std::make_shared< RTGAccelerationSettings >(rtgForceVectorValues, decayScaleFactor, referenceEpoch));

    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToPropagate, centralBodies );

    std::shared_ptr< basic_astrodynamics::AccelerationModel3d > accelerationModel = accelerationModelMap[ "Vehicle"] ["Vehicle"][ 0 ];

    // Dynamic cast acceleration settings to required type and check consistency.
    std::shared_ptr< system_models::RTGAccelerationModel > rtgAccelerationModel =
            std::dynamic_pointer_cast< system_models::RTGAccelerationModel >( accelerationModel );

    // Test Successful Construction
    BOOST_CHECK_EQUAL( rtgAccelerationModel != nullptr, true );

    // Set Keplerian elements for Vehicle.
    Eigen::Vector6d initialStateInKeplerianElements;
    initialStateInKeplerianElements( semiMajorAxisIndex ) = 7200.0E3;
    initialStateInKeplerianElements( eccentricityIndex ) = 0.05;
    initialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    initialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
    initialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
    initialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

    // Set (perturbed) initial state.
    Eigen::Matrix< double, 6, 1 > systemInitialState =
            convertKeplerianToCartesianElements( initialStateInKeplerianElements, earthGravitationalParameter );

    // Create propagator settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, finalEphemerisTime, cowell );

    // Create integrator settings.
    std::shared_ptr< IntegratorSettings< double > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< double > >(
            initialEphemerisTime, 40.0, CoefficientSets::rungeKuttaFehlberg78, 40.0, 40.0, 1.0, 1.0 );

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


    for (int i=0; i<=1; i++)


    {
        // Define parameters to be estimated.
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
                getInitialStateParameterSettings< double >( propagatorSettings, bodies );
        int parameter_size;
        if (i==0){
            parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", rtg_force_vector ) );
            parameter_size=3;
        }
        else if (i==1){
            parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", rtg_force_vector_magnitude ) );
            parameter_size=1;
        }

        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Earth", rotation_pole_position ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Earth", ground_station_position, "Station1" ) );

        // Create parameters
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate< double >( parameterNames, bodies, propagatorSettings );

        printEstimatableParameterEntries( parametersToEstimate );

        std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
        for( std::map< ObservableType, std::vector< LinkDefinition > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
             linkEndIterator != linkEndsPerObservable.end( );
             linkEndIterator++ )
        {
            ObservableType currentObservable = linkEndIterator->first;

            std::vector< LinkDefinition > currentLinkEndsList = linkEndIterator->second;
            for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
            {
                observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                        currentObservable, currentLinkEndsList.at( i ), std::shared_ptr< LightTimeCorrectionSettings >( ) ) );
            }
        }

        // Create orbit determination object.
        OrbitDeterminationManager< double, double > orbitDeterminationManager = OrbitDeterminationManager< double, double >(
                bodies, parametersToEstimate, observationSettingsList, integratorSettings, propagatorSettings );

        // Compute list of observation times.
        std::vector< double > baseTimeList;
        double observationTimeStart = initialEphemerisTime + 1000.0;
        double observationInterval = 60.0;
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
             linkEndIterator != linkEndsPerObservable.end( );
             linkEndIterator++ )
        {
            ObservableType currentObservable = linkEndIterator->first;
            std::vector< LinkDefinition > currentLinkEndsList = linkEndIterator->second;
            for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
            {
                measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings<> >(
                        currentObservable, currentLinkEndsList[ i ], baseTimeList, receiver ) );
            }
        }

        // Simulate observations.
        std::shared_ptr< ObservationCollection<> > observationsAndTimes = simulateObservations< double, double >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

        // Set weights
        std::map< std::shared_ptr< observation_models::ObservationCollectionParser >, double > weightsPerObservationParser;
        weightsPerObservationParser[ observationParser( one_way_range ) ] = 1.0 / ( 1.0 * 1.0 );
        weightsPerObservationParser[ observationParser( angular_position ) ] = 1.0 / ( 1.0E-5 * 1.0E-5 );
        weightsPerObservationParser[ observationParser( one_way_doppler ) ] =
                1.0 / ( 1.0E-11 * 1.0E-11 * physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT );
        observationsAndTimes->setConstantWeightPerObservable( weightsPerObservationParser );

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
        for( unsigned int i = 6; i < 6 + parameter_size; i++ )
        {
            parameterPerturbation[ i ] = 1.0e-6;
        }

        initialParameterEstimate += parameterPerturbation;
        parametersToEstimate->resetParameterValues( initialParameterEstimate );

        // Define estimation input
        std::shared_ptr< EstimationInput< double, double > > estimationInput =
                std::make_shared< EstimationInput< double, double > >( observationsAndTimes );

        estimationInput->defineEstimationSettings( true, true, true, true, false );
        estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 4 ) );

        // Perform estimation
        std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );

        Eigen::VectorXd estimationError = estimationOutput->parameterEstimate_ - truthParameters;
        std::cout << "estimation error: " << ( estimationError ).transpose( ) << std::endl;

        // Check if parameters are correctly estimated
        Eigen::VectorXd estimatedParametervalues = estimationOutput->parameterEstimate_;

        // Initial state.
        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( truthParameters( i ) - estimationOutput->parameterEstimate_( i ) ), 0.1 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( i + 3 ) - estimationOutput->parameterEstimate_( i + 3 ) ), 1.0E-6 );
        }
        // Radiation pressure and drag coefficients.
        //BOOST_CHECK_SMALL( std::fabs( truthParameters( 6 ) - estimationOutput->parameterEstimate_( 6 ) ), 1.0e-4 );
        //BOOST_CHECK_SMALL( std::fabs( truthParameters( 7 ) - estimationOutput->parameterEstimate_( 7 ) ), 1.0e-4 );

        // rtg parameter values.
        for( unsigned int i = 6; i < 6 + parameter_size; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( truthParameters( i ) - estimationOutput->parameterEstimate_( i ) ), 1.0E-9 );
        }

        // Earth pole position.
        for( unsigned int i = 6 + parameter_size; i < 8+parameter_size; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( truthParameters( i ) - estimationOutput->parameterEstimate_( i ) ), 1.0E-12 );
        }

        // Ground station position.
        for( unsigned int i = 8+parameter_size; i < 11+parameter_size; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( truthParameters( i ) - estimationOutput->parameterEstimate_( i ) ), 1.0E-6 );
        }
    }
}

}  // namespace unit_tests

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace tudat
