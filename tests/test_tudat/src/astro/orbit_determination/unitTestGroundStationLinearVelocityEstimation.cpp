/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <limits>
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/estimation_setup/podProcessing.h"

using namespace tudat;
using namespace tudat::observation_models;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::numerical_integrators;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::propagators;
using namespace tudat::basic_astrodynamics;

namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_ground_station_linear_velocity_estimation )

//! Estimate the linear velocity of a ground station, together with the initial state of an Earth orbiter, from ideal range data
BOOST_AUTO_TEST_CASE( test_GroundStationLinearVelocityEstimation )
{
    spice_interface::loadStandardSpiceKernels( );

    double initialEphemerisTime = 1.0e7;
    double finalEphemerisTime = initialEphemerisTime + 86400.0;

    // Create bodies
    BodyListSettings bodySettings = getDefaultBodySettings( { "Earth", "Sun", "Moon" }, "Earth", "ECLIPJ2000" );
    bodySettings.addSettings( "Vehicle" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.at( "Vehicle" )
            ->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris<> >(
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d > >( ), "Earth", "ECLIPJ2000" ) );

    // Create ground stations: Station1 moves linearly (body-fixed velocity in m/s), Station2 has no motion model
    const Eigen::Vector3d stationVelocity = ( Eigen::Vector3d( ) << 0.02, -0.01, 0.03 ).finished( );
    createGroundStation( bodies.at( "Earth" ),
                         "Station1",
                         ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ),
                         coordinate_conversions::geodetic_position,
                         { std::make_shared< LinearGroundStationMotionSettings >( stationVelocity, initialEphemerisTime ) } );
    createGroundStation( bodies.at( "Earth" ),
                         "Station2",
                         ( Eigen::Vector3d( ) << 0.0, -0.55, 2.0 ).finished( ),
                         coordinate_conversions::geodetic_position );

    // Create propagation settings
    SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Vehicle" ][ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Vehicle" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Vehicle" ][ "Moon" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    std::vector< std::string > bodiesToIntegrate = { "Vehicle" };
    std::vector< std::string > centralBodies = { "Earth" };
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    Eigen::Vector6d initialStateInKeplerianElements;
    initialStateInKeplerianElements( semiMajorAxisIndex ) = 7200.0E3;
    initialStateInKeplerianElements( eccentricityIndex ) = 0.05;
    initialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    initialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
    initialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
    initialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );
    Eigen::Matrix< double, 6, 1 > systemInitialState = convertKeplerianToCartesianElements(
            initialStateInKeplerianElements, bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( ) );

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodies,
                    accelerationModelMap,
                    bodiesToIntegrate,
                    systemInitialState,
                    initialEphemerisTime,
                    rungeKuttaFixedStepSettings( 40.0, CoefficientSets::rungeKuttaFehlberg78 ),
                    std::make_shared< PropagationTimeTerminationSettings >( finalEphemerisTime ) );

    // Define parameters: initial state of the vehicle and velocity of Station1
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialStateParameterSettings< double >( propagatorSettings, bodies );
    parameterNames.push_back( groundStationLinearVelocity( "Earth", "Station1" ) );
    std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double >( parameterNames, bodies, propagatorSettings );
    BOOST_CHECK_EQUAL( parametersToEstimate->getParameterSetSize( ), 9 );

    // The velocity cannot be estimated for a station without a linear motion model
    BOOST_CHECK_THROW(
            createParametersToEstimate< double >( { groundStationLinearVelocity( "Earth", "Station2" ) }, bodies, propagatorSettings ),
            std::runtime_error );

    // Define one-way range observations from and to Station1, and to Station2
    std::vector< LinkDefinition > linkEndsList;
    for( const std::string stationName : { "Station1", "Station2" } )
    {
        LinkDefinition linkEnds;
        linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Vehicle", "" );
        linkEnds[ receiver ] = std::make_pair( std::string( "Earth" ), stationName );
        linkEndsList.push_back( linkEnds );
    }
    LinkDefinition uplinkEnds;
    uplinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Earth", "Station1" );
    uplinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Vehicle", "" );
    linkEndsList.push_back( uplinkEnds );

    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    std::vector< double > observationTimes;
    for( unsigned int j = 0; j < 500; j++ )
    {
        observationTimes.push_back( initialEphemerisTime + 1000.0 + static_cast< double >( j ) * 60.0 );
    }
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
    for( const LinkDefinition& linkEnds : linkEndsList )
    {
        observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( one_way_range, linkEnds ) );
        measurementSimulationInput.push_back(
                std::make_shared< TabulatedObservationSimulationSettings<> >( one_way_range, linkEnds, observationTimes, receiver ) );
    }

    // Simulate observations
    OrbitDeterminationManager< double, double > orbitDeterminationManager =
            OrbitDeterminationManager< double, double >( bodies, parametersToEstimate, observationSettingsList, propagatorSettings );
    std::shared_ptr< ObservationCollection<> > observationsAndTimes = simulateObservations< double, double >(
            measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    // Perturb the initial state and the station velocity
    Eigen::VectorXd truthParameters = parametersToEstimate->template getFullParameterValues< double >( );
    BOOST_CHECK_SMALL( ( truthParameters.segment( 6, 3 ) - stationVelocity ).norm( ), 1.0E-15 );
    Eigen::VectorXd initialParameterEstimate = truthParameters;
    initialParameterEstimate.segment( 0, 3 ) += Eigen::Vector3d::Constant( 1.0 );
    initialParameterEstimate.segment( 3, 3 ) += Eigen::Vector3d::Constant( 1.0E-3 );
    initialParameterEstimate.segment( 6, 3 ) += Eigen::Vector3d::Constant( 1.0E-3 );
    parametersToEstimate->resetParameterValues( initialParameterEstimate );

    // Estimate the parameters
    std::shared_ptr< EstimationInput< double, double > > estimationInput =
            std::make_shared< EstimationInput< double, double > >( observationsAndTimes );
    estimationInput->defineEstimationSettings( true, true, true, true, false );
    estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 4 ) );
    std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );

    Eigen::VectorXd estimationError = estimationOutput->parameterEstimate_ - truthParameters;
    std::cout << "estimation error: " << estimationError.transpose( ) << std::endl;
    for( unsigned int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( estimationError( i ) ), 1.0E-6 );
        BOOST_CHECK_SMALL( std::fabs( estimationError( i + 3 ) ), 1.0E-9 );
        BOOST_CHECK_SMALL( std::fabs( estimationError( i + 6 ) ), 1.0E-9 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
