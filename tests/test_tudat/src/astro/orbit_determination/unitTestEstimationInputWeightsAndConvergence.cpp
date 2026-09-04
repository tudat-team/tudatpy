/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <limits>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/createInverseAprioriCovariance.h"
#include "tudat/simulation/estimation_setup/executePlanetaryParameterEstimationTestCase.h"
#include "tudat/simulation/estimation_setup/executeEarthOrbiterParameterEstimationTestCase.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/gravitationalParameter.h"
#include "tudat/simulation/estimation_setup/podProcessing.h"

namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_estimation_input_output )

BOOST_AUTO_TEST_CASE( test_WeightDefinitions )

{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    const double startTime = double( 1.0E7 );
    const int numberOfDaysOfData = 3;

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Mars" );

    // Specify initial time
    double initialEphemerisTime = startTime;
    double finalEphemerisTime = initialEphemerisTime + numberOfDaysOfData * 86400.0;

    // Create bodies needed in simulation
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
            "ECLIPJ2000",
            "IAU_Earth",
            spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
            initialEphemerisTime,
            2.0 * mathematical_constants::PI / ( physical_constants::JULIAN_DAY ) );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( 400.0 );

    bodies.at( "Vehicle" )
            ->setEphemeris( std::make_shared< TabulatedCartesianEphemeris<> >(
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d > >( ), "Earth", "ECLIPJ2000" ) );

    // Create ground stations: same position, but different representation
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "Station1" );
    groundStationNames.push_back( "Station2" );
    groundStationNames.push_back( "Station3" );

    createGroundStation( bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station2", ( Eigen::Vector3d( ) << 0.0, -0.55, 2.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station3", ( Eigen::Vector3d( ) << 0.0, 0.05, 4.0 ).finished( ), geodetic_position );

    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 8, 8 ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;
    bodiesToIntegrate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7200.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

    // Set (perturbed) initial state.
    Eigen::Matrix< double, 6, 1 > systemInitialState =
            convertKeplerianToCartesianElements( asterixInitialStateInKeplerianElements, earthGravitationalParameter );

    // Create propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< double, double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double, double > >(
                    centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, double( finalEphemerisTime ), cowell );

    // Create integrator settings
    std::shared_ptr< IntegratorSettings< double > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< double > >(
            double( initialEphemerisTime ), 40.0, CoefficientSets::rungeKuttaFehlberg78, 40.0, 40.0, 1.0, 1.0 );

    // Define parameters.
    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;

    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = LinkEndId( "Earth", groundStationNames.at( i ) );
        linkEnds[ receiver ] = LinkEndId( "Vehicle", "" );
        stationTransmitterLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = LinkEndId( "Earth", groundStationNames.at( i ) );
        linkEnds[ transmitter ] = LinkEndId( "Vehicle", "" );
        stationReceiverLinkEnds.push_back( linkEnds );
    }

    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 1 ] );

    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 1 ] );
    linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ 2 ] );

    linkEndsPerObservable[ angular_position ].push_back( stationReceiverLinkEnds[ 2 ] );
    linkEndsPerObservable[ angular_position ].push_back( stationTransmitterLinkEnds[ 1 ] );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back(
            std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >( "Vehicle", systemInitialState, "Earth" ) );

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, bodies );

    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;

    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( );
         linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;

        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            observationSettingsList.push_back(
                    std::make_shared< ObservationModelSettings >( currentObservable, currentLinkEndsList.at( i ) ) );
        }
    }

    // Create orbit determination object.
    OrbitDeterminationManager< double, double > orbitDeterminationManager = OrbitDeterminationManager< double, double >(
            bodies, parametersToEstimate, observationSettingsList, integratorSettings, propagatorSettings );

    std::vector< double > baseTimeList;
    double observationTimeStart = initialEphemerisTime + 1000.0;
    double observationInterval = 20.0;
    unsigned int nbObsPerDay = 50;
    for( int i = 0; i < numberOfDaysOfData; i++ )
    {
        for( unsigned int j = 0; j < nbObsPerDay; j++ )
        {
            baseTimeList.push_back( observationTimeStart + static_cast< double >( i ) * 86400.0 +
                                    static_cast< double >( j ) * observationInterval );
        }
    }

    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput =
            getObservationSimulationSettings< double >( linkEndsPerObservable, baseTimeList, receiver );

    // Simulate observations
    std::shared_ptr< ObservationDataset< double, double > > simulatedObservations = simulateObservationDataset< double, double >(
            measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    // Define estimation input
    std::shared_ptr< EstimationInput< double, double > > estimationInput =
            std::make_shared< EstimationInput< double, double > >( simulatedObservations );

    std::map< ObservableType, std::pair< int, int > > observationTypeStartAndSize = simulatedObservations->getObservableTypeStartAndSize( );

    {
        simulatedObservations->setConstantSingleObservationScalarWeight( ObservationSelectionCondition< double, double >::all( ), 0.1 );

        // Define estimation input
        std::shared_ptr< EstimationInput< double, double > > estimationInput =
                std::make_shared< EstimationInput< double, double > >( simulatedObservations );
        std::map< ObservableType, std::pair< int, int > > observationTypeStartAndSize =
                simulatedObservations->getObservableTypeStartAndSize( );

        Eigen::VectorXd totalWeights = estimationInput->getWeightsMatrixDiagonals( );

        for( unsigned int i = 0; i < totalWeights.rows( ); i++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( totalWeights( i ), 0.1, std::numeric_limits< double >::epsilon( ) );
        }
    }

    {
        simulatedObservations->setConstantSingleObservationScalarWeight(
                ObservationSelectionCondition< double, double >::observableType( one_way_range ), 1.0 / ( 3.0 * 3.0 ) );
        simulatedObservations->setConstantSingleObservationScalarWeight(
                ObservationSelectionCondition< double, double >::observableType( angular_position ), 1.0 / ( 1.0E-5 * 1.0E-5 ) );
        simulatedObservations->setConstantSingleObservationScalarWeight(
                ObservationSelectionCondition< double, double >::observableType( one_way_doppler ),
                1.0 / ( 1.0E-11 * 1.0E-11 * SPEED_OF_LIGHT * SPEED_OF_LIGHT ) );
        const std::map< ObservableType, double > expectedWeightPerObservable = {
            { one_way_range, 1.0 / ( 3.0 * 3.0 ) },
            { angular_position, 1.0 / ( 1.0E-5 * 1.0E-5 ) },
            { one_way_doppler, 1.0 / ( 1.0E-11 * 1.0E-11 * SPEED_OF_LIGHT * SPEED_OF_LIGHT ) }
        };

        // Define estimation input
        std::shared_ptr< EstimationInput< double, double > > estimationInput =
                std::make_shared< EstimationInput< double, double > >( simulatedObservations );
        std::map< ObservableType, std::pair< int, int > > observationTypeStartAndSize =
                simulatedObservations->getObservableTypeStartAndSize( );

        Eigen::VectorXd totalWeights = estimationInput->getWeightsMatrixDiagonals( );

        for( const auto& it : expectedWeightPerObservable )
        {
            ObservableType observableType = it.first;
            for( int i = 0; i < observationTypeStartAndSize.at( observableType ).second; i++ )
            {
                BOOST_CHECK_CLOSE_FRACTION( totalWeights( observationTypeStartAndSize.at( observableType ).first + i ),
                                            it.second,
                                            std::numeric_limits< double >::epsilon( ) );
            }
        }
    }

    {
        Eigen::Vector2d angularPositionWeight;
        angularPositionWeight << 0.1, 0.2;
        simulatedObservations->setConstantSingleObservationScalarWeight( ObservationSelectionCondition< double, double >::all( ), 2.0 );
        simulatedObservations->setConstantSingleObservationDiagonalWeight(
                ObservationSelectionCondition< double, double >::observableType( angular_position ), angularPositionWeight );

        // Define estimation input
        std::shared_ptr< EstimationInput< double, double > > estimationInput =
                std::make_shared< EstimationInput< double, double > >( simulatedObservations );
        std::map< ObservableType, std::pair< int, int > > observationTypeStartAndSize =
                simulatedObservations->getObservableTypeStartAndSize( );

        Eigen::VectorXd totalWeights = estimationInput->getWeightsMatrixDiagonals( );

        std::pair< int, int > startEndIndex = observationTypeStartAndSize.at( angular_position );

        for( int i = 0; i < startEndIndex.first; i++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( totalWeights( i ), 2.0, std::numeric_limits< double >::epsilon( ) );
        }

        for( int i = 0; i < startEndIndex.second; i++ )
        {
            if( i % 2 == 0 )
            {
                BOOST_CHECK_CLOSE_FRACTION(
                        totalWeights( startEndIndex.first + i ), angularPositionWeight( 0 ), std::numeric_limits< double >::epsilon( ) );
            }
            else
            {
                BOOST_CHECK_CLOSE_FRACTION(
                        totalWeights( startEndIndex.first + i ), angularPositionWeight( 1 ), std::numeric_limits< double >::epsilon( ) );
            }
        }

        for( unsigned int i = startEndIndex.first + startEndIndex.second; i < totalWeights.rows( ); i++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( totalWeights( i ), 2.0, std::numeric_limits< double >::epsilon( ) );
        }
    }

    // Test tabulated weights
    {
        // Set same tabulated weights to each range observation set
        int sizeRangeObsPerObsSet = nbObsPerDay * numberOfDaysOfData;
        Eigen::VectorXd singleSetRangeWeights =
                Eigen::VectorXd::LinSpaced( sizeRangeObsPerObsSet, 1.0 / ( 3.0 * 3.0 ), 1.0 / ( 4.0 * 4.0 ) );

        // Compute full range weight vector
        unsigned int nbRangeObsSets = simulatedObservations->getObservationSetIdsForObservableType( one_way_range ).size( );
        Eigen::VectorXd rangeWeights = Eigen::VectorXd::Zero( nbRangeObsSets * sizeRangeObsPerObsSet );
        for( unsigned int k = 0; k < nbRangeObsSets; k++ )
        {
            rangeWeights.segment( k * sizeRangeObsPerObsSet, sizeRangeObsPerObsSet ) = singleSetRangeWeights;
        }

        // Set total tabulated weights for all Doppler observation sets
        int totalSizeDopplerObs = static_cast< int >( simulatedObservations->getTotalScalarSizeForObservableType( one_way_doppler ) );
        Eigen::VectorXd dopplerWeights = Eigen::VectorXd::LinSpaced( totalSizeDopplerObs,
                                                                     1.0 / ( 1.0e-11 * SPEED_OF_LIGHT * 1.0e-11 * SPEED_OF_LIGHT ),
                                                                     1.0 / ( 1.5e-11 * SPEED_OF_LIGHT * 1.5e-11 * SPEED_OF_LIGHT ) );

        // Default angular position weights set to 1
        int totalSizeAngularPositionObs =
                static_cast< int >( simulatedObservations->getTotalScalarSizeForObservableType( angular_position ) );
        Eigen::VectorXd angularPositionWeights = Eigen::VectorXd::Ones( totalSizeAngularPositionObs );

        // Concatenate tabulated weights per observable type (default weights for angular_position observables)
        simulatedObservations->setWeightVectorForObservableType( one_way_range, singleSetRangeWeights );
        simulatedObservations->setWeightVectorForObservableType( one_way_doppler, dopplerWeights );
        simulatedObservations->setWeightVectorForObservableType( angular_position, angularPositionWeights );

        // Define estimation input
        std::shared_ptr< EstimationInput< double, double > > estimationInput =
                std::make_shared< EstimationInput< double, double > >( simulatedObservations );
        std::map< ObservableType, std::pair< int, int > > observationTypeStartAndSize =
                simulatedObservations->getObservableTypeStartAndSize( );
        Eigen::VectorXd totalWeights = estimationInput->getWeightsMatrixDiagonals( );

        // Define expected weights per observable
        std::map< ObservableType, Eigen::VectorXd > expectedWeights;
        expectedWeights[ one_way_range ] = rangeWeights;
        expectedWeights[ one_way_doppler ] = dopplerWeights;
        expectedWeights[ angular_position ] = angularPositionWeights;

        for( auto it : expectedWeights )
        {
            for( int i = 0; i < observationTypeStartAndSize.at( it.first ).second; i++ )
            {
                BOOST_CHECK_CLOSE_FRACTION( totalWeights( observationTypeStartAndSize.at( it.first ).first + i ),
                                            it.second( i ),
                                            std::numeric_limits< double >::epsilon( ) );
            }
        }
    }
}

//! Test that best-iteration selection during estimation uses the least-squares cost function
BOOST_AUTO_TEST_CASE( test_CostFunctionBasedBestIterationSelection )
{
    using TimeType = double;
    using StateScalarType = double;

    spice_interface::loadStandardSpiceKernels( );

    const int numberOfDaysOfData = 365;
    const TimeType initialEphemerisTime = TimeType( 1.0E7 );
    const TimeType finalEphemerisTime = initialEphemerisTime + static_cast< TimeType >( numberOfDaysOfData ) * 86400.0;
    const double maximumTimeStep = 3600.0;
    const double buffer = 10.0 * maximumTimeStep;

    std::vector< std::string > bodyNames = { "Earth", "Mars", "Sun", "Moon", "Jupiter", "Saturn" };
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings.at( "Moon" )->ephemerisSettings->resetFrameOrigin( "Sun" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfMars;
    accelerationsOfMars[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfMars[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfMars[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfMars[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfMars[ "Saturn" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Mars" ] = accelerationsOfMars;

    std::vector< std::string > bodiesToIntegrate = { "Mars" };
    std::map< std::string, std::string > centralBodyMap;
    centralBodyMap[ "Mars" ] = "SSB";
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, centralBodyMap );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
            "Mars",
            propagators::getInitialStateOfBody< TimeType, StateScalarType >(
                    "Mars", centralBodyMap.at( "Mars" ), bodies, initialEphemerisTime ),
            centralBodyMap.at( "Mars" ) ) );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate< StateScalarType, TimeType >( parameterNames, bodies );

    std::vector< std::string > centralBodies = { "SSB" };
    std::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
            std::make_shared< IntegratorSettings< TimeType > >( rungeKutta4, initialEphemerisTime, 1800.0 );
    std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType, TimeType > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< StateScalarType, TimeType > >(
                    centralBodies,
                    accelerationModelMap,
                    bodiesToIntegrate,
                    getInitialStateVectorOfBodiesToEstimate( parametersToEstimate ),
                    finalEphemerisTime,
                    cowell );

    LinkEnds linkEnds;
    linkEnds[ receiver ] = LinkEndId( "Earth", "" );
    linkEnds[ transmitter ] = LinkEndId( "Mars", "" );
    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( one_way_range, linkEnds ) );
    observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( angular_position, linkEnds ) );

    OrbitDeterminationManager< StateScalarType, TimeType > orbitDeterminationManager(
            bodies, parametersToEstimate, observationSettingsList, integratorSettings, propagatorSettings );

    std::vector< TimeType > observationTimes;
    observationTimes.reserve( numberOfDaysOfData );
    for( int i = 0; i < numberOfDaysOfData - 1; i++ )
    {
        observationTimes.push_back( initialEphemerisTime + ( static_cast< TimeType >( i ) + 0.5 ) * 86400.0 );
    }

    const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > nominalInitialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< StateScalarType >( );

    simulation_setup::noiseSeed = 0;
    std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > > measurementSimulationInput;
    measurementSimulationInput.push_back(
            std::make_shared< TabulatedObservationSimulationSettings< TimeType > >( one_way_range, linkEnds, observationTimes, receiver ) );
    measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< TimeType > >(
            angular_position, linkEnds, observationTimes, receiver ) );
    addGaussianNoiseFunctionToObservationSimulationSettings( measurementSimulationInput, 0.1, one_way_range );
    addGaussianNoiseFunctionToObservationSimulationSettings( measurementSimulationInput, 3.0E-7, angular_position );

    std::shared_ptr< ObservationDataset< StateScalarType, TimeType > > simulatedObservations =
            simulateObservationDataset< StateScalarType, TimeType >(
                    measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    // Inject deterministic structured biases so range and angular residual improvements compete across iterations.
    Eigen::VectorXd baseRangeObservations = simulatedObservations->getObservationVectorForObservableType( one_way_range );
    Eigen::VectorXd baseAngularObservations = simulatedObservations->getObservationVectorForObservableType( angular_position );
    for( int i = 0; i < baseRangeObservations.size( ); i++ )
    {
        const double cycleArgument = static_cast< double >( i ) / 31.0;
        baseRangeObservations( i ) += 0.25 * std::sin( 2.0 * mathematical_constants::PI * cycleArgument );
    }
    for( int i = 0; i < baseAngularObservations.size( ) / 2; i++ )
    {
        const double cycleArgument = static_cast< double >( i ) / 43.0;
        baseAngularObservations( 2 * i ) += 4.0E-8 * std::cos( 2.0 * mathematical_constants::PI * cycleArgument );
        baseAngularObservations( 2 * i + 1 ) -= 4.0E-8 * std::sin( 2.0 * mathematical_constants::PI * cycleArgument );
    }

    const double rangeBaseWeight = 1.0 / ( 0.1 * 0.1 );
    const double angularBaseWeight = 1.0 / ( 1.0E-8 * 1.0E-8 );
    Eigen::VectorXd rangeWeights = Eigen::VectorXd::Zero( baseRangeObservations.size( ) );
    for( int i = 0; i < rangeWeights.size( ); i++ )
    {
        const double cycleArgument = static_cast< double >( i ) / 29.0;
        const double scaleFactor = 0.25 + 1.75 * ( 0.5 + 0.5 * std::sin( 2.0 * mathematical_constants::PI * cycleArgument ) );
        rangeWeights( i ) = rangeBaseWeight * scaleFactor;
    }
    Eigen::VectorXd angularWeights = Eigen::VectorXd::Zero( baseAngularObservations.size( ) );
    for( int i = 0; i < angularWeights.size( ) / 2; i++ )
    {
        const double cycleArgument = static_cast< double >( i ) / 37.0;
        const double raScaleFactor = 0.2 + 3.8 * ( 0.5 + 0.5 * std::sin( 2.0 * mathematical_constants::PI * cycleArgument ) );
        const double decScaleFactor = 0.2 + 3.8 * ( 0.5 + 0.5 * std::cos( 2.0 * mathematical_constants::PI * cycleArgument ) );
        angularWeights( 2 * i ) = angularBaseWeight * raScaleFactor;
        angularWeights( 2 * i + 1 ) = angularBaseWeight * decScaleFactor;
    }
    simulatedObservations->setWeightVectorForObservableType( one_way_range, rangeWeights );
    simulatedObservations->setWeightVectorForObservableType( angular_position, angularWeights );
    simulatedObservations->setObservationVectorForObservableType( one_way_range, baseRangeObservations );
    simulatedObservations->setObservationVectorForObservableType( angular_position, baseAngularObservations );

    int numberOfDistinctBestIterationCases = 0;

    std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > perturbationCases;
    perturbationCases.push_back( ( Eigen::VectorXd( 6 ) << 5.0E3, -3.0E3, 2.0E3, 0.5, -0.4, 0.3 ).finished( ) );
    perturbationCases.push_back( ( Eigen::VectorXd( 6 ) << 8.0E3, -6.0E3, 4.0E3, 1.2, -0.9, 0.7 ).finished( ) );
    perturbationCases.push_back( ( Eigen::VectorXd( 6 ) << 2.0E4, 1.0E4, -1.5E4, 3.0, -2.0, 1.5 ).finished( ) );
    perturbationCases.push_back( ( Eigen::VectorXd( 6 ) << -1.5E4, 8.0E3, 1.0E4, -2.5, 1.8, -1.2 ).finished( ) );
    perturbationCases.push_back( ( Eigen::VectorXd( 6 ) << 3.0E4, -2.5E4, 1.0E4, 5.0, -4.0, 3.5 ).finished( ) );
    perturbationCases.push_back( ( Eigen::VectorXd( 6 ) << -4.0E4, 2.0E4, -3.0E4, -6.0, 4.5, -3.2 ).finished( ) );
    perturbationCases.push_back( ( Eigen::VectorXd( 6 ) << 5.0E4, 4.0E4, -2.0E4, 7.5, 6.0, -4.5 ).finished( ) );

    for( unsigned int caseIndex = 0; caseIndex < perturbationCases.size( ); caseIndex++ )
    {
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > perturbedInitialParameterEstimate = nominalInitialParameterEstimate;
        perturbedInitialParameterEstimate += perturbationCases.at( caseIndex );
        parametersToEstimate->template resetParameterValues< StateScalarType >( perturbedInitialParameterEstimate );

        std::shared_ptr< EstimationInput< StateScalarType, TimeType > > estimationInput =
                std::make_shared< EstimationInput< StateScalarType, TimeType > >( simulatedObservations );
        estimationInput->defineEstimationSettings( true, true, false, false, true, false );
        estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 10, -1.0, -1.0, 1000 ) );

        std::shared_ptr< EstimationOutput< StateScalarType, TimeType > > estimationOutput =
                orbitDeterminationManager.estimateParameters( estimationInput );

        BOOST_CHECK( estimationOutput->residualHistory_.size( ) >= 2 );
        BOOST_CHECK( estimationOutput->bestIteration_ >= 0 );
        BOOST_CHECK( estimationOutput->bestIteration_ < static_cast< int >( estimationOutput->residualHistory_.size( ) ) );

        const Eigen::VectorXd weights = estimationInput->getWeightsMatrixDiagonals( );
        int minimumCostIteration = -1;
        int minimumRmsIteration = -1;
        double minimumCost = TUDAT_NAN;
        double minimumRms = TUDAT_NAN;
        double rmsFromMinimumCostIteration = TUDAT_NAN;
        std::vector< double > costHistory;
        std::vector< double > rmsHistory;
        costHistory.reserve( estimationOutput->residualHistory_.size( ) );
        rmsHistory.reserve( estimationOutput->residualHistory_.size( ) );
        for( unsigned int i = 0; i < estimationOutput->residualHistory_.size( ); i++ )
        {
            const double currentCost =
                    linear_algebra::computeLeastSquaresCostFunction( weights, estimationOutput->residualHistory_.at( i ) );
            const double currentRms = linear_algebra::getVectorEntryRootMeanSquare( estimationOutput->residualHistory_.at( i ) );
            costHistory.push_back( currentCost );
            rmsHistory.push_back( currentRms );

            if( currentCost < minimumCost || !( minimumCost == minimumCost ) )
            {
                minimumCost = currentCost;
                minimumCostIteration = static_cast< int >( i );
                rmsFromMinimumCostIteration = currentRms;
            }
            if( currentRms < minimumRms || !( minimumRms == minimumRms ) )
            {
                minimumRms = currentRms;
                minimumRmsIteration = static_cast< int >( i );
            }
        }

        BOOST_CHECK_EQUAL( estimationOutput->bestIteration_, minimumCostIteration );
        BOOST_CHECK_CLOSE_FRACTION( estimationOutput->residualStandardDeviation_, rmsFromMinimumCostIteration, 1.0E-15 );

        std::cout << "Perturbation case " << caseIndex << " iteration history:" << std::endl;
        for( unsigned int i = 0; i < costHistory.size( ); i++ )
        {
            std::cout << "  Iteration " << i << ": RMS = " << rmsHistory.at( i ) << ", cost = " << costHistory.at( i ) << std::endl;
        }
        std::cout << "  Estimator best iteration: " << estimationOutput->bestIteration_ << std::endl;
        std::cout << "  Minimum-cost iteration: " << minimumCostIteration << std::endl;
        std::cout << "  Minimum-RMS iteration: " << minimumRmsIteration << std::endl;

        if( minimumCostIteration != minimumRmsIteration )
        {
            numberOfDistinctBestIterationCases++;

            std::cout << "Selected perturbation case for cost-vs-RMS distinction: " << caseIndex << std::endl;
            for( unsigned int i = 0; i < costHistory.size( ); i++ )
            {
                std::cout << "Iteration " << i << ": RMS = " << rmsHistory.at( i ) << ", cost = " << costHistory.at( i ) << std::endl;
            }
            std::cout << "Best iteration from estimator: " << estimationOutput->bestIteration_ << std::endl;
            std::cout << "Minimum-cost iteration: " << minimumCostIteration << std::endl;
            std::cout << "Minimum-RMS iteration: " << minimumRmsIteration << std::endl;

            BOOST_CHECK_EQUAL( estimationOutput->bestIteration_, minimumCostIteration );
            BOOST_CHECK( estimationOutput->bestIteration_ != minimumRmsIteration );
            BOOST_CHECK_CLOSE_FRACTION( estimationOutput->residualStandardDeviation_, rmsFromMinimumCostIteration, 1.0E-15 );
        }
    }

    BOOST_REQUIRE_MESSAGE( numberOfDistinctBestIterationCases > 0,
                           "Could not generate a case where minimum-cost and minimum-RMS iteration differ for mixed range/angular "
                           "estimation in the configured deterministic perturbation cases." );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
