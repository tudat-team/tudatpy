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


#include <limits>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/simulation/estimation_setup/orbitDeterminationTestCases.h"
#include "tudat/simulation/estimation_setup/podProcessing.h"


namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_observations_processing )

std::shared_ptr< ObservationCollection< double, double > > setUpObservationCollectionToTest(
        const double startTime, const int numberOfDaysOfData, const int numberOfObservations, const double observationsInterval,
        std::vector< LinkEnds >& stationReceiverLinkEnds,
        std::vector< LinkEnds >& stationTransmitterLinkEnds,
        std::map< ObservableType, double >& obsStartTimes,
        std::map< ObservableType, std::vector< double > >& baseTimeList,
        SystemOfBodies& bodies )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

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
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
            "ECLIPJ2000", "IAU_Earth", spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
            initialEphemerisTime, 2.0 * mathematical_constants::PI / ( physical_constants::JULIAN_DAY ) );

    bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( 400.0 );


    bodies.at( "Vehicle" )->setEphemeris( std::make_shared< TabulatedCartesianEphemeris< > >(
            std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d > >( ), "Earth", "ECLIPJ2000" ) );


    // Creatre ground stations: same position, but different representation
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
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
            bodies, accelerationMap, bodiesToIntegrate, centralBodies );

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
    Eigen::Matrix< double, 6, 1 > systemInitialState = convertKeplerianToCartesianElements(
            asterixInitialStateInKeplerianElements, earthGravitationalParameter );

    // Create propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< double, double > > propagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< double, double > >
            ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, double( finalEphemerisTime ), cowell );

    // Create integrator settings
    std::shared_ptr< IntegratorSettings< double > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< double > >(
            double( initialEphemerisTime ), 40.0, CoefficientSets::rungeKuttaFehlberg78, 40.0, 40.0, 1.0, 1.0 );

    // Define link ends.
    stationReceiverLinkEnds.clear( );
    stationTransmitterLinkEnds.clear( );

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
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >( "Vehicle", systemInitialState, "Earth" ) );

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, bodies );

    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;

    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;


        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( currentObservable, currentLinkEndsList.at( i ) ) );
        }
    }

    // Create orbit determination object.
    OrbitDeterminationManager< double, double > orbitDeterminationManager = OrbitDeterminationManager< double, double >(
            bodies, parametersToEstimate, observationSettingsList, integratorSettings, propagatorSettings );


    baseTimeList.clear( );
    std::vector< double > rangeObsTimes, dopplerObsTimes, angularPositionObsTimes;
    obsStartTimes = { { one_way_range, initialEphemerisTime + 1000.0 }, { one_way_doppler, initialEphemerisTime + 1000.0 + 86400.0 },
                      { angular_position, initialEphemerisTime + 1000.0 + 2.0 * 86400.0 } };
    for( int j = 0; j < numberOfObservations; j++ )
    {
        rangeObsTimes.push_back( initialEphemerisTime + 1000.0 + static_cast< double >( j ) * observationsInterval );
        dopplerObsTimes.push_back( initialEphemerisTime + 1000.0 + 86400.0 + static_cast< double >( j ) * observationsInterval );

        double angularPositionObsBuffer = ( ( j > int(numberOfObservations/2) ) ? 3000.0 : 1000.0 );
        angularPositionObsTimes.push_back( initialEphemerisTime + angularPositionObsBuffer + 2.0 * 86400.0 + static_cast< double >( j ) * observationsInterval );
    }
    baseTimeList[ one_way_range ] = rangeObsTimes;
    baseTimeList[ one_way_doppler ] = dopplerObsTimes;
    baseTimeList[ angular_position ] = angularPositionObsTimes;


    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
    for ( auto linkEndIterator : linkEndsPerObservable )
    {
        ObservableType currentObservable = linkEndIterator.first;

        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator.second;
        for ( unsigned int i = 0; i < currentLinkEndsList.size( ) ; i++ )
        {
            measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< double > >(
                    currentObservable, currentLinkEndsList.at( i ), baseTimeList.at( currentObservable ), receiver ) );
        }
    }

    // Add elevation angle dependent variables
    std::shared_ptr< ObservationDependentVariableSettings > elevationAngleSettings = elevationAngleDependentVariable( );
    addDependentVariablesToObservationSimulationSettings( measurementSimulationInput,
                                                          std::vector< std::shared_ptr< ObservationDependentVariableSettings > >( { elevationAngleSettings } ), bodies );

    // Simulate observations
    std::shared_ptr< ObservationCollection< double, double > > simulatedObservations =
            simulateObservations< double, double >( measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    // Define estimation input
    std::shared_ptr< EstimationInput< double, double  > > estimationInput = std::make_shared< EstimationInput< double, double > >( simulatedObservations );

    // Retrieve observation sets start index and size
    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > > obsSetsStartAndSize = simulatedObservations->getObservationSetStartAndSize( );

    std::map< std::shared_ptr< observation_models::ObservationCollectionParser >, double > weightPerObservationParser;
    weightPerObservationParser[ observationParser( one_way_range ) ] = 1.0 / ( 3.0 * 3.0 );
    weightPerObservationParser[ observationParser( angular_position ) ] = 1.0 / ( 1.0E-5 * 1.0E-5 );
    weightPerObservationParser[ observationParser( one_way_doppler ) ] = 1.0 / ( 1.0E-11 * 1.0E-11 * SPEED_OF_LIGHT * SPEED_OF_LIGHT );
    simulatedObservations->setConstantWeightPerObservable( weightPerObservationParser );

    return simulatedObservations;

}

BOOST_AUTO_TEST_CASE( test_ObservationCollectionParser )
{
    const double startTime = double( 1.0E7 );
    const int numberOfDaysOfData = 3;
    unsigned int numberOfObservations = 300;
    double obsInterval = 20.0;

    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;
    std::map< ObservableType, double > obsStartTimes;
    std::map< ObservableType, std::vector< double > > baseTimeList;
    SystemOfBodies bodies;
    std::shared_ptr< ObservationCollection< double, double > > simulatedObservations = setUpObservationCollectionToTest(
            startTime, numberOfDaysOfData, numberOfObservations, obsInterval, stationReceiverLinkEnds, stationTransmitterLinkEnds, obsStartTimes, baseTimeList,bodies );

    // Check full list of observable types
    std::vector< ObservableType > observableTypes = simulatedObservations->getObservableTypes( );
    std::vector< ObservableType > trueObservableTypes = { one_way_range, one_way_doppler, angular_position };
    BOOST_CHECK( observableTypes.size( ) == trueObservableTypes.size( ) );
    for ( auto type : trueObservableTypes )
    {
        BOOST_CHECK( std::count( observableTypes.begin( ), observableTypes.end( ), type ) == 1 );
    }

    // Check full list of link ends
    std::vector< LinkEnds > linkEndsList = simulatedObservations->getLinkEnds( );
    std::vector< LinkEnds > trueLinkEndsList = { stationReceiverLinkEnds[ 0 ], stationTransmitterLinkEnds[ 0 ], stationReceiverLinkEnds[ 1 ], stationTransmitterLinkEnds[ 1 ],
                                                 stationReceiverLinkEnds[ 2 ], stationTransmitterLinkEnds[ 2 ] };
    BOOST_CHECK( linkEndsList.size( ) == trueLinkEndsList.size( ) );
    for ( auto linkEnds : trueLinkEndsList )
    {
        BOOST_CHECK( std::count( linkEndsList.begin( ), linkEndsList.end( ), linkEnds ) == 1 );
    }

    // Check full list of body names
    std::vector< std::string > bodyNamesInLinkEnds = simulatedObservations->getBodiesInLinkEnds( );
    std::vector< std::string > trueBodyNamesInLinkEnds = { "Earth", "Vehicle" };
    BOOST_CHECK( bodyNamesInLinkEnds.size( ) == trueBodyNamesInLinkEnds.size( ) );
    for ( auto name : trueBodyNamesInLinkEnds )
    {
        BOOST_CHECK( std::count( bodyNamesInLinkEnds.begin( ), bodyNamesInLinkEnds.end( ), name ) == 1 );
    }

    // Check full list of reference points
    std::vector< std::string > referencePointsInLinkEnds = simulatedObservations->getReferencePointsInLinkEnds( );
    std::vector< std::string > trueReferencePointsInLinkEnds = { "Station1", "Station2", "Station3" };
    BOOST_CHECK( referencePointsInLinkEnds.size( ) == trueReferencePointsInLinkEnds.size( ) );
    for ( auto refPoint : trueReferencePointsInLinkEnds )
    {
        BOOST_CHECK( std::count( referencePointsInLinkEnds.begin( ), referencePointsInLinkEnds.end( ), refPoint ) == 1 );
    }

    // Check full list of time bounds
    std::vector< std::pair< double, double > > timeBoundsList = simulatedObservations->getTimeBoundsList( );
    std::vector< std::pair< double, double > > trueTimeBoundsList;
    for ( auto it : obsStartTimes )
    {
        if ( it.first != angular_position )
        {
            trueTimeBoundsList.push_back( std::make_pair( it.second, it.second + ( numberOfObservations - 1 ) * obsInterval ) );
        }
        else
        {
            trueTimeBoundsList.push_back( std::make_pair( it.second, it.second + ( numberOfObservations - 1 ) * obsInterval + 2000.0 ) );
        }
    }
    BOOST_CHECK( timeBoundsList.size( ) == trueTimeBoundsList.size( ) );
    for ( auto timeBounds : trueTimeBoundsList )
    {
        BOOST_CHECK( std::count( timeBoundsList.begin( ), timeBoundsList.end( ), timeBounds ) == 1 );
    }


    // Check parsing based on observable type
    std::vector< std::shared_ptr< SingleObservationSet< > > > rangeObservationSets = simulatedObservations->getSingleObservationSets( observationParser( one_way_range ) );

    std::vector< std::shared_ptr< SingleObservationSet< > > > manuallyExtractedRangeObservationSets;
    std::map< LinkEnds, std::vector< std::shared_ptr< SingleObservationSet< > > > > manuallyExtractedRangeObservationSetsMap = simulatedObservations->getObservationsSets( ).at( one_way_range );
    for ( auto linkEndsIt : manuallyExtractedRangeObservationSetsMap )
    {
        for ( auto obsSet : linkEndsIt.second )
        {
            manuallyExtractedRangeObservationSets.push_back( obsSet );
        }
    }

    BOOST_CHECK( rangeObservationSets.size( ) == manuallyExtractedRangeObservationSets.size( ) );
    BOOST_CHECK( rangeObservationSets.size( ) == 3 );

    std::vector< LinkEnds > trueRangeLinkEnds = { stationReceiverLinkEnds[ 0 ], stationTransmitterLinkEnds[ 0 ], stationReceiverLinkEnds[ 1 ]  };
    std::vector< std::string > trueRangeBodyNames = { "Earth", "Vehicle" };
    std::vector< std::string > trueRangeReferencePoints = { "Station1", "Station2" };

    std::vector< LinkEnds > rangeLinkEnds;
    for ( auto obsSet : rangeObservationSets )
    {
        rangeLinkEnds.push_back( obsSet->getLinkEnds( ).linkEnds_ );
    }
    std::vector< LinkEnds > rangeLinkEndsFromObsCollection = simulatedObservations->getLinkEnds( observationParser( one_way_range ) );

    BOOST_CHECK( rangeLinkEnds.size( ) == rangeLinkEndsFromObsCollection.size( ) );
    BOOST_CHECK( rangeLinkEnds.size( ) == 3 );

    for ( auto testLinkEnds : trueRangeLinkEnds )
    {
        BOOST_CHECK( std::count( rangeLinkEnds.begin( ), rangeLinkEnds.end( ), testLinkEnds ) == 1 );
        BOOST_CHECK( std::count( rangeLinkEndsFromObsCollection.begin( ), rangeLinkEndsFromObsCollection.end( ), testLinkEnds ) == 1 );
    }

    std::vector< std::string > rangeBodyNamesFromObsCollection = simulatedObservations->getBodiesInLinkEnds( observationParser( one_way_range ) );
    for ( auto testBodyName : trueRangeBodyNames )
    {
        BOOST_CHECK( std::count( rangeBodyNamesFromObsCollection.begin( ), rangeBodyNamesFromObsCollection.end( ), testBodyName ) == 1 );
    }

    std::vector< std::string > rangeRefPointsFromObsCollection = simulatedObservations->getReferencePointsInLinkEnds( observationParser( one_way_range ) );
    for ( auto testRefPoint : trueRangeReferencePoints )
    {
        BOOST_CHECK( std::count( rangeRefPointsFromObsCollection.begin( ), rangeRefPointsFromObsCollection.end( ), testRefPoint ) == 1 );
    }

    // Retrieve observation values
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > rangeValuesFromObservationCollection =
            simulatedObservations->getObservations( observationParser( one_way_range ) );
    for ( unsigned int k = 0 ; k < rangeObservationSets.size( ) ; k++ )
    {
        BOOST_CHECK( rangeObservationSets.at( k )->getObservationsVector( ) == rangeValuesFromObservationCollection.at( k ) );
        BOOST_CHECK( manuallyExtractedRangeObservationSets.at( k )->getObservationsVector( ) == rangeValuesFromObservationCollection.at( k ) );
    }

    // Retrieve observation times
    std::vector< std::vector< double > > rangeTimesFromObservationCollection = simulatedObservations->getObservationTimes( observationParser( one_way_range ) );
    for ( unsigned int k = 0 ; k < rangeObservationSets.size( ) ; k++ )
    {
        BOOST_CHECK( rangeObservationSets.at( k )->getObservationTimes( ) == rangeTimesFromObservationCollection.at( k ) );
        BOOST_CHECK( manuallyExtractedRangeObservationSets.at( k )->getObservationTimes( ) == rangeTimesFromObservationCollection.at( k ) );
    }

    // Retrieve observation sets start and size
    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > > rangeSetsStartAndSize =
            simulatedObservations->getObservationSetStartAndSize( observationParser( one_way_range ) );
    BOOST_CHECK( utilities::createVectorFromMapKeys( rangeSetsStartAndSize ).size( ) == 1 );
    for ( auto linkEndsIt : rangeSetsStartAndSize.at( one_way_range ) )
    {
        for ( unsigned int k = 0 ; k < linkEndsIt.second.size( ) ; k++ )
        {
            BOOST_CHECK( linkEndsIt.second.at( k ) == rangeSetsStartAndSize.at( one_way_range ).at( linkEndsIt.first )[ k ] );
        }
    }

    // Check that pointers to selected observation sets are identical
    for ( unsigned int i = 0 ; i < rangeObservationSets.size( ) ; i++ )
    {
        BOOST_CHECK( rangeObservationSets.at( i ).get( ) == manuallyExtractedRangeObservationSets.at( i ).get( ) );
    }


    // Check parsing based on link ends
    std::vector< std::shared_ptr< SingleObservationSet< > > > linkEnds1ObservationSets = simulatedObservations->getSingleObservationSets(
            observationParser( std::vector< LinkEnds >( { stationReceiverLinkEnds[ 0 ], stationTransmitterLinkEnds[ 0 ] } ) ) );

    // Check parsing based on names
    std::vector< std::shared_ptr< SingleObservationSet< > > > station1ObservationSets = simulatedObservations->getSingleObservationSets(
            observationParser( "Station1", true ) );

    std::vector< std::shared_ptr< SingleObservationSet< > > > manuallyExtractedStation1ObservationSets;
    for ( auto observableIt : simulatedObservations->getObservationsSets( ) )
    {
        for ( auto linkEndsIt : observableIt.second )
        {
            if ( linkEndsIt.first == stationReceiverLinkEnds[ 0 ] || linkEndsIt.first == stationTransmitterLinkEnds[ 0 ]  )
            {
                for ( auto obsSet : linkEndsIt.second )
                {
                    manuallyExtractedStation1ObservationSets.push_back( obsSet );
                }
            }
        }
    }

    // Check that pointers to selected observation sets are identical
    BOOST_CHECK( linkEnds1ObservationSets.size( ) == manuallyExtractedStation1ObservationSets.size( ) );
    BOOST_CHECK( station1ObservationSets.size( ) == manuallyExtractedStation1ObservationSets.size( ) );
    for ( unsigned int i = 0 ; i < linkEnds1ObservationSets.size( ) ; i++ )
    {
        BOOST_CHECK( linkEnds1ObservationSets.at( i ).get( ) == manuallyExtractedStation1ObservationSets.at( i ).get( ) );
        BOOST_CHECK( station1ObservationSets.at( i ).get( ) == manuallyExtractedStation1ObservationSets.at( i ).get( ) );
    }

    // Retrieve observation sets start and size
    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > > station1SetsStartAndSize =
            simulatedObservations->getObservationSetStartAndSize( observationParser( "Station1", true ) );
    BOOST_CHECK( utilities::createVectorFromMapKeys( station1SetsStartAndSize ).size( ) == 1 );
    for ( auto linkEndsIt : rangeSetsStartAndSize.at( one_way_range ) )
    {
        if ( linkEndsIt.first == stationReceiverLinkEnds[ 0 ] || linkEndsIt.first == stationTransmitterLinkEnds[ 0 ]  )
        {
            for ( unsigned int k = 0 ; k < linkEndsIt.second.size( ) ; k++ )
            {
                BOOST_CHECK( linkEndsIt.second.at( k ) == station1SetsStartAndSize.at( one_way_range ).at( linkEndsIt.first )[ k ] );
            }
        }
        else
        {
            BOOST_CHECK( station1SetsStartAndSize.at( one_way_range ).count( linkEndsIt.first ) == 0 );
        }
    }

    // Check parsing based on time bounds
    std::vector< std::shared_ptr< SingleObservationSet< > > > firstDayObservationSets = simulatedObservations->getSingleObservationSets(
            observationParser( std::make_pair( startTime, startTime + 86400.0 ) ) );

    // Test opposite condition
    std::vector< std::shared_ptr< SingleObservationSet< > > > afterFirstDayObservationSets = simulatedObservations->getSingleObservationSets(
            observationParser( std::make_pair( startTime, startTime + 86400.0 ), true ) );

    // Manually retrieve observation sets based on time bounds
    std::vector< std::shared_ptr< SingleObservationSet< > > > manuallyDefinedFirstDayObservationSets, manuallyDefinedAfterFirstDayObservationSets;
    for ( auto observableIt : simulatedObservations->getObservationsSets( ) )
    {
        for ( auto linkEndsIt : observableIt.second )
        {
            for ( auto obsSet : linkEndsIt.second )
            {
                if ( ( obsSet->getTimeBounds( ).first >= startTime ) && ( obsSet->getTimeBounds( ).second <= startTime + 1000.0 + 86400.0 ) )
                {
                    manuallyDefinedFirstDayObservationSets.push_back( obsSet );
                }
                else
                {
                    manuallyDefinedAfterFirstDayObservationSets.push_back( obsSet );
                }
            }
        }
    }

    // Check that pointers to selected observation sets are identical
    BOOST_CHECK( firstDayObservationSets.size( ) == manuallyDefinedFirstDayObservationSets.size( ) );
    BOOST_CHECK( afterFirstDayObservationSets.size( ) == manuallyDefinedAfterFirstDayObservationSets.size( ) );
    for ( unsigned int i = 0 ; i < firstDayObservationSets.size( ) ; i++ )
    {
        BOOST_CHECK( firstDayObservationSets.at( i ).get( ) == manuallyDefinedFirstDayObservationSets.at( i ).get( ) );
        BOOST_CHECK( afterFirstDayObservationSets.at( i ).get( ) == manuallyDefinedAfterFirstDayObservationSets.at( i ).get( ) );
    }

    // Check multi-type parsing when conditions are treated separately
    std::vector< std::shared_ptr< ObservationCollectionParser > > multiTypeParserList;
    multiTypeParserList.push_back( observationParser( one_way_doppler ) );
    multiTypeParserList.push_back( observationParser( "Station1", true ) );
    multiTypeParserList.push_back( observationParser( stationTransmitterLinkEnds[ 1 ] ) );

    std::vector< std::shared_ptr< SingleObservationSet< > > > multiTypeObservationSets = simulatedObservations->getSingleObservationSets(
            observationParser( multiTypeParserList ) );

    std::vector< std::shared_ptr< SingleObservationSet< > > > manuallyDefinedMultiTypeObservationSets;
    for ( auto observableIt : simulatedObservations->getObservationsSets( ) )
    {
        for ( auto linkEndsIt : observableIt.second )
        {
            if ( ( observableIt.first == one_way_doppler ) || ( linkEndsIt.first == stationReceiverLinkEnds[ 0 ] ) || ( linkEndsIt.first == stationTransmitterLinkEnds[ 0 ] ) ||
                    ( linkEndsIt.first == stationTransmitterLinkEnds[ 1 ] ) )
            {
                for ( auto obs : linkEndsIt.second )
                {
                    manuallyDefinedMultiTypeObservationSets.push_back( obs );
                }
            }
        }
    }

    //Check that pointers to selected observation sets are identical
    BOOST_CHECK( multiTypeObservationSets.size( ) ==  manuallyDefinedMultiTypeObservationSets.size( ) );
    for ( unsigned int k = 0 ; k < multiTypeObservationSets.size( ) ; k++ )
    {
        BOOST_CHECK( multiTypeObservationSets.at( k ).get( ) == manuallyDefinedMultiTypeObservationSets.at( k ).get( ) );
    }


    // Check multi-type parsing when all conditions should be met concurrently
    std::vector< std::shared_ptr< ObservationCollectionParser > > multiTypeParserList2;
    multiTypeParserList2.push_back( observationParser( one_way_range ) );
    multiTypeParserList2.push_back( observationParser( std::make_pair( receiver, LinkEndId( "Earth", "Station2" ) ) ) );
    multiTypeParserList2.push_back( observationParser( "Station2", true ) );

    std::vector< std::shared_ptr< SingleObservationSet< > > > multiTypeObservationSets2 = simulatedObservations->getSingleObservationSets(
            observationParser( multiTypeParserList2, true ) );

    std::vector< std::shared_ptr< SingleObservationSet< > > > manuallyDefinedMultiTypeObservationSets2;
    for ( auto observableIt : simulatedObservations->getObservationsSets( ) )
    {
        for ( auto linkEndsIt : observableIt.second )
        {
            if ( ( observableIt.first == one_way_range ) && ( ( linkEndsIt.first == stationReceiverLinkEnds[ 1 ] ) ) )
            {
                for ( auto obs : linkEndsIt.second )
                {
                    manuallyDefinedMultiTypeObservationSets2.push_back( obs );
                }
            }
        }
    }

    //Check that pointers to selected observation sets are identical
    BOOST_CHECK( multiTypeObservationSets2.size( ) ==  manuallyDefinedMultiTypeObservationSets2.size( ) );
    for ( unsigned int k = 0 ; k < multiTypeObservationSets2.size( ) ; k++ )
    {
        BOOST_CHECK( multiTypeObservationSets2.at( k ).get( ) == manuallyDefinedMultiTypeObservationSets2.at( k ).get( ) );
    }

}

BOOST_AUTO_TEST_CASE( test_ObservationsFiltering )
{
    const double startTime = double( 1.0E7 );
    const int numberOfDaysOfData = 3;
    unsigned int numberOfObservations = 300;
    double obsInterval = 20.0;

    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;
    std::map< ObservableType, double > obsStartTimes;
    std::map< ObservableType, std::vector< double > > baseTimeList;
    SystemOfBodies bodies;
    std::shared_ptr< ObservationCollection< double, double > > simulatedObservations = setUpObservationCollectionToTest(
            startTime, numberOfDaysOfData, numberOfObservations, obsInterval, stationReceiverLinkEnds, stationTransmitterLinkEnds, obsStartTimes, baseTimeList, bodies );

    // Retrieve start indices and sizes for observation sets
    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > > obsSetsStartAndSize = simulatedObservations->getObservationSetStartAndSize( );

    // Retrieve range observation sets
    std::vector< std::shared_ptr< SingleObservationSet< > > > rangeObservationSets = simulatedObservations->getSingleObservationSets( observationParser( one_way_range ) );

    // Retrieve range observation values
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > rangeValuesFromObservationCollection =
            simulatedObservations->getObservations( observationParser( one_way_range ) );

    // Retrieve range observation times
    std::vector< std::vector< double > > rangeTimesFromObservationCollection = simulatedObservations->getObservationTimes( observationParser( one_way_range ) );

    // Test filtering based on observation values
    double cutOffValueMean = rangeObservationSets.at( 0 )->getObservationsVector( ).mean( );
    std::vector< double > maximumRangeValues = { rangeObservationSets.at( 0 )->getObservationsVector( ).maxCoeff( ),
                                                 rangeObservationSets.at( 1 )->getObservationsVector( ).maxCoeff( ),
                                                 rangeObservationSets.at( 2 )->getObservationsVector( ).maxCoeff( ) };
    double cutOffValueMax = *std::max_element( maximumRangeValues.begin( ), maximumRangeValues.end( ) );

    std::map< std::shared_ptr< ObservationCollectionParser >, std::shared_ptr< ObservationFilterBase > > filters =
            { { observationParser( one_way_range ), observationFilter( absolute_value_filtering, cutOffValueMean ) } };
    simulatedObservations->filterObservations( filters );

    // Manually compute post-filtering observations
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > manuallyComputedPostFilterRangeValues;
    for ( unsigned int k = 0 ; k < rangeValuesFromObservationCollection.size( ) ; k++ )
    {
        Eigen::Matrix< double, Eigen::Dynamic, 1 > unfilteredObs = rangeValuesFromObservationCollection.at( k );

        std::vector< unsigned int > indicesRemainingObs ;
        for ( unsigned int j = 0 ; j < unfilteredObs.size( ) ; j++ )
        {
            if ( unfilteredObs[ j ] <= cutOffValueMean )
            {
                indicesRemainingObs.push_back( j );
            }
        }

        Eigen::Matrix< double, Eigen::Dynamic, 1 > filteredObs = Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( indicesRemainingObs.size( ) );
        unsigned int filteredObsIndex = 0;
        for ( auto ind : indicesRemainingObs )
        {
            filteredObs[ filteredObsIndex ] = unfilteredObs[ ind ];
            filteredObsIndex += 1;
        }
        manuallyComputedPostFilterRangeValues.push_back( filteredObs );
    }

    // Retrieve range observations post-filtering
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > rangeValuesPostFiltering = simulatedObservations->getObservations( observationParser( one_way_range ) );

    // Check that all range values post filtering are below the cut-off value
    for ( unsigned int k = 0 ; k < rangeValuesPostFiltering.size( ) ; k++ )
    {
        if ( rangeValuesPostFiltering.at( k ).size( ) > 0 )
        {
            BOOST_CHECK( rangeValuesPostFiltering.at( k ).maxCoeff( ) <= cutOffValueMean );
        }
    }

    // Observation sets start and size after filtering
    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > > obsSetsStartAndSizeAfterFilter =
            simulatedObservations->getObservationSetStartAndSize( );
    int startIndex = 0;
    for ( auto observableIt : obsSetsStartAndSizeAfterFilter )
    {
        unsigned int linkEndsIndex = 0;
        for ( auto linkEndsIt : observableIt.second )
        {
            for ( auto indices : linkEndsIt.second )
            {
                BOOST_CHECK( indices.first == startIndex );
                if ( observableIt.first == one_way_range )
                {
                    BOOST_CHECK( indices.second == rangeValuesPostFiltering[ linkEndsIndex ].size( ) );
                }
                if ( observableIt.first == angular_position )
                {
                    BOOST_CHECK( indices.second == 2.0 * numberOfObservations );
                }
                if ( observableIt.first == one_way_doppler )
                {
                    BOOST_CHECK( indices.second == static_cast< int >( numberOfObservations ) );
                }
                startIndex += indices.second;
            }
            linkEndsIndex += 1;
        }
    }

    // Reintroduce observations
    std::map< std::shared_ptr< ObservationCollectionParser >, std::shared_ptr< ObservationFilterBase > > defilter =
            { { observationParser( one_way_range ), observationFilter( absolute_value_filtering, cutOffValueMean, false ) } };
    simulatedObservations->filterObservations( defilter );

    // Retrieve range observations post-defiltering
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > rangeValuesPostDefiltering = simulatedObservations->getObservations( observationParser( one_way_range ) );
    std::vector< std::vector< double > > rangeTimesPostDefiltering = simulatedObservations->getObservationTimes( observationParser( one_way_range ) );

    // Check consistency after re-introducing all observations (important to check observations order)
    for ( unsigned int k = 0 ; k < rangeValuesPostDefiltering.size( ) ; k++ )
    {
        BOOST_CHECK( rangeValuesFromObservationCollection.at( k ) == rangeValuesPostDefiltering.at( k ) );
        BOOST_CHECK( rangeTimesFromObservationCollection.at( k ) == rangeTimesPostDefiltering.at( k ) );
    }

    // Observation sets start and size after de-filtering (should be equal to original ones)
    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > > obsSetsStartAndSizeAfterDefilter =
            simulatedObservations->getObservationSetStartAndSize( );
    for ( auto observableIt : obsSetsStartAndSizeAfterDefilter )
    {
        for ( auto linkEndsIt : observableIt.second )
        {
            for ( unsigned int k = 0 ; k < linkEndsIt.second.size( ) ; k++ )
            {
                BOOST_CHECK( obsSetsStartAndSize.at( observableIt.first ).at( linkEndsIt.first ).at( k ).first == linkEndsIt.second.at( k ).first );
                BOOST_CHECK( obsSetsStartAndSize.at( observableIt.first ).at( linkEndsIt.first ).at( k ).second == linkEndsIt.second.at( k ).second );
            }
        }
    }

    // Test filtering based on epoch values and time bounds
    std::vector< double > rangeObsTimes = baseTimeList.at( one_way_range );
    std::pair< double, double > timeBounds = std::make_pair( rangeObsTimes[ int(numberOfObservations/3) ], rangeObsTimes[ 2*int(numberOfObservations/3) ] );
    std::vector< double > epochsToFilter;
    std::vector< double > epochsLeft;
    for ( auto time : rangeObsTimes )
    {
        if ( time >= timeBounds.first && time <= timeBounds.second )
        {
            epochsToFilter.push_back( time );
        }
        else
        {
            epochsLeft.push_back( time );
        }
    }

    // Filter out observations that are outside the above-defined time bounds, using epochs-based filtering
    simulatedObservations->filterObservations( observationFilter( epochs_filtering, epochsToFilter, true, true ), observationParser( one_way_range ) );

    // Retrieve and check remaining observation times
    std::vector< std::vector< double > > remainingObsTimes = simulatedObservations->getObservationTimes( observationParser( one_way_range ) );
    for ( unsigned int k = 0 ; k < remainingObsTimes.size( ) ; k++ )
    {
        for ( unsigned int j = 0 ; j < remainingObsTimes.at( k ).size( ) ; j++ )
        {
            BOOST_CHECK( remainingObsTimes.at( k ).at( j ) == epochsToFilter.at( j ) );
        }
    }

    // Reintroduce observations
    simulatedObservations->filterObservations( observationFilter( epochs_filtering, epochsToFilter, false, true ), observationParser( one_way_range ) );

    // Retrieve and check remaining observation times
    remainingObsTimes = simulatedObservations->getObservationTimes( observationParser( one_way_range ) );
    for ( unsigned int k = 0 ; k < remainingObsTimes.size( ) ; k++ )
    {
        BOOST_CHECK( remainingObsTimes.at( k ).size( ) == numberOfObservations );
    }


    // Filter out observations that are outside the above defined time bounds, using time bounds filtering directly
    simulatedObservations->filterObservations( observationFilter( time_bounds_filtering, timeBounds, true, true ), observationParser( one_way_range ) );

    // Retrieve and check remaining observation times
    remainingObsTimes = simulatedObservations->getObservationTimes( observationParser( one_way_range ) );
    for ( unsigned int k = 0 ; k < remainingObsTimes.size( ) ; k++ )
    {
        for ( unsigned int j = 0 ; j < remainingObsTimes.at( k ).size( ) ; j++ )
        {
            BOOST_CHECK( remainingObsTimes.at( k ).at( j ) == epochsToFilter.at( j ) );
        }
    }

    // Re-introduce observations
    simulatedObservations->filterObservations( observationFilter( time_bounds_filtering, timeBounds, false, true ), observationParser( one_way_range ) );

    // Retrieve and check remaining observation times
    remainingObsTimes = simulatedObservations->getObservationTimes( observationParser( one_way_range ) );
    for ( unsigned int k = 0 ; k < remainingObsTimes.size( ) ; k++ )
    {
        BOOST_CHECK( remainingObsTimes.at( k ).size( ) == numberOfObservations );
    }

    // Test filtering based on residual values
    double residualCutOffValue = 0.25;
    unsigned int nbRangeStation1 = rangeObservationSets.at( 0 )->getNumberOfObservables( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > residualsStation1 = ( residualCutOffValue + 0.05 ) * Eigen::Matrix< double, Eigen::Dynamic, 1 >::Ones( nbRangeStation1, 1 );
    unsigned int nbRangeStation2 = rangeObservationSets.at( 1 )->getNumberOfObservables( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > residualsStation2 = ( residualCutOffValue - 0.05 ) * Eigen::Matrix< double, Eigen::Dynamic, 1 >::Ones( nbRangeStation2, 1 );
    unsigned int nbRangeStation3 = rangeObservationSets.at( 2 )->getNumberOfObservables( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > residualsStation3 = ( residualCutOffValue - 0.05 ) * Eigen::Matrix< double, Eigen::Dynamic, 1 >::Ones( nbRangeStation3, 1 );
    unsigned int nbLargeResidualsStation3 = 0;
    for ( unsigned int i = 0 ; i < nbRangeStation3 ; i++ )
    {
        if ( i%2 )
        {
            residualsStation3[ i ] = ( residualCutOffValue + 0.05 );
            nbLargeResidualsStation3 += 1;
        }
    }
    rangeObservationSets.at( 0 )->setResiduals( residualsStation1 );
    rangeObservationSets.at( 1 )->setResiduals( residualsStation2 );
    rangeObservationSets.at( 2 )->setResiduals( residualsStation3 );

    std::map< std::shared_ptr< ObservationCollectionParser >, std::shared_ptr< ObservationFilterBase > > residualFilter =
            { { observationParser( one_way_range ), observationFilter( residual_filtering, residualCutOffValue ) } };
    simulatedObservations->filterObservations( residualFilter );

    BOOST_CHECK( rangeObservationSets.at( 0 )->getNumberOfObservables( ) == 0 );
    BOOST_CHECK( rangeObservationSets.at( 0 )->getNumberOfFilteredObservations( ) == nbRangeStation1 );
    BOOST_CHECK( rangeObservationSets.at( 1 )->getNumberOfObservables( ) == nbRangeStation2 );
    BOOST_CHECK( rangeObservationSets.at( 1 )->getNumberOfFilteredObservations( ) == 0 );
    BOOST_CHECK( rangeObservationSets.at( 2 )->getNumberOfObservables( ) == nbLargeResidualsStation3 );
    BOOST_CHECK( rangeObservationSets.at( 2 )->getNumberOfFilteredObservations( ) == nbRangeStation3 - nbLargeResidualsStation3 );

    // Check that all residual values post filtering are below the cut-off value
    std::vector< Eigen::VectorXd > rangeResidualsPostFiltering = simulatedObservations->getResiduals( observationParser( one_way_range ) );
    for ( unsigned int k = 0 ; k < rangeResidualsPostFiltering.size( ) ; k++ )
    {
        if ( rangeResidualsPostFiltering.at( k ).size( ) > 0 )
        {
            BOOST_CHECK( rangeResidualsPostFiltering.at( k ).maxCoeff( ) <= residualCutOffValue );
        }
    }

    // Check that pointers to empty single observation sets are still properly defined post filtering
    std::vector< std::shared_ptr< SingleObservationSet< > > > postFilterRangeObservationSets = simulatedObservations->getSingleObservationSets( observationParser( one_way_range ) );
    for ( unsigned int k= 0 ; k < rangeObservationSets.size( ) ; k++ )
    {
        BOOST_CHECK( rangeObservationSets.at( k ).get( ) == postFilterRangeObservationSets.at( k ).get( ) );
    }

    // Remove empty observation sets post-filtering (one set should be removed)
    simulatedObservations->removeEmptySingleObservationSets( );
    BOOST_CHECK( simulatedObservations->getSingleObservationSets( observationParser( one_way_range ) ).size( ) == 2 );

    // Remove specific observation set
    std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > > indexSetToRemove =
            { { angular_position, { {  stationReceiverLinkEnds[ 2 ], std::vector< unsigned int >( { 0 } ) } } } };
    simulatedObservations->removeSingleObservationSets( indexSetToRemove );
    BOOST_CHECK( simulatedObservations->getSingleObservationSets( observationParser( angular_position ) ).size( ) == 1 );


    // Retrieve elevation angles before filtering
    std::shared_ptr< ObservationDependentVariableSettings > elevationAngleSettings = elevationAngleDependentVariable( );
    std::vector< Eigen::MatrixXd > originalElevationAngles = simulatedObservations->getDependentVariables( elevationAngleSettings, false ).first;
    std::vector< unsigned int > nbObservationsToKeepPerSet;
    for ( auto elevationAngleMatrix : originalElevationAngles )
    {
        unsigned int nbPositiveElevationAngles = 0;
        for ( unsigned int i = 0 ; i < elevationAngleMatrix.rows( ) ; i++ )
        {
            if ( elevationAngleMatrix( i, 0 ) > 0.0 )
            {
                nbPositiveElevationAngles += 1;
            }
        }
        nbObservationsToKeepPerSet.push_back( nbPositiveElevationAngles );
    }

    // Filter out observations with negative elevation angles
    simulatedObservations->filterObservations( observationFilter( elevationAngleSettings, Eigen::Vector1d::Zero( ), true, true ) );

    std::vector< Eigen::MatrixXd > remainingElevationAngles = simulatedObservations->getDependentVariables( elevationAngleSettings, false ).first;
    for ( unsigned int i = 0 ; i < remainingElevationAngles.size( ) ; i++ )
    {
        BOOST_CHECK( remainingElevationAngles[ i ].rows( ) == nbObservationsToKeepPerSet[ i ] );
    }



}

BOOST_AUTO_TEST_CASE( test_ObservationsSplitting )
{

    const double startTime = double( 1.0E7 );
    const int numberOfDaysOfData = 3;
    unsigned int numberOfObservations = 300;
    double obsInterval = 20.0;

    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;
    std::map< ObservableType, double > obsStartTimes;
    std::map< ObservableType, std::vector< double > > baseTimeList;
    SystemOfBodies bodies;
    std::shared_ptr< ObservationCollection< double, double > > simulatedObservations = setUpObservationCollectionToTest(
            startTime, numberOfDaysOfData, numberOfObservations, obsInterval, stationReceiverLinkEnds, stationTransmitterLinkEnds, obsStartTimes, baseTimeList, bodies );


    // Test splitting based on time tags
    unsigned int indexSplitEpoch = static_cast< unsigned int >(numberOfObservations/3);
    std::vector< double > rangeObsTimes = baseTimeList.at( one_way_range );
    double splitEpoch = rangeObsTimes.at( indexSplitEpoch );
    std::vector< double > timeTags = { rangeObsTimes.at( 0 ), splitEpoch, rangeObsTimes.back( ) };
    simulatedObservations->splitObservationSets( observationSetSplitter( time_tags_splitter, timeTags ), observationParser( stationTransmitterLinkEnds[ 0 ] ) );

    // Check that the original observation set got split in two independent sets (before/after splitEpoch)
    std::vector< std::shared_ptr< SingleObservationSet< > > > splitObsSets = simulatedObservations->getSingleObservationSets( observationParser( stationTransmitterLinkEnds[ 0 ] ) );
    BOOST_CHECK( splitObsSets.size( ) == 2 );
    BOOST_CHECK( splitObsSets.at( 0 )->getNumberOfObservables( ) == indexSplitEpoch+1 );
    BOOST_CHECK( splitObsSets.at( 1 )->getNumberOfObservables( ) == numberOfObservations - (indexSplitEpoch+1) );

    // Check that all observations in the first set occur before splitEpoch
    for ( auto time : splitObsSets.at( 0 )->getObservationTimes( ) )
    {
        BOOST_CHECK( time <= splitEpoch );
    }
    // Check that all observations in the first set occur after splitEpoch
    for ( auto time : splitObsSets.at( 1 )->getObservationTimes( ) )
    {
        BOOST_CHECK( time > splitEpoch );
    }


    // Test splitting based on time span. Define max time span such that one observation set would be left with one observation only
    // (to check if the minimum number of observations condition works)
    double maxTimeSpanObs = (indexSplitEpoch - 1) * obsInterval;
    simulatedObservations->splitObservationSets( observationSetSplitter( time_span_splitter, maxTimeSpanObs, 10 ), observationParser( stationTransmitterLinkEnds[ 0 ] ) );

    // Check that we now have three independent sets (out of the former two, only one set could be split while retaining enough observations)
    splitObsSets = simulatedObservations->getSingleObservationSets( observationParser( stationTransmitterLinkEnds[ 0 ] ) );
    BOOST_CHECK( splitObsSets.size( ) == 3 );
    BOOST_CHECK( splitObsSets.at( 0 )->getNumberOfObservables( ) == indexSplitEpoch );
    BOOST_CHECK( splitObsSets.at( 1 )->getNumberOfObservables( ) == indexSplitEpoch );
    BOOST_CHECK( splitObsSets.at( 2 )->getNumberOfObservables( ) == numberOfObservations - (2*indexSplitEpoch+1) );

    // Check that the time span of the new observation set meets the threshold
    for ( unsigned int i = 0 ; i < 3 ; i++ )
    {
        std::vector< double > newObsTimes = splitObsSets.at( i )->getObservationTimes( );
        BOOST_CHECK( ( newObsTimes.back(  ) - newObsTimes.at( 0 ) ) <= maxTimeSpanObs );
    }

    // Test splitting based on max number of observations
    simulatedObservations->splitObservationSets( observationSetSplitter( nb_observations_splitter, 50, 50 ), observationParser( stationTransmitterLinkEnds[ 0 ] ) );

    // Check that we now have five independent sets (only two of the former three sets contain enough observations to be split in two)
    splitObsSets = simulatedObservations->getSingleObservationSets( observationParser( stationTransmitterLinkEnds[ 0 ] ) );
    BOOST_CHECK( splitObsSets.size( ) == 5 );

    // Check the number of observations in the new observation sets
    for ( unsigned int i = 0 ; i < 5 ; i++ )
    {
        BOOST_CHECK( splitObsSets.at( i )->getNumberOfObservables( ) == 50 );
    }

    // Test splitting based on time interval
    simulatedObservations->splitObservationSets( observationSetSplitter( time_interval_splitter, 1500.0 ), observationParser( stationReceiverLinkEnds[ 2 ] ) );

    // Check that the observation set associated with stationReceiverLinkEnds[ 2 ] (of type angular_position) got split into two independent sets
    splitObsSets = simulatedObservations->getSingleObservationSets( observationParser( stationReceiverLinkEnds[ 2 ] ) );
    BOOST_CHECK( splitObsSets.size( ) == 2 );

    // Check that the split is performed properly
    for ( unsigned int k = 0 ; k < 2 ; k++ )
    {
        for ( unsigned int i = 0 ; i < splitObsSets.at( k )->getNumberOfObservables( ) - 1 ; i++ )
        {
            BOOST_CHECK( ( splitObsSets.at( k )->getObservationTimes( ).at( i + 1 ) - splitObsSets.at( k )->getObservationTimes( ).at( i ) ) <= 1500.0 );
        }
    }
    BOOST_CHECK( ( splitObsSets.at( 1 )->getObservationTimes( ).at( 0 ) - splitObsSets.at( 0 )->getObservationTimes( ).back( ) ) > 1500.0 );


}

BOOST_AUTO_TEST_SUITE_END( )

}

}



