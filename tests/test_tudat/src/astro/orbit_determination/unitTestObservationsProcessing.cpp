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
#include <set>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/podProcessing.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::basic_astrodynamics;
using namespace tudat::coordinate_conversions;
using namespace tudat::ephemerides;
using namespace tudat::estimatable_parameters;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::observation_models;
using namespace tudat::orbit_determination;
using namespace tudat::orbital_element_conversions;
using namespace tudat::physical_constants;
using namespace tudat::propagators;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;

BOOST_AUTO_TEST_SUITE( test_observations_processing )

std::shared_ptr< ObservationDataset< double, double > > setUpObservationDatasetToTest(
        const double startTime,
        const int numberOfDaysOfData,
        const int numberOfObservations,
        const double observationsInterval,
        std::vector< LinkEnds >& stationReceiverLinkEnds,
        std::vector< LinkEnds >& stationTransmitterLinkEnds,
        std::map< ObservableType, double >& obsStartTimes,
        std::map< ObservableType, std::vector< double > >& baseTimeList,
        SystemOfBodies& bodies )
{
    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Mars" );

    const double initialEphemerisTime = startTime;
    const double finalEphemerisTime = initialEphemerisTime + numberOfDaysOfData * 86400.0;

    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
            "ECLIPJ2000",
            "IAU_Earth",
            spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
            initialEphemerisTime,
            2.0 * mathematical_constants::PI / physical_constants::JULIAN_DAY );

    bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( 400.0 );

    bodies.at( "Vehicle" )
            ->setEphemeris( std::make_shared< TabulatedCartesianEphemeris<> >(
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d > >( ), "Earth", "ECLIPJ2000" ) );

    const std::vector< std::string > groundStationNames = { "Station1", "Station2", "Station3" };
    createGroundStation( bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station2", ( Eigen::Vector3d( ) << 0.0, -0.55, 2.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station3", ( Eigen::Vector3d( ) << 0.0, 0.05, 4.0 ).finished( ), geodetic_position );

    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 8, 8 ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    std::vector< std::string > bodiesToIntegrate = { "Vehicle" };
    std::vector< std::string > centralBodies = { "Earth" };
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7200.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    const double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Matrix< double, 6, 1 > systemInitialState =
            convertKeplerianToCartesianElements( asterixInitialStateInKeplerianElements, earthGravitationalParameter );

    std::shared_ptr< TranslationalStatePropagatorSettings< double, double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double, double > >(
                    centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, finalEphemerisTime, cowell );

    std::shared_ptr< IntegratorSettings< double > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< double > >(
            initialEphemerisTime, 40.0, CoefficientSets::rungeKuttaFehlberg78, 40.0, 40.0, 1.0, 1.0 );

    stationReceiverLinkEnds.clear( );
    stationTransmitterLinkEnds.clear( );
    for( const std::string& groundStationName : groundStationNames )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = LinkEndId( "Earth", groundStationName );
        linkEnds[ receiver ] = LinkEndId( "Vehicle", "" );
        stationTransmitterLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = LinkEndId( "Earth", groundStationName );
        linkEnds[ transmitter ] = LinkEndId( "Vehicle", "" );
        stationReceiverLinkEnds.push_back( linkEnds );
    }

    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    linkEndsPerObservable[ one_way_range ] = { stationReceiverLinkEnds[ 0 ],
                                               stationTransmitterLinkEnds[ 0 ],
                                               stationReceiverLinkEnds[ 1 ] };
    linkEndsPerObservable[ one_way_doppler ] = { stationReceiverLinkEnds[ 1 ], stationTransmitterLinkEnds[ 2 ] };
    linkEndsPerObservable[ angular_position ] = { stationReceiverLinkEnds[ 2 ], stationTransmitterLinkEnds[ 1 ] };

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back(
            std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >( "Vehicle", systemInitialState, "Earth" ) );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, bodies );

    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    for( const auto& linkEndIterator : linkEndsPerObservable )
    {
        for( const LinkEnds& linkEnds : linkEndIterator.second )
        {
            observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( linkEndIterator.first, linkEnds ) );
        }
    }

    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies, parametersToEstimate, observationSettingsList, integratorSettings, propagatorSettings );

    baseTimeList.clear( );
    std::vector< double > rangeObsTimes, dopplerObsTimes, angularPositionObsTimes;
    obsStartTimes = { { one_way_range, initialEphemerisTime + 1000.0 },
                      { one_way_doppler, initialEphemerisTime + 1000.0 + 86400.0 },
                      { angular_position, initialEphemerisTime + 1000.0 + 2.0 * 86400.0 } };
    for( int j = 0; j < numberOfObservations; ++j )
    {
        rangeObsTimes.push_back( obsStartTimes.at( one_way_range ) + static_cast< double >( j ) * observationsInterval );
        dopplerObsTimes.push_back( obsStartTimes.at( one_way_doppler ) + static_cast< double >( j ) * observationsInterval );

        const double angularPositionObsBuffer = ( ( j > numberOfObservations / 2 ) ? 3000.0 : 1000.0 );
        angularPositionObsTimes.push_back( initialEphemerisTime + angularPositionObsBuffer + 2.0 * 86400.0 +
                                           static_cast< double >( j ) * observationsInterval );
    }
    baseTimeList[ one_way_range ] = rangeObsTimes;
    baseTimeList[ one_way_doppler ] = dopplerObsTimes;
    baseTimeList[ angular_position ] = angularPositionObsTimes;

    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
    for( const auto& linkEndIterator : linkEndsPerObservable )
    {
        for( const LinkEnds& linkEnds : linkEndIterator.second )
        {
            measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< double > >(
                    linkEndIterator.first, linkEnds, baseTimeList.at( linkEndIterator.first ), receiver ) );
        }
    }

    addDependentVariablesToObservationSimulationSettings(
            measurementSimulationInput,
            std::vector< std::shared_ptr< ObservationDependentVariableSettings > >( { elevationAngleDependentVariable( ) } ),
            bodies );

    std::shared_ptr< ObservationDataset< double, double > > dataset = simulateObservationDataset< double, double >(
            measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    dataset->setConstantSingleObservationScalarWeight( ObservationSelectionCondition< double, double >::observableType( one_way_range ),
                                                       1.0 / ( 3.0 * 3.0 ) );
    dataset->setConstantSingleObservationScalarWeight( ObservationSelectionCondition< double, double >::observableType( angular_position ),
                                                       1.0 / ( 1.0E-5 * 1.0E-5 ) );
    dataset->setConstantSingleObservationScalarWeight( ObservationSelectionCondition< double, double >::observableType( one_way_doppler ),
                                                       1.0 / ( 1.0E-11 * 1.0E-11 * SPEED_OF_LIGHT * SPEED_OF_LIGHT ) );

    return dataset;
}

std::vector< unsigned int > getSetIdsForObservable( const ObservationDataset< double, double >& dataset,
                                                    const ObservableType observableType )
{
    std::vector< unsigned int > setIds;
    for( unsigned int setId = 0; setId < dataset.getNumberOfObservationSets( ); ++setId )
    {
        if( dataset.getObservationSetMetadata( setId ).observableType_ == observableType )
        {
            setIds.push_back( setId );
        }
    }
    return setIds;
}

std::vector< unsigned int > getSetIdsForLinkEnd( const ObservationDataset< double, double >& dataset, const LinkEndId& linkEndId )
{
    std::vector< unsigned int > setIds;
    for( unsigned int setId = 0; setId < dataset.getNumberOfObservationSets( ); ++setId )
    {
        const LinkEnds& linkEnds = dataset.getLinkDefinition( dataset.getObservationSetMetadata( setId ).linkDefinitionId_ ).linkEnds_;
        for( const auto& linkEndIterator : linkEnds )
        {
            if( linkEndIterator.second == linkEndId )
            {
                setIds.push_back( setId );
                break;
            }
        }
    }
    return setIds;
}

std::size_t getScalarSizeForSetIds( const ObservationDataset< double, double >& dataset, const std::vector< unsigned int >& setIds )
{
    std::size_t totalSize = 0;
    for( const unsigned int setId : setIds )
    {
        totalSize += dataset.getTotalScalarSizeForSet( setId );
    }
    return totalSize;
}

/*!
 * Verifies dataset metadata and condition-based selection on simulated data.
 *
 * Test outline: creates a representative multi-station range, Doppler and angular
 * dataset. It checks the registered observable/link metadata, selection by
 * observable type and link end, viewer creation and scalar flattened data size.
 */
BOOST_AUTO_TEST_CASE( test_dataset_metadata_and_selection )
{
    const double startTime = 1.0E7;
    const int numberOfDaysOfData = 3;
    const int numberOfObservations = 300;
    const double obsInterval = 20.0;

    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;
    std::map< ObservableType, double > obsStartTimes;
    std::map< ObservableType, std::vector< double > > baseTimeList;
    SystemOfBodies bodies;
    std::shared_ptr< ObservationDataset< double, double > > dataset = setUpObservationDatasetToTest( startTime,
                                                                                                     numberOfDaysOfData,
                                                                                                     numberOfObservations,
                                                                                                     obsInterval,
                                                                                                     stationReceiverLinkEnds,
                                                                                                     stationTransmitterLinkEnds,
                                                                                                     obsStartTimes,
                                                                                                     baseTimeList,
                                                                                                     bodies );

    // The setup must produce the expected number of observable/link combinations and rows.
    BOOST_CHECK_EQUAL( dataset->getNumberOfObservationSets( ), 7 );
    BOOST_CHECK_EQUAL( dataset->getNumberOfObservations( ), 7 * numberOfObservations );

    std::set< ObservableType > observableTypes;
    std::set< LinkEnds > linkEndsList;
    std::set< std::string > bodyNames;
    std::set< std::string > referencePointNames;
    for( unsigned int setId = 0; setId < dataset->getNumberOfObservationSets( ); ++setId )
    {
        const ObservationSetMetadata< double, double >& metadata = dataset->getObservationSetMetadata( setId );
        observableTypes.insert( metadata.observableType_ );
        const LinkEnds& linkEnds = dataset->getLinkDefinition( metadata.linkDefinitionId_ ).linkEnds_;
        linkEndsList.insert( linkEnds );
        for( const auto& linkEndIterator : linkEnds )
        {
            bodyNames.insert( linkEndIterator.second.bodyName_ );
            if( !linkEndIterator.second.stationName_.empty( ) )
            {
                referencePointNames.insert( linkEndIterator.second.stationName_ );
            }
        }
    }

    // Metadata stored in the dataset must recover all observable types, station link definitions and link-end names.
    BOOST_CHECK( observableTypes == std::set< ObservableType >( { one_way_range, one_way_doppler, angular_position } ) );
    BOOST_CHECK_EQUAL( linkEndsList.count( stationReceiverLinkEnds.at( 0 ) ), 1 );
    BOOST_CHECK_EQUAL( linkEndsList.count( stationTransmitterLinkEnds.at( 0 ) ), 1 );
    BOOST_CHECK_EQUAL( linkEndsList.count( stationReceiverLinkEnds.at( 1 ) ), 1 );
    BOOST_CHECK_EQUAL( linkEndsList.count( stationTransmitterLinkEnds.at( 1 ) ), 1 );
    BOOST_CHECK_EQUAL( linkEndsList.count( stationReceiverLinkEnds.at( 2 ) ), 1 );
    BOOST_CHECK_EQUAL( linkEndsList.count( stationTransmitterLinkEnds.at( 2 ) ), 1 );
    BOOST_CHECK( bodyNames == std::set< std::string >( { "Earth", "Vehicle" } ) );
    BOOST_CHECK( referencePointNames == std::set< std::string >( { "Station1", "Station2", "Station3" } ) );

    const std::vector< unsigned int > rangeSetIds = getSetIdsForObservable( *dataset, one_way_range );
    const std::vector< unsigned int > station1SetIds = getSetIdsForLinkEnd( *dataset, LinkEndId( "Earth", "Station1" ) );

    // Dataset metadata helpers must select range sets and Station1 sets across observable types.
    BOOST_CHECK_EQUAL( rangeSetIds.size( ), 3 );
    BOOST_CHECK_EQUAL( station1SetIds.size( ), 2 );

    const std::vector< unsigned int > firstDayRangeIds = dataset->getObservationIdsMatchingCondition(
            ObservationSelectionCondition< double, double >::observableType( one_way_range ) &&
            ObservationSelectionCondition< double, double >::timeBounds( startTime, startTime + 86400.0 ) );

    // Combining observable and time conditions must select one day of range observations for all stations.
    BOOST_CHECK_EQUAL( firstDayRangeIds.size( ), 3 * numberOfObservations );

    const ObservationSelectionCondition< double, double > station2RangeCondition =
            ObservationSelectionCondition< double, double >::observableType( one_way_range ) &&
            ObservationSelectionCondition< double, double >::linkEnd( receiver, LinkEndId( "Earth", "Station2" ) );
    const std::vector< unsigned int > station2RangeIds = dataset->getObservationIdsMatchingCondition( station2RangeCondition );

    // Combining observable and link-end conditions must isolate one station's range set.
    BOOST_CHECK_EQUAL( station2RangeIds.size( ), numberOfObservations );

    const ObservationDatasetViewer< double, double > station2RangeView = dataset->createViewer( station2RangeCondition );

    // A viewer created from the same condition must expose the selected rows and scalar flattened data only.
    BOOST_CHECK_EQUAL( station2RangeView.getNumberOfObservations( ), numberOfObservations );
    BOOST_CHECK_EQUAL( station2RangeView.createEstimationFlattenedObservationData( ).getObservationVector( ).size( ),
                       numberOfObservations );

    const std::size_t expectedScalarSize = getScalarSizeForSetIds( *dataset, getSetIdsForObservable( *dataset, one_way_range ) ) +
            getScalarSizeForSetIds( *dataset, getSetIdsForObservable( *dataset, one_way_doppler ) ) +
            getScalarSizeForSetIds( *dataset, getSetIdsForObservable( *dataset, angular_position ) );

    // The full flattened data size must equal the sum of scalar sizes over all observable types.
    BOOST_CHECK_EQUAL( dataset->createEstimationFlattenedObservationData( ).getObservationVector( ).size( ), expectedScalarSize );
}

/*!
 * Verifies rejection, restoration and reduced dataset creation on simulated data.
 *
 * Test outline: rejects range observations selected by value, confirms active and
 * rejected flattened data/viewers have the expected sizes, restores the data and
 * then creates a reduced dataset for a time-window selection.
 */
BOOST_AUTO_TEST_CASE( test_dataset_rejection_restoration_and_reduced_views )
{
    const double startTime = 1.0E7;
    const int numberOfDaysOfData = 3;
    const int numberOfObservations = 300;
    const double obsInterval = 20.0;

    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;
    std::map< ObservableType, double > obsStartTimes;
    std::map< ObservableType, std::vector< double > > baseTimeList;
    SystemOfBodies bodies;
    std::shared_ptr< ObservationDataset< double, double > > dataset = setUpObservationDatasetToTest( startTime,
                                                                                                     numberOfDaysOfData,
                                                                                                     numberOfObservations,
                                                                                                     obsInterval,
                                                                                                     stationReceiverLinkEnds,
                                                                                                     stationTransmitterLinkEnds,
                                                                                                     obsStartTimes,
                                                                                                     baseTimeList,
                                                                                                     bodies );

    const std::vector< unsigned int > rangeSetIds = getSetIdsForObservable( *dataset, one_way_range );
    double cutOffValueMean = 0.0;
    for( const unsigned int observationId : dataset->getObservationIdsForSet( rangeSetIds.at( 0 ) ) )
    {
        cutOffValueMean += dataset->getObservationValue( observationId )( 0 );
    }
    cutOffValueMean /= static_cast< double >( dataset->getNumberOfObservationsForSet( rangeSetIds.at( 0 ) ) );

    const Eigen::VectorXd rangeLimit = ( Eigen::VectorXd( 1 ) << cutOffValueMean ).finished( );
    const ObservationSelectionCondition< double, double > highRangeValues =
            ObservationSelectionCondition< double, double >::observableType( one_way_range ) &&
            ObservationSelectionCondition< double, double >::observationAbsoluteValueGreaterThan( rangeLimit );

    const std::vector< unsigned int > rejectedRangeIds = dataset->getObservationIdsMatchingCondition( highRangeValues );

    // The value threshold must identify a non-empty subset before mutating row status.
    BOOST_CHECK( !rejectedRangeIds.empty( ) );

    const int originalScalarSize = dataset->createEstimationFlattenedObservationData( ).getObservationVector( ).size( );
    dataset->rejectObservations( highRangeValues, "range value threshold" );

    // Rejection must affect rejected-row queries, rejected viewers and active-only flattened data consistently.
    BOOST_CHECK_EQUAL( dataset->getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::rejected( ) ).size( ),
                       rejectedRangeIds.size( ) );
    BOOST_CHECK_EQUAL( dataset->createViewer( ObservationSelectionCondition< double, double >::rejected( ) ).getNumberOfObservations( ),
                       rejectedRangeIds.size( ) );
    BOOST_CHECK_EQUAL( dataset->createEstimationFlattenedObservationData( true ).getObservationVector( ).size( ), originalScalarSize );
    BOOST_CHECK_EQUAL( dataset->createEstimationFlattenedObservationData( ).getObservationVector( ).size( ),
                       originalScalarSize - static_cast< int >( rejectedRangeIds.size( ) ) );

    const std::shared_ptr< ObservationDataset< double, double > > activeDataset =
            dataset->createNewAndKeep( ObservationSelectionCondition< double, double >::active( ) );

    // A reduced active dataset must contain exactly the observations not rejected above.
    BOOST_CHECK_EQUAL( activeDataset->getNumberOfObservations( ), dataset->getNumberOfObservations( ) - rejectedRangeIds.size( ) );

    dataset->restoreObservations( ObservationSelectionCondition< double, double >::rejected( ) );

    // Restoration must make all observations active again and recover the original flattened data size.
    BOOST_CHECK_EQUAL( dataset->getObservationIdsMatchingCondition( ObservationSelectionCondition< double, double >::active( ) ).size( ),
                       dataset->getNumberOfObservations( ) );
    BOOST_CHECK_EQUAL( dataset->createEstimationFlattenedObservationData( ).getObservationVector( ).size( ), originalScalarSize );

    const std::vector< double > rangeObsTimes = baseTimeList.at( one_way_range );
    const std::pair< double, double > middleRangeWindow =
            std::make_pair( rangeObsTimes.at( numberOfObservations / 3 ), rangeObsTimes.at( 2 * numberOfObservations / 3 ) );
    const ObservationSelectionCondition< double, double > middleRangeValues =
            ObservationSelectionCondition< double, double >::observableType( one_way_range ) &&
            ObservationSelectionCondition< double, double >::timeBounds( middleRangeWindow.first, middleRangeWindow.second );
    const std::shared_ptr< ObservationDataset< double, double > > middleRangeDataset = dataset->createNewAndKeep( middleRangeValues );

    // A reduced range-only time-window dataset must retain all three range sets with the expected inclusive time count.
    BOOST_CHECK_EQUAL( middleRangeDataset->getNumberOfObservationSets( ), 3 );
    BOOST_CHECK_EQUAL( middleRangeDataset->getNumberOfObservations( ),
                       3 * ( 2 * ( numberOfObservations / 3 ) - numberOfObservations / 3 + 1 ) );
}

/*!
 * Verifies sorted dataset flattened data preserves row-associated data.
 *
 * Test outline: inserts unsorted one-way range observations with weights and
 * residuals, requests time sorting and then confirms the sorted rows retain
 * their original observation, weight and residual values in flattened data order.
 */
BOOST_AUTO_TEST_CASE( test_dataset_flattened_data_weights_residuals_and_ordering )
{
    LinkEnds linkEnds;
    linkEnds[ receiver ] = LinkEndId( "A" );
    linkEnds[ transmitter ] = LinkEndId( "A" );

    std::vector< double > observationTimes = { 90.0, 10.0, 70.0, 30.0, 50.0, 20.0, 80.0, 40.0, 60.0, 100.0 };
    std::vector< Eigen::VectorXd > observations;
    std::vector< Eigen::VectorXd > weights;
    std::vector< Eigen::VectorXd > residuals;
    for( std::size_t i = 0; i < observationTimes.size( ); ++i )
    {
        observations.push_back( ( Eigen::VectorXd( 1 ) << static_cast< double >( i ) ).finished( ) );
        weights.push_back( ( Eigen::VectorXd( 1 ) << 100.0 + static_cast< double >( i ) ).finished( ) );
        residuals.push_back( ( Eigen::VectorXd( 1 ) << -static_cast< double >( i ) ).finished( ) );
    }

    ObservationDataset< double, double > dataset;
    const unsigned int setId = dataset.addObservationSet( one_way_range,
                                                          LinkDefinition( linkEnds ),
                                                          observations,
                                                          observationTimes,
                                                          receiver,
                                                          {},
                                                          nullptr,
                                                          nullptr,
                                                          weights,
                                                          residuals,
                                                          true );

    const std::vector< unsigned int >& observationIds = dataset.getObservationIdsForSet( setId );

    // Sorting on insertion must make row times strictly increasing.
    for( std::size_t i = 1; i < observationIds.size( ); ++i )
    {
        BOOST_CHECK( dataset.getObservationTime( observationIds.at( i ) ) > dataset.getObservationTime( observationIds.at( i - 1 ) ) );
    }

    // Sorted rows must retain the weight and residual values associated with their original observation value.
    for( const unsigned int observationId : observationIds )
    {
        const double originalIndex = dataset.getObservationValue( observationId )( 0 );
        BOOST_CHECK_EQUAL( dataset.getWeightValue( observationId )( 0 ), 100.0 + originalIndex );
        BOOST_CHECK_EQUAL( dataset.getResidualValue( observationId )( 0 ), -originalIndex );
    }

    const FlattenedObservationData< double, double > flattenedData = dataset.createEstimationFlattenedObservationData( );

    // Flattened data vectors for observations, weights and residuals must all cover the same sorted scalar rows.
    BOOST_CHECK_EQUAL( flattenedData.getObservationVector( ).size( ), static_cast< int >( observationTimes.size( ) ) );
    BOOST_CHECK_EQUAL( flattenedData.getWeightVector( ).size( ), static_cast< int >( observationTimes.size( ) ) );
    BOOST_CHECK_EQUAL( flattenedData.getResidualVector( ).size( ), static_cast< int >( observationTimes.size( ) ) );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
