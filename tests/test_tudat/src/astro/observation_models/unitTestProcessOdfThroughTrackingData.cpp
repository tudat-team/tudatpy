/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <limits>
#include <string>
#include <vector>

#include "tudat/basics/testMacros.h"
#include "tudat/io/preProcessOdfFile.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/estimation_setup/createObservationDataset.h"
#include "tudat/simulation/estimation_setup/compressDopplerObservationCollection.h"

using namespace tudat::input_output;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;

using namespace tudat;

BOOST_AUTO_TEST_SUITE( test_process_odf_data )

namespace
{

const std::string syntheticDopplerStationName = "DSS-14";

double getSyntheticDopplerStartUtc( )
{
    return basic_astrodynamics::DateTime( 2000, 1, 1, 12, 0, 0.0 ).epoch< double >( );
}

std::vector< double > convertUtcOffsetsToTdbTimes( const std::vector< double >& utcOffsets )
{
    earth_orientation::TerrestrialTimeScaleConverter timeScaleConverter;
    Eigen::Vector3d stationPosition = simulation_setup::getCombinedApproximateGroundStationPositions( ).at( syntheticDopplerStationName );

    std::vector< double > utcTimes;
    std::vector< Eigen::Vector3d > stationPositions;
    for( double utcOffset : utcOffsets )
    {
        utcTimes.push_back( getSyntheticDopplerStartUtc( ) + utcOffset );
        stationPositions.push_back( stationPosition );
    }
    return timeScaleConverter.getCurrentTimes< double >(
            basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, utcTimes, stationPositions );
}

std::vector< double > convertTdbTimesToUtcOffsets( const std::vector< double >& tdbTimes )
{
    if( tdbTimes.empty( ) )
    {
        return {};
    }

    earth_orientation::TerrestrialTimeScaleConverter timeScaleConverter;
    Eigen::Vector3d stationPosition = simulation_setup::getCombinedApproximateGroundStationPositions( ).at( syntheticDopplerStationName );
    std::vector< double > utcTimes = timeScaleConverter.getCurrentTimesFromSinglePosition< double >(
            basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, tdbTimes, stationPosition );

    std::vector< double > utcOffsets;
    for( double utcTime : utcTimes )
    {
        utcOffsets.push_back( utcTime - getSyntheticDopplerStartUtc( ) );
    }
    return utcOffsets;
}

std::shared_ptr< ObservationDataset< double, double > > createSyntheticDopplerDataset( const std::vector< double >& utcOffsets,
                                                                                       const double integrationTime )
{
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > observations;
    for( double utcOffset : utcOffsets )
    {
        Eigen::Matrix< double, Eigen::Dynamic, 1 > observation( 1 );
        observation( 0 ) = utcOffset;
        observations.push_back( observation );
    }

    LinkDefinition linkEnds;
    linkEnds[ transmitter ] = LinkEndId( "Earth", syntheticDopplerStationName );
    linkEnds[ retransmitter ] = LinkEndId( "SyntheticSpacecraft", "" );
    linkEnds[ receiver ] = LinkEndId( "Earth", syntheticDopplerStationName );

    std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings =
            std::make_shared< ObservationAncillarySimulationSettings >( );
    ancillarySettings->setAncillaryDoubleData( doppler_integration_time, integrationTime );

    std::shared_ptr< ObservationDataset< double, double > > dataset = std::make_shared< ObservationDataset< double, double > >( );
    dataset->addObservationSet( dsn_n_way_averaged_doppler,
                                linkEnds,
                                observations,
                                convertUtcOffsetsToTdbTimes( utcOffsets ),
                                receiver,
                                std::vector< Eigen::VectorXd >( ),
                                nullptr,
                                ancillarySettings );
    return dataset;
}

void checkScalarObservations( const std::shared_ptr< ObservationDataset< double, double > >& dataset,
                              const std::vector< double >& expectedObservations )
{
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > observations = dataset->getObservationsForSet( 0 );
    BOOST_CHECK_EQUAL( observations.size( ), expectedObservations.size( ) );
    for( unsigned int i = 0; i < expectedObservations.size( ); i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( observations.at( i )( 0 ), expectedObservations.at( i ), 1.0E-14 );
    }
}

void checkUtcOffsets( const std::shared_ptr< ObservationDataset< double, double > >& dataset,
                      const std::vector< double >& expectedUtcOffsets )
{
    std::vector< double > utcOffsets = convertTdbTimesToUtcOffsets( dataset->getObservationTimesForSet( 0 ) );
    BOOST_CHECK_EQUAL( utcOffsets.size( ), expectedUtcOffsets.size( ) );
    for( unsigned int i = 0; i < expectedUtcOffsets.size( ); i++ )
    {
        BOOST_CHECK_SMALL( utcOffsets.at( i ) - expectedUtcOffsets.at( i ), 1.0E-6 );
    }
}

}  // namespace

BOOST_AUTO_TEST_CASE( testProcessOdfData )
{
    // Define ODF data paths
    std::string odFile = tudat::paths::getTudatTestDataPath( ) + "/odf07155.odf";

    std::string spacecraftName = "MESSENGER";
    const auto trackingDataAndSupplementaryData = loadOdfFile< long double, Time >( { odFile }, spacecraftName, "Earth" );

    spice_interface::loadStandardSpiceKernels( );
    BodyListSettings bodySettings = getDefaultBodySettings( { "Earth", "Sun" }, "SSB", "J2000" );
    bodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );
    SystemOfBodies bodies = createSystemOfBodies< long double, Time >( bodySettings );
    observation_models::setTrackingSupplementaryDataInBodies( bodies, trackingDataAndSupplementaryData.second );

    // Create data structure that handles Observed Data in Tudat
    std::shared_ptr< observation_models::ObservationDataset< long double, Time > > observedObservationDataset =
            observation_models::createObservationDatasetFromTrackingData< long double, Time >( trackingDataAndSupplementaryData.first,
                                                                                               bodies );

    bool checkedDss14DopplerSet = false;

    // Check the observations for NWayRange and DsnNWayAveragedDoppler
    for( unsigned int setId = 0; setId < observedObservationDataset->getNumberOfObservationSets( ); ++setId )
    {
        observation_models::ObservableType observableType = observedObservationDataset->getObservationSetMetadata( setId ).observableType_;
        const observation_models::LinkEnds& linkEnds =
                observedObservationDataset
                        ->getLinkDefinition( observedObservationDataset->getObservationSetMetadata( setId ).linkDefinitionId_ )
                        .linkEnds_;
        auto observations = observedObservationDataset->getObservationsForSet( setId );
        auto observationTimes = observedObservationDataset->getObservationTimesForSet( setId );
        auto ancillarySettings = observedObservationDataset->getAncillarySettings(
                observedObservationDataset->getObservationSetMetadata( setId ).ancillarySettingsId_ );

        if( observations.empty( ) || observationTimes.empty( ) )
        {
            continue;
        }

        // Check NWayRange observable
        if( observableType == observation_models::n_way_range )
        {
            BOOST_CHECK_EQUAL( observation_models::getObservableName( observableType ), "NWayRange" );
            BOOST_CHECK_EQUAL( linkEnds.at( retransmitter ).bodyName_, "MESSENGER" );
            if( linkEnds.at( transmitter ).getReferencePointName( ) != "DSS-14" ||
                linkEnds.at( receiver ).getReferencePointName( ) != "DSS-14" )
            {
                continue;
            }
            BOOST_CHECK_CLOSE_FRACTION( double( observationTimes.at( 0 ) ), 234262616.184812933, 1e-9 );
            BOOST_CHECK_CLOSE_FRACTION( observations[ 0 ].transpose( )( 0 ), 333589.366953747, 1e-9 );

            if( ancillarySettings != nullptr )
            {
                BOOST_CHECK_EQUAL( ancillarySettings->getDoubleData( ).at( observation_models::reception_reference_frequency_band ),
                                   1.000000000 );
                BOOST_CHECK_CLOSE_FRACTION(
                        ancillarySettings->getDoubleData( ).at( observation_models::sequential_range_lowest_ranging_component ),
                        14.000000000,
                        1e-9 );
                BOOST_CHECK_EQUAL( ancillarySettings->getDoubleVectorData( ).at( observation_models::link_ends_delays ).size( ), 3 );
                BOOST_CHECK_EQUAL( ancillarySettings->getDoubleVectorData( ).at( observation_models::frequency_bands ).size( ), 2 );
            }
        }

        // Check DsnNWayAveragedDoppler observable
        if( observableType == observation_models::dsn_n_way_averaged_doppler )
        {
            BOOST_CHECK_EQUAL( observation_models::getObservableName( observableType ), "DsnNWayAveragedDoppler" );
            BOOST_CHECK_EQUAL( linkEnds.at( retransmitter ).bodyName_, "MESSENGER" );
            if( linkEnds.at( transmitter ).getReferencePointName( ) != "DSS-14" ||
                linkEnds.at( receiver ).getReferencePointName( ) != "DSS-14" )
            {
                continue;
            }
            checkedDss14DopplerSet = true;
            BOOST_CHECK_CLOSE_FRACTION( double( observationTimes.at( 0 ) ), 234262457.184812993, 1e-9 );
            BOOST_CHECK_CLOSE_FRACTION( observations[ 0 ].transpose( )( 0 ), 1.563486099, 1e-9 );

            if( ancillarySettings != nullptr )
            {
                BOOST_CHECK_EQUAL( ancillarySettings->getDoubleData( ).at( observation_models::reception_reference_frequency_band ),
                                   1.000000000 );
                BOOST_CHECK_CLOSE_FRACTION(
                        ancillarySettings->getDoubleData( ).at( observation_models::doppler_integration_time ), 60.000000000, 1e-9 );
                BOOST_CHECK_CLOSE_FRACTION( ancillarySettings->getDoubleData( ).at( observation_models::doppler_reference_frequency ),
                                            7177641534.000000000,
                                            1e-9 );
                BOOST_CHECK_EQUAL( ancillarySettings->getDoubleVectorData( ).at( observation_models::link_ends_delays ).size( ), 3 );
                BOOST_CHECK_EQUAL( ancillarySettings->getDoubleVectorData( ).at( observation_models::frequency_bands ).size( ), 2 );
            }
        }
    }
    BOOST_CHECK( checkedDss14DopplerSet );
}

BOOST_AUTO_TEST_CASE( testCompressDopplerDataUsesCadenceRuns )
{
    std::shared_ptr< ObservationDataset< double, double > > originalObservationDataset =
            createSyntheticDopplerDataset( { 0.0, 10.0, 20.0, 30.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0 }, 10.0 );
    std::shared_ptr< ObservationDataset< double, double > > compressedObservationDataset =
            createCompressedDopplerDataset( originalObservationDataset, 3, 1 );

    checkScalarObservations( compressedObservationDataset, { 10.0, 70.0, 100.0 } );
    checkUtcOffsets( compressedObservationDataset, { 10.0, 70.0, 100.0 } );
    BOOST_CHECK_CLOSE_FRACTION(
            compressedObservationDataset
                    ->getAncillarySettings( compressedObservationDataset->getObservationSetMetadata( 0 ).ancillarySettingsId_ )
                    ->getAncillaryDoubleData( doppler_integration_time ),
            30.0,
            1.0E-14 );

    std::shared_ptr< ObservationDataset< double, double > > gapFreeCompressedObservationDataset =
            createCompressedDopplerDataset( createSyntheticDopplerDataset( { 0.0, 10.0, 20.0, 30.0, 40.0, 50.0 }, 10.0 ), 3, 1 );
    checkScalarObservations( gapFreeCompressedObservationDataset, { 10.0, 40.0 } );
    checkUtcOffsets( gapFreeCompressedObservationDataset, { 10.0, 40.0 } );

    std::shared_ptr< ObservationDataset< double, double > > shortRunCompressedObservationSet =
            compressDopplerData( createSyntheticDopplerDataset( { 0.0, 10.0 }, 10.0 ), 0, 3 );
    BOOST_CHECK_EQUAL( shortRunCompressedObservationSet->getObservationsForSet( 0 ).size( ), 0 );
    BOOST_CHECK_EQUAL( shortRunCompressedObservationSet->getObservationTimesForSet( 0 ).size( ), 0 );
}

BOOST_AUTO_TEST_SUITE_END( )
