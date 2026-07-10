/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      If tabs are used as spaces, it doesn't work. The seperator should also be tabs then.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <cstdio>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>

#include <boost/filesystem.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/simulation/estimation_setup/observationCollection.h"

#include "tudat/io/preProcessFdetsFile.h"
#include "tudat/io/readTrackingTxtFile.h"
#include "tudat/simulation/estimation_setup/processTrackingTxtFile.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"

// Some simplifications for shorter lines
namespace tio = tudat::input_output;
namespace tss = tudat::simulation_setup;
namespace tom = tudat::observation_models;

namespace tudat
{
namespace unit_tests
{

using namespace basic_astrodynamics;
using namespace earth_orientation;
using namespace simulation_setup;

namespace
{

std::string createTempPath( const std::string& suffix )
{
    boost::filesystem::path tempPath =
            boost::filesystem::temp_directory_path( ) / boost::filesystem::unique_path( "tudat-tracking-cadence-%%%%%%" + suffix );
    return tempPath.string( );
}

class CoutRedirect
{
public:
    CoutRedirect( ): originalBuffer_( std::cout.rdbuf( buffer_.rdbuf( ) ) ) {}

    ~CoutRedirect( )
    {
        std::cout.rdbuf( originalBuffer_ );
    }

    std::string getOutput( ) const
    {
        return buffer_.str( );
    }

private:
    std::ostringstream buffer_;
    std::streambuf* originalBuffer_;
};

std::shared_ptr< tio::TrackingTxtFileContents > createSyntheticAveragedDopplerTrackingFile(
        const std::vector< double >& observationSeconds,
        const bool addFileNameMetadata = true,
        const double precomputedCadence = std::numeric_limits< double >::quiet_NaN( ) )
{
    const std::string filePath = createTempPath( ".txt" );
    {
        std::ofstream file( filePath.c_str( ) );
        for( unsigned int i = 0; i < observationSeconds.size( ); i++ )
        {
            file << "2000 1 1 12 0 " << observationSeconds.at( i ) << " " << 1000.0 + static_cast< double >( i ) << "\n";
        }
    }

    std::shared_ptr< tio::TrackingTxtFileContents > trackingFile = tio::createTrackingTxtFileContents(
            filePath, { "year", "month", "day", "hour", "minute", "second", "doppler_averaged_frequency_hz" }, '#', " \t" );
    std::remove( filePath.c_str( ) );

    trackingFile->addMetaData( tio::TrackingDataType::receiving_station_name, "TEST_STATION" );
    trackingFile->addMetaData( tio::TrackingDataType::transmitting_station_name, "TEST_STATION" );
    if( addFileNameMetadata )
    {
        trackingFile->addMetaData( tio::TrackingDataType::file_name, "synthetic_ifms_gap_file.tab" );
    }
    if( std::isfinite( precomputedCadence ) )
    {
        trackingFile->addMetaData( tio::TrackingDataType::doppler_integration_time, precomputedCadence );
    }
    return trackingFile;
}

std::shared_ptr< tom::SingleObservationSet< double, double > > createSyntheticAveragedDopplerObservationSet(
        const std::vector< double >& observationSeconds,
        const bool addFileNameMetadata = true,
        const double precomputedCadence = std::numeric_limits< double >::quiet_NaN( ) )
{
    std::shared_ptr< tio::TrackingTxtFileContents > trackingFile =
            createSyntheticAveragedDopplerTrackingFile( observationSeconds, addFileNameMetadata, precomputedCadence );

    std::map< std::string, Eigen::Vector3d > stationPositions;
    stationPositions[ "TEST_STATION" ] = Eigen::Vector3d::Zero( );
    std::shared_ptr< tom::ProcessedTrackingTxtFileContents< double, double > > processedTrackingFile =
            std::make_shared< tom::ProcessedTrackingTxtFileContents< double, double > >(
                    trackingFile, "SyntheticSpacecraft", stationPositions );

    std::shared_ptr< tom::ObservationCollection< double, double > > observationCollection =
            tom::createTrackingTxtFileObservationCollection< double, double >( processedTrackingFile, { tom::dsn_n_way_averaged_doppler } );
    auto observationSets = observationCollection->getObservationsSets( );
    return observationSets.at( tom::dsn_n_way_averaged_doppler ).begin( )->second.at( 0 );
}

}  // namespace

//! Temporary utility function to print arrays to std::cout
template< typename T >
void printArr( const T& arr )
{
    std::cout << "[ ";
    for( const auto& i : arr )
    {
        std::cout << i << ' ';
    }
    std::cout << "]\n";
}

//! Temporary utility function to print link Ends
std::ostream& operator<<( std::ostream& os, const observation_models::LinkEnds& linkEnds )
{
    os << "LinkEnds{\n";
    for( auto& linkEnd : linkEnds )
    {
        os << "\t{type: " << std::to_string( linkEnd.first ) << ", body: " << linkEnd.second.getBodyName( )
           << ", station: " << linkEnd.second.getReferencePointName( ) << "}\n";
    }
    os << "}";
    return os;
}

//! Utility function to get a single block from a datMap that maps keys to vectors
template< typename K, typename V >
std::map< K, V > extractBlockFromVectorMap( const std::map< K, std::vector< V > >& vectorMap, int blockIndex )
{
    std::map< K, V > singleBlock;
    for( const auto& pair : vectorMap )
    {
        if( blockIndex < 0 )
        {
            blockIndex += pair.second.size( );
        }
        singleBlock[ pair.first ] = pair.second.at( blockIndex );
    }
    return singleBlock;
}

//! Function that specifies a standard format for the Viking file. A user can also do this if they often read the same file format
std::shared_ptr< tio::TrackingTxtFileContents > readVikingRangeFile( const std::string& fileName )
{
    std::vector< std::string > columnTypes( { "spacecraft_id",
                                              "dsn_transmitting_station_nr",
                                              "dsn_receiving_station_nr",
                                              "year",
                                              "month_three_letter",
                                              "day",
                                              "hour",
                                              "minute",
                                              "second",
                                              "round_trip_light_time_microseconds",
                                              "light_time_measurement_accuracy_microseconds" } );
    auto vikingFile = tio::createTrackingTxtFileContents( fileName, columnTypes );
    vikingFile->addMetaData( tio::TrackingDataType::file_name, "Viking lander range data" );
    return vikingFile;
}

//! Function that specifies a standard format for a JUICE Fdets file. A user can also do this if they often read the same file format
std::shared_ptr< tio::TrackingTxtFileContents > readJuiceFdetsFile( const std::string& fileName )
{
    std::vector< std::string > columnTypes( {
            "utc_datetime_string",
            "signal_to_noise_ratio",
            "normalised_spectral_max",
            "doppler_measured_frequency_hz",
            "doppler_noise_hz",
    } );

    auto rawFileContents = tio::createTrackingTxtFileContents( fileName, columnTypes, '#', ", \t" );
    rawFileContents->addMetaData( tio::TrackingDataType::file_name, "JUICE Fdets Test File" );
    return rawFileContents;
}

// Setting some path variables for the test files
const std::string vikingRangePath = tudat::paths::getTudatTestDataPath( ) + "vikingrange.txt";
const std::string marsPathfinderRangePath = tudat::paths::getTudatTestDataPath( ) + "mars-pathfinder-range.txt";
const std::string junoRangePath = tudat::paths::getTudatTestDataPath( ) + "juno_range.txt";
const std::string marinerRangePath = tudat::paths::getTudatTestDataPath( ) + "mariner9obs.txt";
const std::string juiceFdetsDopplerPath = tudat::paths::getTudatTestDataPath( ) + "Fdets.jui2023.04.26.Hb.0006.r2i.txt";
const std::string juiceFdetsDopplerWithScanPath = tudat::paths::getTudatTestDataPath( ) + "Fdets.jui2024.08.20.Yg.r2i.txt";

//! Starting the entire test suite
BOOST_AUTO_TEST_SUITE( test_tracking_txt_file_reader );

//! Test raw reading from file with custom `readVikingRangeFile` function for Viking data
BOOST_AUTO_TEST_CASE( VikingRangeDataCustomFunction )
{
    std::shared_ptr< tio::TrackingTxtFileContents > rawVikingFile = readVikingRangeFile( vikingRangePath );

    std::string spacecraftName = "Viking";

    auto rawDataMap = rawVikingFile->getRawDataMap( );
    auto dataMap = rawVikingFile->getDoubleDataMap( );
    auto dataBlock3 = extractBlockFromVectorMap( dataMap, 3 );

    BOOST_CHECK_EQUAL( dataBlock3[ tio::TrackingDataType::spacecraft_id ], 1 );
    BOOST_CHECK_EQUAL( dataBlock3[ tio::TrackingDataType::dsn_transmitting_station_nr ], 43 );
    BOOST_CHECK_EQUAL( dataBlock3[ tio::TrackingDataType::dsn_receiving_station_nr ], 43 );
    BOOST_CHECK_EQUAL( dataBlock3[ tio::TrackingDataType::year ], 1976 );
    BOOST_CHECK_EQUAL( dataBlock3[ tio::TrackingDataType::month ], 7 );
    BOOST_CHECK_EQUAL( dataBlock3[ tio::TrackingDataType::day ], 22 );
    BOOST_CHECK_EQUAL( dataBlock3[ tio::TrackingDataType::hour ], 6 );
    BOOST_CHECK_EQUAL( dataBlock3[ tio::TrackingDataType::minute ], 2 );
    BOOST_CHECK_EQUAL( dataBlock3[ tio::TrackingDataType::second ], 32 );
    BOOST_CHECK_EQUAL( dataBlock3[ tio::TrackingDataType::n_way_light_time ], 2290.150246895 );

    //! Incidentally, check if the correct number of observations are transferred into an observation collection
    auto observationCollection = observation_models::createTrackingTxtFileObservationCollection< double, double >(
            rawVikingFile, spacecraftName, { tom::n_way_range } );
    BOOST_CHECK_EQUAL( observationCollection->getTotalObservableSize( ), 1258 );
}

//! Test raw reading from file without custom reading function for the MarsPathFinder data
BOOST_AUTO_TEST_CASE( marsPathfinderRangeSimpleReading )
{
    std::vector< std::string > fieldTypeVector{
        "spacecraft_id",
        "dsn_transmitting_station_nr",
        "dsn_receiving_station_nr",
        "year",
        "month_three_letter",
        "day",
        "hour",
        "minute",
        "second",
        "round_trip_light_time_microseconds",
        "light_time_measurement_accuracy_microseconds",
    };

    auto rawTrackingFile = tio::createTrackingTxtFileContents( marsPathfinderRangePath, fieldTypeVector, '#', ",: \t" );
    rawTrackingFile->addMetaData( tio::TrackingDataType::spacecraft_transponder_delay, 420.e-6 );
    rawTrackingFile->addMetaData( tio::TrackingDataType::uplink_frequency, 7.2e9 );
    rawTrackingFile->addMetaData( tio::TrackingDataType::downlink_frequency, 8.4e9 );

    auto dataMap = rawTrackingFile->getDoubleDataMap( );
    auto dataBlock4 = extractBlockFromVectorMap( dataMap, 4 );

    BOOST_CHECK_EQUAL( dataBlock4[ tio::TrackingDataType::spacecraft_id ], 3 );
    BOOST_CHECK_EQUAL( dataBlock4[ tio::TrackingDataType::dsn_transmitting_station_nr ], 65 );
    BOOST_CHECK_EQUAL( dataBlock4[ tio::TrackingDataType::dsn_receiving_station_nr ], 65 );
    BOOST_CHECK_EQUAL( dataBlock4[ tio::TrackingDataType::year ], 1997 );
    BOOST_CHECK_EQUAL( dataBlock4[ tio::TrackingDataType::month ], 7 );
    BOOST_CHECK_EQUAL( dataBlock4[ tio::TrackingDataType::day ], 25 );
    BOOST_CHECK_EQUAL( dataBlock4[ tio::TrackingDataType::hour ], 18 );
    BOOST_CHECK_EQUAL( dataBlock4[ tio::TrackingDataType::minute ], 17 );
    BOOST_CHECK_EQUAL( dataBlock4[ tio::TrackingDataType::second ], 02 );
    BOOST_CHECK_EQUAL( dataBlock4[ tio::TrackingDataType::n_way_light_time ], 1420.476556473 );
    BOOST_CHECK_EQUAL( dataBlock4[ tio::TrackingDataType::light_time_measurement_accuracy ], 6.7e-8 );

    BOOST_CHECK_EQUAL( rawTrackingFile->getNumColumns( ), 11 );
    const auto& metaDataDoubleMap = rawTrackingFile->getMetaDataDoubleMap( );
    BOOST_CHECK_EQUAL( metaDataDoubleMap.at( tio::TrackingDataType::spacecraft_transponder_delay ), 420.e-6 );
}

//! Test simple raw reading with Juno data
BOOST_AUTO_TEST_CASE( junoSimpleReading )
{
    std::vector< std::string > fieldTypeVector{ "spacecraft_id",
                                                "dsn_transmitting_station_nr",
                                                "dsn_receiving_station_nr",
                                                "year",
                                                "month",
                                                "day",
                                                "hour",
                                                "minute",
                                                "second",
                                                "round_trip_light_time_seconds",
                                                "light_time_measurement_delay_microseconds",
                                                "planet_nr",
                                                "tdb_spacecraft_seconds_j2000",
                                                "x_planet_frame_km",
                                                "y_planet_frame_km",
                                                "z_planet_frame_km",
                                                "vx_planet_frame_kms",
                                                "vy_planet_frame_kms",
                                                "vz_planet_frame_kms" };

    auto rawTrackingFile = tio::createTrackingTxtFileContents( junoRangePath, fieldTypeVector, '#', ",: \t" );
    auto dataMap = rawTrackingFile->getDoubleDataMap( );
    auto dataBlock0 = extractBlockFromVectorMap( dataMap, 0 );

    BOOST_CHECK_EQUAL( dataBlock0[ tio::TrackingDataType::spacecraft_id ], 61 );
    BOOST_CHECK_EQUAL( dataBlock0[ tio::TrackingDataType::dsn_transmitting_station_nr ], 55 );
    BOOST_CHECK_EQUAL( dataBlock0[ tio::TrackingDataType::dsn_receiving_station_nr ], 55 );
    BOOST_CHECK_EQUAL( dataBlock0[ tio::TrackingDataType::year ], 2016 );
    BOOST_CHECK_EQUAL( dataBlock0[ tio::TrackingDataType::month ], 8 );
    BOOST_CHECK_EQUAL( dataBlock0[ tio::TrackingDataType::day ], 27 );
    BOOST_CHECK_EQUAL( dataBlock0[ tio::TrackingDataType::hour ], 13 );
    BOOST_CHECK_EQUAL( dataBlock0[ tio::TrackingDataType::minute ], 45 );
    BOOST_CHECK_EQUAL( dataBlock0[ tio::TrackingDataType::second ], 6 );
    BOOST_CHECK_EQUAL( dataBlock0[ tio::TrackingDataType::n_way_light_time ], 6355.0487233317 );
    BOOST_CHECK_EQUAL( dataBlock0[ tio::TrackingDataType::light_time_measurement_accuracy ], 0.0 );
    BOOST_CHECK_EQUAL( dataBlock0[ tio::TrackingDataType::planet_nr ], 5 );
    BOOST_CHECK_CLOSE( dataBlock0[ tio::TrackingDataType::tdb_spacecraft_j2000 ], 525574396.542800, 1e-4 );
    BOOST_CHECK_CLOSE( dataBlock0[ tio::TrackingDataType::x_planet_frame ], 976985.733, 1e-4 );
    BOOST_CHECK_CLOSE( dataBlock0[ tio::TrackingDataType::y_planet_frame ], 68435520.227, 1e-4 );
    BOOST_CHECK_CLOSE( dataBlock0[ tio::TrackingDataType::z_planet_frame ], 32772692.214, 1e-4 );
    BOOST_CHECK_CLOSE( dataBlock0[ tio::TrackingDataType::vx_planet_frame ], 0727.110, 1e-4 );
    BOOST_CHECK_CLOSE( dataBlock0[ tio::TrackingDataType::vy_planet_frame ], 26571.899, 1e-4 );
    BOOST_CHECK_CLOSE( dataBlock0[ tio::TrackingDataType::vz_planet_frame ], -51299.726, 1e-4 );

    BOOST_CHECK_EQUAL( rawTrackingFile->getNumColumns( ), 19 );
}

//! Test raw reading and adding metadata with Mariner data
BOOST_AUTO_TEST_CASE( marinerSimpleReading )
{
    std::vector< std::string > fieldTypeVector{ "year",
                                                "month",
                                                "day",
                                                "hour",
                                                "minute",
                                                "second",
                                                "round_trip_light_time_microseconds",
                                                "light_time_measurement_accuracy_microseconds",
                                                "residual_de405_microseconds" };

    auto rawTrackingFile = tio::createTrackingTxtFileContents( marinerRangePath, fieldTypeVector, '#', ",: \t" );
    rawTrackingFile->addMetaData( input_output::TrackingDataType::observation_body, "Earth" );
    rawTrackingFile->addMetaData( input_output::TrackingDataType::observed_body, "Mars" );

    auto dataMap = rawTrackingFile->getDoubleDataMap( );
    auto dataBlockLast = extractBlockFromVectorMap( dataMap, -1 );

    //  1972 10 12 00:06:02  2610383946.989   0.475  -0.226
    BOOST_CHECK_EQUAL( dataBlockLast[ tio::TrackingDataType::year ], 1972 );
    BOOST_CHECK_EQUAL( dataBlockLast[ tio::TrackingDataType::month ], 10 );
    BOOST_CHECK_EQUAL( dataBlockLast[ tio::TrackingDataType::day ], 12 );
    BOOST_CHECK_EQUAL( dataBlockLast[ tio::TrackingDataType::hour ], 0 );
    BOOST_CHECK_EQUAL( dataBlockLast[ tio::TrackingDataType::minute ], 6 );
    BOOST_CHECK_EQUAL( dataBlockLast[ tio::TrackingDataType::second ], 2 );
    BOOST_CHECK_EQUAL( dataBlockLast[ tio::TrackingDataType::n_way_light_time ], 2610.383946989 );
    BOOST_CHECK_CLOSE( dataBlockLast[ tio::TrackingDataType::light_time_measurement_accuracy ], 0.475e-6, 1e-10 );
    BOOST_CHECK_CLOSE( dataBlockLast[ tio::TrackingDataType::residual_de405 ], -0.226e-6, 1e-10 );

    BOOST_CHECK_EQUAL( rawTrackingFile->getNumColumns( ), 9 );

    const auto& metaDataStrMap = rawTrackingFile->getMetaDataStrMap( );
    BOOST_CHECK_EQUAL( metaDataStrMap.at( input_output::TrackingDataType::observation_body ), "Earth" );
    BOOST_CHECK_EQUAL( metaDataStrMap.at( input_output::TrackingDataType::observed_body ), "Mars" );
}

//! Test observation collection and time conversions with Viking Data
BOOST_AUTO_TEST_CASE( TestVikingRangeDataObservationCollection )
{
    // Load the observations from the Viking file
    std::shared_ptr< tio::TrackingTxtFileContents > rawVikingFile = readVikingRangeFile( vikingRangePath );
    auto observationCollection = observation_models::createTrackingTxtFileObservationCollection< double, double >(
            rawVikingFile, "Viking", { tom::n_way_range } );

    // Check size of observations
    BOOST_CHECK_EQUAL( observationCollection->getTotalObservableSize( ), 1258 );

    // Check if n-way-range is present in the collection
    const auto& observationTypeStartSize = observationCollection->getObservationTypeStartAndSize( );
    BOOST_CHECK( observationTypeStartSize.find( tom::n_way_range ) != observationTypeStartSize.end( ) );

    //  Checking the first element with station 63 - 63 in the Viking file
    const tom::LinkDefinition linkDefDsn63( {
            { tom::transmitter, tom::LinkEndId( "Earth", "DSS-63" ) },
            { tom::reflector, tom::LinkEndId( "Viking", "" ) },
            { tom::receiver, tom::LinkEndId( "Earth", "DSS-63" ) },
    } );

    auto observationsAndTimesDsn63 = observationCollection->getSingleLinkObservationsAndTimes( tom::n_way_range, linkDefDsn63 );
    auto observationsDsn63 = observationsAndTimesDsn63.first;
    auto timesDsn63 = observationsAndTimesDsn63.second;

    BOOST_CHECK_CLOSE( observationsDsn63( 0, 0 ),
                       2371564782.809 * 1.0E-6 * physical_constants::SPEED_OF_LIGHT,
                       1e-12 );  // Round trip light time = 2371.564782809 seconds (don't relax from 1e-12)
    DateTime utcObservationTime = DateTime( 1976, 8, 8, 18, 4, 2.0 );
    double testTime = createDefaultTimeConverter( )->getCurrentTime< Time >(
            utc_scale, tdb_scale, utcObservationTime.epoch< Time >( ), getApproximateDsnGroundStationPositions( ).at( "DSS-63" ) );
    BOOST_CHECK_CLOSE( timesDsn63[ 0 ], testTime, ( 10.0 * std::numeric_limits< double >::epsilon( ) ) );  // "1976-08-08T18:04:02.000"
}

//! Test with JUICE data. Using doppler measured frequency observable
BOOST_AUTO_TEST_CASE( TestJuiceFile )
{
    const double dopplerBaseFrequency = 8420.0e6;
    const std::string transmittingStationName = "NNO";
    const std::string receivingStationName = "HOBART12";

    spice_interface::loadStandardSpiceKernels( );

    std::shared_ptr< tio::TrackingTxtFileContents > rawFdetsDopplerFile = readJuiceFdetsFile( juiceFdetsDopplerPath );
    rawFdetsDopplerFile->addMetaData( tio::TrackingDataType::doppler_base_frequency, dopplerBaseFrequency );
    rawFdetsDopplerFile->addMetaData( tio::TrackingDataType::doppler_bandwidth, 2.0e3 );
    rawFdetsDopplerFile->addMetaData( tio::TrackingDataType::receiving_station_name, receivingStationName );
    rawFdetsDopplerFile->addMetaData( tio::TrackingDataType::transmitting_station_name, transmittingStationName );

    // CHECK THE RAW FILE
    BOOST_CHECK_EQUAL( rawFdetsDopplerFile->getNumColumns( ), 5 );

    auto dataMap = rawFdetsDopplerFile->getDoubleDataMap( );
    auto dataBlockLast = extractBlockFromVectorMap( dataMap, -1 );

    auto timeDataBlockLast = dataBlockLast[ tio::TrackingDataType::utc_reception_time_j2000 ];

    DateTime utcObservationTime = DateTime( 2023, 4, 25, 9, 46, 10.0 );

    BOOST_CHECK_CLOSE_FRACTION(
            timeDataBlockLast, utcObservationTime.epoch< double >( ), 10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_EQUAL( dataBlockLast[ tio::TrackingDataType::signal_to_noise ], 6.766652540970647242e+05 );
    BOOST_CHECK_EQUAL( dataBlockLast[ tio::TrackingDataType::spectral_max ], 5.754946258897545704e+03 );
    BOOST_CHECK_EQUAL( dataBlockLast[ tio::TrackingDataType::doppler_measured_frequency ], 5977954.253958693705 );
    BOOST_CHECK_EQUAL( dataBlockLast[ tio::TrackingDataType::doppler_noise ], -1.3803018149815216e-02 );

    // Check the VLBI Station name
    const auto& metaDataStrMap = rawFdetsDopplerFile->getMetaDataStrMap( );
    BOOST_CHECK_EQUAL( metaDataStrMap.at( tio::TrackingDataType::transmitting_station_name ), transmittingStationName );
    BOOST_CHECK_EQUAL( metaDataStrMap.at( tio::TrackingDataType::receiving_station_name ), receivingStationName );

    std::vector< std::shared_ptr< data::TrackingData< double, Time > > > trackingData =
            tio::readFdetsFiles< double, Time >( std::vector< std::string >( { juiceFdetsDopplerPath } ),
                                                 std::vector< double >( { dopplerBaseFrequency } ),
                                                 tio::FdetDateFormat::datetime_string,
                                                 "JUICE",
                                                 std::vector< std::string >( { transmittingStationName } ),
                                                 std::vector< std::string >( { receivingStationName } ) )
                    .first;

    BOOST_REQUIRE_EQUAL( trackingData.size( ), 1 );
    BOOST_CHECK_EQUAL( trackingData.at( 0 )->getObservableType( ), "DopplerMeasuredFrequency" );
    BOOST_CHECK_EQUAL( trackingData.at( 0 )->getReferenceLinkEnd( ), "receiver" );
    BOOST_CHECK_EQUAL( trackingData.at( 0 )->getTimeScale( ), "UTC" );
    BOOST_CHECK_EQUAL( trackingData.at( 0 )->getNumberOfObservations( ), 120 );

    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > observations = trackingData.at( 0 )->getObservations( );
    std::vector< Time > epochs = trackingData.at( 0 )->getObservationEpochs( );

    BOOST_CHECK_EQUAL( observations.at( 0 )( 0 ), dopplerBaseFrequency + 5978760.982806123793 );
    BOOST_CHECK_EQUAL( observations.at( 1 )( 0 ), dopplerBaseFrequency + 5978754.318319843151 );
    BOOST_CHECK_EQUAL( observations.at( 2 )( 0 ), dopplerBaseFrequency + 5978747.672510409728 );
    BOOST_CHECK_EQUAL( observations.at( observations.size( ) - 1 )( 0 ), dopplerBaseFrequency + 5977954.253958693705 );

    BOOST_CHECK_CLOSE_FRACTION( static_cast< double >( utcObservationTime.epoch< Time >( ) ),
                                static_cast< double >( epochs.at( epochs.size( ) - 1 ) ),
                                10.0 * std::numeric_limits< double >::epsilon( ) );

    std::shared_ptr< tom::ObservationCollection< double, Time > > observationCollection =
            tom::createObservationCollection< double, Time >( trackingData );
    std::vector< Time > observationCollectionEpochs = observationCollection->getConcatenatedTimeVector( );
    Time expectedTdbObservationTime = TerrestrialTimeScaleConverter( ).getCurrentTime< Time >(
            utc_scale, tdb_scale, utcObservationTime.epoch< Time >( ), Eigen::Vector3d::Zero( ) );
    BOOST_REQUIRE_EQUAL( observationCollectionEpochs.size( ), epochs.size( ) );
    BOOST_CHECK_SMALL( static_cast< double >( observationCollectionEpochs.back( ) - expectedTdbObservationTime ), 1.0E-12 );
}

BOOST_AUTO_TEST_CASE( TestFdetsFileReaderDateFormatAndScanDetection )
{
    std::shared_ptr< tio::TrackingTxtFileContents > rawFdetsDopplerFile = tio::readFdetsFile( juiceFdetsDopplerPath );
    BOOST_CHECK_EQUAL( rawFdetsDopplerFile->getNumColumns( ), 5 );
    BOOST_CHECK_EQUAL( rawFdetsDopplerFile->getRawColumnTypes( ).at( 0 ), "utc_datetime_string" );

    std::shared_ptr< tio::TrackingTxtFileContents > rawFdetsDopplerFileWithScan = tio::readFdetsFile( juiceFdetsDopplerWithScanPath );

    BOOST_CHECK_EQUAL( rawFdetsDopplerFileWithScan->getNumColumns( ), 6 );
    BOOST_CHECK_EQUAL( rawFdetsDopplerFileWithScan->getRawColumnTypes( ).at( 0 ), "scan_number" );
    BOOST_CHECK_EQUAL( rawFdetsDopplerFileWithScan->getRawColumnTypes( ).at( 1 ), "utc_datetime_string" );

    auto dataMap = rawFdetsDopplerFileWithScan->getDoubleDataMap( );
    auto dataBlockFirst = extractBlockFromVectorMap( dataMap, 0 );

    DateTime utcObservationTime = DateTime( 2024, 8, 20, 17, 29, 51.5 );
    BOOST_CHECK_CLOSE_FRACTION( dataBlockFirst[ tio::TrackingDataType::utc_reception_time_j2000 ],
                                utcObservationTime.epoch< double >( ),
                                10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_EQUAL( dataBlockFirst[ tio::TrackingDataType::scan_nr ], 1 );
    BOOST_CHECK_EQUAL( dataBlockFirst[ tio::TrackingDataType::signal_to_noise ], 2.571405547427670390e+05 );
    BOOST_CHECK_EQUAL( dataBlockFirst[ tio::TrackingDataType::spectral_max ], 6.072471290268301800e+02 );
    BOOST_CHECK_EQUAL( dataBlockFirst[ tio::TrackingDataType::doppler_measured_frequency ], 13682699.425314944237 );
    BOOST_CHECK_EQUAL( dataBlockFirst[ tio::TrackingDataType::doppler_noise ], 5.1043355817910196e-03 );

    BOOST_CHECK_THROW( tio::readFdetsFile( juiceFdetsDopplerPath, tio::FdetDateFormat::pair_of_numbers ), std::runtime_error );
}

//! Test averaged Doppler cadence inference when filtered rows leave middle-of-file gaps
BOOST_AUTO_TEST_CASE( TestAveragedDopplerCadenceGaps )
{
    std::shared_ptr< tom::SingleObservationSet< double, double > > gapObservationSet;
    std::string gapWarning;
    {
        CoutRedirect outputRedirect;
        gapObservationSet = createSyntheticAveragedDopplerObservationSet( { 0.0, 10.0, 20.0, 50.0 }, true, 10.0 );
        gapWarning = outputRedirect.getOutput( );
    }

    BOOST_CHECK_CLOSE_FRACTION(
            gapObservationSet->getAncillarySettings( )->getAncillaryDoubleData( tom::doppler_integration_time ), 10.0, 1.0E-14 );
    BOOST_CHECK( gapWarning.find( "synthetic_ifms_gap_file.tab" ) != std::string::npos );
    BOOST_CHECK( gapWarning.find( "found 1 cadence gap" ) != std::string::npos );
    BOOST_CHECK( gapWarning.find( "nominal cadence 10" ) != std::string::npos );
    BOOST_CHECK( gapWarning.find( "index 3" ) != std::string::npos );
    BOOST_CHECK( gapWarning.find( "observed delta 30" ) != std::string::npos );

    std::shared_ptr< tom::SingleObservationSet< double, double > > gapFreeObservationSet;
    std::string gapFreeWarning;
    {
        CoutRedirect outputRedirect;
        gapFreeObservationSet = createSyntheticAveragedDopplerObservationSet( { 0.0, 10.0, 20.0, 30.0 } );
        gapFreeWarning = outputRedirect.getOutput( );
    }
    BOOST_CHECK_CLOSE_FRACTION(
            gapFreeObservationSet->getAncillarySettings( )->getAncillaryDoubleData( tom::doppler_integration_time ), 10.0, 1.0E-14 );
    BOOST_CHECK( gapFreeWarning.find( "cadence gap" ) == std::string::npos );

    std::string unknownFileWarning;
    {
        CoutRedirect outputRedirect;
        createSyntheticAveragedDopplerObservationSet( { 0.0, 10.0, 30.0 }, false, 10.0 );
        unknownFileWarning = outputRedirect.getOutput( );
    }
    BOOST_CHECK( unknownFileWarning.find( "unknown tracking file" ) != std::string::npos );

    BOOST_CHECK_THROW( createSyntheticAveragedDopplerObservationSet( { 0.0, 10.0, 20.0, 50.0 } ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( TestIfmsCadenceInferredBeforeFiltering )
{
    const std::string filePath = createTempPath( ".TAB" );
    {
        std::ofstream file( filePath.c_str( ) );
        file << "0 2000-01-01T12:00:00.000 1 0.0 1.0 2000-01-01T12:00:00.000 1000.0 0.0 2000.0 2000.0 0.0 1.0\n";
        file << "1 2000-01-01T12:00:10.000 1 10.0 1.0 2000-01-01T12:00:10.000 1000.0 0.0 2001.0 2001.0 0.0 1.0\n";
        file << "2 2000-01-01T12:00:20.000 1 20.0 1.0 2000-01-01T12:00:20.000 1000.0 0.0 -999.999 2002.0 0.0 1.0\n";
        file << "3 2000-01-01T12:00:30.000 1 30.0 1.0 2000-01-01T12:00:30.000 1000.0 0.0 2003.0 2003.0 0.0 1.0\n";
    }

    std::shared_ptr< tio::TrackingTxtFileContents > filteredIfmsFile = tio::readIfmsFile( filePath, false, true );
    std::remove( filePath.c_str( ) );

    BOOST_CHECK_EQUAL( filteredIfmsFile->getNumRows( ), 3 );
    BOOST_CHECK_CLOSE_FRACTION(
            filteredIfmsFile->getMetaDataDoubleMap( ).at( tio::TrackingDataType::doppler_integration_time ), 10.0, 1.0E-14 );

    std::shared_ptr< tom::SingleObservationSet< double, double > > gapObservationSet;
    std::string gapWarning;
    {
        CoutRedirect outputRedirect;
        gapObservationSet = createSyntheticAveragedDopplerObservationSet( { 0.0, 20.0, 40.0 }, true, 10.0 );
        gapWarning = outputRedirect.getOutput( );
    }

    BOOST_CHECK_CLOSE_FRACTION(
            gapObservationSet->getAncillarySettings( )->getAncillaryDoubleData( tom::doppler_integration_time ), 10.0, 1.0E-14 );
    BOOST_CHECK( gapWarning.find( "found 2 cadence gap" ) != std::string::npos );
    BOOST_CHECK( gapWarning.find( "nominal cadence 10" ) != std::string::npos );
}

//! Test reading of ground station locations
// FIXME-DOPPLER: This might need to be moved to another file
BOOST_AUTO_TEST_CASE( GroundStationLocations )
{
    const static std::string pysctrackGroundStationPosFile = tudat::paths::getStationLocationDataPath( ) + "/glo.sit";
    const static std::string pysctrackGroundStationVelFile = tudat::paths::getStationLocationDataPath( ) + "/glo.vel";
    const static std::string pysctrackGroundStationCodesFile = tudat::paths::getStationLocationDataPath( ) + "/ns_codes.dat";

    std::map< std::string, Eigen::Vector3d > stationPosMap =
            utilities::getMapFromFile< std::string, Eigen::Vector3d >( pysctrackGroundStationPosFile, '$', " \t" );
    std::map< std::string, Eigen::Vector3d > stationVelMap =
            utilities::getMapFromFile< std::string, Eigen::Vector3d >( pysctrackGroundStationVelFile, '$', " \t" );
    std::map< std::string, std::string > stationCodeMap =
            utilities::getMapFromFile< std::string, std::string >( pysctrackGroundStationCodesFile, '*', " \t" );
    //  std::map<std::string, Eigen::Vector3d> stationPosMap = simulation_setup::getApproximateGroundStationPositionsFromFile();
    //  std::map<std::string, std::string> stationCodeMap = simulation_setup::getGroundStationCodesFromFile();

    std::map< std::string, Eigen::Vector3d > stationPositionMap = getVlbiStationPositions( );
    std::cout << stationPositionMap.at( "HOBART12" );
    std::map< std::string, Eigen::Vector3d > stationPositionMap2 = getVlbiStationPositions( );
    std::cout << stationPositionMap2.at( "HOBART12" );
    //  Eigen::Vector3d hobartPos = stationPosMap.at("HOBART12");
    //  Eigen::Vector3d hobartVel = stationVelMap.at("HOBART12");
    //  Eigen::Vector3d hobartPos = simulation_setup::getApproximateGroundStationPositionFromFile("HOBART12");
    //  Eigen::Vector3d hobartVel = simulation_setup::getApproximateGroundStationVelocityFromFile("HOBART12");

    //  BOOST_TEST(hobartPos == stationPosMap["HOBART12"], boost::test_tools::per_element());
    //  BOOST_TESTfirstBlockMinMaxResidual(hobartPos == Eigen::Vector3d(-3949990.106, 2522421.118, -4311708.734),
    //  boost::test_tools::per_element()); BOOST_TEST(hobartVel == Eigen::Vector3d(-39.90, 9.00, 39.90), boost::test_tools::per_element());
}

// End test suite
BOOST_AUTO_TEST_SUITE_END( );

}  // namespace unit_tests

}  // namespace tudat
