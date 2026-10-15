/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>

#include <cmath>

#include "tudat/io/basicInputOutput.h"
#include "tudat/io/readCrdFile.h"
#include "tudat/io/readSp3File.h"

namespace tudat
{
namespace unit_tests
{

//! Read actual paired postfit SP3 and CRDv2 products into precise, correction-aware model inputs.
BOOST_AUTO_TEST_CASE( published_lageos_orbit_and_normal_points )
{
    const std::string directory = paths::getTudatTestDataPath( ) + "slr/";

    const auto orbit = input_output::readSp3File( directory + "ilrsa.orb.lageos1.240406.v81.sp3" );

    const auto& samples = orbit.records.at( "L51" );

    // The provider declares a fitted UTC/SLRF2020 orbit, with explicit SI position and velocity.
    BOOST_CHECK_EQUAL( orbit.orbitType, "FIT" );
    BOOST_CHECK_EQUAL( orbit.timeSystem, "UTC" );
    BOOST_CHECK_EQUAL( orbit.coordinateSystem, "SLR20" );
    BOOST_CHECK_EQUAL( samples.size( ), 5040 );
    BOOST_CHECK_EQUAL( orbit.epochInterval, 120.0 );
    BOOST_CHECK_SMALL( samples.front( ).position.x( ) - ( -5064.275346 * 1000.0 ), 1.0e-9 );
    BOOST_CHECK_SMALL( samples.front( ).velocity.x( ) - ( -13983.449624 * 0.1 ), 1.0e-12 );
    BOOST_CHECK( samples.front( ).positionAvailable && samples.front( ).velocityAvailable );
    BOOST_CHECK( std::isnan( samples.front( ).clockCorrection ) );
    BOOST_CHECK( !samples.front( ).orbitPredicted );
    BOOST_CHECK_EQUAL( ( samples.at( 1 ).epoch - samples.front( ).epoch ).getSeconds< double >( ), 120.0 );

    const auto passes = input_output::readCrdFile( directory + "lageos1_20240402.np2" );

    std::size_t grazPoints = 0;

    for( const auto& pass : passes )
    {
        if( pass.configuration_.cdpPadId_ != 7839 )
        {
            continue;
        }

        const auto& configuration = pass.configuration_;

        const auto measurements = input_output::extractNormalPointMeasurements( pass );

        // These real Graz sessions are UTC(GPS) transmit events: atmosphere/CoM absent, calibration already applied.
        BOOST_CHECK_EQUAL( configuration.formatVersion_, 2 );
        BOOST_CHECK_EQUAL( configuration.stationEpochTimeScale_, 4 );
        BOOST_CHECK_EQUAL( configuration.systemNumber_, 34 );
        BOOST_CHECK_EQUAL( configuration.occupancy_, 2 );
        BOOST_CHECK_EQUAL( configuration.troposphericCorrectionApplied_, 0 );
        BOOST_CHECK_EQUAL( configuration.centerOfMassCorrectionApplied_, 0 );
        BOOST_CHECK_EQUAL( configuration.stationSystemDelayApplied_, 1 );
        BOOST_CHECK_EQUAL( configuration.rangeType_, 2 );
        BOOST_CHECK_EQUAL( configuration.systemConfigurations_.at( "0902" ).transmitWavelengthNm_, 532.0 );
        BOOST_CHECK( !pass.data_.calibrationData_.empty( ) );
        BOOST_CHECK( !pass.data_.meteorologicalData_.empty( ) );
        BOOST_CHECK_EQUAL( measurements.size( ), pass.data_.normalPointData_.size( ) );
        for( std::size_t index = 0; index < measurements.size( ); ++index )
        {
            const auto& record = pass.data_.normalPointData_.at( index );

            // Loaded ranges and their event metadata can enter a model without a second calibration correction.
            BOOST_CHECK_EQUAL( record.epochEvent_, 2 );
            BOOST_CHECK_EQUAL( record.systemConfigurationId_, "0902" );
            BOOST_CHECK_EQUAL( measurements.at( index ).second,
                               input_output::convertCrdTwoWayTimeOfFlightToSlrRange( record.twoWayTimeOfFlight_ ) );
            BOOST_CHECK( measurements.at( index ).first > samples.front( ).epoch );
            BOOST_CHECK( measurements.at( index ).first < samples.back( ).epoch );
        }

        grazPoints += measurements.size( );
    }

    // Assert actual station coverage so a changed fixture or selector cannot make the case empty.
    BOOST_CHECK_EQUAL( grazPoints, 8 );
}

}  // namespace unit_tests
}  // namespace tudat
