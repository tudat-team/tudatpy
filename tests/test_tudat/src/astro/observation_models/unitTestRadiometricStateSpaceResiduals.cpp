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

#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/observationAncillarySettings.h"
#include "tudat/astro/observation_models/observationFrequencies.h"
#include "tudat/astro/observation_models/radiometricResidualStateSpaceConversion.h"
#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation_setup/observationCollection.h"
#include "tudat/simulation/estimation_setup/singleObservationSet.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::observation_models;

namespace
{

LinkDefinition getNWayLinkDefinition( const std::string& stationName )
{
    LinkDefinition linkEnds;
    linkEnds[ transmitter ] = LinkEndId( "Earth", stationName );
    linkEnds[ retransmitter ] = LinkEndId( "Spacecraft", "" );
    linkEnds[ receiver ] = LinkEndId( "Earth", stationName );
    return linkEnds;
}

LinkDefinition getOneWayLinkDefinition( const std::string& stationName )
{
    LinkDefinition linkEnds;
    linkEnds[ transmitter ] = LinkEndId( "Spacecraft", "" );
    linkEnds[ receiver ] = LinkEndId( "Earth", stationName );
    return linkEnds;
}

std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > makeScalarSeries( const std::vector< double >& values )
{
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > series;
    for( double value : values )
    {
        Eigen::Matrix< double, Eigen::Dynamic, 1 > observation( 1 );
        observation( 0 ) = value;
        series.push_back( observation );
    }
    return series;
}

std::shared_ptr< SingleObservationSet< double, double > > makeObservationSet(
        const ObservableType observableType,
        const LinkDefinition& linkEnds,
        const std::vector< double >& observations,
        const std::vector< double >& residuals,
        const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings,
        const std::vector< double >& times = std::vector< double >( ) )
{
    std::vector< double > observationTimes = times;
    if( observationTimes.empty( ) )
    {
        observationTimes.resize( observations.size( ) );
        for( unsigned int i = 0; i < observationTimes.size( ); i++ )
        {
            observationTimes[ i ] = static_cast< double >( i + 1 );
        }
    }

    return std::make_shared< SingleObservationSet< double, double > >( observableType,
                                                                       linkEnds,
                                                                       makeScalarSeries( observations ),
                                                                       observationTimes,
                                                                       receiver,
                                                                       std::vector< Eigen::VectorXd >( ),
                                                                       nullptr,
                                                                       ancillarySettings,
                                                                       std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
                                                                       makeScalarSeries( residuals ) );
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_radiometric_state_space_residuals )

BOOST_AUTO_TEST_CASE( testDsnNWayRangeResidualConversion )
{
    const double rangeConversionFactor = 1.5e-1;
    auto ancillarySettings = std::make_shared< ObservationAncillarySimulationSettings >( );
    ancillarySettings->setAncillaryDoubleData( range_conversion_factor, rangeConversionFactor );

    const std::vector< double > residuals = { 2.0, -4.0, 0.5 };
    auto observationSet =
            makeObservationSet( dsn_n_way_range, getNWayLinkDefinition( "DSS-14" ), { 10.0, 11.0, 12.0 }, residuals, ancillarySettings );

    BOOST_CHECK_CLOSE_FRACTION( observationSet->getResidualStateSpaceConversionFactor( ), rangeConversionFactor, 1.0e-15 );

    const auto convertedResiduals = observationSet->getResidualsInStateSpace( );
    const Eigen::VectorXd concatenatedConvertedResiduals = observationSet->getResidualsInStateSpaceVector( );
    BOOST_CHECK_EQUAL( convertedResiduals.size( ), residuals.size( ) );
    BOOST_CHECK_EQUAL( concatenatedConvertedResiduals.size( ), static_cast< int >( residuals.size( ) ) );
    for( unsigned int i = 0; i < residuals.size( ); i++ )
    {
        const double expectedResidual = residuals.at( i ) * rangeConversionFactor;
        BOOST_CHECK_CLOSE_FRACTION( convertedResiduals.at( i )( 0 ), expectedResidual, 1.0e-15 );
        BOOST_CHECK_CLOSE_FRACTION( concatenatedConvertedResiduals( i ), expectedResidual, 1.0e-15 );
    }
}

BOOST_AUTO_TEST_CASE( testDsnNWayRangeMissingConversionFactor )
{
    auto ancillarySettings = std::make_shared< ObservationAncillarySimulationSettings >( );
    auto observationSet = makeObservationSet( dsn_n_way_range, getNWayLinkDefinition( "DSS-14" ), { 10.0 }, { 1.0 }, ancillarySettings );
    BOOST_CHECK_THROW( observationSet->getResidualStateSpaceConversionFactor( ), std::runtime_error );
    BOOST_CHECK_THROW( ( computeResidualStateSpaceConversionFactor( dsn_n_way_range, ancillarySettings ) ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( testDsnNWayAveragedDopplerResidualConversion )
{
    const double referenceFrequency = 7.2e9;
    auto ancillarySettings = getDsnNWayAveragedDopplerAncillarySettings( { x_band, x_band }, x_band, referenceFrequency, 60.0 );

    const double turnaroundRatio = getDsnDefaultTurnaroundRatios( x_band, x_band );
    const double expectedFactor = physical_constants::SPEED_OF_LIGHT / ( 2.0 * turnaroundRatio * referenceFrequency );

    BOOST_CHECK_CLOSE_FRACTION(
            ( computeResidualStateSpaceConversionFactor( dsn_n_way_averaged_doppler, ancillarySettings ) ), expectedFactor, 1.0e-14 );

    const std::vector< double > residuals = { 0.25, -1.5 };
    auto observationSet = makeObservationSet(
            dsn_n_way_averaged_doppler, getNWayLinkDefinition( "DSS-14" ), { 8.4e9, 8.4e9 }, residuals, ancillarySettings );

    BOOST_CHECK_CLOSE_FRACTION( observationSet->getResidualStateSpaceConversionFactor( ), expectedFactor, 1.0e-14 );
    const Eigen::VectorXd convertedResiduals = observationSet->getResidualsInStateSpaceVector( );
    for( unsigned int i = 0; i < residuals.size( ); i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( convertedResiduals( i ), residuals.at( i ) * expectedFactor, 1.0e-14 );
    }
}

BOOST_AUTO_TEST_CASE( testStoredDopplerConversionFactorIsPreferred )
{
    const double referenceFrequency = 7.2e9;
    const double storedFactor = 0.123;
    auto ancillarySettings = getDsnNWayAveragedDopplerAncillarySettings( { x_band, x_band }, x_band, referenceFrequency, 60.0 );
    ancillarySettings->setAncillaryDoubleData( doppler_conversion_factor, storedFactor );

    BOOST_CHECK_CLOSE_FRACTION(
            ( computeDopplerResidualStateSpaceConversionFactor( dsn_n_way_averaged_doppler, ancillarySettings ) ), storedFactor, 1.0e-15 );
}

BOOST_AUTO_TEST_CASE( testOneWayMeasuredFrequencyResidualConversion )
{
    const double referenceCarrierFrequency = 8.4e9;
    auto ancillarySettings = std::make_shared< ObservationAncillarySimulationSettings >( );
    ancillarySettings->setAncillaryDoubleData( doppler_reference_frequency, referenceCarrierFrequency );

    const double expectedFactor = physical_constants::SPEED_OF_LIGHT / referenceCarrierFrequency;
    BOOST_CHECK_CLOSE_FRACTION( ( computeResidualStateSpaceConversionFactor( one_way_doppler_measured_frequency, ancillarySettings ) ),
                                expectedFactor,
                                1.0e-14 );

    auto observationSet = makeObservationSet(
            one_way_doppler_measured_frequency, getOneWayLinkDefinition( "DSS-14" ), { 8.4e9 }, { -3.0 }, ancillarySettings );
    BOOST_CHECK_CLOSE_FRACTION( observationSet->getResidualsInStateSpaceVector( )( 0 ), -3.0 * expectedFactor, 1.0e-14 );
}

BOOST_AUTO_TEST_CASE( testMeasuredFrequencySetLevelConversion )
{
    const double referenceFrequency = 7.2e9;
    auto ancillarySettings = getDsnNWayAveragedDopplerAncillarySettings( { x_band, x_band }, x_band, referenceFrequency, 60.0 );
    const double expectedFactor =
            physical_constants::SPEED_OF_LIGHT / ( 2.0 * getDsnDefaultTurnaroundRatios( x_band, x_band ) * referenceFrequency );

    BOOST_CHECK_CLOSE_FRACTION(
            ( computeResidualStateSpaceConversionFactor( doppler_measured_frequency, ancillarySettings ) ), expectedFactor, 1.0e-14 );
}

BOOST_AUTO_TEST_CASE( testUnsupportedObservableTypesAreRejected )
{
    auto emptyAncillarySettings = std::make_shared< ObservationAncillarySimulationSettings >( );
    auto oneWayDopplerSet =
            makeObservationSet( one_way_doppler, getOneWayLinkDefinition( "DSS-14" ), { 0.1 }, { 0.2 }, emptyAncillarySettings );
    auto twoWayDopplerSet =
            makeObservationSet( two_way_doppler, getNWayLinkDefinition( "DSS-14" ), { 0.1 }, { 0.2 }, emptyAncillarySettings );
    auto oneWayRangeSet =
            makeObservationSet( one_way_range, getOneWayLinkDefinition( "DSS-14" ), { 1.0e6 }, { 3.0 }, emptyAncillarySettings );

    BOOST_CHECK_THROW( oneWayDopplerSet->getResidualStateSpaceConversionFactor( ), std::runtime_error );
    BOOST_CHECK_THROW( oneWayDopplerSet->getResidualsInStateSpace( ), std::runtime_error );
    BOOST_CHECK_THROW( twoWayDopplerSet->getResidualsInStateSpaceVector( ), std::runtime_error );
    BOOST_CHECK_THROW( oneWayRangeSet->getResidualStateSpaceConversionFactor( ), std::runtime_error );
    BOOST_CHECK_THROW( ( computeResidualStateSpaceConversionFactor( one_way_doppler, emptyAncillarySettings ) ), std::runtime_error );
    BOOST_CHECK_THROW( ( computeResidualStateSpaceConversionFactor( two_way_doppler, emptyAncillarySettings ) ), std::runtime_error );
    BOOST_CHECK_THROW( ( computeResidualStateSpaceConversionFactor( dsn_one_way_averaged_doppler, emptyAncillarySettings ) ),
                       std::runtime_error );
}

BOOST_AUTO_TEST_CASE( testCollectionParserOrderAndUnsupportedRejection )
{
    auto rangeAncillarySettings = std::make_shared< ObservationAncillarySimulationSettings >( );
    rangeAncillarySettings->setAncillaryDoubleData( range_conversion_factor, 2.0 );
    auto dopplerAncillarySettings = getDsnNWayAveragedDopplerAncillarySettings( { x_band, x_band }, x_band, 7.2e9, 60.0 );
    auto emptyAncillarySettings = std::make_shared< ObservationAncillarySimulationSettings >( );

    auto rangeSet = makeObservationSet(
            dsn_n_way_range, getNWayLinkDefinition( "DSS-14" ), { 1.0, 2.0 }, { 4.0, 5.0 }, rangeAncillarySettings, { 10.0, 11.0 } );
    auto dopplerSet = makeObservationSet(
            dsn_n_way_averaged_doppler, getNWayLinkDefinition( "DSS-15" ), { 8.4e9 }, { -0.5 }, dopplerAncillarySettings, { 20.0 } );
    auto unsupportedSet =
            makeObservationSet( one_way_doppler, getOneWayLinkDefinition( "DSS-16" ), { 0.1 }, { 0.3 }, emptyAncillarySettings, { 30.0 } );

    std::vector< std::shared_ptr< SingleObservationSet< double, double > > > supportedSets = { rangeSet, dopplerSet };
    ObservationCollection< double, double > supportedCollection( supportedSets );
    const Eigen::VectorXd nativeResiduals = supportedCollection.getConcatenatedResiduals( );
    const Eigen::VectorXd convertedResiduals = supportedCollection.getConcatenatedResidualsInStateSpace( );
    const auto convertedResidualBlocks = supportedCollection.getResidualsInStateSpace( );
    const auto nativeResidualBlocks = supportedCollection.getResiduals( );

    BOOST_CHECK_EQUAL( convertedResidualBlocks.size( ), nativeResidualBlocks.size( ) );
    BOOST_CHECK_EQUAL( convertedResiduals.size( ), nativeResiduals.size( ) );

    Eigen::VectorXd expectedConvertedResiduals = Eigen::VectorXd::Zero( nativeResiduals.size( ) );
    unsigned int currentIndex = 0;
    for( unsigned int i = 0; i < nativeResidualBlocks.size( ); i++ )
    {
        const double conversionFactor = supportedCollection.getSingleObservationSets( ).at( i )->getResidualStateSpaceConversionFactor( );
        BOOST_CHECK_EQUAL( convertedResidualBlocks.at( i ).size( ), nativeResidualBlocks.at( i ).size( ) );
        for( int j = 0; j < nativeResidualBlocks.at( i ).size( ); j++ )
        {
            expectedConvertedResiduals( currentIndex ) = nativeResidualBlocks.at( i )( j ) * conversionFactor;
            BOOST_CHECK_CLOSE_FRACTION(
                    convertedResidualBlocks.at( i )( j ), nativeResidualBlocks.at( i )( j ) * conversionFactor, 1.0e-14 );
            currentIndex++;
        }
    }
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( convertedResiduals, expectedConvertedResiduals, 1.0e-14 );

    const Eigen::VectorXd parsedNativeRangeResiduals = supportedCollection.getConcatenatedResiduals( observationParser( dsn_n_way_range ) );
    const Eigen::VectorXd parsedConvertedRangeResiduals =
            supportedCollection.getConcatenatedResidualsInStateSpace( observationParser( dsn_n_way_range ) );
    Eigen::VectorXd expectedParsedRangeResiduals = parsedNativeRangeResiduals * 2.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( parsedConvertedRangeResiduals, expectedParsedRangeResiduals, 1.0e-15 );

    std::vector< std::shared_ptr< SingleObservationSet< double, double > > > mixedSets = { rangeSet, dopplerSet, unsupportedSet };
    ObservationCollection< double, double > mixedCollection( mixedSets );
    BOOST_CHECK_THROW( mixedCollection.getConcatenatedResidualsInStateSpace( ), std::runtime_error );
    BOOST_CHECK_NO_THROW( mixedCollection.getConcatenatedResiduals( ) );
    BOOST_CHECK_NO_THROW( ( mixedCollection.getConcatenatedResidualsInStateSpace( observationParser( dsn_n_way_range ) ) ) );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
