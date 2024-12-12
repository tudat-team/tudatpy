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

#include "tudat/basics/testMacros.h"

#include "tudat/astro/ground_stations/transmittingFrequencies.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::ground_stations;

BOOST_AUTO_TEST_SUITE( test_transmitting_frequencies )

BOOST_AUTO_TEST_CASE( testPiecewiseLinearFrequencyInterpolator )
{

    std::vector< Time > startTimes = { 1.0, 5.0, 9.0, 10.0 };
    std::vector< Time > endTimes = { 5.0, 9.0, 10.0, 15.0 };
    std::vector< double > rampRates = { 1.0, 2.0, 1.0, 3.0 };
    std::vector< double > startFrequency = { 10.0, 5.0, 10.0, 2.0 };

    PiecewiseLinearFrequencyInterpolator frequencyInterpolator = PiecewiseLinearFrequencyInterpolator(
            startTimes, endTimes, rampRates, startFrequency );

    // Frequency at t = 1.0s
    double frequency = startFrequency[ 0 ];
    double computedFrequency = frequencyInterpolator.template getTemplatedCurrentFrequency< >( 1.0 );
    BOOST_CHECK_EQUAL ( computedFrequency , frequency );

    // Frequency at t = 2.0s
    frequency = startFrequency[ 0 ] + rampRates[ 0 ] * ( 2.0 - startTimes[ 0 ] );
    computedFrequency = frequencyInterpolator.template getTemplatedCurrentFrequency< >( 2.0 );
    BOOST_CHECK_EQUAL ( computedFrequency, frequency );

    // Frequency at t = 5.0s: start frequency of the 1st ramp
    frequency = startFrequency[ 1 ];
    computedFrequency = frequencyInterpolator.template getTemplatedCurrentFrequency< >( 5.0 );
    BOOST_CHECK_EQUAL ( computedFrequency, frequency );

    // Frequency at t = 13.0s
    frequency = startFrequency[ 3 ] + rampRates[ 3 ] * ( 13.0 - startTimes[ 3 ] );
    computedFrequency = frequencyInterpolator.template getTemplatedCurrentFrequency< >( 13.0 );
    BOOST_CHECK_EQUAL ( computedFrequency, frequency );

    // Frequency at t = 15.0s
    frequency = startFrequency[ 3 ] + rampRates[ 3 ] * ( 15.0 - startTimes[ 3 ] );
    computedFrequency = frequencyInterpolator.template getTemplatedCurrentFrequency< >( 15.0 );
    BOOST_CHECK_EQUAL ( computedFrequency, frequency );

    // Integral between t = 2.0s and t = 4.0s
    double integral = 2.0 * ( 11.0 + 13.0 ) / 2.0;
    double computedIntegral = frequencyInterpolator.template getTemplatedFrequencyIntegral< >( 2.0, 4.0 );
    BOOST_CHECK_EQUAL ( computedIntegral, integral );

    // Integral between t = 5.0s and t = 7.0s
    integral = 2.0 * ( 5.0 + 9.0 ) / 2.0;
    computedIntegral = frequencyInterpolator.template getTemplatedFrequencyIntegral< >( 5.0, 7.0 );
    BOOST_CHECK_EQUAL ( computedIntegral, integral );

    // Integral between t = 7.0s and t = 9.0s
    integral = 2.0 * ( 9.0 + 13.0 ) / 2.0;
    computedIntegral = frequencyInterpolator.template getTemplatedFrequencyIntegral< >( 7.0, 9.0 );
    BOOST_CHECK_EQUAL ( computedIntegral, integral );

    // Integral between t = 5.0s and t = 9.0s
    integral = 4.0 * ( 5.0 + 13.0 ) / 2.0;
    computedIntegral = frequencyInterpolator.template getTemplatedFrequencyIntegral< >( 5.0, 9.0 );
    BOOST_CHECK_EQUAL ( integral, computedIntegral );

    // Integral between t = 2.0s and t = 13.0s
    integral = ( 3.0 * ( 11.0 + 14.0 ) + 4.0 * ( 5.0 + 13.0 ) + 1.0 * ( 11.0 + 10.0 ) + 3.0 * ( 2.0 + 11.0 ) ) / 2.0;
    computedIntegral = frequencyInterpolator.template getTemplatedFrequencyIntegral< >( 2.0, 13.0 );
    BOOST_CHECK_EQUAL ( computedIntegral, integral );

    // Check whether error is thrown for vectors with invalid size
    startTimes.push_back( 15.0 );
    bool errorThrown;
    try
    {
        frequencyInterpolator = PiecewiseLinearFrequencyInterpolator(
            startTimes, endTimes, rampRates, startFrequency );
        errorThrown = false;
    }
    catch( std::runtime_error const& )
    {
        errorThrown = true;
    }
    BOOST_CHECK( errorThrown );
//
//    // Check whether error is thrown for discontinuous start/end times
//    startTimes.pop_back( );
//    startTimes.pop_back( );
//    startTimes.push_back( 9.9 );
//    try
//    {
//        frequencyInterpolator = PiecewiseLinearFrequencyInterpolator(
//            startTimes, endTimes, rampRates, startFrequency );
//        errorThrown = false;
//    }
//    catch( std::runtime_error const& )
//    {
//        errorThrown = true;
//    }
//    BOOST_CHECK( errorThrown );
    

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
