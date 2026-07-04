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
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/observation_models/observableTypes.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::observation_models;
using namespace tudat::mathematical_constants;

BOOST_AUTO_TEST_SUITE( test_residual_wrapping )

//! Test that isResidualWrappingRequired returns true for angular types and false for others.
BOOST_AUTO_TEST_CASE( testIsResidualWrappingRequired )
{
    // Types that SHOULD require wrapping
    BOOST_CHECK( isResidualWrappingRequired( angular_position ) );
    BOOST_CHECK( isResidualWrappingRequired( relative_angular_position ) );
    BOOST_CHECK( isResidualWrappingRequired( azimuth_elevation_angle ) );
    BOOST_CHECK( isResidualWrappingRequired( euler_angle_313_observable ) );

    // Types that should NOT require wrapping
    BOOST_CHECK( !isResidualWrappingRequired( one_way_range ) );
    BOOST_CHECK( !isResidualWrappingRequired( one_way_doppler ) );
    BOOST_CHECK( !isResidualWrappingRequired( position_observable ) );
    BOOST_CHECK( !isResidualWrappingRequired( velocity_observable ) );
    BOOST_CHECK( !isResidualWrappingRequired( one_way_differenced_range ) );
    BOOST_CHECK( !isResidualWrappingRequired( n_way_range ) );
    BOOST_CHECK( !isResidualWrappingRequired( two_way_doppler ) );
    BOOST_CHECK( !isResidualWrappingRequired( relative_position_observable ) );
    BOOST_CHECK( !isResidualWrappingRequired( n_way_differenced_range ) );
    BOOST_CHECK( !isResidualWrappingRequired( doppler_measured_frequency ) );
    BOOST_CHECK( !isResidualWrappingRequired( dsn_one_way_averaged_doppler ) );
    BOOST_CHECK( !isResidualWrappingRequired( dsn_n_way_averaged_doppler ) );
    BOOST_CHECK( !isResidualWrappingRequired( dsn_n_way_range ) );
    BOOST_CHECK( !isResidualWrappingRequired( differenced_time_of_arrival ) );
    BOOST_CHECK( !isResidualWrappingRequired( one_way_doppler_measured_frequency ) );
    BOOST_CHECK( !isResidualWrappingRequired( differenced_frequency_of_arrival ) );
    BOOST_CHECK( !isResidualWrappingRequired( pixel_coordinates ) );
}

//! Test the ResidualWrappingRange struct and its helper methods.
BOOST_AUTO_TEST_CASE( testResidualWrappingRangeStruct )
{
    // Default constructor (no wrapping)
    ResidualWrappingRange defaultRange;
    BOOST_CHECK_SMALL( defaultRange.minimumRange, 1.0e-15 );
    BOOST_CHECK_SMALL( defaultRange.maximumRange, 1.0e-15 );
    BOOST_CHECK_SMALL( defaultRange.period( ), 1.0e-15 );
    BOOST_CHECK_SMALL( defaultRange.center( ), 1.0e-15 );

    // [0, 2*pi] range
    ResidualWrappingRange range0Pi( 0.0, 2.0 * PI );
    BOOST_CHECK_CLOSE( range0Pi.minimumRange, 0.0, 1.0e-15 );
    BOOST_CHECK_CLOSE( range0Pi.maximumRange, 2.0 * PI, 1.0e-15 );
    BOOST_CHECK_CLOSE( range0Pi.period( ), 2.0 * PI, 1.0e-15 );
    BOOST_CHECK_CLOSE( range0Pi.center( ), PI, 1.0e-15 );

    // [-pi, pi] range
    ResidualWrappingRange rangeMinusPiPi( -PI, PI );
    BOOST_CHECK_CLOSE( rangeMinusPiPi.minimumRange, -PI, 1.0e-15 );
    BOOST_CHECK_CLOSE( rangeMinusPiPi.maximumRange, PI, 1.0e-15 );
    BOOST_CHECK_CLOSE( rangeMinusPiPi.period( ), 2.0 * PI, 1.0e-15 );
    BOOST_CHECK_SMALL( rangeMinusPiPi.center( ), 1.0e-15 );

    // [-pi/2, pi/2] range
    ResidualWrappingRange rangeHalfPi( -0.5 * PI, 0.5 * PI );
    BOOST_CHECK_CLOSE( rangeHalfPi.minimumRange, -0.5 * PI, 1.0e-15 );
    BOOST_CHECK_CLOSE( rangeHalfPi.maximumRange, 0.5 * PI, 1.0e-15 );
    BOOST_CHECK_CLOSE( rangeHalfPi.period( ), PI, 1.0e-15 );
    BOOST_CHECK_SMALL( rangeHalfPi.center( ), 1.0e-15 );
}

//! Test getResidualWrappingRanges for angular_position.
BOOST_AUTO_TEST_CASE( testAngularPositionWrappingRanges )
{
    std::vector< ResidualWrappingRange > ranges = getResidualWrappingRanges( angular_position );
    BOOST_CHECK_EQUAL( ranges.size( ), 2 );

    // Component 0: RA wraps to [-pi, pi]
    BOOST_CHECK_CLOSE( ranges[ 0 ].minimumRange, -PI, 1.0e-15 );
    BOOST_CHECK_CLOSE( ranges[ 0 ].maximumRange, PI, 1.0e-15 );

    // Component 1: DEC wraps to [-pi/2, pi/2]
    BOOST_CHECK_CLOSE( ranges[ 1 ].minimumRange, -0.5 * PI, 1.0e-15 );
    BOOST_CHECK_CLOSE( ranges[ 1 ].maximumRange, 0.5 * PI, 1.0e-15 );
}

//! Test getResidualWrappingRanges for non-wrapped observable types.
BOOST_AUTO_TEST_CASE( testNonWrappedObservableTypes )
{
    BOOST_CHECK( getResidualWrappingRanges( one_way_range ).empty( ) );
    BOOST_CHECK( getResidualWrappingRanges( one_way_doppler ).empty( ) );
    BOOST_CHECK( getResidualWrappingRanges( position_observable ).empty( ) );
    BOOST_CHECK( getResidualWrappingRanges( velocity_observable ).empty( ) );
    BOOST_CHECK( getResidualWrappingRanges( pixel_coordinates ).empty( ) );
}

//! Test the wrapping function directly for angular_position residuals.
BOOST_AUTO_TEST_CASE( testWrappingOfAngularPositionResiduals )
{
    // We test the wrapping concept directly: a residual with value outside [0, 2*pi]
    // for RA should get wrapped back.

    std::vector< ResidualWrappingRange > ranges = getResidualWrappingRanges( angular_position );

    // Test RA component wrapping: wrap to [-pi, pi]
    double raResidual = 3.0 * PI;  // Should wrap to -PI (since 3*pi - 2*pi*round(1.5) = 3*pi - 4*pi = -pi)
    double period = ranges[ 0 ].period( );
    double center = ranges[ 0 ].center( );
    double wrappedRa = raResidual - period * std::round( ( raResidual - center ) / period );
    BOOST_CHECK_CLOSE( wrappedRa, -PI, 1.0e-15 );
    BOOST_CHECK( wrappedRa >= ranges[ 0 ].minimumRange - 1.0e-12 );
    BOOST_CHECK( wrappedRa <= ranges[ 0 ].maximumRange + 1.0e-12 );

    // Test RA with negative value
    raResidual = -0.5 * PI;  // Already in [-pi, pi], should stay at -0.5*pi
    wrappedRa = raResidual - period * std::round( ( raResidual - center ) / period );
    BOOST_CHECK_CLOSE( wrappedRa, -0.5 * PI, 1.0e-15 );
    BOOST_CHECK( wrappedRa >= ranges[ 0 ].minimumRange - 1.0e-12 );
    BOOST_CHECK( wrappedRa <= ranges[ 0 ].maximumRange + 1.0e-12 );

    // Test RA with value within range (shouldn't change)
    raResidual = 0.5 * PI;
    wrappedRa = raResidual - period * std::round( ( raResidual - center ) / period );
    BOOST_CHECK_CLOSE( wrappedRa, 0.5 * PI, 1.0e-15 );

    // Test DEC component wrapping: wrap to [-pi/2, pi/2]
    double decResidual = 2.0 * PI;  // Should wrap to 0 (since 2*pi, center=0: round((2*pi)/(pi)) = 2, so 2*pi - pi*2 = 0)
    period = ranges[ 1 ].period( );
    center = ranges[ 1 ].center( );
    double wrappedDec = decResidual - period * std::round( ( decResidual - center ) / period );
    BOOST_CHECK_SMALL( wrappedDec, 1.0e-15 );
    BOOST_CHECK( wrappedDec >= ranges[ 1 ].minimumRange - 1.0e-12 );
    BOOST_CHECK( wrappedDec <= ranges[ 1 ].maximumRange + 1.0e-12 );

    // Test DEC with value at range boundary
    decResidual = 0.75 * PI;  // Should wrap to -0.25*pi (since 0.75*pi - pi = -0.25*pi)
    wrappedDec = decResidual - period * std::round( ( decResidual - center ) / period );
    BOOST_CHECK_CLOSE( wrappedDec, -0.25 * PI, 1.0e-15 );
}

//! Test the wrapping function for Euler angle 313 residuals.
BOOST_AUTO_TEST_CASE( testWrappingOfEulerAngleResiduals )
{
    std::vector< ResidualWrappingRange > ranges = getResidualWrappingRanges( euler_angle_313_observable );
    BOOST_CHECK_EQUAL( ranges.size( ), 3 );

    for( int comp = 0; comp < 3; comp++ )
    {
        double period = ranges[ comp ].period( );
        double center = ranges[ comp ].center( );

        // Test wrapping of a value outside [-pi, pi]
        double residual = 1.5 * PI;
        double wrapped = residual - period * std::round( ( residual - center ) / period );
        BOOST_CHECK_CLOSE( wrapped, -0.5 * PI, 1.0e-15 );
        BOOST_CHECK( wrapped >= ranges[ comp ].minimumRange - 1.0e-12 );
        BOOST_CHECK( wrapped <= ranges[ comp ].maximumRange + 1.0e-12 );

        // Test negative wrapping
        residual = -1.5 * PI;
        wrapped = residual - period * std::round( ( residual - center ) / period );
        BOOST_CHECK_CLOSE( wrapped, 0.5 * PI, 1.0e-15 );

        // Test value already within range
        residual = 0.5 * PI;
        wrapped = residual - period * std::round( ( residual - center ) / period );
        BOOST_CHECK_CLOSE( wrapped, 0.5 * PI, 1.0e-15 );

        // Test large offset
        residual = 10.0 * PI;
        wrapped = residual - period * std::round( ( residual - center ) / period );
        BOOST_CHECK_SMALL( wrapped, 1.0e-15 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
