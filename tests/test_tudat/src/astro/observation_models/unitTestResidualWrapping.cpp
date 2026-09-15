/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <utility>
#include <vector>

#include <Eigen/Core>
#include <boost/test/included/unit_test.hpp>

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManagerHelpers.h"

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
    BOOST_CHECK_SMALL( range0Pi.minimumRange, 1.0e-15 );
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

//! Test getResidualWrappingRanges for two-component angular observables.
BOOST_AUTO_TEST_CASE( testAngularObservableWrappingRanges )
{
    const std::vector< ObservableType > observableTypes = { angular_position, relative_angular_position, azimuth_elevation_angle };

    for( const ObservableType observableType : observableTypes )
    {
        const std::vector< ResidualWrappingRange > ranges = getResidualWrappingRanges( observableType );
        BOOST_REQUIRE_EQUAL( ranges.size( ), 2 );

        // Component 0 (RA / azimuth) wraps to [-pi, pi].
        BOOST_CHECK_CLOSE( ranges[ 0 ].minimumRange, -PI, 1.0e-15 );
        BOOST_CHECK_CLOSE( ranges[ 0 ].maximumRange, PI, 1.0e-15 );

        // Component 1 (DEC / elevation) is not periodic.
        BOOST_CHECK_SMALL( ranges[ 1 ].period( ), 1.0e-15 );
    }
}

//! Test getResidualWrappingRanges for the 3-1-3 Euler angle observable.
BOOST_AUTO_TEST_CASE( testEulerAngleWrappingRanges )
{
    const std::vector< ResidualWrappingRange > ranges = getResidualWrappingRanges( euler_angle_313_observable );
    BOOST_REQUIRE_EQUAL( ranges.size( ), 3 );

    BOOST_CHECK_CLOSE( ranges[ 0 ].minimumRange, -PI, 1.0e-15 );
    BOOST_CHECK_CLOSE( ranges[ 0 ].maximumRange, PI, 1.0e-15 );
    BOOST_CHECK_SMALL( ranges[ 1 ].period( ), 1.0e-15 );
    BOOST_CHECK_CLOSE( ranges[ 2 ].minimumRange, -PI, 1.0e-15 );
    BOOST_CHECK_CLOSE( ranges[ 2 ].maximumRange, PI, 1.0e-15 );
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

//! Test the production wrapping function for two-component angular observables.
BOOST_AUTO_TEST_CASE( testWrappingOfAngularObservableResiduals )
{
    const std::vector< ObservableType > observableTypes = { angular_position, relative_angular_position, azimuth_elevation_angle };

    for( const ObservableType observableType : observableTypes )
    {
        Eigen::VectorXd residuals( 10 );
        residuals << 42.0, 3.0 * PI, PI, -3.0 * PI, -PI, 0.5 * PI, 0.75 * PI, 10.0 * PI, -0.75 * PI, -42.0;

        simulation_setup::wrapObservationResiduals< double >( residuals, std::make_pair( 1, 8 ), observableType );

        const Eigen::VectorXd expectedResiduals =
                ( Eigen::VectorXd( 10 ) << 42.0, -PI, PI, PI, -PI, 0.5 * PI, 0.75 * PI, 0.0, -0.75 * PI, -42.0 ).finished( );

        for( int i = 0; i < residuals.rows( ); i++ )
        {
            BOOST_CHECK_SMALL( residuals( i ) - expectedResiduals( i ), 1.0e-14 );
        }
    }
}

//! Test normalized right ascension residual wrapping using the observed declination.
BOOST_AUTO_TEST_CASE( testWrappingOfNormalizedAngularPositionResiduals )
{
    Eigen::VectorXd residuals( 8 );
    residuals << 42.0, 0.75 * PI, 0.2, 1.5 * PI, -0.4, 0.3, 0.1, -42.0;

    Eigen::VectorXd observedObservations( 6 );
    observedObservations << 0.1, PI / 3.0, -0.2, 0.0, 0.0, PI / 2.0;

    ResidualWrappingSettings residualWrappingSettings;
    residualWrappingSettings.normalizeRightAscension = true;
    simulation_setup::wrapObservationResiduals< double >(
            residuals, std::make_pair( 1, 6 ), angular_position, observedObservations, residualWrappingSettings );

    const Eigen::VectorXd expectedResiduals =
            ( Eigen::VectorXd( 8 ) << 42.0, -0.25 * PI, 0.2, -0.5 * PI, -0.4, 0.3, 0.1, -42.0 ).finished( );

    for( int i = 0; i < residuals.rows( ); i++ )
    {
        BOOST_CHECK_SMALL( residuals( i ) - expectedResiduals( i ), 1.0e-14 );
    }
}

//! Test the production wrapping function for Euler angle 313 residuals.
BOOST_AUTO_TEST_CASE( testWrappingOfEulerAngleResiduals )
{
    Eigen::VectorXd residuals( 8 );
    residuals << 42.0, 1.5 * PI, PI, -1.5 * PI, 10.0 * PI, -PI, 3.0 * PI, -42.0;

    simulation_setup::wrapObservationResiduals< double >( residuals, std::make_pair( 1, 6 ), euler_angle_313_observable );

    const Eigen::VectorXd expectedResiduals = ( Eigen::VectorXd( 8 ) << 42.0, -0.5 * PI, PI, 0.5 * PI, 0.0, -PI, -PI, -42.0 ).finished( );

    for( int i = 0; i < residuals.rows( ); i++ )
    {
        BOOST_CHECK_SMALL( residuals( i ) - expectedResiduals( i ), 1.0e-14 );
    }
}

//! Test that the production wrapping function is a no-op for nonperiodic observables.
BOOST_AUTO_TEST_CASE( testWrappingOfNonperiodicObservableResiduals )
{
    Eigen::VectorXd residuals( 3 );
    residuals << 3.0 * PI, -4.0 * PI, 5.0 * PI;
    const Eigen::VectorXd expectedResiduals = residuals;

    simulation_setup::wrapObservationResiduals< double >( residuals, std::make_pair( 0, 3 ), position_observable );

    for( int i = 0; i < residuals.rows( ); i++ )
    {
        BOOST_CHECK_SMALL( residuals( i ) - expectedResiduals( i ), 1.0e-15 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
