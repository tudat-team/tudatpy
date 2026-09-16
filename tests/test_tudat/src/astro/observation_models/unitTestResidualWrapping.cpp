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

#include <limits>
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

//! Test the production wrapping function for two-component angular observables.
BOOST_AUTO_TEST_CASE( testWrappingOfAngularObservableResiduals )
{
    const double fullPeriod = 2.0 * PI;
    const double offset = 10.0 * std::numeric_limits< double >::epsilon( );
    const double nearPositivePeriod = fullPeriod - offset;
    const double nearNegativePeriod = -fullPeriod + offset;
    const std::vector< ObservableType > observableTypes = { angular_position, relative_angular_position, azimuth_elevation_angle };

    for( const ObservableType observableType : observableTypes )
    {
        Eigen::VectorXd residuals( 12 );
        residuals << 0.1, nearPositivePeriod, 0.25, nearNegativePeriod, -0.25, PI + 0.02, 0.1, -PI - 0.02, -0.1, 0.05, -0.2, -0.05;

        simulation_setup::wrapObservationResiduals< double >( residuals, std::make_pair( 1, 10 ), observableType );

        const Eigen::VectorXd expectedResiduals = ( Eigen::VectorXd( 12 ) << 0.1,
                                                    nearPositivePeriod - fullPeriod,
                                                    0.25,
                                                    nearNegativePeriod + fullPeriod,
                                                    -0.25,
                                                    -PI + 0.02,
                                                    0.1,
                                                    PI - 0.02,
                                                    -0.1,
                                                    0.05,
                                                    -0.2,
                                                    -0.05 )
                                                          .finished( );

        for( int i = 0; i < residuals.rows( ); i++ )
        {
            BOOST_CHECK_SMALL( residuals( i ) - expectedResiduals( i ), 1.0e-14 );
        }
        BOOST_CHECK_LT( residuals( 1 ), 0.0 );
        BOOST_CHECK_GT( residuals( 3 ), 0.0 );
    }
}

//! Test normalized right ascension residual wrapping using the observed declination.
BOOST_AUTO_TEST_CASE( testWrappingOfNormalizedAngularPositionResiduals )
{
    Eigen::VectorXd residuals( 6 );
    residuals << PI - 0.02, 0.2, -2.0 * PI + 0.03, -0.4, 0.3, 0.1;

    Eigen::VectorXd observedObservations( 6 );
    observedObservations << 0.1, PI / 3.0, -0.2, 0.0, 0.0, PI / 2.0;

    ResidualWrappingSettings residualWrappingSettings;
    residualWrappingSettings.normalizeRightAscension = true;
    simulation_setup::wrapObservationResiduals< double >(
            residuals, std::make_pair( 0, 6 ), angular_position, observedObservations, residualWrappingSettings );

    const Eigen::VectorXd expectedResiduals = ( Eigen::VectorXd( 6 ) << -0.02, 0.2, 0.03, -0.4, 0.3, 0.1 ).finished( );

    for( int i = 0; i < residuals.rows( ); i++ )
    {
        BOOST_CHECK_SMALL( residuals( i ) - expectedResiduals( i ), 1.0e-14 );
    }
}

//! Test the production wrapping function for Euler angle 313 residuals.
BOOST_AUTO_TEST_CASE( testWrappingOfEulerAngleResiduals )
{
    Eigen::VectorXd residuals( 9 );
    residuals << PI + 0.02, 0.2, -PI - 0.02, -PI - 0.03, -0.1, PI + 0.03, 0.05, 0.15, -0.04;

    simulation_setup::wrapObservationResiduals< double >( residuals, std::make_pair( 0, 9 ), euler_angle_313_observable );

    const Eigen::VectorXd expectedResiduals =
            ( Eigen::VectorXd( 9 ) << -PI + 0.02, 0.2, PI - 0.02, PI - 0.03, -0.1, -PI + 0.03, 0.05, 0.15, -0.04 ).finished( );

    for( int i = 0; i < residuals.rows( ); i++ )
    {
        BOOST_CHECK_SMALL( residuals( i ) - expectedResiduals( i ), 1.0e-14 );
    }
}

//! Test that the production wrapping function is a no-op for nonperiodic observables.
BOOST_AUTO_TEST_CASE( testWrappingOfNonperiodicObservableResiduals )
{
    Eigen::VectorXd residuals( 3 );
    residuals << 2.5, -3.0, 0.75;
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
