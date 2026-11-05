/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_MAIN

#include <algorithm>
#include <limits>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/basics/testMacros.h"

#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/relativity/schwarzschildMetric.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::coordinate_conversions;
using namespace tudat::relativity;
using namespace tudat::physical_constants;

BOOST_AUTO_TEST_SUITE( test_schwarzschild_metric )

BOOST_AUTO_TEST_CASE( testSchwarzschildMetricPerturbation )
{
    const double sunGravitationalParameter = 132712440018.0E9;
    auto ppnParameterSet = std::make_shared< PPNParameterSet >( 1.0, 1.0 );
    auto schwarzschildMetric = std::make_shared< HarmonicSchwarzschildMetric >(
            [ = ]( ) { return sunGravitationalParameter; }, ppnParameterSet, "Sun", false );

    Eigen::Vector6d testState = Eigen::Vector6d::Zero( );
    Eigen::Vector3d sphericalPosition;
    sphericalPosition << 50.0E9, mathematical_constants::PI / 6.0, 3.0 * mathematical_constants::PI / 8.0;
    testState.segment( 0, 3 ) = convertSphericalToCartesian< double >( sphericalPosition );

    schwarzschildMetric->update( testState, 0.0, true, true );

    const Eigen::Matrix4d metricPerturbation = schwarzschildMetric->getCurrentCovariantMetricPeturbation( );

    const double distance = testState.segment( 0, 3 ).norm( );
    const double scaledPotential = sunGravitationalParameter * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT / distance;

    BOOST_CHECK_CLOSE_FRACTION( metricPerturbation( 0, 0 ),
                                2.0 * scaledPotential - 2.0 * ppnParameterSet->getParameterBeta( ) * scaledPotential * scaledPotential,
                                1.0E-14 );

    Eigen::Matrix3d expectedSpatial = 2.0 * ppnParameterSet->getParameterGamma( ) * scaledPotential * Eigen::Matrix3d::Identity( );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( metricPerturbation.block( 1, 1, 3, 3 ), expectedSpatial, 1.0E-14 );

    for( int row = 0; row < 4; ++row )
    {
        for( int col = 0; col < 4; ++col )
        {
            if( row != col )
            {
                BOOST_CHECK_SMALL( metricPerturbation( row, col ), 1.0E-30 );
            }
        }
    }

    const Eigen::Matrix4d covariantMetric = schwarzschildMetric->getCurrentCovariantMetric( );
    const Eigen::Matrix4d contravariantMetric = schwarzschildMetric->getCurrentContravariantMetric( );
    const Eigen::Matrix4d identityCheck = covariantMetric * contravariantMetric;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( identityCheck, Eigen::Matrix4d::Identity( ), 5.0E-13 );
}

// Analytic vs numeric partials: the harmonic-Schwarzschild metric provides analytic Christoffel
// symbols. Validate them against the defining relation
//   Gamma^lambda_{mu nu} = 1/2 g^{lambda sigma} ( d_mu g_{sigma nu} + d_nu g_{sigma mu} - d_sigma g_{mu nu} )
// using central finite differences of the metric. The metric is static, so d_0 g = 0 and only
// spatial derivatives appear. The perturbation h (~1e-9) is differenced rather than the full
// metric (~1) to avoid catastrophic roundoff in the ~1e-16 derivatives.
BOOST_AUTO_TEST_CASE( testSchwarzschildChristoffelAgainstFiniteDifferences )
{
    const double earthGravitationalParameter = 3.986004418E14;
    auto ppnParameterSet = std::make_shared< PPNParameterSet >( 1.0, 1.0 );
    auto metric = std::make_shared< HarmonicSchwarzschildMetric >(
            [ = ]( ) { return earthGravitationalParameter; }, ppnParameterSet, "Earth", false );

    // Generic (non-symmetric, off-axis) evaluation position.
    Eigen::Vector6d state = Eigen::Vector6d::Zero( );
    state.segment( 0, 3 ) = ( Eigen::Vector3d( ) << 7.0E6, -3.0E6, 4.0E6 ).finished( );

    metric->update( state, 0.0, true, true );
    const std::vector< Eigen::Matrix4d > analyticChristoffel = metric->getCurrentChristoffelSymbols( );
    const Eigen::Matrix4d contravariantMetric = metric->getCurrentContravariantMetric( );

    // Spatial central differences of the metric PERTURBATION (d_0 h = 0, static metric).
    auto perturbationAt = [ & ]( const Eigen::Vector3d& position ) {
        Eigen::Vector6d s = state;
        s.segment( 0, 3 ) = position;
        metric->update( s, 0.0, true, false );
        return metric->getCurrentCovariantMetricPeturbation( );
    };
    const double step = 100.0;                                                       // m
    std::vector< Eigen::Matrix4d > metricDerivative( 4, Eigen::Matrix4d::Zero( ) );  // d_sigma g (sigma=0 stays zero)
    for( int axis = 0; axis < 3; ++axis )
    {
        Eigen::Vector3d plus = state.segment( 0, 3 ), minus = state.segment( 0, 3 );
        plus( axis ) += step;
        minus( axis ) -= step;
        metricDerivative[ axis + 1 ] = ( perturbationAt( plus ) - perturbationAt( minus ) ) / ( 2.0 * step );
    }

    // Numeric Christoffel from the defining relation.
    std::vector< Eigen::Matrix4d > numericChristoffel( 4, Eigen::Matrix4d::Zero( ) );
    for( int lambda = 0; lambda < 4; ++lambda )
    {
        for( int mu = 0; mu < 4; ++mu )
        {
            for( int nu = 0; nu < 4; ++nu )
            {
                double value = 0.0;
                for( int sigma = 0; sigma < 4; ++sigma )
                {
                    value += 0.5 * contravariantMetric( lambda, sigma ) *
                            ( metricDerivative[ mu ]( sigma, nu ) + metricDerivative[ nu ]( sigma, mu ) -
                              metricDerivative[ sigma ]( mu, nu ) );
                }
                numericChristoffel[ lambda ]( mu, nu ) = value;
            }
        }
    }

    // Compare against the largest analytic Christoffel magnitude (components span 0 .. ~1e-16,
    // so a per-element fractional tolerance is ill-conditioned on the near-zero entries).
    double maxMagnitude = 0.0, maxDifference = 0.0;
    for( int lambda = 0; lambda < 4; ++lambda )
    {
        maxMagnitude = std::max( maxMagnitude, analyticChristoffel[ lambda ].cwiseAbs( ).maxCoeff( ) );
        maxDifference = std::max( maxDifference, ( analyticChristoffel[ lambda ] - numericChristoffel[ lambda ] ).cwiseAbs( ).maxCoeff( ) );
    }
    BOOST_CHECK_GT( maxMagnitude, 0.0 );
    BOOST_CHECK_SMALL( maxDifference / maxMagnitude, 1.0E-6 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
