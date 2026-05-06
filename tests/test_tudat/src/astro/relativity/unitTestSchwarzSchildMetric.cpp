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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/unit_test.hpp>

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
            [=]( ){ return sunGravitationalParameter; },
            ppnParameterSet,
            "Sun",
            false );

    Eigen::Vector6d testState = Eigen::Vector6d::Zero( );
    Eigen::Vector3d sphericalPosition;
    sphericalPosition << 50.0E9,
            mathematical_constants::PI / 6.0,
            3.0 * mathematical_constants::PI / 8.0;
    testState.segment( 0, 3 ) = convertSphericalToCartesian< double >( sphericalPosition );

    schwarzschildMetric->update( testState, 0.0, true, true );

    const Eigen::Matrix4d metricPerturbation = schwarzschildMetric->getCurrentCovariantMetricPeturbation( );

    const double distance = testState.segment( 0, 3 ).norm( );
    const double scaledPotential = sunGravitationalParameter * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT / distance;

    BOOST_CHECK_CLOSE_FRACTION(
            metricPerturbation( 0, 0 ),
            2.0 * scaledPotential - 2.0 * ppnParameterSet->getParameterBeta( ) * scaledPotential * scaledPotential,
            1.0E-14 );

    Eigen::Matrix3d expectedSpatial = 2.0 * ppnParameterSet->getParameterGamma( ) * scaledPotential * Eigen::Matrix3d::Identity( );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            metricPerturbation.block( 1, 1, 3, 3 ),
            expectedSpatial,
            1.0E-14 );

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

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            identityCheck,
            Eigen::Matrix4d::Identity( ),
            5.0E-13 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
