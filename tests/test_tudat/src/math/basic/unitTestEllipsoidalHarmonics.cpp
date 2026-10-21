/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *          Garmier, R., Barriot, J.-P. (2001), Celestial Mechanics and Dynamical Astronomy, 79(4), 235-275.
 *          Hu, X., Jekeli, C. (2015), Journal of Geodesy, 89(2), 159-177.
 *          See tudat/math/basic/ellipsoidalHarmonics.h for the full citations.
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <vector>

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/included/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/math/basic/ellipsoidalHarmonics.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::basic_mathematics;

BOOST_AUTO_TEST_SUITE( test_ellipsoidal_harmonics )

//! The ellipsoidal coordinates of a point must give back that point through their defining relations.
BOOST_AUTO_TEST_CASE( testEllipsoidalCoordinates )
{
    // Focal parameters of the 67P reference ellipsoid, in metres.
    const double longestSemiAxis = 2.802437330193671E+03;
    const double h = std::sqrt( longestSemiAxis * longestSemiAxis - 2.090307360485069E+03 * 2.090307360485069E+03 );
    const double k = std::sqrt( longestSemiAxis * longestSemiAxis - 1.867442428762624E+03 * 1.867442428762624E+03 );

    std::vector< Eigen::Vector3d > positions;
    positions.push_back( ( Eigen::Vector3d( ) << 3.0E3, 1.0E3, -2.0E3 ).finished( ) );
    positions.push_back( ( Eigen::Vector3d( ) << -4.5E3, 2.2E3, 1.1E3 ).finished( ) );
    positions.push_back( ( Eigen::Vector3d( ) << 0.0, 0.0, 2.81E3 ).finished( ) );
    positions.push_back( ( Eigen::Vector3d( ) << 5.0E3, 1.0, 1.0 ).finished( ) );

    for( unsigned int i = 0; i < positions.size( ); i++ )
    {
        const Eigen::Vector3d ellipsoidal = convertCartesianToEllipsoidalCoordinates( positions.at( i ), h, k );
        const double rho = ellipsoidal( 0 );
        const double phi = ellipsoidal( 1 );
        const double gamma = ellipsoidal( 2 );

        const Eigen::Vector3d recovered = ( Eigen::Vector3d( ) << rho * std::sin( phi ),
                                            std::sqrt( rho * rho - h * h ) * std::cos( phi ) * std::cos( gamma ),
                                            std::sqrt( rho * rho - k * k ) * std::cos( phi ) * std::sin( gamma ) )
                                                  .finished( );

        BOOST_CHECK_SMALL( ( recovered - positions.at( i ) ).norm( ) / positions.at( i ).norm( ), 1.0E-12 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
