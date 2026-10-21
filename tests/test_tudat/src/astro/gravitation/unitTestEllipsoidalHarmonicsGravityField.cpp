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

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/included/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/astro/gravitation/ellipsoidalHarmonicsGravityField.h"
#include "tudat/basics/testMacros.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/ellipsoidalHarmonicsFileReader.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::gravitation;
using namespace tudat::input_output;

BOOST_AUTO_TEST_SUITE( test_ellipsoidal_harmonics_gravity_field )

//! Reference ellipsoid of 67P/Churyumov-Gerasimenko, longest semi-axis first, in metres.
const Eigen::Vector3d semiAxes = ( Eigen::Vector3d( ) << 2.802437330193671E+03, 2.090307360485069E+03, 1.867442428762624E+03 ).finished( );

std::string getDataDirectory( )
{
    return paths::getTudatTestDataPath( ) + "/ellipsoidal_harmonics/";
}

//! Degree 15 field of 67P, with the radial parts either tabulated or integrated at every call.
std::shared_ptr< EllipsoidalHarmonicsGravityField > createTestField( const int radialTableGridPoints = 3200 )
{
    return std::make_shared< EllipsoidalHarmonicsGravityField >(
            semiAxes,
            readEllipsoidalHarmonicCoefficients( getDataDirectory( ) + "67p_ellipsoidal_harmonic_coefficients.txt", 15 ),
            readLameFunctionCoefficients( getDataDirectory( ) + "67p_lame_coefficients.txt", 15 ),
            readEllipsoidalHarmonicTransformFactors( getDataDirectory( ) + "67p_transform_factors.txt", 15 ),
            24,
            radialTableGridPoints,
            0.0,
            "67P_fixed" );
}

//! Test positions outside the reference ellipsoid and away from its long axis.
std::vector< Eigen::Vector3d > getTestPositions( )
{
    return { ( Eigen::Vector3d( ) << 1.0E3, 2.5E3, 2.0E3 ).finished( ),
             ( Eigen::Vector3d( ) << -3.0E3, 2.0E3, -1.5E3 ).finished( ),
             ( Eigen::Vector3d( ) << 0.0, 0.0, 3.5E3 ).finished( ),
             ( Eigen::Vector3d( ) << 8.0E3, 5.0E3, -4.0E3 ).finished( ) };
}

//! The potential must match values from a separate implementation of the same model at 500 points.
BOOST_AUTO_TEST_CASE( testPotentialAgainstReference )
{
    const std::shared_ptr< EllipsoidalHarmonicsGravityField > field = createTestField( );

    std::ifstream referenceFile( getDataDirectory( ) + "67p_reference_potentials.txt" );
    BOOST_REQUIRE( referenceFile.good( ) );

    double maximumRelativeError = 0.0;
    int numberOfPoints = 0;
    std::string line;
    while( std::getline( referenceFile, line ) )
    {
        std::istringstream lineStream( line );
        std::string token;
        std::vector< double > values;
        while( lineStream >> token )
        {
            values.push_back( parseFortranDouble( token ) );
        }
        if( values.size( ) != 4 )
        {
            continue;
        }

        const double potential = field->getGravitationalPotential( Eigen::Vector3d( values.at( 0 ), values.at( 1 ), values.at( 2 ) ) );
        maximumRelativeError = std::max( maximumRelativeError, std::fabs( potential / values.at( 3 ) - 1.0 ) );
        numberOfPoints++;
    }

    BOOST_CHECK_EQUAL( numberOfPoints, 500 );
    BOOST_CHECK_SMALL( maximumRelativeError, 2.5E-9 );
}

//! Far from the body the potential and its gradient must approach those of a point mass.
BOOST_AUTO_TEST_CASE( testFarFieldPointMass )
{
    const std::shared_ptr< EllipsoidalHarmonicsGravityField > field = createTestField( );

    // The gravitational parameter is taken from the degree zero coefficient of the input file.
    const double gravitationalParameter = field->getGravitationalParameter( );
    BOOST_CHECK_CLOSE_FRACTION( gravitationalParameter, 668.78268559825256, 1.0E-10 );

    const double distance = 1.0E4 * semiAxes( 0 );
    const Eigen::Vector3d direction = Eigen::Vector3d( 0.3, 1.0, -0.7 ).normalized( );

    BOOST_CHECK_CLOSE_FRACTION( field->getGravitationalPotential( distance * direction ), gravitationalParameter / distance, 1.0E-5 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( field->getGradientOfPotential( distance * direction ),
                                       ( -gravitationalParameter / ( distance * distance ) * direction ).eval( ),
                                       1.0E-4 );
}

//! The gradient and Hessian must match finite differences of the potential and gradient, also on the long axis.
BOOST_AUTO_TEST_CASE( testPotentialDerivatives )
{
    const std::shared_ptr< EllipsoidalHarmonicsGravityField > field = createTestField( );

    std::vector< Eigen::Vector3d > positions = getTestPositions( );
    positions.push_back( ( Eigen::Vector3d( ) << 1.5 * semiAxes( 0 ), 0.0, 0.0 ).finished( ) );

    for( unsigned int i = 0; i < positions.size( ); i++ )
    {
        const Eigen::Vector3d position = positions.at( i );
        const double step = position.norm( ) * 1.0E-5;

        Eigen::Vector3d numericalGradient;
        Eigen::Matrix3d numericalHessian;
        for( int j = 0; j < 3; j++ )
        {
            Eigen::Vector3d plus = position;
            Eigen::Vector3d minus = position;
            plus( j ) += step;
            minus( j ) -= step;
            numericalGradient( j ) =
                    ( field->getGravitationalPotential( plus ) - field->getGravitationalPotential( minus ) ) / ( 2.0 * step );
            numericalHessian.col( j ) = ( field->getGradientOfPotential( plus ) - field->getGradientOfPotential( minus ) ) / ( 2.0 * step );
        }

        const Eigen::Vector3d gradient = field->getGradientOfPotential( position );
        const Eigen::Matrix3d hessian = field->getHessianOfPotential( position );

        BOOST_CHECK_SMALL( ( gradient - numericalGradient ).norm( ) / numericalGradient.norm( ), 1.0E-6 );
        BOOST_CHECK_SMALL( ( hessian - numericalHessian ).cwiseAbs( ).maxCoeff( ) / numericalHessian.cwiseAbs( ).maxCoeff( ), 1.0E-4 );

        // Outside the body the potential obeys Laplace's equation, so the Hessian has zero trace.
        BOOST_CHECK_SMALL( std::fabs( hessian.trace( ) ) / hessian.cwiseAbs( ).maxCoeff( ), 1.0E-4 );
    }
}

//! Tabulated radial parts must give the same field as integrating them at every call, also beyond the table.
BOOST_AUTO_TEST_CASE( testRadialTable )
{
    const std::shared_ptr< EllipsoidalHarmonicsGravityField > tabulated = createTestField( 3200 );
    const std::shared_ptr< EllipsoidalHarmonicsGravityField > exact = createTestField( 0 );

    // The last position lies beyond the default table, which ends at fifty times the longest semi-axis.
    std::vector< Eigen::Vector3d > positions = getTestPositions( );
    positions.push_back( 100.0 * semiAxes( 0 ) * Eigen::Vector3d( 0.3, 1.0, -0.7 ).normalized( ) );

    for( unsigned int i = 0; i < positions.size( ); i++ )
    {
        const Eigen::Vector3d position = positions.at( i );

        BOOST_CHECK_CLOSE_FRACTION(
                tabulated->getGravitationalPotential( position ), exact->getGravitationalPotential( position ), 1.0E-10 );

        const Eigen::Vector3d exactGradient = exact->getGradientOfPotential( position );
        BOOST_CHECK_SMALL( ( tabulated->getGradientOfPotential( position ) - exactGradient ).norm( ) / exactGradient.norm( ), 1.0E-8 );

        const Eigen::Matrix3d exactHessian = exact->getHessianOfPotential( position );
        BOOST_CHECK_SMALL( ( tabulated->getHessianOfPotential( position ) - exactHessian ).cwiseAbs( ).maxCoeff( ) /
                                   exactHessian.cwiseAbs( ).maxCoeff( ),
                           1.0E-5 );
    }
}

//! Potential, gradient and Hessian share cached values, so the order of the calls must not change the results.
BOOST_AUTO_TEST_CASE( testCachedEvaluation )
{
    const Eigen::Vector3d position = ( Eigen::Vector3d( ) << 5.0E3, 4.0E3, -3.0E3 ).finished( );
    const Eigen::Vector3d otherPosition = ( Eigen::Vector3d( ) << -4.0E3, 6.0E3, 2.0E3 ).finished( );

    // Each reference value comes from a field that has computed nothing else.
    const double potential = createTestField( )->getGravitationalPotential( position );
    const Eigen::Vector3d gradient = createTestField( )->getGradientOfPotential( position );
    const Eigen::Matrix3d hessian = createTestField( )->getHessianOfPotential( position );

    const std::shared_ptr< EllipsoidalHarmonicsGravityField > field = createTestField( );
    field->getHessianOfPotential( position );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( field->getGradientOfPotential( position ), gradient, 1.0E-14 );
    BOOST_CHECK_CLOSE_FRACTION( field->getGravitationalPotential( position ), potential, 1.0E-14 );
    field->getGradientOfPotential( otherPosition );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( field->getHessianOfPotential( position ), hessian, 1.0E-14 );
    field->getGravitationalPotential( otherPosition );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( field->getGradientOfPotential( position ), gradient, 1.0E-14 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
