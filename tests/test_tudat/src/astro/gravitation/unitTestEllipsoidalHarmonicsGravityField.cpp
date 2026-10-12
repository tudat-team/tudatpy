/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
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
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;

BOOST_AUTO_TEST_SUITE( test_ellipsoidal_harmonics_gravity_field )

//! Reference ellipsoid of 67P/Churyumov-Gerasimenko, longest semi-axis first, in metres.
Eigen::Vector3d getReferenceSemiAxes( )
{
    return ( Eigen::Vector3d( ) << 2.802437330193671E+03, 2.090307360485069E+03, 1.867442428762624E+03 ).finished( );
}

//! Gravitational parameter of 67P, in m^3/s^2.
const double gravitationalParameter = 666.2;

std::string getDataDirectory( )
{
    return paths::getTudatTestDataPath( ) + "/ellipsoidal_harmonics/";
}

std::shared_ptr< EllipsoidalHarmonicsGravityField > createTestField( const int maximumDegree = 15 )
{
    return std::make_shared< EllipsoidalHarmonicsGravityField >(
            gravitationalParameter,
            getReferenceSemiAxes( ),
            readEllipsoidalHarmonicCoefficients( getDataDirectory( ) + "67p_ellipsoidal_harmonic_coefficients.txt", maximumDegree ),
            readLameFunctionCoefficients( getDataDirectory( ) + "67p_lame_coefficients.txt", maximumDegree ),
            readEllipsoidalHarmonicTransformFactors( getDataDirectory( ) + "67p_transform_factors.txt", maximumDegree ),
            24,
            "67P_fixed" );
}

//! The focal parameters must follow from the semi-axes, and be ordered.
BOOST_AUTO_TEST_CASE( test_field_is_set_up_from_the_semi_axes )
{
    const std::shared_ptr< EllipsoidalHarmonicsGravityField > field = createTestField( );
    const Eigen::Vector3d semiAxes = getReferenceSemiAxes( );

    BOOST_CHECK_CLOSE_FRACTION(
            field->getFirstFocalParameter( ), std::sqrt( semiAxes( 0 ) * semiAxes( 0 ) - semiAxes( 1 ) * semiAxes( 1 ) ), 1.0E-14 );
    BOOST_CHECK_CLOSE_FRACTION(
            field->getSecondFocalParameter( ), std::sqrt( semiAxes( 0 ) * semiAxes( 0 ) - semiAxes( 2 ) * semiAxes( 2 ) ), 1.0E-14 );

    // The second focal parameter is the larger of the two, which the expansion relies on.
    BOOST_CHECK( field->getSecondFocalParameter( ) > field->getFirstFocalParameter( ) );

    BOOST_CHECK_EQUAL( field->getMaximumDegree( ), 15 );
    BOOST_CHECK_EQUAL( field->getFixedReferenceFrame( ), "67P_fixed" );
    BOOST_CHECK_CLOSE_FRACTION( field->getGravitationalParameter( ), gravitationalParameter, 1.0E-14 );
}

//! Semi-axes that are not strictly decreasing have no valid focal parameters, so are rejected.
BOOST_AUTO_TEST_CASE( test_badly_ordered_semi_axes_are_rejected )
{
    const int maximumDegree = 2;
    const EllipsoidalHarmonicCoefficients coefficients =
            readEllipsoidalHarmonicCoefficients( getDataDirectory( ) + "67p_ellipsoidal_harmonic_coefficients.txt", maximumDegree );
    const EllipsoidalHarmonicMatrices lameCoefficients =
            readLameFunctionCoefficients( getDataDirectory( ) + "67p_lame_coefficients.txt", maximumDegree );
    const EllipsoidalHarmonicMatrices transformFactors =
            readEllipsoidalHarmonicTransformFactors( getDataDirectory( ) + "67p_transform_factors.txt", maximumDegree );

    // Shortest axis given first.
    const Eigen::Vector3d reversed = ( Eigen::Vector3d( ) << 1.0E3, 2.0E3, 3.0E3 ).finished( );
    BOOST_CHECK_THROW(
            EllipsoidalHarmonicsGravityField( gravitationalParameter, reversed, coefficients, lameCoefficients, transformFactors ),
            std::runtime_error );

    // A sphere, where the two focal parameters both collapse to zero.
    const Eigen::Vector3d sphere = ( Eigen::Vector3d( ) << 2.0E3, 2.0E3, 2.0E3 ).finished( );
    BOOST_CHECK_THROW( EllipsoidalHarmonicsGravityField( gravitationalParameter, sphere, coefficients, lameCoefficients, transformFactors ),
                       std::runtime_error );
}

//! The potential must reproduce the independently computed reference values.
BOOST_AUTO_TEST_CASE( test_potential_matches_reference_data )
{
    const std::shared_ptr< EllipsoidalHarmonicsGravityField > field = createTestField( );

    std::ifstream referenceFile( getDataDirectory( ) + "67p_reference_potentials.txt" );
    BOOST_REQUIRE( referenceFile.good( ) );

    double maximumRelativeError = 0.0;
    int numberOfPoints = 0;
    std::string line;

    while( std::getline( referenceFile, line ) )
    {
        if( line.find_first_not_of( " \t\r\n" ) == std::string::npos )
        {
            continue;
        }

        std::istringstream lineStream( line );
        std::string token;
        std::vector< double > values;
        while( lineStream >> token )
        {
            values.push_back( parseFortranDouble( token ) );
        }
        BOOST_REQUIRE_EQUAL( values.size( ), 4u );

        const Eigen::Vector3d position = ( Eigen::Vector3d( ) << values.at( 0 ), values.at( 1 ), values.at( 2 ) ).finished( );
        const double potential = field->getGravitationalPotential( position );

        maximumRelativeError = std::max( maximumRelativeError, std::fabs( ( potential - values.at( 3 ) ) / values.at( 3 ) ) );
        numberOfPoints++;
    }

    BOOST_REQUIRE_EQUAL( numberOfPoints, 500 );
    BOOST_TEST_MESSAGE( "Gravity field potential against reference data: maximum relative error " << maximumRelativeError );
    BOOST_CHECK_SMALL( maximumRelativeError, 2.5E-9 );
}

//! Asking for the same position twice must not change the answer.
BOOST_AUTO_TEST_CASE( test_repeated_evaluation_is_consistent )
{
    const std::shared_ptr< EllipsoidalHarmonicsGravityField > field = createTestField( );

    const Eigen::Vector3d first = ( Eigen::Vector3d( ) << 3.5E3, -1.2E3, 0.9E3 ).finished( );
    const Eigen::Vector3d second = ( Eigen::Vector3d( ) << -4.0E3, 2.0E3, -1.5E3 ).finished( );

    const double firstPotential = field->getGravitationalPotential( first );
    const double secondPotential = field->getGravitationalPotential( second );

    // Going back to the first position must give the first answer again, so the cached values
    // are keyed on position correctly.
    BOOST_CHECK_CLOSE_FRACTION( field->getGravitationalPotential( first ), firstPotential, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( field->getGravitationalPotential( second ), secondPotential, 1.0E-15 );
    BOOST_CHECK( std::fabs( firstPotential - secondPotential ) > 1.0E-6 * std::fabs( firstPotential ) );
}

//! Inside the reference ellipsoid the expansion does not converge, so it must refuse to answer.
BOOST_AUTO_TEST_CASE( test_evaluation_inside_reference_ellipsoid_is_rejected )
{
    const std::shared_ptr< EllipsoidalHarmonicsGravityField > field = createTestField( );
    const Eigen::Vector3d semiAxes = getReferenceSemiAxes( );

    // Well inside.
    BOOST_CHECK_THROW( field->getGravitationalPotential( ( Eigen::Vector3d( ) << 0.0, 0.0, 0.5 * semiAxes( 2 ) ).finished( ) ),
                       std::runtime_error );

    // Just inside the end of the long axis, which is the closest approach the model allows.
    BOOST_CHECK_THROW( field->getGravitationalPotential( ( Eigen::Vector3d( ) << 0.99 * semiAxes( 0 ), 0.0, 0.0 ).finished( ) ),
                       std::runtime_error );

    // Just outside the same point must work.
    BOOST_CHECK_NO_THROW( field->getGravitationalPotential( ( Eigen::Vector3d( ) << 1.01 * semiAxes( 0 ), 0.0, 0.0 ).finished( ) ) );
}

/*! The gradient is not implemented yet and must say so.
 *
 * Without this override the base class would quietly return the gradient of a point mass,
 * which would be wrong but look reasonable.
 */
BOOST_AUTO_TEST_CASE( test_gradient_is_not_silently_wrong )
{
    const std::shared_ptr< EllipsoidalHarmonicsGravityField > field = createTestField( );

    BOOST_CHECK_THROW( field->getGradientOfPotential( ( Eigen::Vector3d( ) << 4.0E3, 1.0E3, 1.0E3 ).finished( ) ), std::runtime_error );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
