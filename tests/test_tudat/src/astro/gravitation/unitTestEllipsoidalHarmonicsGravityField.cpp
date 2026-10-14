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
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;

BOOST_AUTO_TEST_SUITE( test_ellipsoidal_harmonics_gravity_field )

//! Reference ellipsoid of 67P/Churyumov-Gerasimenko, longest semi-axis first, in metres.
Eigen::Vector3d getReferenceSemiAxes( )
{
    return ( Eigen::Vector3d( ) << 2.802437330193671E+03, 2.090307360485069E+03, 1.867442428762624E+03 ).finished( );
}

//! Gravitational parameter the degree zero coefficient of this data set implies, in m^3/s^2.
const double expectedGravitationalParameter = 668.78268559825256;

std::string getDataDirectory( )
{
    return paths::getTudatTestDataPath( ) + "/ellipsoidal_harmonics/";
}

std::shared_ptr< EllipsoidalHarmonicsGravityField > createTestField( const int maximumDegree = 15, const int radialTableGridPoints = 3200 )
{
    return std::make_shared< EllipsoidalHarmonicsGravityField >(
            getReferenceSemiAxes( ),
            readEllipsoidalHarmonicCoefficients( getDataDirectory( ) + "67p_ellipsoidal_harmonic_coefficients.txt", maximumDegree ),
            readLameFunctionCoefficients( getDataDirectory( ) + "67p_lame_coefficients.txt", maximumDegree ),
            readEllipsoidalHarmonicTransformFactors( getDataDirectory( ) + "67p_transform_factors.txt", maximumDegree ),
            24,
            radialTableGridPoints,
            0.0,
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
    // The field takes its gravitational parameter from the degree zero coefficient, so it is not
    // an independent input and must come out at the value that coefficient implies.
    BOOST_CHECK_CLOSE_FRACTION( field->getGravitationalParameter( ), expectedGravitationalParameter, 1.0E-10 );
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
    BOOST_CHECK_THROW( EllipsoidalHarmonicsGravityField( reversed, coefficients, lameCoefficients, transformFactors ), std::runtime_error );

    // A sphere, where the two focal parameters both collapse to zero.
    const Eigen::Vector3d sphere = ( Eigen::Vector3d( ) << 2.0E3, 2.0E3, 2.0E3 ).finished( );
    BOOST_CHECK_THROW( EllipsoidalHarmonicsGravityField( sphere, coefficients, lameCoefficients, transformFactors ), std::runtime_error );
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

//! The gradient must match finite differences of the potential.
BOOST_AUTO_TEST_CASE( test_gradient_matches_finite_differences )
{
    const std::shared_ptr< EllipsoidalHarmonicsGravityField > field = createTestField( );

    // Positions are kept off the long axis, where the gradient cannot be formed.
    std::vector< Eigen::Vector3d > positions;
    positions.push_back( ( Eigen::Vector3d( ) << 1.0E3, 2.5E3, 2.0E3 ).finished( ) );
    positions.push_back( ( Eigen::Vector3d( ) << -3.0E3, 2.0E3, -1.5E3 ).finished( ) );
    positions.push_back( ( Eigen::Vector3d( ) << 0.0, 0.0, 3.5E3 ).finished( ) );
    positions.push_back( ( Eigen::Vector3d( ) << 8.0E3, 5.0E3, -4.0E3 ).finished( ) );

    double largestRelativeError = 0.0;

    for( unsigned int i = 0; i < positions.size( ); i++ )
    {
        const Eigen::Vector3d position = positions.at( i );
        const Eigen::Vector3d analytic = field->getGradientOfPotential( position );

        Eigen::Vector3d numerical;
        for( int j = 0; j < 3; j++ )
        {
            const double step = position.norm( ) * 1.0E-6;
            Eigen::Vector3d plus = position;
            Eigen::Vector3d minus = position;
            plus( j ) += step;
            minus( j ) -= step;
            numerical( j ) = ( field->getGravitationalPotential( plus ) - field->getGravitationalPotential( minus ) ) / ( 2.0 * step );
        }

        largestRelativeError = std::max( largestRelativeError, ( analytic - numerical ).norm( ) / numerical.norm( ) );
    }

    BOOST_TEST_MESSAGE( "Gravity field gradient against finite differences: largest relative error " << largestRelativeError );
    BOOST_CHECK_SMALL( largestRelativeError, 1.0E-7 );
}

//! Mixing potential and gradient calls at the same position must not disturb either.
BOOST_AUTO_TEST_CASE( test_potential_and_gradient_can_be_interleaved )
{
    const std::shared_ptr< EllipsoidalHarmonicsGravityField > field = createTestField( );
    const Eigen::Vector3d position = ( Eigen::Vector3d( ) << 2.0E3, 3.0E3, -1.0E3 ).finished( );

    const double potential = field->getGravitationalPotential( position );
    const Eigen::Vector3d gradient = field->getGradientOfPotential( position );

    BOOST_CHECK_CLOSE_FRACTION( field->getGravitationalPotential( position ), potential, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( field->getGradientOfPotential( position )( 0 ), gradient( 0 ), 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( field->getGravitationalPotential( position ), potential, 1.0E-15 );
}

/*! The gradient must work on the long axis too, where the coordinates are degenerate.
 *
 * The field is perfectly smooth there, but these coordinates are not, so the analytic gradient
 * cannot be formed. Orbits in a plane containing the long axis cross it twice per revolution,
 * so rather than give up, the field differences the potential close to the axis.
 */
BOOST_AUTO_TEST_CASE( test_gradient_works_on_the_long_axis )
{
    const std::shared_ptr< EllipsoidalHarmonicsGravityField > field = createTestField( );
    const double distance = 1.5 * getReferenceSemiAxes( )( 0 );
    const Eigen::Vector3d onAxis = ( Eigen::Vector3d( ) << distance, 0.0, 0.0 ).finished( );

    BOOST_CHECK_NO_THROW( field->getGravitationalPotential( onAxis ) );
    BOOST_CHECK_NO_THROW( field->getGradientOfPotential( onAxis ) );

    const Eigen::Vector3d gradient = field->getGradientOfPotential( onAxis );

    // On the axis, by symmetry of neither the body nor the expansion, the gradient still has to
    // point back towards the body and be of the size the gravitational parameter implies.
    BOOST_CHECK( gradient( 0 ) < 0.0 );
    BOOST_CHECK_CLOSE_FRACTION( gradient.norm( ), field->getGravitationalParameter( ) / ( distance * distance ), 0.5 );

    // Independent check with a different step size, so this tests the derivative rather than
    // one particular way of taking it.
    const double step = distance * 1.0E-4;
    Eigen::Vector3d independent;
    for( int i = 0; i < 3; i++ )
    {
        Eigen::Vector3d forward = onAxis;
        Eigen::Vector3d backward = onAxis;
        forward( i ) += step;
        backward( i ) -= step;
        independent( i ) = ( field->getGravitationalPotential( forward ) - field->getGravitationalPotential( backward ) ) / ( 2.0 * step );
    }

    BOOST_CHECK_SMALL( ( gradient - independent ).norm( ) / independent.norm( ), 1.0E-6 );
}

/*! At the switchover the two ways of getting the gradient must agree.
 *
 * This is what makes the switch invisible: the region is wide enough that the analytic gradient
 * is still accurate at its edge, so nothing jumps as a trajectory crosses it.
 */
BOOST_AUTO_TEST_CASE( test_gradient_agrees_at_the_long_axis_switchover )
{
    const std::shared_ptr< EllipsoidalHarmonicsGravityField > field = createTestField( );
    const double distance = 2.0 * getReferenceSemiAxes( )( 0 );

    // The tolerance is the cosine of the angle from the long axis, so this puts the point right
    // at the edge of the region where the potential is differenced.
    const double angle = field->getLongAxisToleranceForGradient( );
    const Eigen::Vector3d position = ( Eigen::Vector3d( ) << distance * std::cos( angle ), distance * std::sin( angle ), 0.0 ).finished( );

    field->resetLongAxisToleranceForGradient( 0.0 );
    const Eigen::Vector3d analytic = field->getGradientOfPotential( position );

    field->resetLongAxisToleranceForGradient( 1.0 );
    const Eigen::Vector3d differenced = field->getGradientOfPotential( position );

    BOOST_CHECK_SMALL( ( analytic - differenced ).norm( ) / analytic.norm( ), 1.0E-8 );
}

//! The gradient must also refuse to answer inside the reference ellipsoid.
BOOST_AUTO_TEST_CASE( test_gradient_inside_reference_ellipsoid_is_rejected )
{
    const std::shared_ptr< EllipsoidalHarmonicsGravityField > field = createTestField( );

    BOOST_CHECK_THROW( field->getGradientOfPotential( ( Eigen::Vector3d( ) << 0.0, 0.0, 0.5 * getReferenceSemiAxes( )( 2 ) ).finished( ) ),
                       std::runtime_error );
}

/*! The Hessian must match finite differences of the gradient, and obey Laplace's equation.
 *
 * Outside the body the potential is harmonic, so the Hessian has zero trace. Nothing in the
 * implementation enforces that, which makes it a strong independent check.
 */
BOOST_AUTO_TEST_CASE( test_hessian_matches_finite_differences_and_is_traceless )
{
    const std::shared_ptr< EllipsoidalHarmonicsGravityField > field = createTestField( );

    std::vector< Eigen::Vector3d > positions;
    positions.push_back( ( Eigen::Vector3d( ) << 1.0E3, 2.5E3, 2.0E3 ).finished( ) );
    positions.push_back( ( Eigen::Vector3d( ) << -3.0E3, 2.0E3, -1.5E3 ).finished( ) );
    positions.push_back( ( Eigen::Vector3d( ) << 8.0E3, 5.0E3, -4.0E3 ).finished( ) );

    double largestRelativeError = 0.0;
    double largestScaledTrace = 0.0;

    for( unsigned int i = 0; i < positions.size( ); i++ )
    {
        const Eigen::Vector3d position = positions.at( i );
        const Eigen::Matrix3d analytic = field->getHessianOfPotential( position );

        largestScaledTrace = std::max( largestScaledTrace, std::fabs( analytic.trace( ) ) / analytic.cwiseAbs( ).maxCoeff( ) );

        // The Laplacian is just the trace, so it must agree.
        BOOST_CHECK_CLOSE_FRACTION( field->getLaplacianOfPotential( position ), analytic.trace( ), 1.0E-12 );

        Eigen::Matrix3d numerical;
        for( int j = 0; j < 3; j++ )
        {
            const double step = position.norm( ) * 1.0E-6;
            Eigen::Vector3d plus = position;
            Eigen::Vector3d minus = position;
            plus( j ) += step;
            minus( j ) -= step;
            numerical.col( j ) = ( field->getGradientOfPotential( plus ) - field->getGradientOfPotential( minus ) ) / ( 2.0 * step );
        }

        largestRelativeError =
                std::max( largestRelativeError, ( analytic - numerical ).cwiseAbs( ).maxCoeff( ) / numerical.cwiseAbs( ).maxCoeff( ) );
    }

    BOOST_TEST_MESSAGE( "Gravity field Hessian: largest relative error " << largestRelativeError << ", largest scaled trace "
                                                                         << largestScaledTrace );
    BOOST_CHECK_SMALL( largestRelativeError, 1.0E-5 );
    BOOST_CHECK_SMALL( largestScaledTrace, 1.0E-6 );
}

//! The Hessian must work on the long axis too, by the same route the gradient takes there.
BOOST_AUTO_TEST_CASE( test_hessian_works_on_the_long_axis )
{
    const std::shared_ptr< EllipsoidalHarmonicsGravityField > field = createTestField( );
    const Eigen::Vector3d onAxis = ( Eigen::Vector3d( ) << 1.5 * getReferenceSemiAxes( )( 0 ), 0.0, 0.0 ).finished( );

    BOOST_CHECK_NO_THROW( field->getHessianOfPotential( onAxis ) );

    const Eigen::Matrix3d hessian = field->getHessianOfPotential( onAxis );
    BOOST_CHECK_SMALL( ( hessian - hessian.transpose( ) ).cwiseAbs( ).maxCoeff( ) / hessian.cwiseAbs( ).maxCoeff( ), 1.0E-12 );
    BOOST_CHECK_SMALL( std::fabs( hessian.trace( ) ) / hessian.cwiseAbs( ).maxCoeff( ), 1.0E-4 );
}

/*! Asking for the potential, gradient and Hessian in any order must give the same answers.
 *
 * They share cached intermediates that are keyed on position and on how much has been computed,
 * so an order that lets one of them skip its own update would otherwise read values that were
 * never filled in.
 */
BOOST_AUTO_TEST_CASE( test_results_do_not_depend_on_the_order_of_the_calls )
{
    const Eigen::Vector3d position = ( Eigen::Vector3d( ) << 5.0E3, 4.0E3, -3.0E3 ).finished( );

    // Reference values, each from a field that has been asked nothing else.
    const double potential = createTestField( )->getGravitationalPotential( position );
    const Eigen::Vector3d gradient = createTestField( )->getGradientOfPotential( position );
    const Eigen::Matrix3d hessian = createTestField( )->getHessianOfPotential( position );

    BOOST_CHECK( gradient.norm( ) > 0.0 );
    BOOST_CHECK( hessian.cwiseAbs( ).maxCoeff( ) > 0.0 );

    // Hessian first, then the gradient: the gradient must not come back empty.
    {
        const std::shared_ptr< EllipsoidalHarmonicsGravityField > field = createTestField( );
        field->getHessianOfPotential( position );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( field->getGradientOfPotential( position ), gradient, 1.0E-14 );
        BOOST_CHECK_CLOSE_FRACTION( field->getGravitationalPotential( position ), potential, 1.0E-14 );
    }

    // Gradient first, then the Hessian.
    {
        const std::shared_ptr< EllipsoidalHarmonicsGravityField > field = createTestField( );
        field->getGradientOfPotential( position );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( field->getHessianOfPotential( position ), hessian, 1.0E-14 );
    }

    // Everything, twice, in a different order each time.
    {
        const std::shared_ptr< EllipsoidalHarmonicsGravityField > field = createTestField( );
        field->getGravitationalPotential( position );
        field->getHessianOfPotential( position );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( field->getGradientOfPotential( position ), gradient, 1.0E-14 );
        field->getGravitationalPotential( position );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( field->getHessianOfPotential( position ), hessian, 1.0E-14 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( field->getGradientOfPotential( position ), gradient, 1.0E-14 );
    }
}

/*! Tabulating the radial parts must not change the answers.
 *
 * The radial parts depend on the ellipsoidal radius alone, so they are tabulated against it and
 * interpolated, which is what makes repeated evaluation cheap. That is an approximation, so it
 * is checked here against integrating them exactly, for the potential, the gradient and the
 * Hessian alike.
 */
BOOST_AUTO_TEST_CASE( test_tabulated_radial_parts_agree_with_exact_integration )
{
    const std::shared_ptr< EllipsoidalHarmonicsGravityField > tabulated = createTestField( 15, 3200 );
    const std::shared_ptr< EllipsoidalHarmonicsGravityField > exact = createTestField( 15, 0 );

    BOOST_CHECK( tabulated->getRadialTable( ) != nullptr );
    BOOST_CHECK( exact->getRadialTable( ) == nullptr );

    std::vector< Eigen::Vector3d > positions;
    positions.push_back( ( Eigen::Vector3d( ) << 1.0E3, 2.5E3, 2.0E3 ).finished( ) );
    positions.push_back( ( Eigen::Vector3d( ) << -3.0E3, 2.0E3, -1.5E3 ).finished( ) );
    positions.push_back( ( Eigen::Vector3d( ) << 0.0, 0.0, 3.5E3 ).finished( ) );
    positions.push_back( ( Eigen::Vector3d( ) << 8.0E3, 5.0E3, -4.0E3 ).finished( ) );
    positions.push_back( ( Eigen::Vector3d( ) << 4.0E4, -2.0E4, 3.0E4 ).finished( ) );

    double worstPotential = 0.0;
    double worstGradient = 0.0;
    double worstHessian = 0.0;

    for( unsigned int i = 0; i < positions.size( ); i++ )
    {
        const Eigen::Vector3d position = positions.at( i );

        const double exactPotential = exact->getGravitationalPotential( position );
        worstPotential = std::max( worstPotential, std::fabs( tabulated->getGravitationalPotential( position ) / exactPotential - 1.0 ) );

        const Eigen::Vector3d exactGradient = exact->getGradientOfPotential( position );
        worstGradient = std::max( worstGradient,
                                  ( tabulated->getGradientOfPotential( position ) - exactGradient ).norm( ) / exactGradient.norm( ) );

        const Eigen::Matrix3d exactHessian = exact->getHessianOfPotential( position );
        worstHessian = std::max( worstHessian,
                                 ( tabulated->getHessianOfPotential( position ) - exactHessian ).cwiseAbs( ).maxCoeff( ) /
                                         exactHessian.cwiseAbs( ).maxCoeff( ) );
    }

    BOOST_TEST_MESSAGE( "Tabulated against exact: potential " << worstPotential << ", gradient " << worstGradient << ", Hessian "
                                                              << worstHessian );

    // The table is built to be no worse than the precision of the input coefficients, which is
    // about two parts in a billion.
    BOOST_CHECK_SMALL( worstPotential, 1.0E-10 );
    BOOST_CHECK_SMALL( worstGradient, 1.0E-8 );
    BOOST_CHECK_SMALL( worstHessian, 1.0E-5 );
}

//! Beyond the tabulated range the radial parts must fall back on exact integration.
BOOST_AUTO_TEST_CASE( test_radius_beyond_the_table_falls_back_on_exact_integration )
{
    const std::shared_ptr< EllipsoidalHarmonicsGravityField > tabulated = createTestField( 15, 3200 );
    const std::shared_ptr< EllipsoidalHarmonicsGravityField > exact = createTestField( 15, 0 );

    const double beyond = 2.0 * tabulated->getRadialTable( )->getLargestRadius( );
    const Eigen::Vector3d position = beyond * Eigen::Vector3d( 0.3, 1.0, -0.7 ).normalized( );

    BOOST_CHECK( !tabulated->getRadialTable( )->isRadiusInRange( beyond ) );

    // Falling back means the two must agree to the last bit, not merely to the table's accuracy.
    BOOST_CHECK_CLOSE_FRACTION( tabulated->getGravitationalPotential( position ), exact->getGravitationalPotential( position ), 1.0E-14 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
