/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    The reference potentials used here were produced by an independent implementation of the
 *    same model, for comet 67P/Churyumov-Gerasimenko expanded to degree 15. They are given at
 *    500 points on the sphere that encloses the body.
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

#include "tudat/basics/testMacros.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/ellipsoidalHarmonicsFileReader.h"
#include "tudat/math/basic/ellipsoidalHarmonics.h"
#include "tudat/math/basic/legendrePolynomials.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::basic_mathematics;
using namespace tudat::input_output;

BOOST_AUTO_TEST_SUITE( test_ellipsoidal_harmonics )

//! Semi-axes of the reference ellipsoid of 67P/Churyumov-Gerasimenko, in metres.
const double longestSemiAxis = 2.802437330193671E+03;
const double middleSemiAxis = 2.090307360485069E+03;
const double shortestSemiAxis = 1.867442428762624E+03;

//! The two focal parameters derived from them, which is how the model is parametrised.
double getFirstFocalParameter( )
{
    return std::sqrt( ( longestSemiAxis - middleSemiAxis ) * ( longestSemiAxis + middleSemiAxis ) );
}

double getSecondFocalParameter( )
{
    return std::sqrt( ( longestSemiAxis - shortestSemiAxis ) * ( longestSemiAxis + shortestSemiAxis ) );
}

std::string getDataDirectory( )
{
    return paths::getTudatTestDataPath( ) + "/ellipsoidal_harmonics/";
}

//! Everything needed to evaluate the 67P test field, assembled from the three input files.
struct TestField {
    EllipsoidalHarmonicMatrices lameCoefficients;
    EllipsoidalHarmonicMatrices transformFactors;
    EllipsoidalHarmonicCoefficients fileCoefficients;
    EllipsoidalHarmonicCoefficients coefficients;
    EllipsoidalHarmonicCoefficients radialFunctionsAtSurface;
    Eigen::VectorXd quadratureNodes;
    Eigen::VectorXd quadratureWeights;
    int maximumDegree;
};

TestField loadTestField( const int maximumDegree = 15, const int quadratureOrder = 24 )
{
    TestField field;
    field.maximumDegree = maximumDegree;

    field.lameCoefficients = readLameFunctionCoefficients( getDataDirectory( ) + "67p_lame_coefficients.txt", maximumDegree );
    field.transformFactors = readEllipsoidalHarmonicTransformFactors( getDataDirectory( ) + "67p_transform_factors.txt", maximumDegree );
    field.fileCoefficients =
            readEllipsoidalHarmonicCoefficients( getDataDirectory( ) + "67p_ellipsoidal_harmonic_coefficients.txt", maximumDegree );

    computeGaussLegendreNodesAndWeights( quadratureOrder, field.quadratureNodes, field.quadratureWeights );

    // The file coefficients are normalised to the surface of the reference ellipsoid, which is
    // the surface of constant ellipsoidal radius equal to the longest semi-axis.
    field.radialFunctionsAtSurface = createEllipsoidalHarmonicCoefficients( maximumDegree );
    computeRadialFunctions( field.lameCoefficients,
                            longestSemiAxis,
                            getFirstFocalParameter( ),
                            getSecondFocalParameter( ),
                            field.quadratureNodes,
                            field.quadratureWeights,
                            field.radialFunctionsAtSurface );

    field.coefficients = normaliseCoefficientsToReferenceSurface( field.fileCoefficients, field.radialFunctionsAtSurface );

    return field;
}

double evaluatePotential( const TestField& field, const Eigen::Vector3d& position )
{
    const double h = getFirstFocalParameter( );
    const double k = getSecondFocalParameter( );

    const Eigen::Vector3d ellipsoidal = convertCartesianToEllipsoidalCoordinates( position, h, k );

    EllipsoidalHarmonicCoefficients radialFunctions = createEllipsoidalHarmonicCoefficients( field.maximumDegree );
    computeRadialFunctions(
            field.lameCoefficients, ellipsoidal( 0 ), h, k, field.quadratureNodes, field.quadratureWeights, radialFunctions );

    LegendreCache legendreCache( field.maximumDegree, field.maximumDegree, true );
    legendreCache.setComputeFirstDerivatives( false );
    legendreCache.setComputeSecondDerivatives( false );
    legendreCache.update( std::sin( ellipsoidal( 1 ) ) );

    return sumEllipsoidalHarmonicTerms( field.coefficients, radialFunctions, field.transformFactors, legendreCache, ellipsoidal( 2 ) );
}

Eigen::Vector3d evaluateGradient( const TestField& field, const Eigen::Vector3d& position )
{
    const double h = getFirstFocalParameter( );
    const double k = getSecondFocalParameter( );

    const Eigen::Vector3d ellipsoidal = convertCartesianToEllipsoidalCoordinates( position, h, k );

    EllipsoidalHarmonicCoefficients radialFunctions = createEllipsoidalHarmonicCoefficients( field.maximumDegree );
    EllipsoidalHarmonicCoefficients radialFunctionDerivatives = createEllipsoidalHarmonicCoefficients( field.maximumDegree );
    computeRadialFunctionsAndDerivatives( field.lameCoefficients,
                                          ellipsoidal( 0 ),
                                          h,
                                          k,
                                          field.quadratureNodes,
                                          field.quadratureWeights,
                                          radialFunctions,
                                          radialFunctionDerivatives );

    // The gradient needs the Legendre derivatives, unlike the potential.
    LegendreCache legendreCache( field.maximumDegree, field.maximumDegree, true );
    legendreCache.setComputeFirstDerivatives( true );
    legendreCache.setComputeSecondDerivatives( false );
    legendreCache.update( std::sin( ellipsoidal( 1 ) ) );

    return computeEllipsoidalHarmonicPotentialGradient(
            field.coefficients, radialFunctions, radialFunctionDerivatives, field.transformFactors, legendreCache, ellipsoidal, h, k );
}

//! Going to ellipsoidal coordinates and back must return the original position.
BOOST_AUTO_TEST_CASE( test_coordinate_conversion_round_trip )
{
    const double h = getFirstFocalParameter( );
    const double k = getSecondFocalParameter( );

    std::vector< Eigen::Vector3d > positions;
    positions.push_back( ( Eigen::Vector3d( ) << 3.0E3, 1.0E3, -2.0E3 ).finished( ) );
    positions.push_back( ( Eigen::Vector3d( ) << -4.5E3, 2.2E3, 1.1E3 ).finished( ) );
    positions.push_back( ( Eigen::Vector3d( ) << 0.0, 0.0, 2.81E3 ).finished( ) );
    positions.push_back( ( Eigen::Vector3d( ) << 0.0, 3.5E3, 0.0 ).finished( ) );
    // Close to the x axis, where the angle gamma becomes poorly determined.
    positions.push_back( ( Eigen::Vector3d( ) << 5.0E3, 1.0, 1.0 ).finished( ) );

    for( unsigned int i = 0; i < positions.size( ); i++ )
    {
        const Eigen::Vector3d ellipsoidal = convertCartesianToEllipsoidalCoordinates( positions.at( i ), h, k );
        const Eigen::Vector3d recovered = convertEllipsoidalToCartesianCoordinates( ellipsoidal, h, k );

        // The ellipsoidal radius must lie outside the focal parameters for the point to be
        // outside the reference ellipsoid at all.
        BOOST_CHECK( ellipsoidal( 0 ) > k );

        for( int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( recovered( j ) - positions.at( i )( j ), 1.0E-8 );
        }
    }
}

//! The quadrature rule must integrate polynomials of the order it claims to handle exactly.
BOOST_AUTO_TEST_CASE( test_gauss_legendre_rule )
{
    Eigen::VectorXd nodes;
    Eigen::VectorXd weights;
    computeGaussLegendreNodesAndWeights( 5, nodes, weights );

    // Weights of a rule on [-1, 1] always add up to the length of the interval.
    BOOST_CHECK_CLOSE_FRACTION( weights.sum( ), 2.0, 1.0E-14 );

    // A five node rule is exact for polynomials up to degree nine.
    for( int power = 0; power <= 9; power++ )
    {
        double integral = 0.0;
        for( int i = 0; i < nodes.size( ); i++ )
        {
            integral += weights( i ) * std::pow( nodes( i ), power );
        }

        const double expected = ( power % 2 == 0 ) ? 2.0 / ( power + 1 ) : 0.0;
        BOOST_CHECK_SMALL( integral - expected, 1.0E-14 );
    }
}

/*! The expansion needs Legendre functions without the Condon-Shortley sign.
 *
 * Tudat's geodesy-normalised functions have no such sign, which is what this model expects.
 * Getting it wrong flips the L and M groups only, and changes the potential by a few per cent
 * rather than producing anything obviously broken, so it is worth pinning down here.
 */
BOOST_AUTO_TEST_CASE( test_legendre_sign_convention )
{
    const double sineOfPhi = 0.3;

    LegendreCache cache( 2, 2, true );
    cache.update( sineOfPhi );

    BOOST_CHECK_CLOSE_FRACTION( cache.getLegendrePolynomial( 0, 0 ), 1.0, 1.0E-14 );
    BOOST_CHECK_CLOSE_FRACTION( cache.getLegendrePolynomial( 1, 0 ), std::sqrt( 3.0 ) * sineOfPhi, 1.0E-14 );
    BOOST_CHECK_CLOSE_FRACTION( cache.getLegendrePolynomial( 1, 1 ), std::sqrt( 3.0 - 3.0 * sineOfPhi * sineOfPhi ), 1.0E-14 );
}

//! Term counts of the four groups must add up to the 2*degree+1 terms of each degree.
BOOST_AUTO_TEST_CASE( test_term_group_bookkeeping )
{
    for( int degree = 0; degree <= 15; degree++ )
    {
        BOOST_CHECK_EQUAL(
                getNumberOfKTerms( degree ) + getNumberOfLTerms( degree ) + getNumberOfMTerms( degree ) + getNumberOfNTerms( degree ),
                2 * degree + 1 );

        // Every order index must map to exactly one slot, and the slots must be used once each.
        std::vector< int > timesUsed( 2 * degree + 1, 0 );
        for( int orderIndex = 1; orderIndex <= 2 * degree + 1; orderIndex++ )
        {
            EllipsoidalHarmonicTermGroup group;
            int indexInGroup;
            getTermGroupAndIndex( degree, orderIndex, group, indexInGroup );

            BOOST_CHECK( indexInGroup >= 0 );
            BOOST_CHECK( indexInGroup < getNumberOfTermsInGroup( degree, group ) );

            int offset = 0;
            for( int earlier = 0; earlier < static_cast< int >( group ); earlier++ )
            {
                offset += getNumberOfTermsInGroup( degree, static_cast< EllipsoidalHarmonicTermGroup >( earlier ) );
            }
            timesUsed.at( offset + indexInGroup )++;
        }

        for( unsigned int i = 0; i < timesUsed.size( ); i++ )
        {
            BOOST_CHECK_EQUAL( timesUsed.at( i ), 1 );
        }
    }
}

//! Values taken from the three input files by hand, to pin down the file layouts.
BOOST_AUTO_TEST_CASE( test_input_file_reading )
{
    // Fortran style exponents must read the same as normal ones.
    BOOST_CHECK_CLOSE_FRACTION( parseFortranDouble( "0.297750446271375D+00" ), 0.297750446271375, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( parseFortranDouble( "-0.254717206845257D-02" ), -0.254717206845257E-2, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( parseFortranDouble( "0.100000000000000E+01" ), 1.0, 1.0E-15 );

    const EllipsoidalHarmonicCoefficients coefficients =
            readEllipsoidalHarmonicCoefficients( getDataDirectory( ) + "67p_ellipsoidal_harmonic_coefficients.txt", 15 );

    BOOST_CHECK_EQUAL( coefficients.getMaximumDegree( ), 15 );

    // Degree zero has a single term, the one that carries most of the field.
    BOOST_CHECK_CLOSE_FRACTION( coefficients.kTerms.at( 0 )( 0 ), 0.297750446271375, 1.0E-14 );

    // Degree one splits its three values across the K, L and M groups in that order.
    BOOST_CHECK_CLOSE_FRACTION( coefficients.kTerms.at( 1 )( 0 ), 0.185269531533431E-2, 1.0E-14 );
    BOOST_CHECK_CLOSE_FRACTION( coefficients.lTerms.at( 1 )( 0 ), -0.254717206845257E-2, 1.0E-14 );
    BOOST_CHECK_CLOSE_FRACTION( coefficients.mTerms.at( 1 )( 0 ), -0.331405394064142E-2, 1.0E-14 );

    // Converting a polynomial from powers of (1-x) to powers of x:
    // 2 + 4(1-x) + 3(1-x)^2 is the same polynomial as 9 - 10x + 3x^2.
    Eigen::VectorXd inPowersOfOneMinusX( 3 );
    inPowersOfOneMinusX << 2.0, 4.0, 3.0;
    const Eigen::VectorXd inPowersOfX = convertPolynomialToPowersOfArgument( inPowersOfOneMinusX );
    BOOST_CHECK_CLOSE_FRACTION( inPowersOfX( 0 ), 9.0, 1.0E-14 );
    BOOST_CHECK_CLOSE_FRACTION( inPowersOfX( 1 ), -10.0, 1.0E-14 );
    BOOST_CHECK_CLOSE_FRACTION( inPowersOfX( 2 ), 3.0, 1.0E-14 );

    const EllipsoidalHarmonicMatrices lameCoefficients =
            readLameFunctionCoefficients( getDataDirectory( ) + "67p_lame_coefficients.txt", 15 );

    // Degree two has two K terms, so a two by two matrix of polynomials. The values below are
    // the file's two lines after conversion to powers of (rho/h)^2.
    BOOST_CHECK_EQUAL( lameCoefficients.kTerms.at( 2 ).rows( ), 2 );
    BOOST_CHECK_EQUAL( lameCoefficients.kTerms.at( 2 ).cols( ), 2 );
    BOOST_CHECK_CLOSE_FRACTION( lameCoefficients.kTerms.at( 2 )( 0, 0 ), 0.73547651534661, 1.0E-12 );
    BOOST_CHECK_CLOSE_FRACTION( lameCoefficients.kTerms.at( 2 )( 1, 0 ), -1.99597131372873, 1.0E-12 );
    BOOST_CHECK_CLOSE_FRACTION( lameCoefficients.kTerms.at( 2 )( 0, 1 ), -4.15492427553361, 1.0E-12 );
    BOOST_CHECK_CLOSE_FRACTION( lameCoefficients.kTerms.at( 2 )( 1, 1 ), 3.66518380400427, 1.0E-12 );

    const EllipsoidalHarmonicMatrices transformFactors =
            readEllipsoidalHarmonicTransformFactors( getDataDirectory( ) + "67p_transform_factors.txt", 15 );

    // Transform factors are stored a row per term, unconverted.
    BOOST_CHECK_CLOSE_FRACTION( transformFactors.kTerms.at( 2 )( 0, 0 ), 0.995385926022259, 1.0E-14 );
    BOOST_CHECK_CLOSE_FRACTION( transformFactors.kTerms.at( 2 )( 0, 1 ), -0.959523750452055E-1, 1.0E-14 );
    BOOST_CHECK_CLOSE_FRACTION( transformFactors.kTerms.at( 2 )( 1, 0 ), -0.959523750452093E-1, 1.0E-14 );
    BOOST_CHECK_CLOSE_FRACTION( transformFactors.kTerms.at( 2 )( 1, 1 ), -0.995385926022254, 1.0E-14 );

    // Degree one has no N terms at all.
    BOOST_CHECK_EQUAL( transformFactors.nTerms.at( 1 ).rows( ), 0 );
}

//! The potential must reproduce the independently computed reference values.
BOOST_AUTO_TEST_CASE( test_potential_matches_reference_data )
{
    const TestField field = loadTestField( );
    const int maximumDegree = field.maximumDegree;
    const double h = getFirstFocalParameter( );
    const double k = getSecondFocalParameter( );

    std::ifstream referenceFile( getDataDirectory( ) + "67p_reference_potentials.txt" );
    BOOST_REQUIRE( referenceFile.good( ) );

    EllipsoidalHarmonicCoefficients radialFunctions = createEllipsoidalHarmonicCoefficients( maximumDegree );
    LegendreCache legendreCache( maximumDegree, maximumDegree, true );
    legendreCache.setComputeFirstDerivatives( false );
    legendreCache.setComputeSecondDerivatives( false );

    double maximumRelativeError = 0.0;
    double sumOfSquaredRelativeErrors = 0.0;
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
        const double referencePotential = values.at( 3 );

        const Eigen::Vector3d ellipsoidal = convertCartesianToEllipsoidalCoordinates( position, h, k );
        computeRadialFunctions(
                field.lameCoefficients, ellipsoidal( 0 ), h, k, field.quadratureNodes, field.quadratureWeights, radialFunctions );
        legendreCache.update( std::sin( ellipsoidal( 1 ) ) );

        const double potential =
                sumEllipsoidalHarmonicTerms( field.coefficients, radialFunctions, field.transformFactors, legendreCache, ellipsoidal( 2 ) );

        const double relativeError = ( potential - referencePotential ) / referencePotential;
        maximumRelativeError = std::max( maximumRelativeError, std::fabs( relativeError ) );
        sumOfSquaredRelativeErrors += relativeError * relativeError;
        numberOfPoints++;
    }

    BOOST_REQUIRE_EQUAL( numberOfPoints, 500 );

    const double rootMeanSquareRelativeError = std::sqrt( sumOfSquaredRelativeErrors / numberOfPoints );
    BOOST_TEST_MESSAGE( "Ellipsoidal harmonic potential against reference data: maximum relative error "
                        << maximumRelativeError << ", root mean square " << rootMeanSquareRelativeError );

    // The remaining error comes from the precision of the input coefficients, not from the
    // quadrature: integrating more finely does not reduce it.
    BOOST_CHECK_SMALL( maximumRelativeError, 2.5E-9 );
    BOOST_CHECK_SMALL( rootMeanSquareRelativeError, 2.0E-10 );

    // The convenience function must agree with the step by step evaluation above.
    const Eigen::Vector3d testPosition = ( Eigen::Vector3d( ) << 2.0E3, -1.5E3, 1.2E3 ).finished( );
    const double stepByStep = evaluatePotential( field, testPosition );
    const double inOneCall = computeEllipsoidalHarmonicGravitationalPotential(
            testPosition, h, k, field.coefficients, field.lameCoefficients, field.transformFactors, 24 );

    BOOST_CHECK_CLOSE_FRACTION( inOneCall, stepByStep, 1.0E-14 );
}

//! Raising the quadrature order must not change the potential, since order 24 already converged.
BOOST_AUTO_TEST_CASE( test_quadrature_order_is_sufficient )
{
    const int maximumDegree = 15;
    const double h = getFirstFocalParameter( );
    const double k = getSecondFocalParameter( );

    const EllipsoidalHarmonicMatrices lameCoefficients =
            readLameFunctionCoefficients( getDataDirectory( ) + "67p_lame_coefficients.txt", maximumDegree );

    EllipsoidalHarmonicCoefficients atOrder24 = createEllipsoidalHarmonicCoefficients( maximumDegree );
    EllipsoidalHarmonicCoefficients atOrder48 = createEllipsoidalHarmonicCoefficients( maximumDegree );

    Eigen::VectorXd nodes24, weights24, nodes48, weights48;
    computeGaussLegendreNodesAndWeights( 24, nodes24, weights24 );
    computeGaussLegendreNodesAndWeights( 48, nodes48, weights48 );

    computeRadialFunctions( lameCoefficients, longestSemiAxis, h, k, nodes24, weights24, atOrder24 );
    computeRadialFunctions( lameCoefficients, longestSemiAxis, h, k, nodes48, weights48, atOrder48 );

    for( int degree = 0; degree <= maximumDegree; degree++ )
    {
        for( int groupIndex = 0; groupIndex < 4; groupIndex++ )
        {
            const EllipsoidalHarmonicTermGroup group = static_cast< EllipsoidalHarmonicTermGroup >( groupIndex );
            const Eigen::VectorXd& coarse = atOrder24.getGroup( degree, group );
            const Eigen::VectorXd& fine = atOrder48.getGroup( degree, group );

            for( int term = 0; term < coarse.size( ); term++ )
            {
                BOOST_CHECK_CLOSE_FRACTION( coarse( term ), fine( term ), 1.0E-10 );
            }
        }
    }
}

/*! Far from the body the field must look like a point mass.
 *
 * This checks that the coefficients were un-normalised correctly: the potential times the
 * distance has to settle on the degree zero coefficient divided by its own radial part at the
 * reference surface. For this field that works out at about 669 m^3/s^2, the expected mass
 * parameter of 67P.
 *
 * It settles slowly, as one over distance, because this expansion is not centred on the centre
 * of mass and so has a non-zero degree one term. That is what the decay check below pins down.
 */
BOOST_AUTO_TEST_CASE( test_far_field_approaches_point_mass )
{
    const TestField field = loadTestField( );

    const double expectedGravitationalParameter =
            field.fileCoefficients.kTerms.at( 0 )( 0 ) / field.radialFunctionsAtSurface.kTerms.at( 0 )( 0 );

    BOOST_TEST_MESSAGE( "Gravitational parameter implied by the degree zero coefficient: " << expectedGravitationalParameter
                                                                                           << " m^3/s^2" );

    // Loose sanity bound, so an accidental change of normalisation would be caught.
    BOOST_CHECK_CLOSE_FRACTION( expectedGravitationalParameter, 669.0, 1.0E-2 );

    std::vector< Eigen::Vector3d > directions;
    directions.push_back( ( Eigen::Vector3d( ) << 1.0, 0.0, 0.0 ).finished( ) );
    directions.push_back( ( Eigen::Vector3d( ) << 0.0, 1.0, 0.0 ).finished( ) );
    directions.push_back( ( Eigen::Vector3d( ) << 0.0, 0.0, 1.0 ).finished( ) );
    directions.push_back( ( Eigen::Vector3d( ) << 1.0, 1.0, 1.0 ).finished( ) );
    directions.push_back( ( Eigen::Vector3d( ) << -2.0, 1.0, -3.0 ).finished( ) );

    // Largest deviation of potential times distance from the point mass value, over all
    // directions, at a given distance.
    auto getLargestDeviation = [ & ]( const double distance ) {
        double largest = 0.0;
        for( unsigned int i = 0; i < directions.size( ); i++ )
        {
            const Eigen::Vector3d position = distance * directions.at( i ).normalized( );
            const double potential = evaluatePotential( field, position );
            largest = std::max( largest, std::fabs( distance * potential / expectedGravitationalParameter - 1.0 ) );
        }
        return largest;
    };

    const double deviationAtThousandSemiAxes = getLargestDeviation( 1.0E3 * longestSemiAxis );
    const double deviationAtTenThousandSemiAxes = getLargestDeviation( 1.0E4 * longestSemiAxis );

    BOOST_TEST_MESSAGE( "Deviation from a point mass: " << deviationAtThousandSemiAxes << " at 1e3 semi-axes, "
                                                        << deviationAtTenThousandSemiAxes << " at 1e4" );

    BOOST_CHECK_SMALL( deviationAtThousandSemiAxes, 1.0E-4 );
    BOOST_CHECK_SMALL( deviationAtTenThousandSemiAxes, 1.0E-5 );

    // Ten times further out must be about ten times closer to a point mass.
    BOOST_CHECK( deviationAtTenThousandSemiAxes < 0.2 * deviationAtThousandSemiAxes );

    // Close in, the higher degree terms still matter, so the field is not simply a point mass.
    BOOST_CHECK( getLargestDeviation( 1.2 * longestSemiAxis ) > 1.0E-3 );
}

//! The radial parts computed alongside their derivatives must match the value-only routine.
BOOST_AUTO_TEST_CASE( test_radial_functions_agree_with_value_only_routine )
{
    const TestField field = loadTestField( );
    const double h = getFirstFocalParameter( );
    const double k = getSecondFocalParameter( );
    const double rho = 1.4 * longestSemiAxis;

    EllipsoidalHarmonicCoefficients valuesOnly = createEllipsoidalHarmonicCoefficients( field.maximumDegree );
    computeRadialFunctions( field.lameCoefficients, rho, h, k, field.quadratureNodes, field.quadratureWeights, valuesOnly );

    EllipsoidalHarmonicCoefficients values = createEllipsoidalHarmonicCoefficients( field.maximumDegree );
    EllipsoidalHarmonicCoefficients derivatives = createEllipsoidalHarmonicCoefficients( field.maximumDegree );
    computeRadialFunctionsAndDerivatives(
            field.lameCoefficients, rho, h, k, field.quadratureNodes, field.quadratureWeights, values, derivatives );

    for( int degree = 0; degree <= field.maximumDegree; degree++ )
    {
        for( int groupIndex = 0; groupIndex < 4; groupIndex++ )
        {
            const EllipsoidalHarmonicTermGroup group = static_cast< EllipsoidalHarmonicTermGroup >( groupIndex );
            const Eigen::VectorXd& a = valuesOnly.getGroup( degree, group );
            const Eigen::VectorXd& b = values.getGroup( degree, group );
            for( int term = 0; term < a.size( ); term++ )
            {
                BOOST_CHECK_CLOSE_FRACTION( a( term ), b( term ), 1.0E-14 );
            }
        }
    }
}

//! The derivative of each radial part must match a finite difference of the radial part.
BOOST_AUTO_TEST_CASE( test_radial_function_derivatives_match_finite_differences )
{
    const TestField field = loadTestField( );
    const double h = getFirstFocalParameter( );
    const double k = getSecondFocalParameter( );

    const double rho = 1.5 * longestSemiAxis;
    const double step = rho * 1.0E-6;

    EllipsoidalHarmonicCoefficients values = createEllipsoidalHarmonicCoefficients( field.maximumDegree );
    EllipsoidalHarmonicCoefficients derivatives = createEllipsoidalHarmonicCoefficients( field.maximumDegree );
    computeRadialFunctionsAndDerivatives(
            field.lameCoefficients, rho, h, k, field.quadratureNodes, field.quadratureWeights, values, derivatives );

    EllipsoidalHarmonicCoefficients forward = createEllipsoidalHarmonicCoefficients( field.maximumDegree );
    EllipsoidalHarmonicCoefficients backward = createEllipsoidalHarmonicCoefficients( field.maximumDegree );
    computeRadialFunctions( field.lameCoefficients, rho + step, h, k, field.quadratureNodes, field.quadratureWeights, forward );
    computeRadialFunctions( field.lameCoefficients, rho - step, h, k, field.quadratureNodes, field.quadratureWeights, backward );

    double largestRelativeError = 0.0;
    for( int degree = 0; degree <= field.maximumDegree; degree++ )
    {
        for( int groupIndex = 0; groupIndex < 4; groupIndex++ )
        {
            const EllipsoidalHarmonicTermGroup group = static_cast< EllipsoidalHarmonicTermGroup >( groupIndex );
            const Eigen::VectorXd& analytic = derivatives.getGroup( degree, group );
            const Eigen::VectorXd& plus = forward.getGroup( degree, group );
            const Eigen::VectorXd& minus = backward.getGroup( degree, group );

            for( int term = 0; term < analytic.size( ); term++ )
            {
                const double numerical = ( plus( term ) - minus( term ) ) / ( 2.0 * step );
                if( std::fabs( numerical ) > 0.0 )
                {
                    largestRelativeError = std::max( largestRelativeError, std::fabs( analytic( term ) / numerical - 1.0 ) );
                }
            }
        }
    }

    BOOST_TEST_MESSAGE( "Radial part derivatives against finite differences: largest relative error " << largestRelativeError );
    BOOST_CHECK_SMALL( largestRelativeError, 1.0E-6 );
}

//! The coordinate Jacobian must match finite differences of the forward coordinate map.
BOOST_AUTO_TEST_CASE( test_coordinate_jacobian_matches_finite_differences )
{
    const double h = getFirstFocalParameter( );
    const double k = getSecondFocalParameter( );

    const Eigen::Vector3d ellipsoidal = ( Eigen::Vector3d( ) << 1.6 * longestSemiAxis, 0.4, -1.1 ).finished( );
    const Eigen::Matrix3d analytic = computeEllipsoidalCoordinateJacobian( ellipsoidal, h, k );

    Eigen::Matrix3d numerical;
    for( int i = 0; i < 3; i++ )
    {
        const double step = ( i == 0 ) ? ellipsoidal( 0 ) * 1.0E-7 : 1.0E-7;
        Eigen::Vector3d plus = ellipsoidal;
        Eigen::Vector3d minus = ellipsoidal;
        plus( i ) += step;
        minus( i ) -= step;
        numerical.col( i ) =
                ( convertEllipsoidalToCartesianCoordinates( plus, h, k ) - convertEllipsoidalToCartesianCoordinates( minus, h, k ) ) /
                ( 2.0 * step );
    }

    for( int i = 0; i < 3; i++ )
    {
        for( int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( ( analytic( i, j ) - numerical( i, j ) ) / analytic.cwiseAbs( ).maxCoeff( ), 1.0E-6 );
        }
    }
}

//! The analytic gradient must match finite differences of the potential.
BOOST_AUTO_TEST_CASE( test_gradient_matches_finite_differences )
{
    const TestField field = loadTestField( );

    // Directions are kept away from the long axis, where the gradient is singular.
    std::vector< Eigen::Vector3d > positions;
    positions.push_back( ( Eigen::Vector3d( ) << 1.0E3, 2.5E3, 2.0E3 ).finished( ) );
    positions.push_back( ( Eigen::Vector3d( ) << -3.0E3, 2.0E3, -1.5E3 ).finished( ) );
    positions.push_back( ( Eigen::Vector3d( ) << 0.0, 4.0E3, 0.0 ).finished( ) );
    positions.push_back( ( Eigen::Vector3d( ) << 0.0, 0.0, 3.5E3 ).finished( ) );
    positions.push_back( ( Eigen::Vector3d( ) << 6.0E3, 5.0E3, -4.0E3 ).finished( ) );
    positions.push_back( ( Eigen::Vector3d( ) << 1.5E4, -1.0E4, 8.0E3 ).finished( ) );

    double largestRelativeError = 0.0;

    for( unsigned int i = 0; i < positions.size( ); i++ )
    {
        const Eigen::Vector3d position = positions.at( i );
        const Eigen::Vector3d analytic = evaluateGradient( field, position );

        Eigen::Vector3d numerical;
        for( int j = 0; j < 3; j++ )
        {
            const double step = position.norm( ) * 1.0E-6;
            Eigen::Vector3d plus = position;
            Eigen::Vector3d minus = position;
            plus( j ) += step;
            minus( j ) -= step;
            numerical( j ) = ( evaluatePotential( field, plus ) - evaluatePotential( field, minus ) ) / ( 2.0 * step );
        }

        largestRelativeError = std::max( largestRelativeError, ( analytic - numerical ).norm( ) / numerical.norm( ) );
    }

    BOOST_TEST_MESSAGE( "Gradient against finite differences: largest relative error " << largestRelativeError );
    BOOST_CHECK_SMALL( largestRelativeError, 1.0E-7 );
}

//! Far from the body the gradient must become that of a point mass.
BOOST_AUTO_TEST_CASE( test_gradient_approaches_point_mass )
{
    const TestField field = loadTestField( );

    const double expectedGravitationalParameter =
            field.fileCoefficients.kTerms.at( 0 )( 0 ) / field.radialFunctionsAtSurface.kTerms.at( 0 )( 0 );

    const double distance = 1.0E4 * longestSemiAxis;
    const Eigen::Vector3d direction = Eigen::Vector3d( 0.3, 1.0, -0.7 ).normalized( );
    const Eigen::Vector3d position = distance * direction;

    const Eigen::Vector3d gradient = evaluateGradient( field, position );

    // The potential is positive and falls off as one over distance, so its gradient points at
    // the body and has size mu over distance squared.
    BOOST_CHECK_CLOSE_FRACTION( gradient.norm( ), expectedGravitationalParameter / ( distance * distance ), 1.0E-4 );
    BOOST_CHECK_CLOSE_FRACTION( gradient.normalized( ).dot( direction ), -1.0, 1.0E-5 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
