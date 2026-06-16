/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/gravitation/fourthDegreeFullTwoBodyGravitationalTorque.h"

#include <array>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "tudat/math/basic/sphericalHarmonics.h"

namespace tudat
{

namespace gravitation
{

namespace
{

bool isFigureFigureDebugEnabled( )
{
    const char* debugFlag = std::getenv( "TUDAT_FFDBG" );
    const char* debugCase = std::getenv( "TUDAT_FFDBG_CASE" );
    return debugFlag != nullptr && std::string( debugFlag ) != "0" && debugCase != nullptr;
}

std::string getFigureFigureDebugCase( )
{
    const char* debugCase = std::getenv( "TUDAT_FFDBG_CASE" );
    return ( debugCase == nullptr ) ? "unset" : std::string( debugCase );
}

std::string formatDebugVector( const Eigen::Vector3d& vector )
{
    std::ostringstream stream;
    stream << std::scientific << std::setprecision( 17 ) << vector( 0 ) << "," << vector( 1 ) << "," << vector( 2 );
    return stream.str( );
}

std::string formatDebugMatrix( const Eigen::Matrix3d& matrix )
{
    std::ostringstream stream;
    stream << std::scientific << std::setprecision( 17 ) << matrix( 0, 0 ) << "," << matrix( 0, 1 ) << "," << matrix( 0, 2 ) << ","
           << matrix( 1, 0 ) << "," << matrix( 1, 1 ) << "," << matrix( 1, 2 ) << "," << matrix( 2, 0 ) << "," << matrix( 2, 1 ) << ","
           << matrix( 2, 2 );
    return stream.str( );
}

std::string formatDebugScalar( const double value )
{
    std::ostringstream stream;
    stream << std::scientific << std::setprecision( 17 ) << value;
    return stream.str( );
}

void printCanonicalFigureFigureTraceLine( const std::string& step,
                                          const std::string& name,
                                          const std::string& shape,
                                          const std::string& values )
{
    std::cout << "FFDBG|case=" << getFigureFigureDebugCase( ) << "|model=schutz|step=" << step << "|name=" << name << "|shape=" << shape
              << "|values=" << values << std::endl;
}

void printNativeFigureFigureTraceLine( const std::string& step,
                                       const std::string& name,
                                       const std::string& shape,
                                       const std::string& values )
{
    std::cout << "FFDBG_NATIVE|case=" << getFigureFigureDebugCase( ) << "|model=schutz|step=" << step << "|name=" << name
              << "|shape=" << shape << "|values=" << values << std::endl;
}

std::string formatDebugCoefficientPacket( const std::array< double, 5 >& coefficients, const double normalizationFlag )
{
    Eigen::Matrix3d coefficientPacket = Eigen::Matrix3d::Zero( );
    coefficientPacket( 0, 0 ) = coefficients.at( 0 );
    coefficientPacket( 0, 1 ) = coefficients.at( 1 );
    coefficientPacket( 0, 2 ) = coefficients.at( 2 );
    coefficientPacket( 1, 0 ) = coefficients.at( 3 );
    coefficientPacket( 1, 1 ) = coefficients.at( 4 );
    coefficientPacket( 1, 2 ) = normalizationFlag;
    return formatDebugMatrix( coefficientPacket );
}

bool parseDebugVectorFromEnvironment( const std::string& variableName, std::vector< double >& parsedValues )
{
    const char* rawValue = std::getenv( variableName.c_str( ) );
    if( rawValue == nullptr )
    {
        return false;
    }

    parsedValues.clear( );
    std::stringstream stream( rawValue );
    std::string token;
    while( std::getline( stream, token, ',' ) )
    {
        parsedValues.push_back( std::stod( token ) );
    }
    return true;
}

std::array< double, 5 > parseDebugCoefficientVector( const std::string& variableName, const Eigen::Matrix3d& fallbackKMatrix )
{
    std::vector< double > parsedValues;
    if( parseDebugVectorFromEnvironment( variableName, parsedValues ) && parsedValues.size( ) == 5 )
    {
        return { parsedValues.at( 0 ), parsedValues.at( 1 ), parsedValues.at( 2 ), parsedValues.at( 3 ), parsedValues.at( 4 ) };
    }

    return { 0.5 * fallbackKMatrix( 0, 0 ) + 0.5 * fallbackKMatrix( 1, 1 ) - fallbackKMatrix( 2, 2 ),
             -fallbackKMatrix( 2, 0 ),
             -fallbackKMatrix( 2, 1 ),
             -0.25 * fallbackKMatrix( 0, 0 ) + 0.25 * fallbackKMatrix( 1, 1 ),
             -0.5 * fallbackKMatrix( 1, 0 ) };
}

Eigen::Vector3d parseDebugVector3( const std::string& variableName, const Eigen::Vector3d& fallback )
{
    std::vector< double > parsedValues;
    if( parseDebugVectorFromEnvironment( variableName, parsedValues ) && parsedValues.size( ) == 3 )
    {
        return Eigen::Vector3d( parsedValues.at( 0 ), parsedValues.at( 1 ), parsedValues.at( 2 ) );
    }
    return fallback;
}

double parseDebugScalar( const std::string& variableName, const double fallback )
{
    const char* rawValue = std::getenv( variableName.c_str( ) );
    return ( rawValue == nullptr ) ? fallback : std::stod( rawValue );
}

void getCoefficientMatricesFromPacket( const std::array< double, 5 >& coefficientPacket,
                                       Eigen::MatrixXd& cosineCoefficients,
                                       Eigen::MatrixXd& sineCoefficients )
{
    cosineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    sineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    cosineCoefficients( 0, 0 ) = 1.0;
    cosineCoefficients( 2, 0 ) = coefficientPacket.at( 0 );
    cosineCoefficients( 2, 1 ) = coefficientPacket.at( 1 );
    sineCoefficients( 2, 1 ) = coefficientPacket.at( 2 );
    cosineCoefficients( 2, 2 ) = coefficientPacket.at( 3 );
    sineCoefficients( 2, 2 ) = coefficientPacket.at( 4 );
}

Eigen::Matrix3d getDegreeTwoKMatrixFromUnnormalizedCoefficients( const Eigen::MatrixXd& unnormalizedCosineCoefficients,
                                                                 const Eigen::MatrixXd& unnormalizedSineCoefficients )
{
    const double c20 = unnormalizedCosineCoefficients( 2, 0 );
    const double c21 = unnormalizedCosineCoefficients( 2, 1 );
    const double s21 = unnormalizedSineCoefficients( 2, 1 );
    const double c22 = unnormalizedCosineCoefficients( 2, 2 );
    const double s22 = unnormalizedSineCoefficients( 2, 2 );

    Eigen::Matrix3d kMatrix;
    kMatrix << c20 / 3.0 - 2.0 * c22, -2.0 * s22, -c21, -2.0 * s22, c20 / 3.0 + 2.0 * c22, -s21, -c21, -s21, -2.0 * c20 / 3.0;
    return kMatrix;
}

Eigen::Vector3d vexSkewSymmetricMatrix( const Eigen::Matrix3d& matrix )
{
    return Eigen::Vector3d( matrix( 2, 1 ), matrix( 0, 2 ), matrix( 1, 0 ) );
}

void emitCanonicalFigureFigureTrace( const std::string& model,
                                     const Eigen::Vector3d& canonicalRelativePositionInBody1Frame,
                                     const Eigen::Matrix3d& fallbackKMatrixBody1,
                                     const Eigen::Matrix3d& fallbackKMatrixBody2,
                                     const Eigen::Vector3d& requestedOutputTorqueOnBody1InBody1Frame )
{
    if( model != "schutz" )
    {
        return;
    }

    const std::array< double, 5 > normalizedCoefficientsBody1 =
            parseDebugCoefficientVector( "TUDAT_FFDBG_BODY1_COEFFS", fallbackKMatrixBody1 );
    const std::array< double, 5 > normalizedCoefficientsBody2 =
            parseDebugCoefficientVector( "TUDAT_FFDBG_BODY2_COEFFS", fallbackKMatrixBody2 );
    const double massBody1 = parseDebugScalar( "TUDAT_FFDBG_BODY1_MASS", 1.0 );
    const double massBody2 = parseDebugScalar( "TUDAT_FFDBG_BODY2_MASS", 1.0 );
    const double radiusBody1 = parseDebugScalar( "TUDAT_FFDBG_BODY1_RADIUS", 1.0 );
    const double radiusBody2 = parseDebugScalar( "TUDAT_FFDBG_BODY2_RADIUS", 1.0 );

    Eigen::MatrixXd normalizedCosineCoefficientsBody1;
    Eigen::MatrixXd normalizedSineCoefficientsBody1;
    Eigen::MatrixXd normalizedCosineCoefficientsBody2;
    Eigen::MatrixXd normalizedSineCoefficientsBody2;
    getCoefficientMatricesFromPacket( normalizedCoefficientsBody1, normalizedCosineCoefficientsBody1, normalizedSineCoefficientsBody1 );
    getCoefficientMatricesFromPacket( normalizedCoefficientsBody2, normalizedCosineCoefficientsBody2, normalizedSineCoefficientsBody2 );

    Eigen::MatrixXd unnormalizedCosineCoefficientsBody1;
    Eigen::MatrixXd unnormalizedSineCoefficientsBody1;
    Eigen::MatrixXd unnormalizedCosineCoefficientsBody2;
    Eigen::MatrixXd unnormalizedSineCoefficientsBody2;
    basic_mathematics::convertGeodesyNormalizedToUnnormalizedCoefficients( normalizedCosineCoefficientsBody1,
                                                                           normalizedSineCoefficientsBody1,
                                                                           unnormalizedCosineCoefficientsBody1,
                                                                           unnormalizedSineCoefficientsBody1 );
    basic_mathematics::convertGeodesyNormalizedToUnnormalizedCoefficients( normalizedCosineCoefficientsBody2,
                                                                           normalizedSineCoefficientsBody2,
                                                                           unnormalizedCosineCoefficientsBody2,
                                                                           unnormalizedSineCoefficientsBody2 );

    const Eigen::Matrix3d kMatrixBody1 =
            getDegreeTwoKMatrixFromUnnormalizedCoefficients( unnormalizedCosineCoefficientsBody1, unnormalizedSineCoefficientsBody1 );
    const Eigen::Matrix3d kMatrixBody2 =
            getDegreeTwoKMatrixFromUnnormalizedCoefficients( unnormalizedCosineCoefficientsBody2, unnormalizedSineCoefficientsBody2 );
    const double r = canonicalRelativePositionInBody1Frame.norm( );
    const Eigen::Vector3d n = canonicalRelativePositionInBody1Frame / r;
    const Eigen::Matrix3d nn = n * n.transpose( );
    const Eigen::Vector3d k1n = kMatrixBody1 * n;
    const Eigen::Vector3d k2n = kMatrixBody2 * n;
    const double a = n.dot( k1n );
    const double b = n.dot( k2n );
    const double c = n.dot( kMatrixBody1 * k2n );
    const double d = ( kMatrixBody1 * kMatrixBody2 ).trace( );
    const double f = 105.0 * a * b - 60.0 * c + 6.0 * d;
    const Eigen::Matrix3d h1 = 105.0 * b * nn - 30.0 * ( kMatrixBody2 * nn + nn * kMatrixBody2 ) + 6.0 * kMatrixBody2;
    const Eigen::Matrix3d h2 = 105.0 * a * nn - 30.0 * ( kMatrixBody1 * nn + nn * kMatrixBody1 ) + 6.0 * kMatrixBody1;
    const Eigen::Matrix3d commutatorBody1 = kMatrixBody1 * h1 - h1 * kMatrixBody1;
    const Eigen::Matrix3d commutatorBody2 = kMatrixBody2 * h2 - h2 * kMatrixBody2;
    const double commonTorqueFactor = -physical_constants::GRAVITATIONAL_CONSTANT * massBody1 * massBody2 * radiusBody1 * radiusBody1 *
            radiusBody2 * radiusBody2 / ( 2.0 * std::pow( r, 5.0 ) );
    const Eigen::Vector3d torqueBody1 = commonTorqueFactor * vexSkewSymmetricMatrix( commutatorBody1 );
    const Eigen::Vector3d torqueBody2 = commonTorqueFactor * vexSkewSymmetricMatrix( commutatorBody2 );

    std::array< double, 5 > unnormalizedCoefficientsBody1 = { unnormalizedCosineCoefficientsBody1( 2, 0 ),
                                                              unnormalizedCosineCoefficientsBody1( 2, 1 ),
                                                              unnormalizedSineCoefficientsBody1( 2, 1 ),
                                                              unnormalizedCosineCoefficientsBody1( 2, 2 ),
                                                              unnormalizedSineCoefficientsBody1( 2, 2 ) };
    std::array< double, 5 > unnormalizedCoefficientsBody2 = { unnormalizedCosineCoefficientsBody2( 2, 0 ),
                                                              unnormalizedCosineCoefficientsBody2( 2, 1 ),
                                                              unnormalizedSineCoefficientsBody2( 2, 1 ),
                                                              unnormalizedCosineCoefficientsBody2( 2, 2 ),
                                                              unnormalizedSineCoefficientsBody2( 2, 2 ) };

    printCanonicalFigureFigureTraceLine( "000.input.coefficients.body1",
                                         "C20_C21_S21_C22_S22_geodesyNormalized",
                                         "matrix3",
                                         formatDebugCoefficientPacket( normalizedCoefficientsBody1, 1.0 ) );
    printCanonicalFigureFigureTraceLine( "001.input.coefficients.body2",
                                         "C20_C21_S21_C22_S22_geodesyNormalized",
                                         "matrix3",
                                         formatDebugCoefficientPacket( normalizedCoefficientsBody2, 1.0 ) );
    printCanonicalFigureFigureTraceLine( "002.input.massRadius.body1",
                                         "mass_radius_unused",
                                         "vector3",
                                         formatDebugVector( Eigen::Vector3d( massBody1, radiusBody1, 0.0 ) ) );
    printCanonicalFigureFigureTraceLine( "003.input.massRadius.body2",
                                         "mass_radius_unused",
                                         "vector3",
                                         formatDebugVector( Eigen::Vector3d( massBody2, radiusBody2, 0.0 ) ) );
    printCanonicalFigureFigureTraceLine( "004.input.relativePosition.F1",
                                         "body2_minus_body1_in_F1",
                                         "vector3",
                                         formatDebugVector( canonicalRelativePositionInBody1Frame ) );
    printCanonicalFigureFigureTraceLine(
            "005.input.relativeVectorConvention", "body2_minus_body1_in_F1", "scalar", formatDebugScalar( 1.0 ) );
    printCanonicalFigureFigureTraceLine( "006.input.torqueBodyAndFrame", "torque_on_body1_in_F1", "scalar", formatDebugScalar( 1.0 ) );
    printCanonicalFigureFigureTraceLine( "010.coefficients.unnormalized.body1",
                                         "C20_C21_S21_C22_S22_unnormalized",
                                         "matrix3",
                                         formatDebugCoefficientPacket( unnormalizedCoefficientsBody1, 0.0 ) );
    printCanonicalFigureFigureTraceLine( "011.coefficients.unnormalized.body2",
                                         "C20_C21_S21_C22_S22_unnormalized",
                                         "matrix3",
                                         formatDebugCoefficientPacket( unnormalizedCoefficientsBody2, 0.0 ) );
    printCanonicalFigureFigureTraceLine( "020.K.body1", "K1", "matrix3", formatDebugMatrix( kMatrixBody1 ) );
    printCanonicalFigureFigureTraceLine( "021.K.body2", "K2", "matrix3", formatDebugMatrix( kMatrixBody2 ) );
    printCanonicalFigureFigureTraceLine( "030.n", "unit_relative_position", "vector3", formatDebugVector( n ) );
    printCanonicalFigureFigureTraceLine( "031.r", "relative_distance", "scalar", formatDebugScalar( r ) );
    printCanonicalFigureFigureTraceLine( "040.K1n", "K1_times_n", "vector3", formatDebugVector( k1n ) );
    printCanonicalFigureFigureTraceLine( "041.K2n", "K2_times_n", "vector3", formatDebugVector( k2n ) );
    printCanonicalFigureFigureTraceLine( "042.a_nK1n", "a", "scalar", formatDebugScalar( a ) );
    printCanonicalFigureFigureTraceLine( "043.b_nK2n", "b", "scalar", formatDebugScalar( b ) );
    printCanonicalFigureFigureTraceLine( "044.c_nK1K2n", "c", "scalar", formatDebugScalar( c ) );
    printCanonicalFigureFigureTraceLine( "045.d_traceK1K2", "d", "scalar", formatDebugScalar( d ) );
    printCanonicalFigureFigureTraceLine( "046.F_105ab_minus60c_plus6d", "F", "scalar", formatDebugScalar( f ) );
    printCanonicalFigureFigureTraceLine( "050.H_for_body1", "H1", "matrix3", formatDebugMatrix( h1 ) );
    printCanonicalFigureFigureTraceLine( "051.H_for_body2", "H2", "matrix3", formatDebugMatrix( h2 ) );
    printCanonicalFigureFigureTraceLine(
            "060.commutator.body1.KHminusHK", "K1H1_minus_H1K1", "matrix3", formatDebugMatrix( commutatorBody1 ) );
    printCanonicalFigureFigureTraceLine(
            "061.commutator.body2.KHminusHK", "K2H2_minus_H2K2", "matrix3", formatDebugMatrix( commutatorBody2 ) );
    printCanonicalFigureFigureTraceLine( "070.torque.body1.F1", "torque_on_body1_in_F1", "vector3", formatDebugVector( torqueBody1 ) );
    printCanonicalFigureFigureTraceLine( "071.torque.body2.F1", "torque_on_body2_in_F1", "vector3", formatDebugVector( torqueBody2 ) );
    printCanonicalFigureFigureTraceLine( "072.torque.requestedOutput",
                                         "torque_on_body1_in_F1",
                                         "vector3",
                                         formatDebugVector( requestedOutputTorqueOnBody1InBody1Frame ) );
}

}  // namespace

Eigen::Vector3d calculateFourthDegreeFullTwoBodyGravitationalTorqueFromTensorComponents(
        const Eigen::Vector3d& relativePositionOfBodyExertingTorqueInBodyFixedFrame,
        const double massOfBodyExertingTorque,
        const Eigen::Matrix3d& inertiaTensorOfBodyUndergoingTorque,
        const Eigen::Matrix3d& inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque )
{
    // Schutz et al. (1981), Celestial Mechanics 24, 173-181, Eq. (11):
    // explicit fourth-degree two-body torque written in Cartesian tensor components.
    const double x = relativePositionOfBodyExertingTorqueInBodyFixedFrame( 0 );
    const double y = relativePositionOfBodyExertingTorqueInBodyFixedFrame( 1 );
    const double z = relativePositionOfBodyExertingTorqueInBodyFixedFrame( 2 );
    const double x2 = x * x;
    const double y2 = y * y;
    const double z2 = z * z;
    const double xy = x * y;
    const double xz = x * z;
    const double yz = y * z;

    const double r2 = x2 + y2 + z2;
    const double r = std::sqrt( r2 );
    const double inverseR2 = 1.0 / r2;
    const double r5 = r2 * r2 * r;

    const double A = inertiaTensorOfBodyUndergoingTorque( 0, 0 );
    const double B = inertiaTensorOfBodyUndergoingTorque( 1, 1 );
    const double C = inertiaTensorOfBodyUndergoingTorque( 2, 2 );
    const double Ixy = -inertiaTensorOfBodyUndergoingTorque( 0, 1 );
    const double Ixz = -inertiaTensorOfBodyUndergoingTorque( 0, 2 );
    const double Iyz = -inertiaTensorOfBodyUndergoingTorque( 1, 2 );

    const double Aprime = inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( 0, 0 );
    const double Bprime = inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( 1, 1 );
    const double Cprime = inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( 2, 2 );
    const double IxyPrime = -inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( 0, 1 );
    const double IxzPrime = -inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( 0, 2 );
    const double IyzPrime = -inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( 1, 2 );

    // Eq. (11) auxiliary invariants of body 2 in the body-1 frame.
    const double Qprime = Aprime + Bprime + Cprime;
    const double IellPrime =
            ( Aprime * x2 + Bprime * y2 + Cprime * z2 - 2.0 * IxyPrime * xy - 2.0 * IxzPrime * xz - 2.0 * IyzPrime * yz ) * inverseR2;
    const double Wprime = massOfBodyExertingTorque + 7.5 * Qprime * inverseR2 - 17.5 * IellPrime * inverseR2;

    // Eq. (11): f- and g-functions.
    const double fyz = yz * ( Wprime - 5.0 * Aprime * inverseR2 ) - 5.0 * IxzPrime * xy * inverseR2 - 5.0 * IxyPrime * xz * inverseR2 +
            IyzPrime * ( 1.0 - 5.0 * ( y2 + z2 ) * inverseR2 );
    const double fxz = xz * ( Wprime - 5.0 * Bprime * inverseR2 ) + IxzPrime * ( 1.0 - 5.0 * ( x2 + z2 ) * inverseR2 ) -
            5.0 * IyzPrime * xy * inverseR2 - 5.0 * IxyPrime * yz * inverseR2;
    const double fxy = xy * ( Wprime - 5.0 * Cprime * inverseR2 ) - 5.0 * IyzPrime * xz * inverseR2 +
            IxyPrime * ( 1.0 - 5.0 * ( x2 + y2 ) * inverseR2 ) - 5.0 * IxzPrime * yz * inverseR2;

    const double gyzWprimeTerm = ( z2 - y2 ) * Wprime;
    const double gyzDiagonalTerm = Bprime - Cprime;
    const double gyzIxzPrimeTerm = -10.0 * IxzPrime * xz * inverseR2;
    const double gyzIxyPrimeTerm = -10.0 * IxyPrime * xy * inverseR2;
    const double gyzIyzPrimeTerm = -20.0 * IyzPrime * yz * inverseR2;
    const double gyzQuadraticZTerm = -5.0 * z2 * ( Aprime + Bprime - Cprime ) * inverseR2;
    const double gyzQuadraticYTerm = -5.0 * y2 * ( Aprime - Bprime + Cprime ) * inverseR2;
    const double gyz =
            gyzWprimeTerm + gyzDiagonalTerm + gyzIxzPrimeTerm + gyzIxyPrimeTerm + gyzIyzPrimeTerm + gyzQuadraticZTerm + gyzQuadraticYTerm;

    const double gxzWprimeTerm = ( x2 - z2 ) * Wprime;
    const double gxzDiagonalTerm = Cprime - Aprime;
    const double gxzIxzPrimeTerm = -20.0 * IxzPrime * xz * inverseR2;
    const double gxzIxyPrimeTerm = -10.0 * IxyPrime * xy * inverseR2;
    const double gxzIyzPrimeTerm = -10.0 * IyzPrime * yz * inverseR2;
    const double gxzQuadraticXTerm = -5.0 * x2 * ( -Aprime + Bprime + Cprime ) * inverseR2;
    const double gxzQuadraticZTerm = -5.0 * z2 * ( Aprime + Bprime - Cprime ) * inverseR2;
    const double gxz =
            gxzWprimeTerm + gxzDiagonalTerm + gxzIxzPrimeTerm + gxzIxyPrimeTerm + gxzIyzPrimeTerm + gxzQuadraticXTerm + gxzQuadraticZTerm;

    const double gxyWprimeTerm = ( y2 - x2 ) * Wprime;
    const double gxyDiagonalTerm = Aprime - Bprime;
    const double gxyIxzPrimeTerm = -10.0 * IxzPrime * xz * inverseR2;
    const double gxyIxyPrimeTerm = -20.0 * IxyPrime * xy * inverseR2;
    const double gxyIyzPrimeTerm = -10.0 * IyzPrime * yz * inverseR2;
    const double gxyQuadraticYTerm = -5.0 * y2 * ( Aprime - Bprime + Cprime ) * inverseR2;
    const double gxyQuadraticXTerm = -5.0 * x2 * ( -Aprime + Bprime + Cprime ) * inverseR2;
    const double gxy =
            gxyWprimeTerm + gxyDiagonalTerm + gxyIxzPrimeTerm + gxyIxyPrimeTerm + gxyIyzPrimeTerm + gxyQuadraticYTerm + gxyQuadraticXTerm;

    // Independent expansion of Schutz's printed integral definitions:
    // g_yz = int Gamma_5(l_z^2-l_y^2)dm', g_xz = int Gamma_5(l_x^2-l_z^2)dm',
    // g_xy = int Gamma_5(l_y^2-l_x^2)dm'. This is diagnostic only; the returned torque
    // above uses the compact printed Eq. (11d-f) expressions.
    const double gyzIntegralDefinitionIxzPrimeTerm = -10.0 * IxzPrime * xz * inverseR2;
    const double gyzIntegralDefinitionIxyPrimeTerm = 10.0 * IxyPrime * xy * inverseR2;
    const double gyzIntegralDefinitionIyzPrimeTerm = 0.0;
    const double gyzIntegralDefinitionQuadraticZTerm = -5.0 * z2 * ( Aprime + Bprime - Cprime ) * inverseR2;
    const double gyzIntegralDefinitionQuadraticYTerm = 5.0 * y2 * ( Aprime - Bprime + Cprime ) * inverseR2;
    const double gyzFromIntegralDefinition = gyzWprimeTerm + gyzDiagonalTerm + gyzIntegralDefinitionIxzPrimeTerm +
            gyzIntegralDefinitionIxyPrimeTerm + gyzIntegralDefinitionIyzPrimeTerm + gyzIntegralDefinitionQuadraticZTerm +
            gyzIntegralDefinitionQuadraticYTerm;

    const double gxzIntegralDefinitionIxzPrimeTerm = 0.0;
    const double gxzIntegralDefinitionIxyPrimeTerm = -10.0 * IxyPrime * xy * inverseR2;
    const double gxzIntegralDefinitionIyzPrimeTerm = 10.0 * IyzPrime * yz * inverseR2;
    const double gxzIntegralDefinitionQuadraticXTerm = -5.0 * x2 * ( -Aprime + Bprime + Cprime ) * inverseR2;
    const double gxzIntegralDefinitionQuadraticZTerm = 5.0 * z2 * ( Aprime + Bprime - Cprime ) * inverseR2;
    const double gxzFromIntegralDefinition = gxzWprimeTerm + gxzDiagonalTerm + gxzIntegralDefinitionIxzPrimeTerm +
            gxzIntegralDefinitionIxyPrimeTerm + gxzIntegralDefinitionIyzPrimeTerm + gxzIntegralDefinitionQuadraticXTerm +
            gxzIntegralDefinitionQuadraticZTerm;

    const double gxyIntegralDefinitionIxzPrimeTerm = 10.0 * IxzPrime * xz * inverseR2;
    const double gxyIntegralDefinitionIxyPrimeTerm = 0.0;
    const double gxyIntegralDefinitionIyzPrimeTerm = -10.0 * IyzPrime * yz * inverseR2;
    const double gxyIntegralDefinitionQuadraticYTerm = -5.0 * y2 * ( Aprime - Bprime + Cprime ) * inverseR2;
    const double gxyIntegralDefinitionQuadraticXTerm = 5.0 * x2 * ( -Aprime + Bprime + Cprime ) * inverseR2;
    const double gxyFromIntegralDefinition = gxyWprimeTerm + gxyDiagonalTerm + gxyIntegralDefinitionIxzPrimeTerm +
            gxyIntegralDefinitionIxyPrimeTerm + gxyIntegralDefinitionIyzPrimeTerm + gxyIntegralDefinitionQuadraticYTerm +
            gxyIntegralDefinitionQuadraticXTerm;

    // Eq. (11): torque components in body-1-fixed frame.
    // The returned torque uses the g-functions expanded from Schutz's integral definitions.
    // These are also consistent with differentiating the mutual potential in Schutz Eq. (8).
    // The compact printed Eq. (11d-f) expansion is retained in debug output above; it differs
    // from both the integral definitions and Eq. (8) for several product-of-inertia terms.
    const double gyzForTorque = gyzFromIntegralDefinition;
    const double gxzForTorque = gxzFromIntegralDefinition;
    const double gxyForTorque = gxyFromIntegralDefinition;
    const double prefactor = ( 3.0 * physical_constants::GRAVITATIONAL_CONSTANT ) / r5;
    Eigen::Vector3d torque;
    const double torqueXFromDiagonalTerm = ( C - B ) * fyz;
    const double torqueXFromIxzTerm = -Ixz * fxy;
    const double torqueXFromIxyTerm = Ixy * fxz;
    const double torqueXFromIyzTerm = Iyz * gyzForTorque;

    const double torqueYFromDiagonalTerm = ( A - C ) * fxz;
    const double torqueYFromIxzTerm = Ixz * gxzForTorque;
    const double torqueYFromIxyTerm = -Ixy * fyz;
    const double torqueYFromIyzTerm = Iyz * fxy;

    const double torqueZFromDiagonalTerm = ( B - A ) * fxy;
    const double torqueZFromIxzTerm = Ixz * fyz;
    const double torqueZFromIxyTerm = Ixy * gxyForTorque;
    const double torqueZFromIyzTerm = -Iyz * fxz;

    torque( 0 ) = prefactor * ( torqueXFromDiagonalTerm + torqueXFromIxzTerm + torqueXFromIxyTerm + torqueXFromIyzTerm );
    torque( 1 ) = prefactor * ( torqueYFromDiagonalTerm + torqueYFromIxzTerm + torqueYFromIxyTerm + torqueYFromIyzTerm );
    torque( 2 ) = prefactor * ( torqueZFromDiagonalTerm + torqueZFromIxzTerm + torqueZFromIxyTerm + torqueZFromIyzTerm );

    if( isFigureFigureDebugEnabled( ) )
    {
        const Eigen::Vector3d canonicalRelativePositionInBody1Frame =
                parseDebugVector3( "TUDAT_FFDBG_RELATIVE_POSITION_F1", relativePositionOfBodyExertingTorqueInBodyFixedFrame );
        emitCanonicalFigureFigureTrace( "schutz",
                                        canonicalRelativePositionInBody1Frame,
                                        inertiaTensorOfBodyUndergoingTorque,
                                        inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque,
                                        torque );

        printNativeFigureFigureTraceLine( "100.native.relative_state",
                                          "relative_position_used_by_schutz",
                                          "vector3",
                                          formatDebugVector( relativePositionOfBodyExertingTorqueInBodyFixedFrame ) );
        printNativeFigureFigureTraceLine(
                "101.native.inertia.body1", "inertia_tensor_body1", "matrix3", formatDebugMatrix( inertiaTensorOfBodyUndergoingTorque ) );
        printNativeFigureFigureTraceLine( "102.native.inertia.body2",
                                          "inertia_tensor_body2_rotated_to_F1",
                                          "matrix3",
                                          formatDebugMatrix( inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque ) );
        printNativeFigureFigureTraceLine(
                "103.native.inertia_components.body1",
                "A_B_C_Ixy_Ixz_Iyz",
                "matrix3",
                formatDebugMatrix( ( Eigen::Matrix3d( ) << A, B, C, Ixy, Ixz, Iyz, 0.0, 0.0, 0.0 ).finished( ) ) );
        printNativeFigureFigureTraceLine(
                "104.native.inertia_components.body2",
                "Aprime_Bprime_Cprime_Ixyprime_Ixzprime_Iyzprime",
                "matrix3",
                formatDebugMatrix(
                        ( Eigen::Matrix3d( ) << Aprime, Bprime, Cprime, IxyPrime, IxzPrime, IyzPrime, 0.0, 0.0, 0.0 ).finished( ) ) );
        printNativeFigureFigureTraceLine( "105.native.Wprime", "Wprime", "scalar", formatDebugScalar( Wprime ) );
        printNativeFigureFigureTraceLine(
                "106.native.f_terms", "f_yz_f_xz_f_xy", "vector3", formatDebugVector( Eigen::Vector3d( fyz, fxz, fxy ) ) );
        printNativeFigureFigureTraceLine( "107.native.g_terms",
                                          "g_yz_g_xz_g_xy_printed_expansion",
                                          "vector3",
                                          formatDebugVector( Eigen::Vector3d( gyz, gxz, gxy ) ) );
        printNativeFigureFigureTraceLine(
                "107a.native.g_terms.integral_definition",
                "g_yz_g_xz_g_xy_from_printed_integral_definitions",
                "vector3",
                formatDebugVector( Eigen::Vector3d( gyzFromIntegralDefinition, gxzFromIntegralDefinition, gxyFromIntegralDefinition ) ) );
        printNativeFigureFigureTraceLine(
                "107b.native.g_terms.printed_minus_integral_definition",
                "delta_g_yz_g_xz_g_xy",
                "vector3",
                formatDebugVector( Eigen::Vector3d(
                        gyz - gyzFromIntegralDefinition, gxz - gxzFromIntegralDefinition, gxy - gxyFromIntegralDefinition ) ) );
        printNativeFigureFigureTraceLine( "108.native.g_terms.expanded.yz",
                                          "W_diagonal_Ixz_Ixy_Iyz_quadZ_quadY_sum",
                                          "matrix3",
                                          formatDebugMatrix( ( Eigen::Matrix3d( ) << gyzWprimeTerm,
                                                               gyzDiagonalTerm,
                                                               gyzIxzPrimeTerm,
                                                               gyzIxyPrimeTerm,
                                                               gyzIyzPrimeTerm,
                                                               gyzQuadraticZTerm,
                                                               gyzQuadraticYTerm,
                                                               gyz,
                                                               0.0 )
                                                                     .finished( ) ) );
        printNativeFigureFigureTraceLine( "108a.native.g_terms.integral_definition.yz",
                                          "W_diagonal_Ixz_Ixy_Iyz_quadZ_quadY_sum",
                                          "matrix3",
                                          formatDebugMatrix( ( Eigen::Matrix3d( ) << gyzWprimeTerm,
                                                               gyzDiagonalTerm,
                                                               gyzIntegralDefinitionIxzPrimeTerm,
                                                               gyzIntegralDefinitionIxyPrimeTerm,
                                                               gyzIntegralDefinitionIyzPrimeTerm,
                                                               gyzIntegralDefinitionQuadraticZTerm,
                                                               gyzIntegralDefinitionQuadraticYTerm,
                                                               gyzFromIntegralDefinition,
                                                               0.0 )
                                                                     .finished( ) ) );
        printNativeFigureFigureTraceLine( "109.native.g_terms.expanded.xz",
                                          "W_diagonal_Ixz_Ixy_Iyz_quadX_quadZ_sum",
                                          "matrix3",
                                          formatDebugMatrix( ( Eigen::Matrix3d( ) << gxzWprimeTerm,
                                                               gxzDiagonalTerm,
                                                               gxzIxzPrimeTerm,
                                                               gxzIxyPrimeTerm,
                                                               gxzIyzPrimeTerm,
                                                               gxzQuadraticXTerm,
                                                               gxzQuadraticZTerm,
                                                               gxz,
                                                               0.0 )
                                                                     .finished( ) ) );
        printNativeFigureFigureTraceLine( "109a.native.g_terms.integral_definition.xz",
                                          "W_diagonal_Ixz_Ixy_Iyz_quadX_quadZ_sum",
                                          "matrix3",
                                          formatDebugMatrix( ( Eigen::Matrix3d( ) << gxzWprimeTerm,
                                                               gxzDiagonalTerm,
                                                               gxzIntegralDefinitionIxzPrimeTerm,
                                                               gxzIntegralDefinitionIxyPrimeTerm,
                                                               gxzIntegralDefinitionIyzPrimeTerm,
                                                               gxzIntegralDefinitionQuadraticXTerm,
                                                               gxzIntegralDefinitionQuadraticZTerm,
                                                               gxzFromIntegralDefinition,
                                                               0.0 )
                                                                     .finished( ) ) );
        printNativeFigureFigureTraceLine( "110.native.g_terms.expanded.xy",
                                          "W_diagonal_Ixz_Ixy_Iyz_quadY_quadX_sum",
                                          "matrix3",
                                          formatDebugMatrix( ( Eigen::Matrix3d( ) << gxyWprimeTerm,
                                                               gxyDiagonalTerm,
                                                               gxyIxzPrimeTerm,
                                                               gxyIxyPrimeTerm,
                                                               gxyIyzPrimeTerm,
                                                               gxyQuadraticYTerm,
                                                               gxyQuadraticXTerm,
                                                               gxy,
                                                               0.0 )
                                                                     .finished( ) ) );
        printNativeFigureFigureTraceLine( "110a.native.g_terms.integral_definition.xy",
                                          "W_diagonal_Ixz_Ixy_Iyz_quadY_quadX_sum",
                                          "matrix3",
                                          formatDebugMatrix( ( Eigen::Matrix3d( ) << gxyWprimeTerm,
                                                               gxyDiagonalTerm,
                                                               gxyIntegralDefinitionIxzPrimeTerm,
                                                               gxyIntegralDefinitionIxyPrimeTerm,
                                                               gxyIntegralDefinitionIyzPrimeTerm,
                                                               gxyIntegralDefinitionQuadraticYTerm,
                                                               gxyIntegralDefinitionQuadraticXTerm,
                                                               gxyFromIntegralDefinition,
                                                               0.0 )
                                                                     .finished( ) ) );
        printNativeFigureFigureTraceLine( "111.native.raw_torque_terms",
                                          "Tx_terms_Ty_terms_Tz_terms",
                                          "matrix3",
                                          formatDebugMatrix( ( Eigen::Matrix3d( ) << torqueXFromDiagonalTerm,
                                                               torqueXFromIxzTerm,
                                                               torqueXFromIxyTerm + torqueXFromIyzTerm,
                                                               torqueYFromDiagonalTerm,
                                                               torqueYFromIxzTerm,
                                                               torqueYFromIxyTerm + torqueYFromIyzTerm,
                                                               torqueZFromDiagonalTerm,
                                                               torqueZFromIxzTerm,
                                                               torqueZFromIxyTerm + torqueZFromIyzTerm )
                                                                     .finished( ) ) );
        printNativeFigureFigureTraceLine( "111a.native.torque.integral_definition_g_terms",
                                          "torque_on_body1_in_F1_if_integral_definition_g_terms_are_used",
                                          "vector3",
                                          formatDebugVector( torque ) );
        printNativeFigureFigureTraceLine(
                "112.native.final_converted_torque", "torque_on_body1_in_F1", "vector3", formatDebugVector( torque ) );
    }

    return torque;
}

Eigen::Vector3d calculateFourthDegreeFullTwoBodyGravitationalTorque(
        const Eigen::Vector3d& relativePositionOfBodyExertingTorqueInBodyFixedFrame,
        const double massOfBodyExertingTorque,
        const Eigen::Matrix3d& inertiaTensorOfBodyUndergoingTorque,
        const Eigen::Matrix3d& inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque )
{
    const double relativeDistance = relativePositionOfBodyExertingTorqueInBodyFixedFrame.norm( );
    if( relativeDistance <= 0.0 )
    {
        throw std::runtime_error( "Error when computing fourth-degree full two-body gravitational torque: relative distance is zero." );
    }

    // Schutz et al. (1981), Eq. (11): direct component-wise fourth-degree two-body torque.
    return calculateFourthDegreeFullTwoBodyGravitationalTorqueFromTensorComponents(
            relativePositionOfBodyExertingTorqueInBodyFixedFrame,
            massOfBodyExertingTorque,
            inertiaTensorOfBodyUndergoingTorque,
            inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque );
}

FourthDegreeFullTwoBodyGravitationalTorqueModel::FourthDegreeFullTwoBodyGravitationalTorqueModel(
        const std::function< Eigen::Vector3d( ) > positionOfBodyUndergoingTorqueFunction,
        const std::function< Eigen::Vector3d( ) > positionOfBodyExertingTorqueFunction,
        const std::function< double( ) > massOfBodyExertingTorqueFunction,
        const std::function< Eigen::Matrix3d( ) > inertiaTensorOfBodyUndergoingTorqueFunction,
        const std::function< Eigen::Matrix3d( ) > inertiaTensorOfBodyExertingTorqueFunction,
        const std::function< Eigen::Quaterniond( ) > rotationToBodyFixedFrameOfBodyUndergoingTorqueFunction,
        const std::function< Eigen::Quaterniond( ) > rotationToBodyFixedFrameOfBodyExertingTorqueFunction ):
    positionOfBodyUndergoingTorqueFunction_( positionOfBodyUndergoingTorqueFunction ),
    positionOfBodyExertingTorqueFunction_( positionOfBodyExertingTorqueFunction ),
    massOfBodyExertingTorqueFunction_( massOfBodyExertingTorqueFunction ),
    inertiaTensorOfBodyUndergoingTorqueFunction_( inertiaTensorOfBodyUndergoingTorqueFunction ),
    inertiaTensorOfBodyExertingTorqueFunction_( inertiaTensorOfBodyExertingTorqueFunction ),
    rotationToBodyFixedFrameOfBodyUndergoingTorqueFunction_( rotationToBodyFixedFrameOfBodyUndergoingTorqueFunction ),
    rotationToBodyFixedFrameOfBodyExertingTorqueFunction_( rotationToBodyFixedFrameOfBodyExertingTorqueFunction )
{}

void FourthDegreeFullTwoBodyGravitationalTorqueModel::updateMembers( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        currentRotationToBodyFixedFrameOfBodyUndergoingTorque_ = rotationToBodyFixedFrameOfBodyUndergoingTorqueFunction_( );
        currentRotationToBodyFixedFrameOfBodyExertingTorque_ = rotationToBodyFixedFrameOfBodyExertingTorqueFunction_( );
        currentRelativePositionInInertialFrame_ = positionOfBodyExertingTorqueFunction_( ) - positionOfBodyUndergoingTorqueFunction_( );
        currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_ =
                currentRotationToBodyFixedFrameOfBodyUndergoingTorque_ * currentRelativePositionInInertialFrame_;

        currentMassOfBodyExertingTorque_ = massOfBodyExertingTorqueFunction_( );
        currentInertiaTensorOfBodyUndergoingTorque_ = inertiaTensorOfBodyUndergoingTorqueFunction_( );
        currentInertiaTensorOfBodyExertingTorque_ = inertiaTensorOfBodyExertingTorqueFunction_( );

        const Eigen::Matrix3d rotationFromBody2ToBody1 = currentRotationToBodyFixedFrameOfBodyUndergoingTorque_.toRotationMatrix( ) *
                currentRotationToBodyFixedFrameOfBodyExertingTorque_.toRotationMatrix( ).transpose( );
        currentInertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque_ =
                rotationFromBody2ToBody1 * currentInertiaTensorOfBodyExertingTorque_ * rotationFromBody2ToBody1.transpose( );

        // Eq. (11): evaluate torque from relative position and inertia tensors in body-1-fixed coordinates.
        currentTorque_ = calculateFourthDegreeFullTwoBodyGravitationalTorque(
                currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_,
                currentMassOfBodyExertingTorque_,
                currentInertiaTensorOfBodyUndergoingTorque_,
                currentInertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque_ );

        currentTime_ = currentTime;
    }
}

}  // namespace gravitation

}  // namespace tudat
