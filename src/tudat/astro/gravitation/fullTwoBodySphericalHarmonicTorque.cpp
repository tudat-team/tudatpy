#include "tudat/astro/gravitation/fullTwoBodySphericalHarmonicTorque.h"

#include <boost/math/special_functions/factorials.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <set>
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
    std::cout << "FFDBG|case=" << getFigureFigureDebugCase( ) << "|model=dmr|step=" << step << "|name=" << name << "|shape=" << shape
              << "|values=" << values << std::endl;
}

void printNativeFigureFigureTraceLine( const std::string& step,
                                       const std::string& name,
                                       const std::string& shape,
                                       const std::string& values )
{
    std::cout << "FFDBG_NATIVE|case=" << getFigureFigureDebugCase( ) << "|model=dmr|step=" << step << "|name=" << name << "|shape=" << shape
              << "|values=" << values << std::endl;
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

std::array< double, 5 > parseDebugCoefficientVector( const std::string& variableName,
                                                     const Eigen::MatrixXd& cosineCoefficients,
                                                     const Eigen::MatrixXd& sineCoefficients )
{
    std::vector< double > parsedValues;
    if( parseDebugVectorFromEnvironment( variableName, parsedValues ) && parsedValues.size( ) == 5 )
    {
        return { parsedValues.at( 0 ), parsedValues.at( 1 ), parsedValues.at( 2 ), parsedValues.at( 3 ), parsedValues.at( 4 ) };
    }

    return { cosineCoefficients( 2, 0 ),
             cosineCoefficients( 2, 1 ),
             sineCoefficients( 2, 1 ),
             cosineCoefficients( 2, 2 ),
             sineCoefficients( 2, 2 ) };
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

void emitCanonicalFigureFigureTrace( const Eigen::Vector3d& canonicalRelativePositionInBody1Frame,
                                     const Eigen::MatrixXd& cosineCoefficientsOfBody1,
                                     const Eigen::MatrixXd& sineCoefficientsOfBody1,
                                     const Eigen::MatrixXd& cosineCoefficientsOfBody2,
                                     const Eigen::MatrixXd& sineCoefficientsOfBody2,
                                     const double radiusBody1,
                                     const double radiusBody2,
                                     const Eigen::Vector3d& requestedOutputTorqueOnBody1InBody1Frame )
{
    const std::array< double, 5 > normalizedCoefficientsBody1 =
            parseDebugCoefficientVector( "TUDAT_FFDBG_BODY1_COEFFS", cosineCoefficientsOfBody1, sineCoefficientsOfBody1 );
    const std::array< double, 5 > normalizedCoefficientsBody2 =
            parseDebugCoefficientVector( "TUDAT_FFDBG_BODY2_COEFFS", cosineCoefficientsOfBody2, sineCoefficientsOfBody2 );
    const double massBody1 = parseDebugScalar( "TUDAT_FFDBG_BODY1_MASS", 1.0 );
    const double massBody2 = parseDebugScalar( "TUDAT_FFDBG_BODY2_MASS", 1.0 );

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

//! Apply the angular momentum operator to one Wigner-D coefficient entry.
/*!
 * Implements the Cartesian form of \hat{J}(D^l_{m,k}) used in the torque evaluation
 * of Dirkx et al. (2019), Eq. (67), with torque definition from Eq. (60).
 */
Eigen::Vector3cd computeAngularMomentumOperatorOnWignerCoefficientFromWignerMatrices( const std::vector< Eigen::MatrixXcd >& wignerMatrices,
                                                                                      const int degree,
                                                                                      const int originalOrder,
                                                                                      const int newOrder )
{
    const std::complex< double > imaginaryUnit( 0.0, 1.0 );
    const double inverseSquareRootTwo = 1.0 / std::sqrt( 2.0 );

    const auto getWignerCoefficient = [ & ]( const int requestedOrderM, const int requestedOrderK ) {
        if( std::abs( requestedOrderM ) > degree || std::abs( requestedOrderK ) > degree )
        {
            return std::complex< double >( 0.0, 0.0 );
        }
        return wignerMatrices.at( degree )( requestedOrderM + degree, requestedOrderK + degree );
    };

    const double plusScaling =
            std::sqrt( std::max( 0.0, static_cast< double >( degree * ( degree + 1 ) - originalOrder * ( originalOrder - 1 ) ) ) ) /
            std::sqrt( 2.0 );
    const double minusScaling =
            std::sqrt( std::max( 0.0, static_cast< double >( degree * ( degree + 1 ) - originalOrder * ( originalOrder + 1 ) ) ) ) /
            std::sqrt( 2.0 );

    // Eq. (67): Cartesian components of \hat{J}(D^l_{m,k}) assembled from ladder-operator contributions.
    const std::complex< double > angularMomentumPlus = imaginaryUnit * plusScaling * getWignerCoefficient( originalOrder - 1, newOrder );
    const std::complex< double > angularMomentumMinus =
            imaginaryUnit * ( -minusScaling ) * getWignerCoefficient( originalOrder + 1, newOrder );
    const std::complex< double > angularMomentumZero =
            imaginaryUnit * static_cast< double >( -originalOrder ) * getWignerCoefficient( originalOrder, newOrder );

    Eigen::Vector3cd angularMomentumInCartesianBasis;
    angularMomentumInCartesianBasis( 0 ) = ( angularMomentumMinus - angularMomentumPlus ) * inverseSquareRootTwo;
    angularMomentumInCartesianBasis( 1 ) = imaginaryUnit * ( angularMomentumMinus + angularMomentumPlus ) * inverseSquareRootTwo;
    angularMomentumInCartesianBasis( 2 ) = angularMomentumZero;
    // Returned \hat{J}-mapped term is used directly in the Eq. (60) torque evaluation.
    return angularMomentumInCartesianBasis;
}

}  // namespace

std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > getBody2TorqueCombinationsToUse(
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse )
{
    // Build unique combinations required for body-2 spin torque accumulation (Eq. (67) term set).
    std::set< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > body2TorqueCombinationSet;
    for( const auto& combination : coefficientCombinationsToUse )
    {
        const unsigned int degreeOfBody1ForCombination = std::get< 0 >( combination );
        const unsigned int orderOfBody1ForCombination = std::get< 1 >( combination );
        const unsigned int degreeOfBody2ForCombination = std::get< 2 >( combination );
        for( unsigned int orderOfBody2ForTorque = 0; orderOfBody2ForTorque <= degreeOfBody2ForCombination; orderOfBody2ForTorque++ )
        {
            body2TorqueCombinationSet.insert( std::make_tuple(
                    degreeOfBody1ForCombination, orderOfBody1ForCombination, degreeOfBody2ForCombination, orderOfBody2ForTorque ) );
        }
    }

    return std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >( body2TorqueCombinationSet.begin( ),
                                                                                                body2TorqueCombinationSet.end( ) );
}

void computeTransformedAngularMomentumCoefficientsFromWignerMatrices(
        const Eigen::MatrixXd& cosineCoefficientsBody2,
        const Eigen::MatrixXd& sineCoefficientsBody2,
        const std::vector< Eigen::MatrixXcd >& wignerMatrices,
        const bool areCoefficientsNormalized,
        std::array< Eigen::MatrixXd, 3 >& transformedCosineCoefficientsBody2AngularMomentum,
        std::array< Eigen::MatrixXd, 3 >& transformedSineCoefficientsBody2AngularMomentum )
{
    // Explicit-matrix path (used by quaternion-derivative partials).
    // Eq. (67): computes \hat{J}-mapped coefficient fields; these feed the Eq. (60) torque relation.
    for( int i = 0; i < 3; i++ )
    {
        transformedCosineCoefficientsBody2AngularMomentum.at( i ).setZero( cosineCoefficientsBody2.rows( ),
                                                                           cosineCoefficientsBody2.cols( ) );
        transformedSineCoefficientsBody2AngularMomentum.at( i ).setZero( sineCoefficientsBody2.rows( ), sineCoefficientsBody2.cols( ) );
    }

    for( int degree = 0; degree < cosineCoefficientsBody2.rows( ); degree++ )
    {
        for( int orderM = 0; ( orderM <= degree && orderM < cosineCoefficientsBody2.cols( ) ); orderM++ )
        {
            double orderMMultiplier;
            if( !areCoefficientsNormalized )
            {
                orderMMultiplier = std::sqrt( boost::math::factorial< double >( degree - orderM ) /
                                              boost::math::factorial< double >( degree + orderM ) );
            }
            else
            {
                orderMMultiplier = ( orderM == 0 ? 1.0 : 1.0 / std::sqrt( 2.0 ) );
            }

            const Eigen::Vector3cd orderZeroAngularMomentumD =
                    computeAngularMomentumOperatorOnWignerCoefficientFromWignerMatrices( wignerMatrices, degree, orderM, 0 );
            for( int i = 0; i < 3; i++ )
            {
                transformedCosineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) +=
                        orderMMultiplier * orderZeroAngularMomentumD( i ).real( ) * cosineCoefficientsBody2( degree, 0 );
                transformedSineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) +=
                        orderMMultiplier * orderZeroAngularMomentumD( i ).imag( ) * cosineCoefficientsBody2( degree, 0 );
            }

            const int maximumOrderAtDegree = std::min( degree, static_cast< int >( cosineCoefficientsBody2.cols( ) - 1 ) );
            for( int orderK = 1; orderK <= maximumOrderAtDegree; orderK++ )
            {
                double currentMultiplier;
                if( !areCoefficientsNormalized )
                {
                    currentMultiplier = std::sqrt( boost::math::factorial< double >( degree + orderK ) /
                                                   boost::math::factorial< double >( degree - orderK ) ) *
                            orderMMultiplier;
                }
                else
                {
                    currentMultiplier = std::sqrt( 2.0 ) * orderMMultiplier;
                }

                const double signMultiplier = ( ( ( orderK % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );
                const Eigen::Vector3cd positiveOrderAngularMomentumD =
                        computeAngularMomentumOperatorOnWignerCoefficientFromWignerMatrices( wignerMatrices, degree, orderM, orderK );
                const Eigen::Vector3cd negativeOrderAngularMomentumD =
                        computeAngularMomentumOperatorOnWignerCoefficientFromWignerMatrices( wignerMatrices, degree, orderM, -orderK );

                for( int i = 0; i < 3; i++ )
                {
                    transformedCosineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) += 0.5 * currentMultiplier *
                            ( ( signMultiplier * positiveOrderAngularMomentumD( i ).real( ) + negativeOrderAngularMomentumD( i ).real( ) ) *
                                      cosineCoefficientsBody2( degree, orderK ) +
                              ( signMultiplier * positiveOrderAngularMomentumD( i ).imag( ) - negativeOrderAngularMomentumD( i ).imag( ) ) *
                                      sineCoefficientsBody2( degree, orderK ) );

                    transformedSineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) += 0.5 * currentMultiplier *
                            ( ( signMultiplier * positiveOrderAngularMomentumD( i ).imag( ) + negativeOrderAngularMomentumD( i ).imag( ) ) *
                                      cosineCoefficientsBody2( degree, orderK ) +
                              ( -signMultiplier * positiveOrderAngularMomentumD( i ).real( ) +
                                negativeOrderAngularMomentumD( i ).real( ) ) *
                                      sineCoefficientsBody2( degree, orderK ) );
                }
            }

            double cosineFinalScaling = ( ( ( orderM % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );
            double sineFinalScaling = ( ( ( ( orderM + 1 ) % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );
            if( orderM > 0 )
            {
                cosineFinalScaling *= 2.0;
                sineFinalScaling *= 2.0;
            }

            for( int i = 0; i < 3; i++ )
            {
                transformedCosineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) *= cosineFinalScaling;
                transformedSineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) *= sineFinalScaling;
            }
        }
    }
}

void computeTransformedAngularMomentumCoefficientsFromWignerCache(
        const Eigen::MatrixXd& cosineCoefficientsBody2,
        const Eigen::MatrixXd& sineCoefficientsBody2,
        const std::shared_ptr< basic_mathematics::WignerDMatricesCache >& wignerCache,
        const bool areCoefficientsNormalized,
        std::array< Eigen::MatrixXd, 3 >& transformedCosineCoefficientsBody2AngularMomentum,
        std::array< Eigen::MatrixXd, 3 >& transformedSineCoefficientsBody2AngularMomentum )
{
    // Cached path (used during nominal model evaluation).
    // Eq. (67): computes \hat{J}-mapped coefficient fields; these feed the Eq. (60) torque relation.
    for( int i = 0; i < 3; i++ )
    {
        transformedCosineCoefficientsBody2AngularMomentum.at( i ).setZero( cosineCoefficientsBody2.rows( ),
                                                                           cosineCoefficientsBody2.cols( ) );
        transformedSineCoefficientsBody2AngularMomentum.at( i ).setZero( sineCoefficientsBody2.rows( ), sineCoefficientsBody2.cols( ) );
    }

    for( int degree = 0; degree < cosineCoefficientsBody2.rows( ); degree++ )
    {
        for( int orderM = 0; ( orderM <= degree && orderM < cosineCoefficientsBody2.cols( ) ); orderM++ )
        {
            double orderMMultiplier;
            if( !areCoefficientsNormalized )
            {
                orderMMultiplier = std::sqrt( boost::math::factorial< double >( degree - orderM ) /
                                              boost::math::factorial< double >( degree + orderM ) );
            }
            else
            {
                orderMMultiplier = ( orderM == 0 ? 1.0 : 1.0 / std::sqrt( 2.0 ) );
            }

            const Eigen::Vector3cd orderZeroAngularMomentumD =
                    wignerCache->getAngularMomentumOperatorOnWignerCoefficient( degree, orderM, 0 );
            for( int i = 0; i < 3; i++ )
            {
                transformedCosineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) +=
                        orderMMultiplier * orderZeroAngularMomentumD( i ).real( ) * cosineCoefficientsBody2( degree, 0 );
                transformedSineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) +=
                        orderMMultiplier * orderZeroAngularMomentumD( i ).imag( ) * cosineCoefficientsBody2( degree, 0 );
            }

            const int maximumOrderAtDegree = std::min( degree, static_cast< int >( cosineCoefficientsBody2.cols( ) - 1 ) );
            for( int orderK = 1; orderK <= maximumOrderAtDegree; orderK++ )
            {
                double currentMultiplier;
                if( !areCoefficientsNormalized )
                {
                    currentMultiplier = std::sqrt( boost::math::factorial< double >( degree + orderK ) /
                                                   boost::math::factorial< double >( degree - orderK ) ) *
                            orderMMultiplier;
                }
                else
                {
                    currentMultiplier = std::sqrt( 2.0 ) * orderMMultiplier;
                }

                const double signMultiplier = ( ( ( orderK % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );
                const Eigen::Vector3cd positiveOrderAngularMomentumD =
                        wignerCache->getAngularMomentumOperatorOnWignerCoefficient( degree, orderM, orderK );
                const Eigen::Vector3cd negativeOrderAngularMomentumD =
                        wignerCache->getAngularMomentumOperatorOnWignerCoefficient( degree, orderM, -orderK );

                for( int i = 0; i < 3; i++ )
                {
                    transformedCosineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) += 0.5 * currentMultiplier *
                            ( ( signMultiplier * positiveOrderAngularMomentumD( i ).real( ) + negativeOrderAngularMomentumD( i ).real( ) ) *
                                      cosineCoefficientsBody2( degree, orderK ) +
                              ( signMultiplier * positiveOrderAngularMomentumD( i ).imag( ) - negativeOrderAngularMomentumD( i ).imag( ) ) *
                                      sineCoefficientsBody2( degree, orderK ) );

                    transformedSineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) += 0.5 * currentMultiplier *
                            ( ( signMultiplier * positiveOrderAngularMomentumD( i ).imag( ) + negativeOrderAngularMomentumD( i ).imag( ) ) *
                                      cosineCoefficientsBody2( degree, orderK ) +
                              ( -signMultiplier * positiveOrderAngularMomentumD( i ).real( ) +
                                negativeOrderAngularMomentumD( i ).real( ) ) *
                                      sineCoefficientsBody2( degree, orderK ) );
                }
            }

            double cosineFinalScaling = ( ( ( orderM % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );
            double sineFinalScaling = ( ( ( ( orderM + 1 ) % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );
            if( orderM > 0 )
            {
                cosineFinalScaling *= 2.0;
                sineFinalScaling *= 2.0;
            }

            for( int i = 0; i < 3; i++ )
            {
                transformedCosineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) *= cosineFinalScaling;
                transformedSineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) *= sineFinalScaling;
            }
        }
    }
}

void FullTwoBodySphericalHarmonicTorque::computeTransformedAngularMomentumCoefficients(
        const Eigen::MatrixXd& cosineCoefficientsBody2,
        const Eigen::MatrixXd& sineCoefficientsBody2,
        const std::shared_ptr< basic_mathematics::WignerDMatricesCache >& wignerCache,
        const bool areCoefficientsNormalized,
        std::array< Eigen::MatrixXd, 3 >& transformedCosineCoefficientsBody2AngularMomentum,
        std::array< Eigen::MatrixXd, 3 >& transformedSineCoefficientsBody2AngularMomentum )
{
    // Forward to cache-based utility; kept as class entry point for reuse by analytical partials.
    // Eq. (67): this function provides the \hat{J}-mapped coefficients used in Eq. (60).
    computeTransformedAngularMomentumCoefficientsFromWignerCache( cosineCoefficientsBody2,
                                                                  sineCoefficientsBody2,
                                                                  wignerCache,
                                                                  areCoefficientsNormalized,
                                                                  transformedCosineCoefficientsBody2AngularMomentum,
                                                                  transformedSineCoefficientsBody2AngularMomentum );
}

void FullTwoBodySphericalHarmonicTorque::updateMembers( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        const bool printFigureFigureDebugOutput = isFigureFigureDebugEnabled( );

        // Step 1: synchronize the acceleration model state (effective coefficients, SH cache, relative geometry).
        accelerationBetweenBodies_->updateMembers( currentTime );

        // Step 2: retrieve current fields and precompute \hat{J}-transformed body-2 coefficients
        // used in torque Eq. (67) / Eq. (60) evaluation.
        const std::shared_ptr< EffectiveMutualSphericalHarmonicsField > effectiveMutualPotentialField =
                accelerationBetweenBodies_->getEffectiveMutualPotentialField( );
        const std::shared_ptr< basic_mathematics::WignerDMatricesCache > wignerCache =
                effectiveMutualPotentialField->getTransformationCache( )->getWignerDMatricesCache( );

        const Eigen::MatrixXd& cosineCoefficientsOfBody1 = effectiveMutualPotentialField->getCosineCoefficientsOfBody1( );
        const Eigen::MatrixXd& sineCoefficientsOfBody1 = effectiveMutualPotentialField->getSineCoefficientsOfBody1( );
        const Eigen::MatrixXd& cosineCoefficientsOfBody2 = effectiveMutualPotentialField->getCosineCoefficientsOfBody2( );
        const Eigen::MatrixXd& sineCoefficientsOfBody2 = effectiveMutualPotentialField->getSineCoefficientsOfBody2( );
        computeTransformedAngularMomentumCoefficients( cosineCoefficientsOfBody2,
                                                       sineCoefficientsOfBody2,
                                                       wignerCache,
                                                       accelerationBetweenBodies_->getAreCoefficientsNormalized( ),
                                                       transformedCosineCoefficientsBody2AngularMomentum_,
                                                       transformedSineCoefficientsBody2AngularMomentum_ );

        const Eigen::Vector3d bodyFixedRelativePosition = accelerationBetweenBodies_->getCurrentBodyFixedRelativePosition( );
        const double currentDistance = bodyFixedRelativePosition.norm( );
        const double preMultiplier = accelerationBetweenBodies_->getCurrentGravitationalParameter( ) / currentDistance;

        if( printFigureFigureDebugOutput )
        {
            printNativeFigureFigureTraceLine(
                    "100.native.input_scalars",
                    "mu_R1_R2_normalized",
                    "vector3",
                    formatDebugVector( Eigen::Vector3d( accelerationBetweenBodies_->getCurrentGravitationalParameter( ),
                                                        accelerationBetweenBodies_->getEquatorialRadiusOfBody1( ),
                                                        accelerationBetweenBodies_->getEquatorialRadiusOfBody2( ) ) ) );
            printNativeFigureFigureTraceLine( "101.native.relative_state",
                                              "relative_position_used_by_dmr",
                                              "vector3",
                                              formatDebugVector( bodyFixedRelativePosition ) );
            printNativeFigureFigureTraceLine(
                    "102.native.rotation",
                    "R_F2_to_F1",
                    "matrix3",
                    formatDebugMatrix( accelerationBetweenBodies_->getCurrentRotationFromBody2ToBody1( ).toRotationMatrix( ) ) );
            printNativeFigureFigureTraceLine( "103.native.body1_coefficients",
                                              "C20_C21_S21_C22_S22",
                                              "matrix3",
                                              formatDebugMatrix( ( Eigen::Matrix3d( ) << cosineCoefficientsOfBody1( 2, 0 ),
                                                                   cosineCoefficientsOfBody1( 2, 1 ),
                                                                   sineCoefficientsOfBody1( 2, 1 ),
                                                                   cosineCoefficientsOfBody1( 2, 2 ),
                                                                   sineCoefficientsOfBody1( 2, 2 ),
                                                                   0.0,
                                                                   0.0,
                                                                   0.0,
                                                                   0.0 )
                                                                         .finished( ) ) );
            printNativeFigureFigureTraceLine( "104.native.body2_coefficients",
                                              "C20_C21_S21_C22_S22",
                                              "matrix3",
                                              formatDebugMatrix( ( Eigen::Matrix3d( ) << cosineCoefficientsOfBody2( 2, 0 ),
                                                                   cosineCoefficientsOfBody2( 2, 1 ),
                                                                   sineCoefficientsOfBody2( 2, 1 ),
                                                                   cosineCoefficientsOfBody2( 2, 2 ),
                                                                   sineCoefficientsOfBody2( 2, 2 ),
                                                                   0.0,
                                                                   0.0,
                                                                   0.0,
                                                                   0.0 )
                                                                         .finished( ) ) );
            printNativeFigureFigureTraceLine(
                    "105.native.transformed_body2_coefficients",
                    "C20_C21_S21_C22_S22",
                    "matrix3",
                    formatDebugMatrix(
                            ( Eigen::Matrix3d( ) << effectiveMutualPotentialField->getTransformedCosineCoefficientsOfBody2( )( 2, 0 ),
                              effectiveMutualPotentialField->getTransformedCosineCoefficientsOfBody2( )( 2, 1 ),
                              effectiveMutualPotentialField->getTransformedSineCoefficientsOfBody2( )( 2, 1 ),
                              effectiveMutualPotentialField->getTransformedCosineCoefficientsOfBody2( )( 2, 2 ),
                              effectiveMutualPotentialField->getTransformedSineCoefficientsOfBody2( )( 2, 2 ),
                              0.0,
                              0.0,
                              0.0,
                              0.0 )
                                    .finished( ) ) );
            printNativeFigureFigureTraceLine(
                    "106.native.J_body2_coefficients.cosine",
                    "J_C20_C21_C22_rows_xyz",
                    "matrix3",
                    formatDebugMatrix( ( Eigen::Matrix3d( ) << transformedCosineCoefficientsBody2AngularMomentum_.at( 0 )( 2, 0 ),
                                         transformedCosineCoefficientsBody2AngularMomentum_.at( 0 )( 2, 1 ),
                                         transformedCosineCoefficientsBody2AngularMomentum_.at( 0 )( 2, 2 ),
                                         transformedCosineCoefficientsBody2AngularMomentum_.at( 1 )( 2, 0 ),
                                         transformedCosineCoefficientsBody2AngularMomentum_.at( 1 )( 2, 1 ),
                                         transformedCosineCoefficientsBody2AngularMomentum_.at( 1 )( 2, 2 ),
                                         transformedCosineCoefficientsBody2AngularMomentum_.at( 2 )( 2, 0 ),
                                         transformedCosineCoefficientsBody2AngularMomentum_.at( 2 )( 2, 1 ),
                                         transformedCosineCoefficientsBody2AngularMomentum_.at( 2 )( 2, 2 ) )
                                               .finished( ) ) );
            printNativeFigureFigureTraceLine(
                    "107.native.J_body2_coefficients.sine",
                    "J_S20_S21_S22_rows_xyz",
                    "matrix3",
                    formatDebugMatrix( ( Eigen::Matrix3d( ) << transformedSineCoefficientsBody2AngularMomentum_.at( 0 )( 2, 0 ),
                                         transformedSineCoefficientsBody2AngularMomentum_.at( 0 )( 2, 1 ),
                                         transformedSineCoefficientsBody2AngularMomentum_.at( 0 )( 2, 2 ),
                                         transformedSineCoefficientsBody2AngularMomentum_.at( 1 )( 2, 0 ),
                                         transformedSineCoefficientsBody2AngularMomentum_.at( 1 )( 2, 1 ),
                                         transformedSineCoefficientsBody2AngularMomentum_.at( 1 )( 2, 2 ),
                                         transformedSineCoefficientsBody2AngularMomentum_.at( 2 )( 2, 0 ),
                                         transformedSineCoefficientsBody2AngularMomentum_.at( 2 )( 2, 1 ),
                                         transformedSineCoefficientsBody2AngularMomentum_.at( 2 )( 2, 2 ) )
                                               .finished( ) ) );
        }

        const std::vector< double >& radius1Powers = accelerationBetweenBodies_->getRadius1Powers( );
        const std::vector< double >& radius2Powers = accelerationBetweenBodies_->getRadius2Powers( );
        const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache =
                accelerationBetweenBodies_->getSphericalHarmonicsCache( );

        // Step 3: evaluate body-2 spin torque contribution in frame F1 from Eq. (67), then apply Eq. (60).
        Eigen::Vector3d body2TorqueInBodyFixedFrameOfBody1 = Eigen::Vector3d::Zero( );
        for( unsigned int i = 0; i < body2TorqueCombinationsToUse_.size( ); i++ )
        {
            const int degreeOfBody1 = std::get< 0 >( body2TorqueCombinationsToUse_.at( i ) );
            const int orderOfBody1 = std::get< 1 >( body2TorqueCombinationsToUse_.at( i ) );
            const int degreeOfBody2 = std::get< 2 >( body2TorqueCombinationsToUse_.at( i ) );
            const int orderOfBody2 = std::get< 3 >( body2TorqueCombinationsToUse_.at( i ) );

            const double equatorialRadiusRatioPower = radius1Powers.at( degreeOfBody1 ) * radius2Powers.at( degreeOfBody2 );
            const int totalDegree = degreeOfBody1 + degreeOfBody2;

            // Expand each stored non-negative-order tuple into all signed-order combinations
            // consistent with Eq. (49)/(67) real-coefficient convention.
            for( int j = 0; j < 4; j++ )
            {
                int signedOrderOfBody1 = 0;
                int signedOrderOfBody2 = 0;
                const bool computeTerm =
                        getSignedOrdersForCombinationCase( j, orderOfBody1, orderOfBody2, signedOrderOfBody1, signedOrderOfBody2 );
                if( !computeTerm )
                {
                    continue;
                }

                const int totalOrder = std::abs( signedOrderOfBody1 + signedOrderOfBody2 );
                const double legendrePolynomial =
                        sphericalHarmonicsCache->getLegendreCache( ).getLegendrePolynomial( totalDegree, totalOrder );
                const double cosineOfMultipleLongitude = sphericalHarmonicsCache->getCosineOfMultipleLongitude( totalOrder );
                const double sineOfMultipleLongitude = sphericalHarmonicsCache->getSineOfMultipleLongitude( totalOrder );

                const double signOrderBody1 = ( signedOrderOfBody1 < 0 ) ? -1.0 : 1.0;
                const double signOrderBody2 = ( signedOrderOfBody2 < 0 ) ? -1.0 : 1.0;
                const double signTotalOrder = ( ( signedOrderOfBody1 + signedOrderOfBody2 ) < 0 ) ? -1.0 : 1.0;

                const double body1CosineCoefficient = cosineCoefficientsOfBody1( degreeOfBody1, std::abs( signedOrderOfBody1 ) );
                const double body1SineCoefficient = sineCoefficientsOfBody1( degreeOfBody1, std::abs( signedOrderOfBody1 ) );
                const double multiplier =
                        getMutualPotentialEffectiveCoefficientMultiplier( degreeOfBody1,
                                                                          signedOrderOfBody1,
                                                                          degreeOfBody2,
                                                                          signedOrderOfBody2,
                                                                          accelerationBetweenBodies_->getAreCoefficientsNormalized( ) );

                Eigen::Vector3d angularMomentumTransformedCosineCoefficientsBody2;
                Eigen::Vector3d angularMomentumTransformedSineCoefficientsBody2;
                for( int k = 0; k < 3; k++ )
                {
                    angularMomentumTransformedCosineCoefficientsBody2( k ) =
                            transformedCosineCoefficientsBody2AngularMomentum_.at( k )( degreeOfBody2, std::abs( signedOrderOfBody2 ) );
                    angularMomentumTransformedSineCoefficientsBody2( k ) =
                            transformedSineCoefficientsBody2AngularMomentum_.at( k )( degreeOfBody2, std::abs( signedOrderOfBody2 ) );
                }

                const Eigen::Vector3d effectiveAngularMomentumCosineCoefficients =
                        ( body1CosineCoefficient * angularMomentumTransformedCosineCoefficientsBody2 -
                          signOrderBody1 * signOrderBody2 * body1SineCoefficient * angularMomentumTransformedSineCoefficientsBody2 ) *
                        multiplier;
                const Eigen::Vector3d effectiveAngularMomentumSineCoefficients =
                        ( signOrderBody2 * body1CosineCoefficient * angularMomentumTransformedSineCoefficientsBody2 +
                          signOrderBody1 * body1SineCoefficient * angularMomentumTransformedCosineCoefficientsBody2 ) *
                        signTotalOrder * multiplier;

                // Scalar basis and effective angular-momentum coefficients form one Eq. (67) term.
                const Eigen::Vector3d currentEq67Contribution = equatorialRadiusRatioPower * legendrePolynomial *
                        ( effectiveAngularMomentumCosineCoefficients * cosineOfMultipleLongitude +
                          effectiveAngularMomentumSineCoefficients * sineOfMultipleLongitude );
                body2TorqueInBodyFixedFrameOfBody1 += currentEq67Contribution;

                if( printFigureFigureDebugOutput )
                {
                    printNativeFigureFigureTraceLine( "200.native.eq67.term." + std::to_string( i ) + "." + std::to_string( j ),
                                                      "l1_m1_l2_m2_msum_radiusPower_P_cos_sin",
                                                      "matrix3",
                                                      formatDebugMatrix( ( Eigen::Matrix3d( ) << static_cast< double >( degreeOfBody1 ),
                                                                           static_cast< double >( signedOrderOfBody1 ),
                                                                           static_cast< double >( degreeOfBody2 ),
                                                                           static_cast< double >( signedOrderOfBody2 ),
                                                                           static_cast< double >( signedOrderOfBody1 + signedOrderOfBody2 ),
                                                                           equatorialRadiusRatioPower,
                                                                           legendrePolynomial,
                                                                           cosineOfMultipleLongitude,
                                                                           sineOfMultipleLongitude )
                                                                                 .finished( ) ) );
                    printNativeFigureFigureTraceLine(
                            "201.native.eq67.scaling." + std::to_string( i ) + "." + std::to_string( j ),
                            "multiplier_body1C_body1S",
                            "vector3",
                            formatDebugVector( Eigen::Vector3d( multiplier, body1CosineCoefficient, body1SineCoefficient ) ) );
                    printNativeFigureFigureTraceLine(
                            "202.native.eq67.angularMomentum." + std::to_string( i ) + "." + std::to_string( j ),
                            "Jcos_Jsin_contribution",
                            "matrix3",
                            formatDebugMatrix( ( Eigen::Matrix3d( ) << angularMomentumTransformedCosineCoefficientsBody2( 0 ),
                                                 angularMomentumTransformedCosineCoefficientsBody2( 1 ),
                                                 angularMomentumTransformedCosineCoefficientsBody2( 2 ),
                                                 angularMomentumTransformedSineCoefficientsBody2( 0 ),
                                                 angularMomentumTransformedSineCoefficientsBody2( 1 ),
                                                 angularMomentumTransformedSineCoefficientsBody2( 2 ),
                                                 currentEq67Contribution( 0 ),
                                                 currentEq67Contribution( 1 ),
                                                 currentEq67Contribution( 2 ) )
                                                       .finished( ) ) );
                    printNativeFigureFigureTraceLine(
                            "203.native.eq67.effectiveAndRunning." + std::to_string( i ) + "." + std::to_string( j ),
                            "Ceff_Seff_running",
                            "matrix3",
                            formatDebugMatrix( ( Eigen::Matrix3d( ) << effectiveAngularMomentumCosineCoefficients( 0 ),
                                                 effectiveAngularMomentumCosineCoefficients( 1 ),
                                                 effectiveAngularMomentumCosineCoefficients( 2 ),
                                                 effectiveAngularMomentumSineCoefficients( 0 ),
                                                 effectiveAngularMomentumSineCoefficients( 1 ),
                                                 effectiveAngularMomentumSineCoefficients( 2 ),
                                                 body2TorqueInBodyFixedFrameOfBody1( 0 ),
                                                 body2TorqueInBodyFixedFrameOfBody1( 1 ),
                                                 body2TorqueInBodyFixedFrameOfBody1( 2 ) )
                                                       .finished( ) ) );
                }
            }
        }
        const Eigen::Vector3d body2TorqueEq67SumBeforePremultiplier = body2TorqueInBodyFixedFrameOfBody1;
        body2TorqueInBodyFixedFrameOfBody1 *= -preMultiplier;
        // Eq. (60): M_2 = -\hat{J}(V_1-2), with preMultiplier carrying the common -GM/r factor.

        // Step 4: compute total relative-frame torque from translational side using Eq. (68),
        // then isolate body-1 torque by subtracting body-2 contribution.
        const Eigen::Vector3d totalTorqueInBodyFixedFrameOfBody1 =
                bodyFixedRelativePosition.cross( accelerationBetweenBodies_->getMutualPotentialGradient( ) );
        const Eigen::Vector3d body1TorqueInBodyFixedFrameOfBody1 = totalTorqueInBodyFixedFrameOfBody1 - body2TorqueInBodyFixedFrameOfBody1;

        if( printFigureFigureDebugOutput )
        {
            printNativeFigureFigureTraceLine( "300.native.eq67_sum",
                                              "sum_before_premultiplier",
                                              "vector3",
                                              formatDebugVector( body2TorqueEq67SumBeforePremultiplier ) );
            printNativeFigureFigureTraceLine( "301.native.eq67_premultiplied_body2",
                                              "body2_torque_in_F1",
                                              "vector3",
                                              formatDebugVector( body2TorqueInBodyFixedFrameOfBody1 ) );
            printNativeFigureFigureTraceLine(
                    "302.native.eq68_balance",
                    "gradient_total_body1",
                    "matrix3",
                    formatDebugMatrix( ( Eigen::Matrix3d( ) << accelerationBetweenBodies_->getMutualPotentialGradient( )( 0 ),
                                         accelerationBetweenBodies_->getMutualPotentialGradient( )( 1 ),
                                         accelerationBetweenBodies_->getMutualPotentialGradient( )( 2 ),
                                         totalTorqueInBodyFixedFrameOfBody1( 0 ),
                                         totalTorqueInBodyFixedFrameOfBody1( 1 ),
                                         totalTorqueInBodyFixedFrameOfBody1( 2 ),
                                         body1TorqueInBodyFixedFrameOfBody1( 0 ),
                                         body1TorqueInBodyFixedFrameOfBody1( 1 ),
                                         body1TorqueInBodyFixedFrameOfBody1( 2 ) )
                                               .finished( ) ) );
        }

        // Step 5: return requested body's torque, applying frame mapping for body 2 using Eq. (69).
        if( acceleratedBodyIsBody1_ )
        {
            currentTorque_ = -body1TorqueInBodyFixedFrameOfBody1;
        }
        else
        {
            currentTorque_ =
                    -( accelerationBetweenBodies_->getCurrentRotationFromBody2ToBody1( ).inverse( ) * body2TorqueInBodyFixedFrameOfBody1 );
        }

        if( printFigureFigureDebugOutput )
        {
            const Eigen::Vector3d canonicalRelativePositionInBody1Frame =
                    parseDebugVector3( "TUDAT_FFDBG_RELATIVE_POSITION_F1", -bodyFixedRelativePosition );
            const Eigen::Vector3d requestedOutputTorqueOnBody1InBody1Frame = acceleratedBodyIsBody1_
                    ? currentTorque_
                    : accelerationBetweenBodies_->getCurrentRotationFromBody2ToBody1( ) * ( -currentTorque_ );
            emitCanonicalFigureFigureTrace( canonicalRelativePositionInBody1Frame,
                                            cosineCoefficientsOfBody1,
                                            sineCoefficientsOfBody1,
                                            cosineCoefficientsOfBody2,
                                            sineCoefficientsOfBody2,
                                            accelerationBetweenBodies_->getEquatorialRadiusOfBody1( ),
                                            accelerationBetweenBodies_->getEquatorialRadiusOfBody2( ),
                                            requestedOutputTorqueOnBody1InBody1Frame );

            printNativeFigureFigureTraceLine( "400.native.final_converted_torque",
                                              "torque_on_body1_in_F1",
                                              "vector3",
                                              formatDebugVector( requestedOutputTorqueOnBody1InBody1Frame ) );
        }

        currentTime_ = currentTime;
    }
}

}  // namespace gravitation

}  // namespace tudat
