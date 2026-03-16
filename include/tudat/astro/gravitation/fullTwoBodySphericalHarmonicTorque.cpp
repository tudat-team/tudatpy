#include "tudat/astro/gravitation/fullTwoBodySphericalHarmonicTorque.h"

#include <boost/math/special_functions/factorials.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace tudat
{

namespace gravitation
{

namespace
{

bool isFullTwoBodyTorqueDebugEnabled( )
{
    static const bool isEnabled = [ ]( )
    {
        const char* flag = std::getenv( "TUDAT_DEBUG_FULL_TWO_BODY_TORQUE" );
        return ( flag != nullptr && std::string( flag ) != "0" );
    }( );
    return isEnabled;
}

std::string debugStatusFromVectorDifference( const Eigen::Vector3d& difference, const double tolerance )
{
    return ( difference.norm( ) <= tolerance ? "OK" : "MISMATCH" );
}

bool isWignerCacheIdentityAtDegree(
        const std::shared_ptr< basic_mathematics::WignerDMatricesCache >& wignerCache,
        const int degree,
        const double tolerance )
{
    for( int orderM = -degree; orderM <= degree; orderM++ )
    {
        for( int orderK = -degree; orderK <= degree; orderK++ )
        {
            const std::complex< double > coefficient = wignerCache->getWignerDCoefficient( degree, orderM, orderK );
            const double expectedReal = ( orderM == orderK ? 1.0 : 0.0 );
            if( std::fabs( coefficient.real( ) - expectedReal ) > tolerance ||
                std::fabs( coefficient.imag( ) ) > tolerance )
            {
                return false;
            }
        }
    }
    return true;
}

bool isBody1C20Only(
        const Eigen::MatrixXd& cosineCoefficientsOfBody1,
        const Eigen::MatrixXd& sineCoefficientsOfBody1,
        const double tolerance = 1.0E-30 )
{
    if( cosineCoefficientsOfBody1.rows( ) <= 2 || cosineCoefficientsOfBody1.cols( ) <= 0 ||
        sineCoefficientsOfBody1.rows( ) <= 2 || sineCoefficientsOfBody1.cols( ) <= 0 )
    {
        return false;
    }

    for( int row = 0; row < cosineCoefficientsOfBody1.rows( ); row++ )
    {
        for( int col = 0; col < cosineCoefficientsOfBody1.cols( ); col++ )
        {
            if( row == 0 && col == 0 )
            {
                continue;
            }
            if( row == 2 && col == 0 )
            {
                continue;
            }
            if( std::fabs( cosineCoefficientsOfBody1( row, col ) ) > tolerance )
            {
                return false;
            }
        }
    }

    for( int row = 0; row < sineCoefficientsOfBody1.rows( ); row++ )
    {
        for( int col = 0; col < sineCoefficientsOfBody1.cols( ); col++ )
        {
            if( std::fabs( sineCoefficientsOfBody1( row, col ) ) > tolerance )
            {
                return false;
            }
        }
    }

    return true;
}

enum class DegreeTwoCoefficientCase
{
    none,
    c20,
    c21,
    s21,
    c22,
    s22,
    multiple
};

std::string coefficientCaseToString( const DegreeTwoCoefficientCase coefficientCase )
{
    switch( coefficientCase )
    {
    case DegreeTwoCoefficientCase::none:
        return "none";
    case DegreeTwoCoefficientCase::c20:
        return "C20";
    case DegreeTwoCoefficientCase::c21:
        return "C21";
    case DegreeTwoCoefficientCase::s21:
        return "S21";
    case DegreeTwoCoefficientCase::c22:
        return "C22";
    case DegreeTwoCoefficientCase::s22:
        return "S22";
    case DegreeTwoCoefficientCase::multiple:
        return "multiple";
    default:
        return "unknown";
    }
}

DegreeTwoCoefficientCase detectSingleActiveDegreeTwoCoefficientCase(
        const Eigen::MatrixXd& transformedCosineCoefficientsOfBody2,
        const Eigen::MatrixXd& transformedSineCoefficientsOfBody2,
        double& coefficientValue,
        const double tolerance = 1.0E-30 )
{
    coefficientValue = 0.0;
    std::vector< std::pair< DegreeTwoCoefficientCase, double > > nonZeroCoefficients;

    if( transformedCosineCoefficientsOfBody2.rows( ) > 2 && transformedCosineCoefficientsOfBody2.cols( ) > 0 )
    {
        if( std::fabs( transformedCosineCoefficientsOfBody2( 2, 0 ) ) > tolerance )
        {
            nonZeroCoefficients.emplace_back( DegreeTwoCoefficientCase::c20, transformedCosineCoefficientsOfBody2( 2, 0 ) );
        }
    }

    if( transformedCosineCoefficientsOfBody2.rows( ) > 2 && transformedCosineCoefficientsOfBody2.cols( ) > 1 &&
        transformedSineCoefficientsOfBody2.rows( ) > 2 && transformedSineCoefficientsOfBody2.cols( ) > 1 )
    {
        if( std::fabs( transformedCosineCoefficientsOfBody2( 2, 1 ) ) > tolerance )
        {
            nonZeroCoefficients.emplace_back( DegreeTwoCoefficientCase::c21, transformedCosineCoefficientsOfBody2( 2, 1 ) );
        }
        if( std::fabs( transformedSineCoefficientsOfBody2( 2, 1 ) ) > tolerance )
        {
            nonZeroCoefficients.emplace_back( DegreeTwoCoefficientCase::s21, transformedSineCoefficientsOfBody2( 2, 1 ) );
        }
    }

    if( transformedCosineCoefficientsOfBody2.rows( ) > 2 && transformedCosineCoefficientsOfBody2.cols( ) > 2 &&
        transformedSineCoefficientsOfBody2.rows( ) > 2 && transformedSineCoefficientsOfBody2.cols( ) > 2 )
    {
        if( std::fabs( transformedCosineCoefficientsOfBody2( 2, 2 ) ) > tolerance )
        {
            nonZeroCoefficients.emplace_back( DegreeTwoCoefficientCase::c22, transformedCosineCoefficientsOfBody2( 2, 2 ) );
        }
        if( std::fabs( transformedSineCoefficientsOfBody2( 2, 2 ) ) > tolerance )
        {
            nonZeroCoefficients.emplace_back( DegreeTwoCoefficientCase::s22, transformedSineCoefficientsOfBody2( 2, 2 ) );
        }
    }

    if( nonZeroCoefficients.empty( ) )
    {
        return DegreeTwoCoefficientCase::none;
    }
    if( nonZeroCoefficients.size( ) > 1 )
    {
        return DegreeTwoCoefficientCase::multiple;
    }
    coefficientValue = nonZeroCoefficients.at( 0 ).second;
    return nonZeroCoefficients.at( 0 ).first;
}

Eigen::Vector3d computeExpectedSpecificC20DegreeTwoTorqueFromDocument(
        const Eigen::Vector3d& relativePositionInBody1Frame,
        const double gravitationalParameterExertingBody,
        const double referenceRadiusBody1,
        const double referenceRadiusBody2,
        const double normalizedC20OfBody1,
        const DegreeTwoCoefficientCase coefficientCase,
        const double normalizedCoefficientValueOfBody2 )
{
    const double x = relativePositionInBody1Frame( 0 );
    const double y = relativePositionInBody1Frame( 1 );
    const double z = relativePositionInBody1Frame( 2 );
    const double r2 = relativePositionInBody1Frame.squaredNorm( );
    const double r = std::sqrt( r2 );
    const double x2 = x * x;
    const double y2 = y * y;
    const double z2 = z * z;
    const double x4 = x2 * x2;
    const double y4 = y2 * y2;
    const double z4 = z2 * z2;
    const double sqrtThree = std::sqrt( 3.0 );

    const double commonPrefactor =
            gravitationalParameterExertingBody *
            referenceRadiusBody1 * referenceRadiusBody1 * referenceRadiusBody2 * referenceRadiusBody2 /
            std::pow( r, 9.0 ) * normalizedC20OfBody1;

    Eigen::Vector3d analyticalSpecificTorque = Eigen::Vector3d::Zero( );
    switch( coefficientCase )
    {
    case DegreeTwoCoefficientCase::c20:
        analyticalSpecificTorque << y * z * ( 3.0 * r2 - 7.0 * z2 ),
                -x * z * ( 3.0 * r2 - 7.0 * z2 ),
                0.0;
        analyticalSpecificTorque *= 75.0 / 2.0 * commonPrefactor * normalizedCoefficientValueOfBody2;
        break;
    case DegreeTwoCoefficientCase::c21:
        analyticalSpecificTorque << x * y * ( x2 + y2 - 6.0 * z2 ),
                -0.2 * ( 4.0 * x4 + 3.0 * x2 * y2 - 27.0 * x2 * z2 - y4 + 3.0 * y2 * z2 + 4.0 * z4 ),
                0.0;
        analyticalSpecificTorque *= 25.0 * sqrtThree * commonPrefactor * normalizedCoefficientValueOfBody2;
        break;
    case DegreeTwoCoefficientCase::s21:
        analyticalSpecificTorque << -0.2 * ( x4 - 3.0 * x2 * y2 - 3.0 * x2 * z2 - 4.0 * y4 + 27.0 * y2 * z2 - 4.0 * z4 ),
                -x * y * ( x2 + y2 - 6.0 * z2 ),
                0.0;
        analyticalSpecificTorque *= 25.0 * sqrtThree * commonPrefactor * normalizedCoefficientValueOfBody2;
        break;
    case DegreeTwoCoefficientCase::c22:
        analyticalSpecificTorque << -y * z * ( 9.0 * x2 - 5.0 * y2 + 2.0 * z2 ),
                x * z * ( 5.0 * x2 - 9.0 * y2 - 2.0 * z2 ),
                0.0;
        analyticalSpecificTorque *= 25.0 * sqrtThree / 2.0 * commonPrefactor * normalizedCoefficientValueOfBody2;
        break;
    case DegreeTwoCoefficientCase::s22:
        analyticalSpecificTorque << x * z * ( x2 - 6.0 * y2 + z2 ),
                y * z * ( 6.0 * x2 - y2 - z2 ),
                0.0;
        analyticalSpecificTorque *= 25.0 * sqrtThree * commonPrefactor * normalizedCoefficientValueOfBody2;
        break;
    default:
        break;
    }

    return analyticalSpecificTorque;
}

}

Eigen::Vector3cd FullTwoBodySphericalHarmonicTorque::computeAngularMomentumOperatorOnWignerCoefficient(
        const std::shared_ptr< basic_mathematics::WignerDMatricesCache >& wignerCache,
        const int degree,
        const int orderM,
        const int orderK )
{
    const std::complex< double > imaginaryUnit( 0.0, 1.0 );
    const double inverseSquareRootTwo = 1.0 / std::sqrt( 2.0 );

    const auto getWignerCoefficient = [ & ]( const int requestedOrderM, const int requestedOrderK )
    {
        if( std::abs( requestedOrderM ) > degree || std::abs( requestedOrderK ) > degree )
        {
            return std::complex< double >( 0.0, 0.0 );
        }
        return wignerCache->getWignerDCoefficient( degree, requestedOrderM, requestedOrderK );
    };

    const double plusScaling = std::sqrt(
                std::max( 0.0, static_cast< double >(
                              degree * ( degree + 1 ) - orderM * ( orderM - 1 ) ) ) ) / std::sqrt( 2.0 );
    const double minusScaling = std::sqrt(
                std::max( 0.0, static_cast< double >(
                              degree * ( degree + 1 ) - orderM * ( orderM + 1 ) ) ) ) / std::sqrt( 2.0 );

    const std::complex< double > angularMomentumPlus =
            imaginaryUnit * plusScaling * getWignerCoefficient( orderM - 1, orderK );
    const std::complex< double > angularMomentumMinus =
            imaginaryUnit * ( -minusScaling ) * getWignerCoefficient( orderM + 1, orderK );
    const std::complex< double > angularMomentumZero =
            imaginaryUnit * static_cast< double >( -orderM ) * getWignerCoefficient( orderM, orderK );

    Eigen::Vector3cd angularMomentumInCartesianBasis;
    angularMomentumInCartesianBasis( 0 ) = ( angularMomentumMinus - angularMomentumPlus ) * inverseSquareRootTwo;
    angularMomentumInCartesianBasis( 1 ) = imaginaryUnit * ( angularMomentumMinus + angularMomentumPlus ) * inverseSquareRootTwo;
    angularMomentumInCartesianBasis( 2 ) = angularMomentumZero;

    if( isFullTwoBodyTorqueDebugEnabled( ) && degree == 2 && std::abs( orderM ) <= 2 && std::abs( orderK ) <= 2 &&
        isWignerCacheIdentityAtDegree( wignerCache, degree, 1.0E-14 ) )
    {
        const double plusScalingOrderM = std::sqrt(
                    std::max( 0.0, static_cast< double >(
                                  degree * ( degree + 1 ) - orderM * ( orderM - 1 ) ) ) ) / std::sqrt( 2.0 );
        const double minusScalingOrderM = std::sqrt(
                    std::max( 0.0, static_cast< double >(
                                  degree * ( degree + 1 ) - orderM * ( orderM + 1 ) ) ) ) / std::sqrt( 2.0 );
        const std::complex< double > angularMomentumPlusOrderM =
                imaginaryUnit * plusScalingOrderM * getWignerCoefficient( orderM - 1, orderK );
        const std::complex< double > angularMomentumMinusOrderM =
                imaginaryUnit * ( -minusScalingOrderM ) * getWignerCoefficient( orderM + 1, orderK );
        const std::complex< double > angularMomentumZeroOrderM =
                imaginaryUnit * static_cast< double >( -orderM ) * getWignerCoefficient( orderM, orderK );

        Eigen::Vector3cd expectedAngularMomentumInCartesianBasis;
        expectedAngularMomentumInCartesianBasis( 0 ) =
                ( angularMomentumMinusOrderM - angularMomentumPlusOrderM ) * inverseSquareRootTwo;
        expectedAngularMomentumInCartesianBasis( 1 ) =
                imaginaryUnit * ( angularMomentumMinusOrderM + angularMomentumPlusOrderM ) * inverseSquareRootTwo;
        expectedAngularMomentumInCartesianBasis( 2 ) = angularMomentumZeroOrderM;

        const Eigen::Vector3d conventionDifferenceReal =
                ( angularMomentumInCartesianBasis - expectedAngularMomentumInCartesianBasis ).real( );
        const Eigen::Vector3d conventionDifferenceImaginary =
                ( angularMomentumInCartesianBasis - expectedAngularMomentumInCartesianBasis ).imag( );
        std::cout << "[FTB-DBG][STEP J_operator]"
                  << " degree=" << degree
                  << " m=" << orderM
                  << " k=" << orderK
                  << " current_real=" << angularMomentumInCartesianBasis.real( ).transpose( )
                  << " current_imag=" << angularMomentumInCartesianBasis.imag( ).transpose( )
                  << " expected_m_index_real=" << expectedAngularMomentumInCartesianBasis.real( ).transpose( )
                  << " expected_m_index_imag=" << expectedAngularMomentumInCartesianBasis.imag( ).transpose( )
                  << " diff_real=" << conventionDifferenceReal.transpose( )
                  << " diff_imag=" << conventionDifferenceImaginary.transpose( )
                  << " status_real=" << debugStatusFromVectorDifference( conventionDifferenceReal, 1.0E-15 )
                  << " status_imag=" << debugStatusFromVectorDifference( conventionDifferenceImaginary, 1.0E-15 )
                  << std::endl;
    }

    return angularMomentumInCartesianBasis;
}

void FullTwoBodySphericalHarmonicTorque::computeTransformedAngularMomentumCoefficients(
        const Eigen::MatrixXd& cosineCoefficientsBody2,
        const Eigen::MatrixXd& sineCoefficientsBody2,
        const std::shared_ptr< basic_mathematics::WignerDMatricesCache >& wignerCache,
        const bool areCoefficientsNormalized,
        std::array< Eigen::MatrixXd, 3 >& transformedCosineCoefficientsBody2AngularMomentum,
        std::array< Eigen::MatrixXd, 3 >& transformedSineCoefficientsBody2AngularMomentum )
{
    for( int i = 0; i < 3; i++ )
    {
        transformedCosineCoefficientsBody2AngularMomentum.at( i ).setZero(
                    cosineCoefficientsBody2.rows( ), cosineCoefficientsBody2.cols( ) );
        transformedSineCoefficientsBody2AngularMomentum.at( i ).setZero(
                    sineCoefficientsBody2.rows( ), sineCoefficientsBody2.cols( ) );
    }

    for( int degree = 0; degree < cosineCoefficientsBody2.rows( ); degree++ )
    {
        for( int orderM = 0; ( orderM <= degree && orderM < cosineCoefficientsBody2.cols( ) ); orderM++ )
        {
            double orderMMultiplier;
            if( !areCoefficientsNormalized )
            {
                orderMMultiplier = std::sqrt(
                            boost::math::factorial< double >( degree - orderM ) /
                            boost::math::factorial< double >( degree + orderM ) );
            }
            else
            {
                orderMMultiplier = ( orderM == 0 ? 1.0 : 1.0 / std::sqrt( 2.0 ) );
            }

            const Eigen::Vector3cd orderZeroAngularMomentumD = computeAngularMomentumOperatorOnWignerCoefficient(
                        wignerCache, degree, orderM, 0 );
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
                    currentMultiplier = std::sqrt(
                                boost::math::factorial< double >( degree + orderK ) /
                                boost::math::factorial< double >( degree - orderK ) ) * orderMMultiplier;
                }
                else
                {
                    currentMultiplier = std::sqrt( 2.0 ) * orderMMultiplier;
                }

                const double signMultiplier = ( ( ( orderK % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );
                const Eigen::Vector3cd positiveOrderAngularMomentumD = computeAngularMomentumOperatorOnWignerCoefficient(
                            wignerCache, degree, orderM, orderK );
                const Eigen::Vector3cd negativeOrderAngularMomentumD = computeAngularMomentumOperatorOnWignerCoefficient(
                            wignerCache, degree, orderM, -orderK );

                for( int i = 0; i < 3; i++ )
                {
                    transformedCosineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) +=
                            0.5 * currentMultiplier *
                            ( ( signMultiplier * positiveOrderAngularMomentumD( i ).real( ) + negativeOrderAngularMomentumD( i ).real( ) ) *
                              cosineCoefficientsBody2( degree, orderK ) +
                              ( signMultiplier * positiveOrderAngularMomentumD( i ).imag( ) - negativeOrderAngularMomentumD( i ).imag( ) ) *
                              sineCoefficientsBody2( degree, orderK ) );

                    transformedSineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) +=
                            0.5 * currentMultiplier *
                            ( ( signMultiplier * positiveOrderAngularMomentumD( i ).imag( ) + negativeOrderAngularMomentumD( i ).imag( ) ) *
                              cosineCoefficientsBody2( degree, orderK ) +
                              ( -signMultiplier * positiveOrderAngularMomentumD( i ).real( ) + negativeOrderAngularMomentumD( i ).real( ) ) *
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

    if( isFullTwoBodyTorqueDebugEnabled( ) &&
        cosineCoefficientsBody2.rows( ) > 2 && cosineCoefficientsBody2.cols( ) > 2 )
    {
        const double finiteDifferenceStep = 1.0E-8;
        const std::array< Eigen::Vector3d, 3 > cartesianAxes = {
            Eigen::Vector3d::UnitX( ), Eigen::Vector3d::UnitY( ), Eigen::Vector3d::UnitZ( ) };

        for( int component = 0; component < 3; component++ )
        {
            basic_mathematics::SphericalHarmonicTransformationCache finiteDifferenceCache(
                        static_cast< int >( cosineCoefficientsBody2.rows( ) - 1 ),
                        static_cast< int >( cosineCoefficientsBody2.cols( ) - 1 ) );
            Eigen::MatrixXd cosinePlus;
            Eigen::MatrixXd sinePlus;
            Eigen::MatrixXd cosineMinus;
            Eigen::MatrixXd sineMinus;

            const Eigen::Quaterniond plusRotation( Eigen::AngleAxisd(
                    finiteDifferenceStep, cartesianAxes.at( component ) ) );
            finiteDifferenceCache.updateFromQuaternion( plusRotation );
            finiteDifferenceCache.transformCoefficientsAtDegree(
                        cosineCoefficientsBody2, sineCoefficientsBody2,
                        cosinePlus, sinePlus, areCoefficientsNormalized );

            const Eigen::Quaterniond minusRotation( Eigen::AngleAxisd(
                    -finiteDifferenceStep, cartesianAxes.at( component ) ) );
            finiteDifferenceCache.updateFromQuaternion( minusRotation );
            finiteDifferenceCache.transformCoefficientsAtDegree(
                        cosineCoefficientsBody2, sineCoefficientsBody2,
                        cosineMinus, sineMinus, areCoefficientsNormalized );

            for( int order = 0; order <= 2; order++ )
            {
                const double finiteDifferenceCosine =
                        ( cosinePlus( 2, order ) - cosineMinus( 2, order ) ) /
                        ( 2.0 * finiteDifferenceStep );
                const double finiteDifferenceSine =
                        ( sinePlus( 2, order ) - sineMinus( 2, order ) ) /
                        ( 2.0 * finiteDifferenceStep );
                const double cosineDifference =
                        transformedCosineCoefficientsBody2AngularMomentum.at( component )( 2, order ) - finiteDifferenceCosine;
                const double sineDifference =
                        transformedSineCoefficientsBody2AngularMomentum.at( component )( 2, order ) - finiteDifferenceSine;

                std::cout << "[FTB-DBG][STEP J_fd_check]"
                          << " component=" << component
                          << " (l,m)=(2," << order << ")"
                          << " JC=" << transformedCosineCoefficientsBody2AngularMomentum.at( component )( 2, order )
                          << " dC_dtheta_fd=" << finiteDifferenceCosine
                          << " JC_minus_fd=" << cosineDifference
                          << " JS=" << transformedSineCoefficientsBody2AngularMomentum.at( component )( 2, order )
                          << " dS_dtheta_fd=" << finiteDifferenceSine
                          << " JS_minus_fd=" << sineDifference
                          << std::endl;
            }
        }

        for( int component = 0; component < 3; component++ )
        {
            std::cout << "[FTB-DBG][STEP J_coefficients]"
                      << " component=" << component
                      << " JC(2,0..2)="
                      << transformedCosineCoefficientsBody2AngularMomentum.at( component )( 2, 0 ) << " "
                      << transformedCosineCoefficientsBody2AngularMomentum.at( component )( 2, 1 ) << " "
                      << transformedCosineCoefficientsBody2AngularMomentum.at( component )( 2, 2 )
                      << " JS(2,0..2)="
                      << transformedSineCoefficientsBody2AngularMomentum.at( component )( 2, 0 ) << " "
                      << transformedSineCoefficientsBody2AngularMomentum.at( component )( 2, 1 ) << " "
                      << transformedSineCoefficientsBody2AngularMomentum.at( component )( 2, 2 )
                      << std::endl;
        }
    }
}

void FullTwoBodySphericalHarmonicTorque::updateMembers( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        accelerationBetweenBodies_->updateMembers( currentTime );

        const std::shared_ptr< EffectiveMutualSphericalHarmonicsField > effectiveMutualPotentialField =
                accelerationBetweenBodies_->getEffectiveMutualPotentialField( );
        const std::shared_ptr< basic_mathematics::WignerDMatricesCache > wignerCache =
                effectiveMutualPotentialField->getTransformationCache( )->getWignerDMatricesCache( );

        const Eigen::MatrixXd cosineCoefficientsOfBody1 = effectiveMutualPotentialField->getCosineCoefficientsOfBody1( );
        const Eigen::MatrixXd sineCoefficientsOfBody1 = effectiveMutualPotentialField->getSineCoefficientsOfBody1( );
        const Eigen::MatrixXd cosineCoefficientsOfBody2 = effectiveMutualPotentialField->getCosineCoefficientsOfBody2( );
        const Eigen::MatrixXd sineCoefficientsOfBody2 = effectiveMutualPotentialField->getSineCoefficientsOfBody2( );
        const Eigen::MatrixXd transformedCosineCoefficientsOfBody2 =
                effectiveMutualPotentialField->getTransformedCosineCoefficientsOfBody2( );
        const Eigen::MatrixXd transformedSineCoefficientsOfBody2 =
                effectiveMutualPotentialField->getTransformedSineCoefficientsOfBody2( );

        std::array< Eigen::MatrixXd, 3 > transformedCosineCoefficientsBody2AngularMomentum;
        std::array< Eigen::MatrixXd, 3 > transformedSineCoefficientsBody2AngularMomentum;
        computeTransformedAngularMomentumCoefficients(
                    cosineCoefficientsOfBody2,
                    sineCoefficientsOfBody2,
                    wignerCache,
                    accelerationBetweenBodies_->getAreCoefficientsNormalized( ),
                    transformedCosineCoefficientsBody2AngularMomentum,
                    transformedSineCoefficientsBody2AngularMomentum );

        const Eigen::Vector3d bodyFixedRelativePosition = accelerationBetweenBodies_->getCurrentBodyFixedRelativePosition( );
        const double currentDistance = bodyFixedRelativePosition.norm( );
        const double preMultiplier = accelerationBetweenBodies_->getCurrentGravitationalParameter( ) / currentDistance;

        const std::vector< double > radius1Powers = accelerationBetweenBodies_->getRadius1Powers( );
        const std::vector< double > radius2Powers = accelerationBetweenBodies_->getRadius2Powers( );
        const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache =
                accelerationBetweenBodies_->getSphericalHarmonicsCache( );

        std::set< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > body2TorqueCombinationSet;
        for( const auto& combination : coefficientCombinationsToUse_ )
        {
            const unsigned int degreeOfBody1ForCombination = std::get< 0 >( combination );
            const unsigned int orderOfBody1ForCombination = std::get< 1 >( combination );
            const unsigned int degreeOfBody2ForCombination = std::get< 2 >( combination );
            for( unsigned int orderOfBody2ForTorque = 0; orderOfBody2ForTorque <= degreeOfBody2ForCombination; orderOfBody2ForTorque++ )
            {
                body2TorqueCombinationSet.insert( std::make_tuple(
                        degreeOfBody1ForCombination,
                        orderOfBody1ForCombination,
                        degreeOfBody2ForCombination,
                        orderOfBody2ForTorque ) );
            }
        }
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >
                body2TorqueCombinationsToUse(
                    body2TorqueCombinationSet.begin( ),
                    body2TorqueCombinationSet.end( ) );

        Eigen::Vector3d body2TorqueInBodyFixedFrameOfBody1 = Eigen::Vector3d::Zero( );
        for( unsigned int i = 0; i < body2TorqueCombinationsToUse.size( ); i++ )
        {
            const int degreeOfBody1 = std::get< 0 >( body2TorqueCombinationsToUse.at( i ) );
            const int orderOfBody1 = std::get< 1 >( body2TorqueCombinationsToUse.at( i ) );
            const int degreeOfBody2 = std::get< 2 >( body2TorqueCombinationsToUse.at( i ) );
            const int orderOfBody2 = std::get< 3 >( body2TorqueCombinationsToUse.at( i ) );

            const double equatorialRadiusRatioPower =
                    radius1Powers.at( degreeOfBody1 ) * radius2Powers.at( degreeOfBody2 );
            const int totalDegree = degreeOfBody1 + degreeOfBody2;

            for( unsigned int j = 0; j < 4; j++ )
            {
                int signedOrderOfBody1 = orderOfBody1;
                int signedOrderOfBody2 = orderOfBody2;
                bool computeTerm = true;

                switch( j )
                {
                case 0:
                    signedOrderOfBody1 = orderOfBody1;
                    signedOrderOfBody2 = orderOfBody2;
                    break;
                case 1:
                    signedOrderOfBody1 = -orderOfBody1;
                    signedOrderOfBody2 = orderOfBody2;
                    computeTerm = ( orderOfBody1 != 0 );
                    break;
                case 2:
                    signedOrderOfBody1 = orderOfBody1;
                    signedOrderOfBody2 = -orderOfBody2;
                    computeTerm = ( orderOfBody2 != 0 );
                    break;
                case 3:
                    signedOrderOfBody1 = -orderOfBody1;
                    signedOrderOfBody2 = -orderOfBody2;
                    computeTerm = ( orderOfBody1 != 0 && orderOfBody2 != 0 );
                    break;
                default:
                    break;
                }

                if( !computeTerm )
                {
                    continue;
                }

                const int totalOrder = std::abs( signedOrderOfBody1 + signedOrderOfBody2 );
                const double legendrePolynomial = sphericalHarmonicsCache->getLegendreCache( ).getLegendrePolynomial(
                            totalDegree, totalOrder );
                const double cosineOfMultipleLongitude = sphericalHarmonicsCache->getCosineOfMultipleLongitude( totalOrder );
                const double sineOfMultipleLongitude = sphericalHarmonicsCache->getSineOfMultipleLongitude( totalOrder );

                const double signOrderBody1 = ( signedOrderOfBody1 < 0 ) ? -1.0 : 1.0;
                const double signOrderBody2 = ( signedOrderOfBody2 < 0 ) ? -1.0 : 1.0;
                const double signTotalOrder = ( ( signedOrderOfBody1 + signedOrderOfBody2 ) < 0 ) ? -1.0 : 1.0;

                const double body1CosineCoefficient = cosineCoefficientsOfBody1(
                            degreeOfBody1, std::abs( signedOrderOfBody1 ) );
                const double body1SineCoefficient = sineCoefficientsOfBody1(
                            degreeOfBody1, std::abs( signedOrderOfBody1 ) );
                const double multiplier = getMutualPotentialEffectiveCoefficientMultiplier(
                            degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2,
                            accelerationBetweenBodies_->getAreCoefficientsNormalized( ) );

                Eigen::Vector3d angularMomentumTransformedCosineCoefficientsBody2;
                Eigen::Vector3d angularMomentumTransformedSineCoefficientsBody2;
                for( int k = 0; k < 3; k++ )
                {
                    angularMomentumTransformedCosineCoefficientsBody2( k ) =
                            transformedCosineCoefficientsBody2AngularMomentum.at( k )( degreeOfBody2, std::abs( signedOrderOfBody2 ) );
                    angularMomentumTransformedSineCoefficientsBody2( k ) =
                            transformedSineCoefficientsBody2AngularMomentum.at( k )( degreeOfBody2, std::abs( signedOrderOfBody2 ) );
                }

                const Eigen::Vector3d effectiveAngularMomentumCosineCoefficients =
                        ( body1CosineCoefficient * angularMomentumTransformedCosineCoefficientsBody2 -
                          signOrderBody1 * signOrderBody2 * body1SineCoefficient * angularMomentumTransformedSineCoefficientsBody2 ) *
                        multiplier;
                const Eigen::Vector3d effectiveAngularMomentumSineCoefficients =
                        ( signOrderBody2 * body1CosineCoefficient * angularMomentumTransformedSineCoefficientsBody2 +
                          signOrderBody1 * body1SineCoefficient * angularMomentumTransformedCosineCoefficientsBody2 ) *
                        signTotalOrder * multiplier;

                body2TorqueInBodyFixedFrameOfBody1 +=
                        equatorialRadiusRatioPower * legendrePolynomial *
                        ( effectiveAngularMomentumCosineCoefficients * cosineOfMultipleLongitude +
                          effectiveAngularMomentumSineCoefficients * sineOfMultipleLongitude );

                if( isFullTwoBodyTorqueDebugEnabled( ) &&
                    degreeOfBody1 == 2 && orderOfBody1 == 0 && degreeOfBody2 == 2 && std::abs( orderOfBody2 ) <= 2 )
                {
                    std::cout << "[FTB-DBG][STEP body2_torque_term]"
                              << " base=(l1,m1,l2,m2)=(" << degreeOfBody1 << "," << orderOfBody1 << ","
                              << degreeOfBody2 << "," << orderOfBody2 << ")"
                              << " signedOrders=(" << signedOrderOfBody1 << "," << signedOrderOfBody2 << ")"
                              << " j_case=" << j
                              << " multiplier=" << multiplier
                              << " C1=" << body1CosineCoefficient
                              << " S1=" << body1SineCoefficient
                              << " JCeff=" << effectiveAngularMomentumCosineCoefficients.transpose( )
                              << " JSeff=" << effectiveAngularMomentumSineCoefficients.transpose( )
                              << " term_contribution="
                              << ( equatorialRadiusRatioPower * legendrePolynomial *
                                   ( effectiveAngularMomentumCosineCoefficients * cosineOfMultipleLongitude +
                                     effectiveAngularMomentumSineCoefficients * sineOfMultipleLongitude ) ).transpose( )
                              << std::endl;
                }
            }
        }
        body2TorqueInBodyFixedFrameOfBody1 *= -preMultiplier;

        if( isFullTwoBodyTorqueDebugEnabled( ) )
        {
            const double finiteDifferenceStep = 1.0E-8;
            const std::array< Eigen::Vector3d, 3 > cartesianAxes = {
                Eigen::Vector3d::UnitX( ), Eigen::Vector3d::UnitY( ), Eigen::Vector3d::UnitZ( ) };
            Eigen::Vector3d finiteDifferenceBody2Torque = Eigen::Vector3d::Zero( );
            const Eigen::Quaterniond referenceBody2ToBody1Rotation =
                    accelerationBetweenBodies_->getCurrentRotationFromBody2ToBody1( );
            std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > potentialCache =
                    std::make_shared< basic_mathematics::SphericalHarmonicsCache >(
                        cosineCoefficientsOfBody1.rows( ) + cosineCoefficientsOfBody2.rows( ),
                        cosineCoefficientsOfBody1.cols( ) + cosineCoefficientsOfBody2.cols( ) );

            for( int component = 0; component < 3; component++ )
            {
                const Eigen::Quaterniond positiveRotation =
                        Eigen::Quaterniond( Eigen::AngleAxisd( finiteDifferenceStep, cartesianAxes.at( component ) ) ) *
                        referenceBody2ToBody1Rotation;
                effectiveMutualPotentialField->computeCurrentEffectiveCoefficients( positiveRotation );
                const double positivePotential =
                        effectiveMutualPotentialField->getGravitationalPotential(
                            bodyFixedRelativePosition, potentialCache );

                const Eigen::Quaterniond negativeRotation =
                        Eigen::Quaterniond( Eigen::AngleAxisd( -finiteDifferenceStep, cartesianAxes.at( component ) ) ) *
                        referenceBody2ToBody1Rotation;
                effectiveMutualPotentialField->computeCurrentEffectiveCoefficients( negativeRotation );
                const double negativePotential =
                        effectiveMutualPotentialField->getGravitationalPotential(
                            bodyFixedRelativePosition, potentialCache );

                finiteDifferenceBody2Torque( component ) =
                        -( positivePotential - negativePotential ) / ( 2.0 * finiteDifferenceStep );
            }

            effectiveMutualPotentialField->computeCurrentEffectiveCoefficients( referenceBody2ToBody1Rotation );

            const Eigen::Vector3d body2TorqueDifference =
                    body2TorqueInBodyFixedFrameOfBody1 - finiteDifferenceBody2Torque;
            std::cout << "[FTB-DBG][STEP body2_torque_fd]"
                      << " analytic_body2_specific=" << body2TorqueInBodyFixedFrameOfBody1.transpose( )
                      << " fd_body2_specific=" << finiteDifferenceBody2Torque.transpose( )
                      << " diff=" << body2TorqueDifference.transpose( )
                      << " status=" << debugStatusFromVectorDifference( body2TorqueDifference, 1.0E-11 )
                      << std::endl;
        }

        const Eigen::Vector3d totalTorqueInBodyFixedFrameOfBody1 =
                bodyFixedRelativePosition.cross( accelerationBetweenBodies_->getMutualPotentialGradient( ) );
        const Eigen::Vector3d body1TorqueInBodyFixedFrameOfBody1 =
                totalTorqueInBodyFixedFrameOfBody1 - body2TorqueInBodyFixedFrameOfBody1;

        if( acceleratedBodyIsBody1_ )
        {
            currentTorque_ = -body1TorqueInBodyFixedFrameOfBody1;
        }
        else
        {
            currentTorque_ = -( accelerationBetweenBodies_->getCurrentRotationFromBody2ToBody1( ).inverse( ) *
                    body2TorqueInBodyFixedFrameOfBody1 );
        }

        if( isFullTwoBodyTorqueDebugEnabled( ) )
        {
            std::cout << "[FTB-DBG][STEP torque_assembly]"
                      << " total_specific=" << totalTorqueInBodyFixedFrameOfBody1.transpose( )
                      << " body2_specific=" << body2TorqueInBodyFixedFrameOfBody1.transpose( )
                      << " body1_specific=" << body1TorqueInBodyFixedFrameOfBody1.transpose( )
                      << " exported_specific=" << currentTorque_.transpose( )
                      << std::endl;
        }

        if( isFullTwoBodyTorqueDebugEnabled( ) &&
            isBody1C20Only( cosineCoefficientsOfBody1, sineCoefficientsOfBody1 ) )
        {
            double activeCoefficientValue = 0.0;
            const DegreeTwoCoefficientCase activeCoefficientCase =
                    detectSingleActiveDegreeTwoCoefficientCase(
                        transformedCosineCoefficientsOfBody2,
                        transformedSineCoefficientsOfBody2,
                        activeCoefficientValue );

            if( activeCoefficientCase != DegreeTwoCoefficientCase::none &&
                activeCoefficientCase != DegreeTwoCoefficientCase::multiple )
            {
                const Eigen::Vector3d expectedSpecificTorque =
                        computeExpectedSpecificC20DegreeTwoTorqueFromDocument(
                            bodyFixedRelativePosition,
                            accelerationBetweenBodies_->getCurrentGravitationalParameter( ),
                            accelerationBetweenBodies_->getEquatorialRadiusOfBody1( ),
                            accelerationBetweenBodies_->getEquatorialRadiusOfBody2( ),
                            cosineCoefficientsOfBody1( 2, 0 ),
                            activeCoefficientCase,
                            activeCoefficientValue );
                const Eigen::Vector3d differenceFromExpected = currentTorque_ - expectedSpecificTorque;
                const Eigen::Vector3d totalDifferenceFromExpected =
                        totalTorqueInBodyFixedFrameOfBody1 - expectedSpecificTorque;
                std::cout << "[FTB-DBG][STEP torque_expected]"
                          << " coeff_case=" << coefficientCaseToString( activeCoefficientCase )
                          << " coeff_value=" << activeCoefficientValue
                          << " expected_specific=" << expectedSpecificTorque.transpose( )
                          << " exported_specific=" << currentTorque_.transpose( )
                          << " exported_minus_expected=" << differenceFromExpected.transpose( )
                          << " exported_status="
                          << debugStatusFromVectorDifference( differenceFromExpected, 1.0E-15 )
                          << " total_specific=" << totalTorqueInBodyFixedFrameOfBody1.transpose( )
                          << " total_minus_expected=" << totalDifferenceFromExpected.transpose( )
                          << " total_status="
                          << debugStatusFromVectorDifference( totalDifferenceFromExpected, 1.0E-15 )
                          << std::endl;
            }
            else
            {
                std::cout << "[FTB-DBG][STEP torque_expected] skipped_case="
                          << coefficientCaseToString( activeCoefficientCase ) << std::endl;
            }
        }

        currentTime_ = currentTime;
    }
}

}

}
