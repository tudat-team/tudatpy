/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/rotational_dynamics_partials/fullTwoBodySphericalHarmonicGravitationalTorquePartial.h"

#include <complex>

#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicCosineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicSineCoefficients.h"
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/basic/sphericalHarmonics.h"

namespace tudat
{

namespace acceleration_partials
{

namespace detail
{

std::array< Eigen::Matrix3d, 4 > getDerivativeOfBodyFixedToInertialRotationMatrixWrtQuaternionForFullTwoBodyTorque(
        const Eigen::Quaterniond& rotationFromInertialToBodyFixedFrame )
{
    std::vector< Eigen::Matrix3d > derivativeList( 4, Eigen::Matrix3d::Zero( ) );
    linear_algebra::computePartialDerivativeOfRotationMatrixWrtQuaternion(
            linear_algebra::convertQuaternionToVectorFormat( rotationFromInertialToBodyFixedFrame.inverse( ) ),
            derivativeList );

    std::array< Eigen::Matrix3d, 4 > derivativeArray;
    for( int i = 0; i < 4; i++ )
    {
        derivativeArray.at( i ) = derivativeList.at( i );
    }
    return derivativeArray;
}

Eigen::Matrix4d getLeftQuaternionMultiplicationMatrix( const Eigen::Vector4d& quaternion )
{
    Eigen::Matrix4d multiplicationMatrix = Eigen::Matrix4d::Zero( );
    multiplicationMatrix <<
            quaternion( 0 ), -quaternion( 1 ), -quaternion( 2 ), -quaternion( 3 ),
            quaternion( 1 ),  quaternion( 0 ), -quaternion( 3 ),  quaternion( 2 ),
            quaternion( 2 ),  quaternion( 3 ),  quaternion( 0 ), -quaternion( 1 ),
            quaternion( 3 ), -quaternion( 2 ),  quaternion( 1 ),  quaternion( 0 );
    return multiplicationMatrix;
}

Eigen::Matrix4d getRightQuaternionMultiplicationMatrix( const Eigen::Vector4d& quaternion )
{
    Eigen::Matrix4d multiplicationMatrix = Eigen::Matrix4d::Zero( );
    multiplicationMatrix <<
            quaternion( 0 ), -quaternion( 1 ), -quaternion( 2 ), -quaternion( 3 ),
            quaternion( 1 ),  quaternion( 0 ),  quaternion( 3 ), -quaternion( 2 ),
            quaternion( 2 ), -quaternion( 3 ),  quaternion( 0 ),  quaternion( 1 ),
            quaternion( 3 ),  quaternion( 2 ), -quaternion( 1 ),  quaternion( 0 );
    return multiplicationMatrix;
}

std::array< Eigen::MatrixXcd, 4 > computeDerivativeOfDegreeOneWignerDMatrixWrtRelativeQuaternion(
        const Eigen::Vector4d& relativeQuaternionVector )
{
    const std::complex< double > imaginaryUnit( 0.0, 1.0 );
    const std::complex< double > cayleyKleinA =
            std::complex< double >( relativeQuaternionVector( 0 ), -relativeQuaternionVector( 3 ) );
    const std::complex< double > cayleyKleinB =
            std::complex< double >( relativeQuaternionVector( 1 ), -relativeQuaternionVector( 2 ) );
    const std::complex< double > cayleyKleinAConjugate = std::conj( cayleyKleinA );
    const std::complex< double > cayleyKleinBConjugate = std::conj( cayleyKleinB );

    std::array< std::complex< double >, 4 > partialOfCayleyKleinA;
    std::array< std::complex< double >, 4 > partialOfCayleyKleinB;
    std::array< std::complex< double >, 4 > partialOfCayleyKleinAConjugate;
    std::array< std::complex< double >, 4 > partialOfCayleyKleinBConjugate;

    partialOfCayleyKleinA.at( 0 ) = 1.0;
    partialOfCayleyKleinA.at( 1 ) = 0.0;
    partialOfCayleyKleinA.at( 2 ) = 0.0;
    partialOfCayleyKleinA.at( 3 ) = -imaginaryUnit;

    partialOfCayleyKleinB.at( 0 ) = 0.0;
    partialOfCayleyKleinB.at( 1 ) = 1.0;
    partialOfCayleyKleinB.at( 2 ) = -imaginaryUnit;
    partialOfCayleyKleinB.at( 3 ) = 0.0;

    for( int i = 0; i < 4; i++ )
    {
        partialOfCayleyKleinAConjugate.at( i ) = std::conj( partialOfCayleyKleinA.at( i ) );
        partialOfCayleyKleinBConjugate.at( i ) = std::conj( partialOfCayleyKleinB.at( i ) );
    }

    std::array< Eigen::MatrixXcd, 4 > derivatives;
    for( int i = 0; i < 4; i++ )
    {
        derivatives.at( i ) = Eigen::MatrixXcd::Zero( 3, 3 );

        derivatives.at( i )( 2, 2 ) = 2.0 * cayleyKleinA * partialOfCayleyKleinA.at( i );
        derivatives.at( i )( 2, 1 ) =
                -std::sqrt( 2.0 ) *
                ( partialOfCayleyKleinA.at( i ) * cayleyKleinBConjugate +
                  cayleyKleinA * partialOfCayleyKleinBConjugate.at( i ) );
        derivatives.at( i )( 2, 0 ) = 2.0 * cayleyKleinBConjugate * partialOfCayleyKleinBConjugate.at( i );

        derivatives.at( i )( 1, 2 ) =
                std::sqrt( 2.0 ) *
                ( partialOfCayleyKleinA.at( i ) * cayleyKleinB +
                  cayleyKleinA * partialOfCayleyKleinB.at( i ) );
        derivatives.at( i )( 1, 1 ) =
                ( i == 0 ? 2.0 * relativeQuaternionVector( 0 ) :
                  ( i == 1 ? -2.0 * relativeQuaternionVector( 1 ) :
                    ( i == 2 ? -2.0 * relativeQuaternionVector( 2 ) :
                      2.0 * relativeQuaternionVector( 3 ) ) ) );
        derivatives.at( i )( 1, 0 ) =
                -std::sqrt( 2.0 ) *
                ( partialOfCayleyKleinAConjugate.at( i ) * cayleyKleinBConjugate +
                  cayleyKleinAConjugate * partialOfCayleyKleinBConjugate.at( i ) );

        derivatives.at( i )( 0, 2 ) = 2.0 * cayleyKleinB * partialOfCayleyKleinB.at( i );
        derivatives.at( i )( 0, 1 ) =
                std::sqrt( 2.0 ) *
                ( partialOfCayleyKleinAConjugate.at( i ) * cayleyKleinB +
                  cayleyKleinAConjugate * partialOfCayleyKleinB.at( i ) );
        derivatives.at( i )( 0, 0 ) = 2.0 * cayleyKleinAConjugate * partialOfCayleyKleinAConjugate.at( i );
    }
    return derivatives;
}

double getWignerRecursionCoefficientMinusOne( const int degree, const int rowIndex, const int columnIndex )
{
    const int orderM = rowIndex - degree;
    const int orderK = columnIndex - degree;
    return std::sqrt(
            static_cast< double >( ( degree + orderK ) * ( degree + orderK - 1 ) ) /
            static_cast< double >( ( degree + orderM ) * ( degree + orderM - 1 ) ) );
}

double getWignerRecursionCoefficientZero( const int degree, const int rowIndex, const int columnIndex )
{
    const int orderM = rowIndex - degree;
    const int orderK = columnIndex - degree;
    return std::sqrt(
            static_cast< double >( 2 * ( degree + orderK ) * ( degree - orderK ) ) /
            static_cast< double >( ( degree + orderM ) * ( degree + orderM - 1 ) ) );
}

double getWignerRecursionCoefficientOne( const int degree, const int rowIndex, const int columnIndex )
{
    const int orderM = rowIndex - degree;
    const int orderK = columnIndex - degree;
    return std::sqrt(
            static_cast< double >( ( degree - orderK ) * ( degree - orderK - 1 ) ) /
            static_cast< double >( ( degree + orderM ) * ( degree + orderM - 1 ) ) );
}

std::array< std::vector< Eigen::MatrixXcd >, 4 > computeDerivativeOfWignerDMatricesWrtRelativeQuaternion(
        const Eigen::Quaterniond& relativeRotationFromBody2ToBody1,
        const std::shared_ptr< basic_mathematics::WignerDMatricesCache >& wignerCache )
{
    const std::vector< Eigen::MatrixXcd >& wignerMatrices = wignerCache->getWignerDMatrices( );
    const int maximumDegree = static_cast< int >( wignerMatrices.size( ) ) - 1;

    std::array< std::vector< Eigen::MatrixXcd >, 4 > derivatives;
    for( int i = 0; i < 4; i++ )
    {
        derivatives.at( i ).resize( maximumDegree + 1 );
        derivatives.at( i ).at( 0 ) = Eigen::MatrixXcd::Zero( 1, 1 );
    }

    if( maximumDegree == 0 )
    {
        return derivatives;
    }

    const Eigen::Vector4d relativeQuaternionVector =
            linear_algebra::convertQuaternionToVectorFormat( relativeRotationFromBody2ToBody1 );
    const std::array< Eigen::MatrixXcd, 4 > degreeOneDerivatives =
            computeDerivativeOfDegreeOneWignerDMatrixWrtRelativeQuaternion( relativeQuaternionVector );

    for( int i = 0; i < 4; i++ )
    {
        derivatives.at( i ).at( 1 ) = degreeOneDerivatives.at( i );
    }

    if( maximumDegree == 1 )
    {
        return derivatives;
    }

    const Eigen::MatrixXcd& degreeOneWignerMatrix = wignerMatrices.at( 1 );

    for( int degree = 2; degree <= maximumDegree; degree++ )
    {
        for( int i = 0; i < 4; i++ )
        {
            derivatives.at( i ).at( degree ) = Eigen::MatrixXcd::Zero( 2 * degree + 1, 2 * degree + 1 );
        }

        for( int rowIndex = degree; rowIndex <= 2 * degree; rowIndex++ )
        {
            for( int columnIndex = 0; columnIndex <= 2 * degree; columnIndex++ )
            {
                if( rowIndex - 2 < 0 )
                {
                    continue;
                }

                for( int i = 0; i < 4; i++ )
                {
                    std::complex< double > derivativeEntry = std::complex< double >( 0.0, 0.0 );

                    if( columnIndex > 1 )
                    {
                        const double coefficient = getWignerRecursionCoefficientMinusOne( degree, rowIndex, columnIndex );
                        derivativeEntry += coefficient *
                                ( degreeOneDerivatives.at( i )( 2, 2 ) *
                                          wignerMatrices.at( degree - 1 )( rowIndex - 2, columnIndex - 2 ) +
                                  degreeOneWignerMatrix( 2, 2 ) *
                                          derivatives.at( i ).at( degree - 1 )( rowIndex - 2, columnIndex - 2 ) );
                    }
                    if( columnIndex > 0 && columnIndex <= 2 * degree - 1 )
                    {
                        const double coefficient = getWignerRecursionCoefficientZero( degree, rowIndex, columnIndex );
                        derivativeEntry += coefficient *
                                ( degreeOneDerivatives.at( i )( 2, 1 ) *
                                          wignerMatrices.at( degree - 1 )( rowIndex - 2, columnIndex - 1 ) +
                                  degreeOneWignerMatrix( 2, 1 ) *
                                          derivatives.at( i ).at( degree - 1 )( rowIndex - 2, columnIndex - 1 ) );
                    }
                    if( columnIndex < 2 * degree - 1 )
                    {
                        const double coefficient = getWignerRecursionCoefficientOne( degree, rowIndex, columnIndex );
                        derivativeEntry += coefficient *
                                ( degreeOneDerivatives.at( i )( 2, 0 ) *
                                          wignerMatrices.at( degree - 1 )( rowIndex - 2, columnIndex ) +
                                  degreeOneWignerMatrix( 2, 0 ) *
                                          derivatives.at( i ).at( degree - 1 )( rowIndex - 2, columnIndex ) );
                    }

                    derivatives.at( i ).at( degree )( rowIndex, columnIndex ) = derivativeEntry;
                }
            }
        }

        for( int rowIndex = 0; rowIndex < degree; rowIndex++ )
        {
            const int orderM = rowIndex - degree;
            for( int columnIndex = 0; columnIndex <= 2 * degree; columnIndex++ )
            {
                const int orderK = columnIndex - degree;
                const double signMultiplier = ( ( ( orderM - orderK ) % 2 ) == 0 ) ? 1.0 : -1.0;
                for( int i = 0; i < 4; i++ )
                {
                    derivatives.at( i ).at( degree )( rowIndex, columnIndex ) =
                            signMultiplier *
                            std::conj( derivatives.at( i ).at( degree )( -orderM + degree, -orderK + degree ) );
                }
            }
        }
    }

    return derivatives;
}

void transformCoefficientsWithWignerMatrices(
        const Eigen::MatrixXd& originalCosineCoefficients,
        const Eigen::MatrixXd& originalSineCoefficients,
        const std::vector< Eigen::MatrixXcd >& wignerMatrices,
        Eigen::MatrixXd& transformedCosineCoefficients,
        Eigen::MatrixXd& transformedSineCoefficients,
        const bool areCoefficientsNormalized )
{
    transformedCosineCoefficients.setZero( originalCosineCoefficients.rows( ), originalCosineCoefficients.cols( ) );
    transformedSineCoefficients.setZero( originalSineCoefficients.rows( ), originalSineCoefficients.cols( ) );

    for( int degree = 0; degree < originalCosineCoefficients.rows( ); degree++ )
    {
        for( int orderM = 0; ( orderM <= degree && orderM < originalCosineCoefficients.cols( ) ); orderM++ )
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

            const std::complex< double > orderZeroD = wignerMatrices.at( degree )( orderM + degree, degree );
            transformedCosineCoefficients( degree, orderM ) +=
                    orderMMultiplier * orderZeroD.real( ) * originalCosineCoefficients( degree, 0 );
            transformedSineCoefficients( degree, orderM ) +=
                    orderMMultiplier * orderZeroD.imag( ) * originalCosineCoefficients( degree, 0 );

            for( int orderK = 1; orderK <= degree; orderK++ )
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

                const double signMultiplier = ( ( ( orderK % 2 ) == 0 ) ? 1.0 : -1.0 );
                const std::complex< double > positiveOrderD =
                        wignerMatrices.at( degree )( orderM + degree, orderK + degree );
                const std::complex< double > negativeOrderD =
                        wignerMatrices.at( degree )( orderM + degree, -orderK + degree );

                transformedCosineCoefficients( degree, orderM ) +=
                        0.5 * currentMultiplier *
                        ( ( signMultiplier * positiveOrderD.real( ) + negativeOrderD.real( ) ) *
                                  originalCosineCoefficients( degree, orderK ) +
                          ( signMultiplier * positiveOrderD.imag( ) - negativeOrderD.imag( ) ) *
                                  originalSineCoefficients( degree, orderK ) );

                transformedSineCoefficients( degree, orderM ) +=
                        0.5 * currentMultiplier *
                        ( ( signMultiplier * positiveOrderD.imag( ) + negativeOrderD.imag( ) ) *
                                  originalCosineCoefficients( degree, orderK ) +
                          ( -signMultiplier * positiveOrderD.real( ) + negativeOrderD.real( ) ) *
                                  originalSineCoefficients( degree, orderK ) );
            }

            double cosineScaling = ( ( ( orderM % 2 ) == 0 ) ? 1.0 : -1.0 );
            double sineScaling = ( ( ( ( orderM + 1 ) % 2 ) == 0 ) ? 1.0 : -1.0 );
            if( orderM > 0 )
            {
                cosineScaling *= 2.0;
                sineScaling *= 2.0;
            }

            transformedCosineCoefficients( degree, orderM ) *= cosineScaling;
            transformedSineCoefficients( degree, orderM ) *= sineScaling;
        }
    }
}

Eigen::Vector3cd computeAngularMomentumOperatorOnWignerCoefficientFromMatrices(
        const std::vector< Eigen::MatrixXcd >& wignerMatrices,
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
        return wignerMatrices.at( degree )( requestedOrderM + degree, requestedOrderK + degree );
    };

    const double plusScaling = std::sqrt(
            std::max( 0.0, static_cast< double >( degree * ( degree + 1 ) - orderM * ( orderM - 1 ) ) ) ) / 2.0;
    const double minusScaling = std::sqrt(
            std::max( 0.0, static_cast< double >( degree * ( degree + 1 ) - orderM * ( orderM + 1 ) ) ) ) / 2.0;

    const std::complex< double > angularMomentumPlus =
            imaginaryUnit * plusScaling * getWignerCoefficient( orderM - 1, orderK );
    const std::complex< double > angularMomentumMinus =
            imaginaryUnit * ( -minusScaling ) * getWignerCoefficient( orderM + 1, orderK );
    const std::complex< double > angularMomentumZero =
            imaginaryUnit * static_cast< double >( -orderM ) * getWignerCoefficient( orderM, orderK );

    Eigen::Vector3cd angularMomentumInCartesianBasis;
    angularMomentumInCartesianBasis( 0 ) = ( angularMomentumMinus - angularMomentumPlus ) * inverseSquareRootTwo;
    angularMomentumInCartesianBasis( 1 ) =
            imaginaryUnit * ( angularMomentumMinus + angularMomentumPlus ) * inverseSquareRootTwo;
    angularMomentumInCartesianBasis( 2 ) = angularMomentumZero;
    return angularMomentumInCartesianBasis;
}

void computeTransformedAngularMomentumCoefficientsWithWignerMatrices(
        const Eigen::MatrixXd& cosineCoefficientsBody2,
        const Eigen::MatrixXd& sineCoefficientsBody2,
        const std::vector< Eigen::MatrixXcd >& wignerMatrices,
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

            const Eigen::Vector3cd orderZeroAngularMomentumD = computeAngularMomentumOperatorOnWignerCoefficientFromMatrices(
                    wignerMatrices, degree, orderM, 0 );
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

                const double signMultiplier = ( ( ( orderK % 2 ) == 0 ) ? 1.0 : -1.0 );
                const Eigen::Vector3cd positiveOrderAngularMomentumD =
                        computeAngularMomentumOperatorOnWignerCoefficientFromMatrices(
                                wignerMatrices, degree, orderM, orderK );
                const Eigen::Vector3cd negativeOrderAngularMomentumD =
                        computeAngularMomentumOperatorOnWignerCoefficientFromMatrices(
                                wignerMatrices, degree, orderM, -orderK );

                for( int i = 0; i < 3; i++ )
                {
                    transformedCosineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) +=
                            0.5 * currentMultiplier *
                            ( ( signMultiplier * positiveOrderAngularMomentumD( i ).real( ) +
                                negativeOrderAngularMomentumD( i ).real( ) ) *
                                      cosineCoefficientsBody2( degree, orderK ) +
                              ( signMultiplier * positiveOrderAngularMomentumD( i ).imag( ) -
                                negativeOrderAngularMomentumD( i ).imag( ) ) *
                                      sineCoefficientsBody2( degree, orderK ) );

                    transformedSineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) +=
                            0.5 * currentMultiplier *
                            ( ( signMultiplier * positiveOrderAngularMomentumD( i ).imag( ) +
                                negativeOrderAngularMomentumD( i ).imag( ) ) *
                                      cosineCoefficientsBody2( degree, orderK ) +
                              ( -signMultiplier * positiveOrderAngularMomentumD( i ).real( ) +
                                negativeOrderAngularMomentumD( i ).real( ) ) *
                                      sineCoefficientsBody2( degree, orderK ) );
                }
            }

            double cosineScaling = ( ( ( orderM % 2 ) == 0 ) ? 1.0 : -1.0 );
            double sineScaling = ( ( ( ( orderM + 1 ) % 2 ) == 0 ) ? 1.0 : -1.0 );
            if( orderM > 0 )
            {
                cosineScaling *= 2.0;
                sineScaling *= 2.0;
            }

            for( int i = 0; i < 3; i++ )
            {
                transformedCosineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) *= cosineScaling;
                transformedSineCoefficientsBody2AngularMomentum.at( i )( degree, orderM ) *= sineScaling;
            }
        }
    }
}

bool getSignedOrdersForVariant(
        const int variant,
        const unsigned int orderOfBody1,
        const unsigned int orderOfBody2,
        int& signedOrderOfBody1,
        int& signedOrderOfBody2 )
{
    bool computeTerm = true;

    switch( variant )
    {
    case 0:
        signedOrderOfBody1 = static_cast< int >( orderOfBody1 );
        signedOrderOfBody2 = static_cast< int >( orderOfBody2 );
        break;
    case 1:
        signedOrderOfBody1 = -static_cast< int >( orderOfBody1 );
        signedOrderOfBody2 = static_cast< int >( orderOfBody2 );
        computeTerm = ( orderOfBody1 != 0 );
        break;
    case 2:
        signedOrderOfBody1 = static_cast< int >( orderOfBody1 );
        signedOrderOfBody2 = -static_cast< int >( orderOfBody2 );
        computeTerm = ( orderOfBody2 != 0 );
        break;
    case 3:
        signedOrderOfBody1 = -static_cast< int >( orderOfBody1 );
        signedOrderOfBody2 = -static_cast< int >( orderOfBody2 );
        computeTerm = ( orderOfBody1 != 0 && orderOfBody2 != 0 );
        break;
    default:
        computeTerm = false;
        break;
    }

    return computeTerm;
}

Eigen::Vector2d computeCurrentScalarBasisFunctionValues(
        const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache >& sphericalHarmonicsCache,
        const double preMultiplier,
        const double equatorialRadiusRatioPower,
        const int totalDegree,
        const int totalOrder )
{
    Eigen::Vector2d scalarBasisFunctionValues = Eigen::Vector2d::Zero( );
    const double legendrePolynomial =
            sphericalHarmonicsCache->getLegendreCache( ).getLegendrePolynomial( totalDegree, totalOrder );
    scalarBasisFunctionValues( 0 ) =
            preMultiplier * equatorialRadiusRatioPower * legendrePolynomial *
            sphericalHarmonicsCache->getCosineOfMultipleLongitude( totalOrder );
    scalarBasisFunctionValues( 1 ) =
            preMultiplier * equatorialRadiusRatioPower * legendrePolynomial *
            sphericalHarmonicsCache->getSineOfMultipleLongitude( totalOrder );
    return scalarBasisFunctionValues;
}

Eigen::Matrix< double, 3, 2 > computeCurrentBodyFixedBasisFunctionGradients(
        const Eigen::Vector3d& bodyFixedRelativePosition,
        const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache >& sphericalHarmonicsCache,
        const double cosineOfLatitude,
        const double preMultiplier,
        const double equatorialRadiusRatioPower,
        const int totalDegree,
        const int totalOrder )
{
    const double legendrePolynomial =
            sphericalHarmonicsCache->getLegendreCache( ).getLegendrePolynomial( totalDegree, totalOrder );
    const double legendrePolynomialDerivative =
            sphericalHarmonicsCache->getLegendreCache( ).getLegendrePolynomialDerivative( totalDegree, totalOrder );

    Eigen::Matrix< double, 3, 2 > bodyFixedBasisFunctionGradients = Eigen::Matrix< double, 3, 2 >::Zero( );
    bodyFixedBasisFunctionGradients.col( 0 ) =
            coordinate_conversions::convertSphericalToCartesianGradient(
                    basic_mathematics::computePotentialGradient(
                            bodyFixedRelativePosition.norm( ),
                            equatorialRadiusRatioPower,
                            sphericalHarmonicsCache->getCosineOfMultipleLongitude( totalOrder ),
                            sphericalHarmonicsCache->getSineOfMultipleLongitude( totalOrder ),
                            cosineOfLatitude,
                            preMultiplier,
                            totalDegree,
                            totalOrder,
                            1.0,
                            0.0,
                            legendrePolynomial,
                            legendrePolynomialDerivative ),
                    bodyFixedRelativePosition );
    bodyFixedBasisFunctionGradients.col( 1 ) =
            coordinate_conversions::convertSphericalToCartesianGradient(
                    basic_mathematics::computePotentialGradient(
                            bodyFixedRelativePosition.norm( ),
                            equatorialRadiusRatioPower,
                            sphericalHarmonicsCache->getCosineOfMultipleLongitude( totalOrder ),
                            sphericalHarmonicsCache->getSineOfMultipleLongitude( totalOrder ),
                            cosineOfLatitude,
                            preMultiplier,
                            totalDegree,
                            totalOrder,
                            0.0,
                            1.0,
                            legendrePolynomial,
                            legendrePolynomialDerivative ),
                    bodyFixedRelativePosition );
    return bodyFixedBasisFunctionGradients;
}

}  // namespace detail

FullTwoBodySphericalHarmonicGravitationalTorquePartial::FullTwoBodySphericalHarmonicGravitationalTorquePartial(
        const std::shared_ptr< gravitation::FullTwoBodySphericalHarmonicTorque > torqueModel,
        const std::shared_ptr< FullTwoBodySphericalHarmonicsGravityPartial > accelerationPartial,
        const std::string& acceleratedBody,
        const std::string& acceleratingBody ):
    TorquePartial( acceleratedBody, acceleratingBody, basic_astrodynamics::full_two_body_spherical_harmonic_gravitational_torque ),
    torqueModel_( torqueModel ),
    accelerationPartial_( accelerationPartial ),
    accelerationModel_( torqueModel == nullptr ? nullptr : torqueModel->getAccelerationBetweenBodies( ) ),
    effectiveMutualPotentialField_(
            accelerationModel_ == nullptr ? nullptr : accelerationModel_->getEffectiveMutualPotentialField( ) ),
    coefficientCombinationsToUse_(
            torqueModel == nullptr ? std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >( ) :
                                     torqueModel->getCoefficientCombinationsToUse( ) ),
    currentPartialWrtQuaternionOfBodyUndergoingTorque_( Eigen::Matrix< double, 3, 4 >::Zero( ) ),
    currentPartialWrtQuaternionOfBodyExertingTorque_( Eigen::Matrix< double, 3, 4 >::Zero( ) ),
    currentPartialWrtPositionOfBodyUndergoingTorque_( Eigen::Matrix3d::Zero( ) ),
    currentPartialWrtPositionOfBodyExertingTorque_( Eigen::Matrix3d::Zero( ) ),
    currentRotationFromInertialToBodyFixedFrameOfBodyUndergoingTorque_( Eigen::Matrix3d::Identity( ) ),
    currentRotationFromInertialToBodyFixedFrameOfBodyExertingTorque_( Eigen::Matrix3d::Identity( ) ),
    currentRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorque_( Eigen::Matrix3d::Identity( ) ),
    currentRelativePositionInInertialFrame_( Eigen::Vector3d::Zero( ) ),
    currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_( Eigen::Vector3d::Zero( ) ),
    currentMutualPotentialGradientInBodyFixedFrameOfBodyUndergoingTorque_( Eigen::Vector3d::Zero( ) ),
    currentBody2SpinTorquePartialWrtBodyFixedRelativePosition_( Eigen::Matrix3d::Zero( ) ),
    currentDistance_( TUDAT_NAN ),
    currentCosineOfLatitude_( TUDAT_NAN ),
    currentPreMultiplier_( TUDAT_NAN )
{
    if( torqueModel_ == nullptr )
    {
        throw std::runtime_error(
                "Error when creating FullTwoBodySphericalHarmonicGravitationalTorquePartial, torque model is nullptr." );
    }
    if( accelerationPartial_ == nullptr || accelerationModel_ == nullptr || effectiveMutualPotentialField_ == nullptr )
    {
        throw std::runtime_error(
                "Error when creating FullTwoBodySphericalHarmonicGravitationalTorquePartial, analytical dependencies are nullptr." );
    }
    if( !torqueModel_->getAcceleratedBodyIsBody1( ) )
    {
        throw std::runtime_error(
                "Error when creating FullTwoBodySphericalHarmonicGravitationalTorquePartial, only the body-1 torque path is "
                "implemented analytically." );
    }

    const std::vector< int >& accelerationEffectiveIndices =
            accelerationPartial_->getEffectiveIndicesForCoefficientCombinations( );
    effectiveIndicesForCoefficientCombinations_ = accelerationEffectiveIndices;
}

std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
FullTwoBodySphericalHarmonicGravitationalTorquePartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > )
{
    return std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );
}

std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
FullTwoBodySphericalHarmonicGravitationalTorquePartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    using namespace estimatable_parameters;

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunction =
            std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );

    if( parameter->getParameterName( ).second.first == bodyUndergoingTorque_ )
    {
        switch( parameter->getParameterName( ).first )
        {
        case spherical_harmonics_cosine_coefficient_block:
        {
            std::shared_ptr< SphericalHarmonicsCosineCoefficients > coefficientsParameter =
                    std::dynamic_pointer_cast< SphericalHarmonicsCosineCoefficients >( parameter );
            partialFunction = std::make_pair(
                    std::bind(
                            &FullTwoBodySphericalHarmonicGravitationalTorquePartial::
                                    wrtCosineSphericalHarmonicCoefficientsOfBodyUndergoingTorque,
                            this,
                            std::placeholders::_1,
                            coefficientsParameter->getBlockIndices( ) ),
                    coefficientsParameter->getParameterSize( ) );
            break;
        }
        case spherical_harmonics_sine_coefficient_block:
        {
            std::shared_ptr< SphericalHarmonicsSineCoefficients > coefficientsParameter =
                    std::dynamic_pointer_cast< SphericalHarmonicsSineCoefficients >( parameter );
            partialFunction = std::make_pair(
                    std::bind(
                            &FullTwoBodySphericalHarmonicGravitationalTorquePartial::
                                    wrtSineSphericalHarmonicCoefficientsOfBodyUndergoingTorque,
                            this,
                            std::placeholders::_1,
                            coefficientsParameter->getBlockIndices( ) ),
                    coefficientsParameter->getParameterSize( ) );
            break;
        }
        default:
            break;
        }
    }
    else if( parameter->getParameterName( ).second.first == bodyExertingTorque_ )
    {
        switch( parameter->getParameterName( ).first )
        {
        case spherical_harmonics_cosine_coefficient_block:
        {
            std::shared_ptr< SphericalHarmonicsCosineCoefficients > coefficientsParameter =
                    std::dynamic_pointer_cast< SphericalHarmonicsCosineCoefficients >( parameter );
            partialFunction = std::make_pair(
                    std::bind(
                            &FullTwoBodySphericalHarmonicGravitationalTorquePartial::
                                    wrtCosineSphericalHarmonicCoefficientsOfBodyExertingTorque,
                            this,
                            std::placeholders::_1,
                            coefficientsParameter->getBlockIndices( ) ),
                    coefficientsParameter->getParameterSize( ) );
            break;
        }
        case spherical_harmonics_sine_coefficient_block:
        {
            std::shared_ptr< SphericalHarmonicsSineCoefficients > coefficientsParameter =
                    std::dynamic_pointer_cast< SphericalHarmonicsSineCoefficients >( parameter );
            partialFunction = std::make_pair(
                    std::bind(
                            &FullTwoBodySphericalHarmonicGravitationalTorquePartial::
                                    wrtSineSphericalHarmonicCoefficientsOfBodyExertingTorque,
                            this,
                            std::placeholders::_1,
                            coefficientsParameter->getBlockIndices( ) ),
                    coefficientsParameter->getParameterSize( ) );
            break;
        }
        default:
            break;
        }
    }

    return partialFunction;
}

void FullTwoBodySphericalHarmonicGravitationalTorquePartial::wrtOrientationOfAcceleratedBody(
        Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const bool addContribution,
        const int startRow,
        const int startColumn )
{
    if( addContribution )
    {
        partialMatrix.block( startRow, startColumn, 3, 4 ) += currentPartialWrtQuaternionOfBodyUndergoingTorque_;
    }
    else
    {
        partialMatrix.block( startRow, startColumn, 3, 4 ) -= currentPartialWrtQuaternionOfBodyUndergoingTorque_;
    }
}

void FullTwoBodySphericalHarmonicGravitationalTorquePartial::wrtOrientationOfAcceleratingBody(
        Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const bool addContribution,
        const int startRow,
        const int startColumn )
{
    if( addContribution )
    {
        partialMatrix.block( startRow, startColumn, 3, 4 ) += currentPartialWrtQuaternionOfBodyExertingTorque_;
    }
    else
    {
        partialMatrix.block( startRow, startColumn, 3, 4 ) -= currentPartialWrtQuaternionOfBodyExertingTorque_;
    }
}

bool FullTwoBodySphericalHarmonicGravitationalTorquePartial::isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
        const std::pair< std::string, std::string >& stateReferencePoint,
        const propagators::IntegratedStateType integratedStateType )
{
    return ( integratedStateType == propagators::translational_state &&
             ( stateReferencePoint.first == bodyUndergoingTorque_ || stateReferencePoint.first == bodyExertingTorque_ ) );
}

void FullTwoBodySphericalHarmonicGravitationalTorquePartial::wrtNonRotationalStateOfAdditionalBody(
        Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const std::pair< std::string, std::string >& stateReferencePoint,
        const propagators::IntegratedStateType integratedStateType )
{
    if( integratedStateType == propagators::translational_state )
    {
        if( stateReferencePoint.first == bodyUndergoingTorque_ )
        {
            partialMatrix.block( 0, 0, 3, 3 ) += currentPartialWrtPositionOfBodyUndergoingTorque_;
        }
        else if( stateReferencePoint.first == bodyExertingTorque_ )
        {
            partialMatrix.block( 0, 0, 3, 3 ) += currentPartialWrtPositionOfBodyExertingTorque_;
        }
    }
}

void FullTwoBodySphericalHarmonicGravitationalTorquePartial::addBody2SpinTorquePartialWrtBody1Coefficient(
        Eigen::Vector3d& partial,
        const int degree,
        const int order,
        const bool wrtCosineCoefficient ) const
{
    const Eigen::MatrixXd cosineCoefficientsOfBody1 = effectiveMutualPotentialField_->getCosineCoefficientsOfBody1( );
    const Eigen::MatrixXd sineCoefficientsOfBody1 = effectiveMutualPotentialField_->getSineCoefficientsOfBody1( );
    const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache =
            accelerationModel_->getSphericalHarmonicsCache( );

    for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
    {
        const int degreeOfBody1 = std::get< 0 >( coefficientCombinationsToUse_.at( i ) );
        const int orderOfBody1 = std::get< 1 >( coefficientCombinationsToUse_.at( i ) );
        const int degreeOfBody2 = std::get< 2 >( coefficientCombinationsToUse_.at( i ) );
        const int orderOfBody2 = std::get< 3 >( coefficientCombinationsToUse_.at( i ) );

        if( degreeOfBody1 != degree || orderOfBody1 != order )
        {
            continue;
        }

        const double equatorialRadiusRatioPower =
                currentRadius1Powers_.at( degreeOfBody1 ) * currentRadius2Powers_.at( degreeOfBody2 );
        const int totalDegree = degreeOfBody1 + degreeOfBody2;

        for( int variant = 0; variant < 4; variant++ )
        {
            int signedOrderOfBody1 = 0;
            int signedOrderOfBody2 = 0;
            if( !detail::getSignedOrdersForVariant(
                        variant, orderOfBody1, orderOfBody2, signedOrderOfBody1, signedOrderOfBody2 ) )
            {
                continue;
            }

            const int totalOrder = std::abs( signedOrderOfBody1 + signedOrderOfBody2 );
            const Eigen::Vector2d scalarBasisFunctionValues = detail::computeCurrentScalarBasisFunctionValues(
                    sphericalHarmonicsCache, currentPreMultiplier_, equatorialRadiusRatioPower, totalDegree, totalOrder );

            const double signOrderBody1 = ( signedOrderOfBody1 < 0 ) ? -1.0 : 1.0;
            const double signOrderBody2 = ( signedOrderOfBody2 < 0 ) ? -1.0 : 1.0;
            const double signTotalOrder = ( ( signedOrderOfBody1 + signedOrderOfBody2 ) < 0 ) ? -1.0 : 1.0;
            const double multiplier = effectiveMutualPotentialField_->getMultiplier(
                    degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );

            Eigen::Vector3d transformedAngularMomentumCosineCoefficientsBody2;
            Eigen::Vector3d transformedAngularMomentumSineCoefficientsBody2;
            for( int component = 0; component < 3; component++ )
            {
                transformedAngularMomentumCosineCoefficientsBody2( component ) =
                        currentTransformedCosineCoefficientsBody2AngularMomentum_.at( component )(
                                degreeOfBody2, std::abs( signedOrderOfBody2 ) );
                transformedAngularMomentumSineCoefficientsBody2( component ) =
                        currentTransformedSineCoefficientsBody2AngularMomentum_.at( component )(
                                degreeOfBody2, std::abs( signedOrderOfBody2 ) );
            }

            Eigen::Vector3d partialOfEffectiveAngularMomentumCosineCoefficients = Eigen::Vector3d::Zero( );
            Eigen::Vector3d partialOfEffectiveAngularMomentumSineCoefficients = Eigen::Vector3d::Zero( );
            if( wrtCosineCoefficient )
            {
                partialOfEffectiveAngularMomentumCosineCoefficients =
                        multiplier * transformedAngularMomentumCosineCoefficientsBody2;
                partialOfEffectiveAngularMomentumSineCoefficients =
                        signOrderBody2 * signTotalOrder * multiplier *
                        transformedAngularMomentumSineCoefficientsBody2;
            }
            else
            {
                partialOfEffectiveAngularMomentumCosineCoefficients =
                        -signOrderBody1 * signOrderBody2 * multiplier *
                        transformedAngularMomentumSineCoefficientsBody2;
                partialOfEffectiveAngularMomentumSineCoefficients =
                        signOrderBody1 * signTotalOrder * multiplier *
                        transformedAngularMomentumCosineCoefficientsBody2;
            }

            partial += partialOfEffectiveAngularMomentumCosineCoefficients * scalarBasisFunctionValues( 0 ) +
                    partialOfEffectiveAngularMomentumSineCoefficients * scalarBasisFunctionValues( 1 );
        }
    }
}

void FullTwoBodySphericalHarmonicGravitationalTorquePartial::addBody2SpinTorquePartialWrtBody2Coefficient(
        Eigen::Vector3d& partial,
        const int degree,
        const int order,
        const bool wrtCosineCoefficient ) const
{
    Eigen::MatrixXd basisCosineCoefficients = Eigen::MatrixXd::Zero(
            currentTransformedCosineCoefficientsBody2_.rows( ), currentTransformedCosineCoefficientsBody2_.cols( ) );
    Eigen::MatrixXd basisSineCoefficients = Eigen::MatrixXd::Zero(
            currentTransformedSineCoefficientsBody2_.rows( ), currentTransformedSineCoefficientsBody2_.cols( ) );
    if( wrtCosineCoefficient )
    {
        basisCosineCoefficients( degree, order ) = 1.0;
    }
    else
    {
        basisSineCoefficients( degree, order ) = 1.0;
    }

    std::array< Eigen::MatrixXd, 3 > transformedCosineCoefficientsBody2AngularMomentumPartials;
    std::array< Eigen::MatrixXd, 3 > transformedSineCoefficientsBody2AngularMomentumPartials;
    torqueModel_->computeTransformedAngularMomentumCoefficients(
            basisCosineCoefficients,
            basisSineCoefficients,
            effectiveMutualPotentialField_->getTransformationCache( )->getWignerDMatricesCache( ),
            accelerationModel_->getAreCoefficientsNormalized( ),
            transformedCosineCoefficientsBody2AngularMomentumPartials,
            transformedSineCoefficientsBody2AngularMomentumPartials );

    const Eigen::MatrixXd cosineCoefficientsOfBody1 = effectiveMutualPotentialField_->getCosineCoefficientsOfBody1( );
    const Eigen::MatrixXd sineCoefficientsOfBody1 = effectiveMutualPotentialField_->getSineCoefficientsOfBody1( );
    const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache =
            accelerationModel_->getSphericalHarmonicsCache( );

    for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
    {
        const int degreeOfBody1 = std::get< 0 >( coefficientCombinationsToUse_.at( i ) );
        const int orderOfBody1 = std::get< 1 >( coefficientCombinationsToUse_.at( i ) );
        const int degreeOfBody2 = std::get< 2 >( coefficientCombinationsToUse_.at( i ) );
        const int orderOfBody2 = std::get< 3 >( coefficientCombinationsToUse_.at( i ) );

        if( degreeOfBody2 != degree )
        {
            continue;
        }

        const double body1CosineCoefficient = cosineCoefficientsOfBody1( degreeOfBody1, orderOfBody1 );
        const double body1SineCoefficient = sineCoefficientsOfBody1( degreeOfBody1, orderOfBody1 );
        const double equatorialRadiusRatioPower =
                currentRadius1Powers_.at( degreeOfBody1 ) * currentRadius2Powers_.at( degreeOfBody2 );
        const int totalDegree = degreeOfBody1 + degreeOfBody2;

        for( int variant = 0; variant < 4; variant++ )
        {
            int signedOrderOfBody1 = 0;
            int signedOrderOfBody2 = 0;
            if( !detail::getSignedOrdersForVariant(
                        variant, orderOfBody1, orderOfBody2, signedOrderOfBody1, signedOrderOfBody2 ) )
            {
                continue;
            }

            const int totalOrder = std::abs( signedOrderOfBody1 + signedOrderOfBody2 );
            const Eigen::Vector2d scalarBasisFunctionValues = detail::computeCurrentScalarBasisFunctionValues(
                    sphericalHarmonicsCache, currentPreMultiplier_, equatorialRadiusRatioPower, totalDegree, totalOrder );

            const double signOrderBody1 = ( signedOrderOfBody1 < 0 ) ? -1.0 : 1.0;
            const double signOrderBody2 = ( signedOrderOfBody2 < 0 ) ? -1.0 : 1.0;
            const double signTotalOrder = ( ( signedOrderOfBody1 + signedOrderOfBody2 ) < 0 ) ? -1.0 : 1.0;
            const double multiplier = effectiveMutualPotentialField_->getMultiplier(
                    degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );

            Eigen::Vector3d partialOfTransformedAngularMomentumCosineCoefficientsBody2 = Eigen::Vector3d::Zero( );
            Eigen::Vector3d partialOfTransformedAngularMomentumSineCoefficientsBody2 = Eigen::Vector3d::Zero( );
            for( int component = 0; component < 3; component++ )
            {
                partialOfTransformedAngularMomentumCosineCoefficientsBody2( component ) =
                        transformedCosineCoefficientsBody2AngularMomentumPartials.at( component )(
                                degreeOfBody2, std::abs( signedOrderOfBody2 ) );
                partialOfTransformedAngularMomentumSineCoefficientsBody2( component ) =
                        transformedSineCoefficientsBody2AngularMomentumPartials.at( component )(
                                degreeOfBody2, std::abs( signedOrderOfBody2 ) );
            }

            const Eigen::Vector3d partialOfEffectiveAngularMomentumCosineCoefficients =
                    ( body1CosineCoefficient * partialOfTransformedAngularMomentumCosineCoefficientsBody2 -
                      signOrderBody1 * signOrderBody2 * body1SineCoefficient *
                              partialOfTransformedAngularMomentumSineCoefficientsBody2 ) *
                    multiplier;
            const Eigen::Vector3d partialOfEffectiveAngularMomentumSineCoefficients =
                    ( signOrderBody2 * body1CosineCoefficient *
                              partialOfTransformedAngularMomentumSineCoefficientsBody2 +
                      signOrderBody1 * body1SineCoefficient *
                              partialOfTransformedAngularMomentumCosineCoefficientsBody2 ) *
                    signTotalOrder * multiplier;

            partial += partialOfEffectiveAngularMomentumCosineCoefficients * scalarBasisFunctionValues( 0 ) +
                    partialOfEffectiveAngularMomentumSineCoefficients * scalarBasisFunctionValues( 1 );
        }
    }
}

void FullTwoBodySphericalHarmonicGravitationalTorquePartial::wrtCosineSphericalHarmonicCoefficientsOfBodyUndergoingTorque(
        Eigen::MatrixXd& partialMatrix,
        const std::vector< std::pair< int, int > >& blockIndices )
{
    partialMatrix = Eigen::MatrixXd::Zero( 3, blockIndices.size( ) );

    const Eigen::Matrix3d bodyFixedCrossProductMatrix =
            linear_algebra::getCrossProductMatrix( currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_ );
    const std::vector< Eigen::Matrix< double, 3, 2 > >& bodyFixedPartialsWrtEffectiveCoefficients =
            accelerationPartial_->getCurrentBodyFixedPartialsWrtEffectiveCoefficients( );

    for( unsigned int i = 0; i < blockIndices.size( ); i++ )
    {
        Eigen::Vector3d partial = Eigen::Vector3d::Zero( );
        addBody2SpinTorquePartialWrtBody1Coefficient(
                partial, blockIndices.at( i ).first, blockIndices.at( i ).second, true );

        for( unsigned int j = 0; j < coefficientCombinationsToUse_.size( ); j++ )
        {
            const int degreeOfBody1 = std::get< 0 >( coefficientCombinationsToUse_.at( j ) );
            const int orderOfBody1 = std::get< 1 >( coefficientCombinationsToUse_.at( j ) );
            const int degreeOfBody2 = std::get< 2 >( coefficientCombinationsToUse_.at( j ) );
            const int orderOfBody2 = std::get< 3 >( coefficientCombinationsToUse_.at( j ) );

            if( blockIndices.at( i ).first == degreeOfBody1 && blockIndices.at( i ).second == orderOfBody1 )
            {
                const int effectiveIndex = effectiveIndicesForCoefficientCombinations_.at( j );
                const double coefficientMultiplier = effectiveMutualPotentialField_->getMultiplier(
                        degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );
                if( degreeOfBody2 < currentTransformedCosineCoefficientsBody2_.rows( ) &&
                    orderOfBody2 < currentTransformedCosineCoefficientsBody2_.cols( ) )
                {
                    const Eigen::Vector2d effectiveCoefficientsPartial =
                            coefficientMultiplier *
                            ( Eigen::Vector2d( ) << currentTransformedCosineCoefficientsBody2_( degreeOfBody2, orderOfBody2 ),
                              currentTransformedSineCoefficientsBody2_( degreeOfBody2, orderOfBody2 ) )
                                    .finished( );
                    partial -= bodyFixedCrossProductMatrix *
                            bodyFixedPartialsWrtEffectiveCoefficients.at( effectiveIndex ) *
                            effectiveCoefficientsPartial;
                }
            }
        }

        partialMatrix.col( i ) = partial;
    }
}

void FullTwoBodySphericalHarmonicGravitationalTorquePartial::wrtSineSphericalHarmonicCoefficientsOfBodyUndergoingTorque(
        Eigen::MatrixXd& partialMatrix,
        const std::vector< std::pair< int, int > >& blockIndices )
{
    partialMatrix = Eigen::MatrixXd::Zero( 3, blockIndices.size( ) );

    const Eigen::Matrix3d bodyFixedCrossProductMatrix =
            linear_algebra::getCrossProductMatrix( currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_ );
    const std::vector< Eigen::Matrix< double, 3, 2 > >& bodyFixedPartialsWrtEffectiveCoefficients =
            accelerationPartial_->getCurrentBodyFixedPartialsWrtEffectiveCoefficients( );

    for( unsigned int i = 0; i < blockIndices.size( ); i++ )
    {
        Eigen::Vector3d partial = Eigen::Vector3d::Zero( );
        addBody2SpinTorquePartialWrtBody1Coefficient(
                partial, blockIndices.at( i ).first, blockIndices.at( i ).second, false );

        for( unsigned int j = 0; j < coefficientCombinationsToUse_.size( ); j++ )
        {
            const int degreeOfBody1 = std::get< 0 >( coefficientCombinationsToUse_.at( j ) );
            const int orderOfBody1 = std::get< 1 >( coefficientCombinationsToUse_.at( j ) );
            const int degreeOfBody2 = std::get< 2 >( coefficientCombinationsToUse_.at( j ) );
            const int orderOfBody2 = std::get< 3 >( coefficientCombinationsToUse_.at( j ) );

            if( blockIndices.at( i ).first == degreeOfBody1 && blockIndices.at( i ).second == orderOfBody1 )
            {
                const int effectiveIndex = effectiveIndicesForCoefficientCombinations_.at( j );
                const double coefficientMultiplier = effectiveMutualPotentialField_->getMultiplier(
                        degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );
                if( degreeOfBody2 < currentTransformedCosineCoefficientsBody2_.rows( ) &&
                    orderOfBody2 < currentTransformedCosineCoefficientsBody2_.cols( ) )
                {
                    const Eigen::Vector2d effectiveCoefficientsPartial =
                            coefficientMultiplier *
                            ( Eigen::Vector2d( ) << -currentTransformedSineCoefficientsBody2_( degreeOfBody2, orderOfBody2 ),
                              currentTransformedCosineCoefficientsBody2_( degreeOfBody2, orderOfBody2 ) )
                                    .finished( );
                    partial -= bodyFixedCrossProductMatrix *
                            bodyFixedPartialsWrtEffectiveCoefficients.at( effectiveIndex ) *
                            effectiveCoefficientsPartial;
                }
            }
        }

        partialMatrix.col( i ) = partial;
    }
}

void FullTwoBodySphericalHarmonicGravitationalTorquePartial::wrtCosineSphericalHarmonicCoefficientsOfBodyExertingTorque(
        Eigen::MatrixXd& partialMatrix,
        const std::vector< std::pair< int, int > >& blockIndices )
{
    partialMatrix = Eigen::MatrixXd::Zero( 3, blockIndices.size( ) );

    const Eigen::Matrix3d bodyFixedCrossProductMatrix =
            linear_algebra::getCrossProductMatrix( currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_ );
    const std::vector< Eigen::Matrix< double, 3, 2 > >& bodyFixedPartialsWrtEffectiveCoefficients =
            accelerationPartial_->getCurrentBodyFixedPartialsWrtEffectiveCoefficients( );
    const std::vector< Eigen::Matrix2d >& effectiveCoefficientsWrtTransformedBody2Coefficients =
            accelerationPartial_->getCurrentEffectiveCoefficientsWrtTransformedBody2Coefficients( );

    for( unsigned int i = 0; i < blockIndices.size( ); i++ )
    {
        Eigen::Vector3d partial = Eigen::Vector3d::Zero( );
        addBody2SpinTorquePartialWrtBody2Coefficient(
                partial, blockIndices.at( i ).first, blockIndices.at( i ).second, true );

        Eigen::MatrixXd transformedCosinePartials;
        Eigen::MatrixXd transformedSinePartials;
        accelerationPartial_->calculateCurrentTransformedBody2CoefficientPartials(
                blockIndices.at( i ).first, blockIndices.at( i ).second, true, transformedCosinePartials, transformedSinePartials );

        for( unsigned int j = 0; j < coefficientCombinationsToUse_.size( ); j++ )
        {
            const int degreeOfBody2 = std::get< 2 >( coefficientCombinationsToUse_.at( j ) );
            const int orderOfBody2 = std::get< 3 >( coefficientCombinationsToUse_.at( j ) );

            if( degreeOfBody2 == blockIndices.at( i ).first &&
                degreeOfBody2 < transformedCosinePartials.rows( ) &&
                orderOfBody2 < transformedCosinePartials.cols( ) )
            {
                const int effectiveIndex = effectiveIndicesForCoefficientCombinations_.at( j );
                const Eigen::Vector2d transformedCoefficientPartials =
                        ( Eigen::Vector2d( ) << transformedCosinePartials( degreeOfBody2, orderOfBody2 ),
                          transformedSinePartials( degreeOfBody2, orderOfBody2 ) )
                                .finished( );
                partial -= bodyFixedCrossProductMatrix *
                        bodyFixedPartialsWrtEffectiveCoefficients.at( effectiveIndex ) *
                        ( effectiveCoefficientsWrtTransformedBody2Coefficients.at( effectiveIndex ) *
                          transformedCoefficientPartials );
            }
        }

        partialMatrix.col( i ) = partial;
    }
}

void FullTwoBodySphericalHarmonicGravitationalTorquePartial::wrtSineSphericalHarmonicCoefficientsOfBodyExertingTorque(
        Eigen::MatrixXd& partialMatrix,
        const std::vector< std::pair< int, int > >& blockIndices )
{
    partialMatrix = Eigen::MatrixXd::Zero( 3, blockIndices.size( ) );

    const Eigen::Matrix3d bodyFixedCrossProductMatrix =
            linear_algebra::getCrossProductMatrix( currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_ );
    const std::vector< Eigen::Matrix< double, 3, 2 > >& bodyFixedPartialsWrtEffectiveCoefficients =
            accelerationPartial_->getCurrentBodyFixedPartialsWrtEffectiveCoefficients( );
    const std::vector< Eigen::Matrix2d >& effectiveCoefficientsWrtTransformedBody2Coefficients =
            accelerationPartial_->getCurrentEffectiveCoefficientsWrtTransformedBody2Coefficients( );

    for( unsigned int i = 0; i < blockIndices.size( ); i++ )
    {
        Eigen::Vector3d partial = Eigen::Vector3d::Zero( );
        addBody2SpinTorquePartialWrtBody2Coefficient(
                partial, blockIndices.at( i ).first, blockIndices.at( i ).second, false );

        Eigen::MatrixXd transformedCosinePartials;
        Eigen::MatrixXd transformedSinePartials;
        accelerationPartial_->calculateCurrentTransformedBody2CoefficientPartials(
                blockIndices.at( i ).first, blockIndices.at( i ).second, false, transformedCosinePartials, transformedSinePartials );

        for( unsigned int j = 0; j < coefficientCombinationsToUse_.size( ); j++ )
        {
            const int degreeOfBody2 = std::get< 2 >( coefficientCombinationsToUse_.at( j ) );
            const int orderOfBody2 = std::get< 3 >( coefficientCombinationsToUse_.at( j ) );

            if( degreeOfBody2 == blockIndices.at( i ).first &&
                degreeOfBody2 < transformedCosinePartials.rows( ) &&
                orderOfBody2 < transformedCosinePartials.cols( ) )
            {
                const int effectiveIndex = effectiveIndicesForCoefficientCombinations_.at( j );
                const Eigen::Vector2d transformedCoefficientPartials =
                        ( Eigen::Vector2d( ) << transformedCosinePartials( degreeOfBody2, orderOfBody2 ),
                          transformedSinePartials( degreeOfBody2, orderOfBody2 ) )
                                .finished( );
                partial -= bodyFixedCrossProductMatrix *
                        bodyFixedPartialsWrtEffectiveCoefficients.at( effectiveIndex ) *
                        ( effectiveCoefficientsWrtTransformedBody2Coefficients.at( effectiveIndex ) *
                          transformedCoefficientPartials );
            }
        }

        partialMatrix.col( i ) = partial;
    }
}

void FullTwoBodySphericalHarmonicGravitationalTorquePartial::update( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        torqueModel_->updateMembers( currentTime );
        accelerationPartial_->update( currentTime );

        currentRotationFromInertialToBodyFixedFrameOfBodyUndergoingTorque_ =
                accelerationModel_->getCurrentRotationFromInertialToBody1( ).toRotationMatrix( );
        currentRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorque_ =
                accelerationModel_->getCurrentRotationFromBody2ToBody1( ).toRotationMatrix( );
        currentRotationFromInertialToBodyFixedFrameOfBodyExertingTorque_ =
                currentRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorque_.transpose( ) *
                currentRotationFromInertialToBodyFixedFrameOfBodyUndergoingTorque_;

        currentRelativePositionInInertialFrame_ = accelerationModel_->getCurrentRelativePosition( );
        currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_ =
                accelerationModel_->getCurrentBodyFixedRelativePosition( );
        currentMutualPotentialGradientInBodyFixedFrameOfBodyUndergoingTorque_ =
                accelerationModel_->getMutualPotentialGradient( );

        currentTransformedCosineCoefficientsBody2_ = effectiveMutualPotentialField_->getTransformedCosineCoefficientsOfBody2( );
        currentTransformedSineCoefficientsBody2_ = effectiveMutualPotentialField_->getTransformedSineCoefficientsOfBody2( );
        torqueModel_->computeTransformedAngularMomentumCoefficients(
                effectiveMutualPotentialField_->getCosineCoefficientsOfBody2( ),
                effectiveMutualPotentialField_->getSineCoefficientsOfBody2( ),
                effectiveMutualPotentialField_->getTransformationCache( )->getWignerDMatricesCache( ),
                accelerationModel_->getAreCoefficientsNormalized( ),
                currentTransformedCosineCoefficientsBody2AngularMomentum_,
                currentTransformedSineCoefficientsBody2AngularMomentum_ );

        currentDistance_ = currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_.norm( );
        currentCosineOfLatitude_ = std::sqrt(
                currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_( 0 ) *
                        currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_( 0 ) +
                currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_( 1 ) *
                        currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_( 1 ) ) / currentDistance_;
        currentPreMultiplier_ = accelerationModel_->getCurrentGravitationalParameter( ) / currentDistance_;
        currentRadius1Powers_ = accelerationModel_->getRadius1Powers( );
        currentRadius2Powers_ = accelerationModel_->getRadius2Powers( );

        const Eigen::MatrixXd cosineCoefficientsOfBody1 = effectiveMutualPotentialField_->getCosineCoefficientsOfBody1( );
        const Eigen::MatrixXd sineCoefficientsOfBody1 = effectiveMutualPotentialField_->getSineCoefficientsOfBody1( );
        const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache =
                accelerationModel_->getSphericalHarmonicsCache( );

        currentBody2SpinTorquePartialWrtBodyFixedRelativePosition_.setZero( );

        for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
        {
            const int degreeOfBody1 = std::get< 0 >( coefficientCombinationsToUse_.at( i ) );
            const int orderOfBody1 = std::get< 1 >( coefficientCombinationsToUse_.at( i ) );
            const int degreeOfBody2 = std::get< 2 >( coefficientCombinationsToUse_.at( i ) );
            const int orderOfBody2 = std::get< 3 >( coefficientCombinationsToUse_.at( i ) );

            const double body1CosineCoefficient = cosineCoefficientsOfBody1( degreeOfBody1, orderOfBody1 );
            const double body1SineCoefficient = sineCoefficientsOfBody1( degreeOfBody1, orderOfBody1 );
            const double equatorialRadiusRatioPower =
                    currentRadius1Powers_.at( degreeOfBody1 ) * currentRadius2Powers_.at( degreeOfBody2 );
            const int totalDegree = degreeOfBody1 + degreeOfBody2;

            for( int variant = 0; variant < 4; variant++ )
            {
                int signedOrderOfBody1 = 0;
                int signedOrderOfBody2 = 0;
                if( !detail::getSignedOrdersForVariant(
                            variant, orderOfBody1, orderOfBody2, signedOrderOfBody1, signedOrderOfBody2 ) )
                {
                    continue;
                }

                const int totalOrder = std::abs( signedOrderOfBody1 + signedOrderOfBody2 );
                const Eigen::Matrix< double, 3, 2 > bodyFixedBasisFunctionGradients =
                        detail::computeCurrentBodyFixedBasisFunctionGradients(
                                currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_,
                                sphericalHarmonicsCache,
                                currentCosineOfLatitude_,
                                currentPreMultiplier_,
                                equatorialRadiusRatioPower,
                                totalDegree,
                                totalOrder );

                const double signOrderBody1 = ( signedOrderOfBody1 < 0 ) ? -1.0 : 1.0;
                const double signOrderBody2 = ( signedOrderOfBody2 < 0 ) ? -1.0 : 1.0;
                const double signTotalOrder = ( ( signedOrderOfBody1 + signedOrderOfBody2 ) < 0 ) ? -1.0 : 1.0;
                const double multiplier = effectiveMutualPotentialField_->getMultiplier(
                        degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );

                Eigen::Vector3d transformedAngularMomentumCosineCoefficientsBody2 = Eigen::Vector3d::Zero( );
                Eigen::Vector3d transformedAngularMomentumSineCoefficientsBody2 = Eigen::Vector3d::Zero( );
                for( int component = 0; component < 3; component++ )
                {
                    transformedAngularMomentumCosineCoefficientsBody2( component ) =
                            currentTransformedCosineCoefficientsBody2AngularMomentum_.at( component )(
                                    degreeOfBody2, std::abs( signedOrderOfBody2 ) );
                    transformedAngularMomentumSineCoefficientsBody2( component ) =
                            currentTransformedSineCoefficientsBody2AngularMomentum_.at( component )(
                                    degreeOfBody2, std::abs( signedOrderOfBody2 ) );
                }

                const Eigen::Vector3d effectiveAngularMomentumCosineCoefficients =
                        ( body1CosineCoefficient * transformedAngularMomentumCosineCoefficientsBody2 -
                          signOrderBody1 * signOrderBody2 * body1SineCoefficient *
                                  transformedAngularMomentumSineCoefficientsBody2 ) *
                        multiplier;
                const Eigen::Vector3d effectiveAngularMomentumSineCoefficients =
                        ( signOrderBody2 * body1CosineCoefficient *
                                  transformedAngularMomentumSineCoefficientsBody2 +
                          signOrderBody1 * body1SineCoefficient *
                                  transformedAngularMomentumCosineCoefficientsBody2 ) *
                        signTotalOrder * multiplier;

                currentBody2SpinTorquePartialWrtBodyFixedRelativePosition_ +=
                        effectiveAngularMomentumCosineCoefficients * bodyFixedBasisFunctionGradients.col( 0 ).transpose( ) +
                        effectiveAngularMomentumSineCoefficients * bodyFixedBasisFunctionGradients.col( 1 ).transpose( );
            }
        }

        const Eigen::Matrix3d bodyFixedMutualPotentialTorquePartialWrtBodyFixedRelativePosition =
                linear_algebra::getCrossProductMatrix( currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_ ) *
                        accelerationPartial_->getCurrentBodyFixedPartialWrtPosition( ) -
                linear_algebra::getCrossProductMatrix( currentMutualPotentialGradientInBodyFixedFrameOfBodyUndergoingTorque_ );
        const Eigen::Matrix3d totalTorquePartialWrtBodyFixedRelativePosition =
                currentBody2SpinTorquePartialWrtBodyFixedRelativePosition_ -
                bodyFixedMutualPotentialTorquePartialWrtBodyFixedRelativePosition;

        currentPartialWrtPositionOfBodyUndergoingTorque_ =
                totalTorquePartialWrtBodyFixedRelativePosition *
                currentRotationFromInertialToBodyFixedFrameOfBodyUndergoingTorque_;
        currentPartialWrtPositionOfBodyExertingTorque_ = -currentPartialWrtPositionOfBodyUndergoingTorque_;

        const std::array< Eigen::Matrix3d, 4 > derivativeOfBodyUndergoingTorqueRotationFromBodyFixedToInertialWrtQuaternion =
                detail::getDerivativeOfBodyFixedToInertialRotationMatrixWrtQuaternionForFullTwoBodyTorque(
                        accelerationModel_->getCurrentRotationFromInertialToBody1( ) );

        const Eigen::Vector4d quaternionVectorOfBodyUndergoingTorque =
                linear_algebra::convertQuaternionToVectorFormat(
                        accelerationModel_->getCurrentRotationFromInertialToBody1( ) );
        const Eigen::Vector4d quaternionVectorOfBodyExertingTorque =
                linear_algebra::convertQuaternionToVectorFormat(
                        Eigen::Quaterniond( currentRotationFromInertialToBodyFixedFrameOfBodyExertingTorque_ ) );
        const Eigen::Vector4d conjugatedQuaternionVectorOfBodyExertingTorque =
                ( Eigen::Vector4d( ) << quaternionVectorOfBodyExertingTorque( 0 ),
                  -quaternionVectorOfBodyExertingTorque( 1 ),
                  -quaternionVectorOfBodyExertingTorque( 2 ),
                  -quaternionVectorOfBodyExertingTorque( 3 ) )
                        .finished( );

        const Eigen::Matrix4d partialOfRelativeQuaternionWrtQuaternionOfBodyUndergoingTorque =
                detail::getRightQuaternionMultiplicationMatrix( conjugatedQuaternionVectorOfBodyExertingTorque );
        const Eigen::Matrix4d partialOfRelativeQuaternionWrtQuaternionOfBodyExertingTorque =
                detail::getLeftQuaternionMultiplicationMatrix( quaternionVectorOfBodyUndergoingTorque ) *
                ( Eigen::Vector4d( 1.0, -1.0, -1.0, -1.0 ).asDiagonal( ) );

        const std::array< std::vector< Eigen::MatrixXcd >, 4 > derivativeOfWignerDMatricesWrtRelativeQuaternion =
                detail::computeDerivativeOfWignerDMatricesWrtRelativeQuaternion(
                        accelerationModel_->getCurrentRotationFromBody2ToBody1( ),
                        effectiveMutualPotentialField_->getTransformationCache( )->getWignerDMatricesCache( ) );

        std::array< Eigen::Vector3d, 4 > partialOfBody2SpinTorqueWrtRelativeQuaternion;
        std::array< Eigen::Vector3d, 4 > partialOfMutualPotentialGradientWrtRelativeQuaternion;
        for( int quaternionIndex = 0; quaternionIndex < 4; quaternionIndex++ )
        {
            partialOfBody2SpinTorqueWrtRelativeQuaternion.at( quaternionIndex ).setZero( );
            partialOfMutualPotentialGradientWrtRelativeQuaternion.at( quaternionIndex ).setZero( );

            Eigen::MatrixXd partialOfTransformedCosineCoefficientsBody2;
            Eigen::MatrixXd partialOfTransformedSineCoefficientsBody2;
            detail::transformCoefficientsWithWignerMatrices(
                    effectiveMutualPotentialField_->getCosineCoefficientsOfBody2( ),
                    effectiveMutualPotentialField_->getSineCoefficientsOfBody2( ),
                    derivativeOfWignerDMatricesWrtRelativeQuaternion.at( quaternionIndex ),
                    partialOfTransformedCosineCoefficientsBody2,
                    partialOfTransformedSineCoefficientsBody2,
                    accelerationModel_->getAreCoefficientsNormalized( ) );

            std::array< Eigen::MatrixXd, 3 > partialOfTransformedCosineCoefficientsBody2AngularMomentum;
            std::array< Eigen::MatrixXd, 3 > partialOfTransformedSineCoefficientsBody2AngularMomentum;
            detail::computeTransformedAngularMomentumCoefficientsWithWignerMatrices(
                    effectiveMutualPotentialField_->getCosineCoefficientsOfBody2( ),
                    effectiveMutualPotentialField_->getSineCoefficientsOfBody2( ),
                    derivativeOfWignerDMatricesWrtRelativeQuaternion.at( quaternionIndex ),
                    accelerationModel_->getAreCoefficientsNormalized( ),
                    partialOfTransformedCosineCoefficientsBody2AngularMomentum,
                    partialOfTransformedSineCoefficientsBody2AngularMomentum );

            for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
            {
                const int degreeOfBody1 = std::get< 0 >( coefficientCombinationsToUse_.at( i ) );
                const int orderOfBody1 = std::get< 1 >( coefficientCombinationsToUse_.at( i ) );
                const int degreeOfBody2 = std::get< 2 >( coefficientCombinationsToUse_.at( i ) );
                const int orderOfBody2 = std::get< 3 >( coefficientCombinationsToUse_.at( i ) );

                const double body1CosineCoefficient = cosineCoefficientsOfBody1( degreeOfBody1, orderOfBody1 );
                const double body1SineCoefficient = sineCoefficientsOfBody1( degreeOfBody1, orderOfBody1 );
                const double equatorialRadiusRatioPower =
                        currentRadius1Powers_.at( degreeOfBody1 ) * currentRadius2Powers_.at( degreeOfBody2 );
                const int totalDegree = degreeOfBody1 + degreeOfBody2;

                for( int variant = 0; variant < 4; variant++ )
                {
                    int signedOrderOfBody1 = 0;
                    int signedOrderOfBody2 = 0;
                    if( !detail::getSignedOrdersForVariant(
                                variant, orderOfBody1, orderOfBody2, signedOrderOfBody1, signedOrderOfBody2 ) )
                    {
                        continue;
                    }

                    const int totalOrder = std::abs( signedOrderOfBody1 + signedOrderOfBody2 );
                    const Eigen::Vector2d scalarBasisFunctionValues = detail::computeCurrentScalarBasisFunctionValues(
                            sphericalHarmonicsCache, currentPreMultiplier_, equatorialRadiusRatioPower, totalDegree, totalOrder );
                    const Eigen::Matrix< double, 3, 2 > bodyFixedBasisFunctionGradients =
                            detail::computeCurrentBodyFixedBasisFunctionGradients(
                                    currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_,
                                    sphericalHarmonicsCache,
                                    currentCosineOfLatitude_,
                                    currentPreMultiplier_,
                                    equatorialRadiusRatioPower,
                                    totalDegree,
                                    totalOrder );

                    const double signOrderBody1 = ( signedOrderOfBody1 < 0 ) ? -1.0 : 1.0;
                    const double signOrderBody2 = ( signedOrderOfBody2 < 0 ) ? -1.0 : 1.0;
                    const double signTotalOrder = ( ( signedOrderOfBody1 + signedOrderOfBody2 ) < 0 ) ? -1.0 : 1.0;
                    const double multiplier = effectiveMutualPotentialField_->getMultiplier(
                            degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );

                    const int signedEffectiveIndex = effectiveMutualPotentialField_->getEffectiveIndex(
                            degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );
                    const Eigen::Vector2d partialOfTransformedCoefficients =
                            ( Eigen::Vector2d( ) <<
                              partialOfTransformedCosineCoefficientsBody2( degreeOfBody2, std::abs( signedOrderOfBody2 ) ),
                              partialOfTransformedSineCoefficientsBody2( degreeOfBody2, std::abs( signedOrderOfBody2 ) ) )
                                    .finished( );
                    partialOfMutualPotentialGradientWrtRelativeQuaternion.at( quaternionIndex ) +=
                            bodyFixedBasisFunctionGradients *
                            ( accelerationPartial_->getCurrentEffectiveCoefficientsWrtTransformedBody2Coefficients( ).at(
                                      signedEffectiveIndex ) *
                              partialOfTransformedCoefficients );

                    Eigen::Vector3d partialOfTransformedAngularMomentumCosineCoefficientsBody2 = Eigen::Vector3d::Zero( );
                    Eigen::Vector3d partialOfTransformedAngularMomentumSineCoefficientsBody2 = Eigen::Vector3d::Zero( );
                    for( int component = 0; component < 3; component++ )
                    {
                        partialOfTransformedAngularMomentumCosineCoefficientsBody2( component ) =
                                partialOfTransformedCosineCoefficientsBody2AngularMomentum.at( component )(
                                        degreeOfBody2, std::abs( signedOrderOfBody2 ) );
                        partialOfTransformedAngularMomentumSineCoefficientsBody2( component ) =
                                partialOfTransformedSineCoefficientsBody2AngularMomentum.at( component )(
                                        degreeOfBody2, std::abs( signedOrderOfBody2 ) );
                    }

                    const Eigen::Vector3d partialOfEffectiveAngularMomentumCosineCoefficients =
                            ( body1CosineCoefficient * partialOfTransformedAngularMomentumCosineCoefficientsBody2 -
                              signOrderBody1 * signOrderBody2 * body1SineCoefficient *
                                      partialOfTransformedAngularMomentumSineCoefficientsBody2 ) *
                            multiplier;
                    const Eigen::Vector3d partialOfEffectiveAngularMomentumSineCoefficients =
                            ( signOrderBody2 * body1CosineCoefficient *
                                      partialOfTransformedAngularMomentumSineCoefficientsBody2 +
                              signOrderBody1 * body1SineCoefficient *
                                      partialOfTransformedAngularMomentumCosineCoefficientsBody2 ) *
                            signTotalOrder * multiplier;

                    partialOfBody2SpinTorqueWrtRelativeQuaternion.at( quaternionIndex ) +=
                            partialOfEffectiveAngularMomentumCosineCoefficients * scalarBasisFunctionValues( 0 ) +
                            partialOfEffectiveAngularMomentumSineCoefficients * scalarBasisFunctionValues( 1 );
                }
            }
        }

        const Eigen::Matrix3d bodyFixedCrossProductMatrix =
                linear_algebra::getCrossProductMatrix( currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_ );

        for( int i = 0; i < 4; i++ )
        {
            const Eigen::Vector3d partialOfBodyFixedRelativePositionWrtQuaternionOfBodyUndergoingTorque =
                    derivativeOfBodyUndergoingTorqueRotationFromBodyFixedToInertialWrtQuaternion.at( i ).transpose( ) *
                    currentRelativePositionInInertialFrame_;

            Eigen::Vector3d coefficientContributionWrtQuaternionOfBodyUndergoingTorque = Eigen::Vector3d::Zero( );
            Eigen::Vector3d coefficientContributionWrtQuaternionOfBodyExertingTorque = Eigen::Vector3d::Zero( );
            for( int relativeQuaternionIndex = 0; relativeQuaternionIndex < 4; relativeQuaternionIndex++ )
            {
                coefficientContributionWrtQuaternionOfBodyUndergoingTorque +=
                        ( partialOfBody2SpinTorqueWrtRelativeQuaternion.at( relativeQuaternionIndex ) -
                          bodyFixedCrossProductMatrix *
                                  partialOfMutualPotentialGradientWrtRelativeQuaternion.at( relativeQuaternionIndex ) ) *
                        partialOfRelativeQuaternionWrtQuaternionOfBodyUndergoingTorque( relativeQuaternionIndex, i );
                coefficientContributionWrtQuaternionOfBodyExertingTorque +=
                        ( partialOfBody2SpinTorqueWrtRelativeQuaternion.at( relativeQuaternionIndex ) -
                          bodyFixedCrossProductMatrix *
                                  partialOfMutualPotentialGradientWrtRelativeQuaternion.at( relativeQuaternionIndex ) ) *
                        partialOfRelativeQuaternionWrtQuaternionOfBodyExertingTorque( relativeQuaternionIndex, i );
            }

            currentPartialWrtQuaternionOfBodyUndergoingTorque_.col( i ) =
                    totalTorquePartialWrtBodyFixedRelativePosition *
                            partialOfBodyFixedRelativePositionWrtQuaternionOfBodyUndergoingTorque +
                    coefficientContributionWrtQuaternionOfBodyUndergoingTorque;
            currentPartialWrtQuaternionOfBodyExertingTorque_.col( i ) =
                    coefficientContributionWrtQuaternionOfBodyExertingTorque;
        }

        currentTime_ = currentTime;
    }
}

}  // namespace acceleration_partials

}  // namespace tudat
