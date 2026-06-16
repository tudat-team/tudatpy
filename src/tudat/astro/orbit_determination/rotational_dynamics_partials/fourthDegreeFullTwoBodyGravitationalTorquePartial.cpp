/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/rotational_dynamics_partials/fourthDegreeFullTwoBodyGravitationalTorquePartial.h"

#include <array>

#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicCosineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicSineCoefficients.h"
#include "tudat/math/basic/legendrePolynomials.h"
#include "tudat/math/basic/linearAlgebra.h"

namespace tudat
{

namespace acceleration_partials
{

namespace detail
{

enum FourthDegreeIndependentInertiaTensorComponentIndices {
    aComponentIndex = 0,
    bComponentIndex = 1,
    cComponentIndex = 2,
    ixyComponentIndex = 3,
    ixzComponentIndex = 4,
    iyzComponentIndex = 5
};

enum FourthDegreeAuxiliaryFunctionIndices {
    fyzFunctionIndex = 0,
    fxzFunctionIndex = 1,
    fxyFunctionIndex = 2,
    gyzFunctionIndex = 3,
    gxzFunctionIndex = 4,
    gxyFunctionIndex = 5
};

Eigen::Matrix< double, 6, 1 > getIndependentInertiaTensorComponentsFromMatrix( const Eigen::Matrix3d& inertiaTensor )
{
    Eigen::Matrix< double, 6, 1 > independentInertiaTensorComponents;
    independentInertiaTensorComponents << inertiaTensor( 0, 0 ), inertiaTensor( 1, 1 ), inertiaTensor( 2, 2 ), -inertiaTensor( 0, 1 ),
            -inertiaTensor( 0, 2 ), -inertiaTensor( 1, 2 );
    return independentInertiaTensorComponents;
}

std::array< Eigen::Matrix3d, 4 > getDerivativeOfBodyFixedToInertialRotationMatrixWrtQuaternion(
        const Eigen::Quaterniond& rotationFromInertialToBodyFixedFrame )
{
    std::vector< Eigen::Matrix3d > derivativeList( 4, Eigen::Matrix3d::Zero( ) );
    linear_algebra::computePartialDerivativeOfRotationMatrixWrtQuaternion(
            linear_algebra::convertQuaternionToVectorFormat( rotationFromInertialToBodyFixedFrame.inverse( ) ), derivativeList );

    std::array< Eigen::Matrix3d, 4 > derivativeArray;
    for( int i = 0; i < 4; i++ )
    {
        derivativeArray.at( i ) = derivativeList.at( i );
    }
    return derivativeArray;
}

FourthDegreeTorqueAuxiliaryQuantities computeFourthDegreeTorqueAuxiliaryQuantities(
        const Eigen::Vector3d& relativePositionOfBodyExertingTorqueInBodyFixedFrameOfBodyUndergoingTorque,
        const double massOfBodyExertingTorque,
        const Eigen::Matrix< double, 6, 1 >& independentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque )
{
    FourthDegreeTorqueAuxiliaryQuantities auxiliaryQuantities;

    auxiliaryQuantities.xCoordinate = relativePositionOfBodyExertingTorqueInBodyFixedFrameOfBodyUndergoingTorque( 0 );
    auxiliaryQuantities.yCoordinate = relativePositionOfBodyExertingTorqueInBodyFixedFrameOfBodyUndergoingTorque( 1 );
    auxiliaryQuantities.zCoordinate = relativePositionOfBodyExertingTorqueInBodyFixedFrameOfBodyUndergoingTorque( 2 );

    auxiliaryQuantities.xCoordinateSquared = auxiliaryQuantities.xCoordinate * auxiliaryQuantities.xCoordinate;
    auxiliaryQuantities.yCoordinateSquared = auxiliaryQuantities.yCoordinate * auxiliaryQuantities.yCoordinate;
    auxiliaryQuantities.zCoordinateSquared = auxiliaryQuantities.zCoordinate * auxiliaryQuantities.zCoordinate;
    auxiliaryQuantities.xyTerm = auxiliaryQuantities.xCoordinate * auxiliaryQuantities.yCoordinate;
    auxiliaryQuantities.xzTerm = auxiliaryQuantities.xCoordinate * auxiliaryQuantities.zCoordinate;
    auxiliaryQuantities.yzTerm = auxiliaryQuantities.yCoordinate * auxiliaryQuantities.zCoordinate;

    auxiliaryQuantities.relativeDistanceSquared =
            auxiliaryQuantities.xCoordinateSquared + auxiliaryQuantities.yCoordinateSquared + auxiliaryQuantities.zCoordinateSquared;
    auxiliaryQuantities.inverseRelativeDistanceSquared = 1.0 / auxiliaryQuantities.relativeDistanceSquared;
    auxiliaryQuantities.inverseRelativeDistanceToFourthPower =
            auxiliaryQuantities.inverseRelativeDistanceSquared * auxiliaryQuantities.inverseRelativeDistanceSquared;
    auxiliaryQuantities.relativeDistanceToFifthPower = auxiliaryQuantities.relativeDistanceSquared *
            auxiliaryQuantities.relativeDistanceSquared * std::sqrt( auxiliaryQuantities.relativeDistanceSquared );
    auxiliaryQuantities.torquePrefactor =
            3.0 * physical_constants::GRAVITATIONAL_CONSTANT / auxiliaryQuantities.relativeDistanceToFifthPower;

    auxiliaryQuantities.aComponentOfBodyExertingTorque =
            independentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( aComponentIndex );
    auxiliaryQuantities.bComponentOfBodyExertingTorque =
            independentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( bComponentIndex );
    auxiliaryQuantities.cComponentOfBodyExertingTorque =
            independentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( cComponentIndex );
    auxiliaryQuantities.ixyComponentOfBodyExertingTorque =
            independentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( ixyComponentIndex );
    auxiliaryQuantities.ixzComponentOfBodyExertingTorque =
            independentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( ixzComponentIndex );
    auxiliaryQuantities.iyzComponentOfBodyExertingTorque =
            independentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( iyzComponentIndex );

    auxiliaryQuantities.traceOfInertiaTensorOfBodyExertingTorque = auxiliaryQuantities.aComponentOfBodyExertingTorque +
            auxiliaryQuantities.bComponentOfBodyExertingTorque + auxiliaryQuantities.cComponentOfBodyExertingTorque;

    auxiliaryQuantities.contractedInertiaTensorOfBodyExertingTorque =
            auxiliaryQuantities.aComponentOfBodyExertingTorque * auxiliaryQuantities.xCoordinateSquared +
            auxiliaryQuantities.bComponentOfBodyExertingTorque * auxiliaryQuantities.yCoordinateSquared +
            auxiliaryQuantities.cComponentOfBodyExertingTorque * auxiliaryQuantities.zCoordinateSquared -
            2.0 * auxiliaryQuantities.ixyComponentOfBodyExertingTorque * auxiliaryQuantities.xyTerm -
            2.0 * auxiliaryQuantities.ixzComponentOfBodyExertingTorque * auxiliaryQuantities.xzTerm -
            2.0 * auxiliaryQuantities.iyzComponentOfBodyExertingTorque * auxiliaryQuantities.yzTerm;

    auxiliaryQuantities.wPrimeQuantity = massOfBodyExertingTorque +
            7.5 * auxiliaryQuantities.traceOfInertiaTensorOfBodyExertingTorque * auxiliaryQuantities.inverseRelativeDistanceSquared -
            17.5 * auxiliaryQuantities.contractedInertiaTensorOfBodyExertingTorque *
                    auxiliaryQuantities.inverseRelativeDistanceToFourthPower;

    auxiliaryQuantities.fyzFunction = auxiliaryQuantities.yzTerm *
                    ( auxiliaryQuantities.wPrimeQuantity -
                      5.0 * auxiliaryQuantities.aComponentOfBodyExertingTorque * auxiliaryQuantities.inverseRelativeDistanceSquared ) -
            5.0 * auxiliaryQuantities.ixzComponentOfBodyExertingTorque * auxiliaryQuantities.xyTerm *
                    auxiliaryQuantities.inverseRelativeDistanceSquared -
            5.0 * auxiliaryQuantities.ixyComponentOfBodyExertingTorque * auxiliaryQuantities.xzTerm *
                    auxiliaryQuantities.inverseRelativeDistanceSquared +
            auxiliaryQuantities.iyzComponentOfBodyExertingTorque *
                    ( 1.0 -
                      5.0 * ( auxiliaryQuantities.yCoordinateSquared + auxiliaryQuantities.zCoordinateSquared ) *
                              auxiliaryQuantities.inverseRelativeDistanceSquared );

    auxiliaryQuantities.fxzFunction = auxiliaryQuantities.xzTerm *
                    ( auxiliaryQuantities.wPrimeQuantity -
                      5.0 * auxiliaryQuantities.bComponentOfBodyExertingTorque * auxiliaryQuantities.inverseRelativeDistanceSquared ) +
            auxiliaryQuantities.ixzComponentOfBodyExertingTorque *
                    ( 1.0 -
                      5.0 * ( auxiliaryQuantities.xCoordinateSquared + auxiliaryQuantities.zCoordinateSquared ) *
                              auxiliaryQuantities.inverseRelativeDistanceSquared ) -
            5.0 * auxiliaryQuantities.iyzComponentOfBodyExertingTorque * auxiliaryQuantities.xyTerm *
                    auxiliaryQuantities.inverseRelativeDistanceSquared -
            5.0 * auxiliaryQuantities.ixyComponentOfBodyExertingTorque * auxiliaryQuantities.yzTerm *
                    auxiliaryQuantities.inverseRelativeDistanceSquared;

    auxiliaryQuantities.fxyFunction = auxiliaryQuantities.xyTerm *
                    ( auxiliaryQuantities.wPrimeQuantity -
                      5.0 * auxiliaryQuantities.cComponentOfBodyExertingTorque * auxiliaryQuantities.inverseRelativeDistanceSquared ) -
            5.0 * auxiliaryQuantities.iyzComponentOfBodyExertingTorque * auxiliaryQuantities.xzTerm *
                    auxiliaryQuantities.inverseRelativeDistanceSquared +
            auxiliaryQuantities.ixyComponentOfBodyExertingTorque *
                    ( 1.0 -
                      5.0 * ( auxiliaryQuantities.xCoordinateSquared + auxiliaryQuantities.yCoordinateSquared ) *
                              auxiliaryQuantities.inverseRelativeDistanceSquared ) -
            5.0 * auxiliaryQuantities.ixzComponentOfBodyExertingTorque * auxiliaryQuantities.yzTerm *
                    auxiliaryQuantities.inverseRelativeDistanceSquared;

    auxiliaryQuantities.gyzFunction =
            ( auxiliaryQuantities.zCoordinateSquared - auxiliaryQuantities.yCoordinateSquared ) * auxiliaryQuantities.wPrimeQuantity +
            auxiliaryQuantities.bComponentOfBodyExertingTorque - auxiliaryQuantities.cComponentOfBodyExertingTorque -
            10.0 * auxiliaryQuantities.ixzComponentOfBodyExertingTorque * auxiliaryQuantities.xzTerm *
                    auxiliaryQuantities.inverseRelativeDistanceSquared -
            10.0 * auxiliaryQuantities.ixyComponentOfBodyExertingTorque * auxiliaryQuantities.xyTerm *
                    auxiliaryQuantities.inverseRelativeDistanceSquared -
            20.0 * auxiliaryQuantities.iyzComponentOfBodyExertingTorque * auxiliaryQuantities.yzTerm *
                    auxiliaryQuantities.inverseRelativeDistanceSquared -
            5.0 * auxiliaryQuantities.zCoordinateSquared *
                    ( auxiliaryQuantities.aComponentOfBodyExertingTorque + auxiliaryQuantities.bComponentOfBodyExertingTorque -
                      auxiliaryQuantities.cComponentOfBodyExertingTorque ) *
                    auxiliaryQuantities.inverseRelativeDistanceSquared -
            5.0 * auxiliaryQuantities.yCoordinateSquared *
                    ( auxiliaryQuantities.aComponentOfBodyExertingTorque - auxiliaryQuantities.bComponentOfBodyExertingTorque +
                      auxiliaryQuantities.cComponentOfBodyExertingTorque ) *
                    auxiliaryQuantities.inverseRelativeDistanceSquared;

    auxiliaryQuantities.gxzFunction =
            ( auxiliaryQuantities.xCoordinateSquared - auxiliaryQuantities.zCoordinateSquared ) * auxiliaryQuantities.wPrimeQuantity +
            auxiliaryQuantities.cComponentOfBodyExertingTorque - auxiliaryQuantities.aComponentOfBodyExertingTorque -
            20.0 * auxiliaryQuantities.ixzComponentOfBodyExertingTorque * auxiliaryQuantities.xzTerm *
                    auxiliaryQuantities.inverseRelativeDistanceSquared -
            10.0 * auxiliaryQuantities.ixyComponentOfBodyExertingTorque * auxiliaryQuantities.xyTerm *
                    auxiliaryQuantities.inverseRelativeDistanceSquared -
            10.0 * auxiliaryQuantities.iyzComponentOfBodyExertingTorque * auxiliaryQuantities.yzTerm *
                    auxiliaryQuantities.inverseRelativeDistanceSquared -
            5.0 * auxiliaryQuantities.xCoordinateSquared *
                    ( -auxiliaryQuantities.aComponentOfBodyExertingTorque + auxiliaryQuantities.bComponentOfBodyExertingTorque +
                      auxiliaryQuantities.cComponentOfBodyExertingTorque ) *
                    auxiliaryQuantities.inverseRelativeDistanceSquared -
            5.0 * auxiliaryQuantities.zCoordinateSquared *
                    ( auxiliaryQuantities.aComponentOfBodyExertingTorque + auxiliaryQuantities.bComponentOfBodyExertingTorque -
                      auxiliaryQuantities.cComponentOfBodyExertingTorque ) *
                    auxiliaryQuantities.inverseRelativeDistanceSquared;

    auxiliaryQuantities.gxyFunction =
            ( auxiliaryQuantities.yCoordinateSquared - auxiliaryQuantities.xCoordinateSquared ) * auxiliaryQuantities.wPrimeQuantity +
            auxiliaryQuantities.aComponentOfBodyExertingTorque - auxiliaryQuantities.bComponentOfBodyExertingTorque -
            10.0 * auxiliaryQuantities.ixzComponentOfBodyExertingTorque * auxiliaryQuantities.xzTerm *
                    auxiliaryQuantities.inverseRelativeDistanceSquared -
            20.0 * auxiliaryQuantities.ixyComponentOfBodyExertingTorque * auxiliaryQuantities.xyTerm *
                    auxiliaryQuantities.inverseRelativeDistanceSquared -
            10.0 * auxiliaryQuantities.iyzComponentOfBodyExertingTorque * auxiliaryQuantities.yzTerm *
                    auxiliaryQuantities.inverseRelativeDistanceSquared -
            5.0 * auxiliaryQuantities.yCoordinateSquared *
                    ( auxiliaryQuantities.aComponentOfBodyExertingTorque - auxiliaryQuantities.bComponentOfBodyExertingTorque +
                      auxiliaryQuantities.cComponentOfBodyExertingTorque ) *
                    auxiliaryQuantities.inverseRelativeDistanceSquared -
            5.0 * auxiliaryQuantities.xCoordinateSquared *
                    ( -auxiliaryQuantities.aComponentOfBodyExertingTorque + auxiliaryQuantities.bComponentOfBodyExertingTorque +
                      auxiliaryQuantities.cComponentOfBodyExertingTorque ) *
                    auxiliaryQuantities.inverseRelativeDistanceSquared;

    return auxiliaryQuantities;
}

Eigen::Vector3d computeFourthDegreeTorqueFunctionVector(
        const FourthDegreeTorqueAuxiliaryQuantities& auxiliaryQuantities,
        const Eigen::Matrix< double, 6, 1 >& independentInertiaTensorComponentsOfBodyUndergoingTorque )
{
    const double aComponentOfBodyUndergoingTorque = independentInertiaTensorComponentsOfBodyUndergoingTorque( aComponentIndex );
    const double bComponentOfBodyUndergoingTorque = independentInertiaTensorComponentsOfBodyUndergoingTorque( bComponentIndex );
    const double cComponentOfBodyUndergoingTorque = independentInertiaTensorComponentsOfBodyUndergoingTorque( cComponentIndex );
    const double ixyComponentOfBodyUndergoingTorque = independentInertiaTensorComponentsOfBodyUndergoingTorque( ixyComponentIndex );
    const double ixzComponentOfBodyUndergoingTorque = independentInertiaTensorComponentsOfBodyUndergoingTorque( ixzComponentIndex );
    const double iyzComponentOfBodyUndergoingTorque = independentInertiaTensorComponentsOfBodyUndergoingTorque( iyzComponentIndex );

    Eigen::Vector3d torqueFunctionVector;
    torqueFunctionVector( 0 ) = ( cComponentOfBodyUndergoingTorque - bComponentOfBodyUndergoingTorque ) * auxiliaryQuantities.fyzFunction -
            ixzComponentOfBodyUndergoingTorque * auxiliaryQuantities.fxyFunction +
            ixyComponentOfBodyUndergoingTorque * auxiliaryQuantities.fxzFunction +
            iyzComponentOfBodyUndergoingTorque * auxiliaryQuantities.gyzFunction;
    torqueFunctionVector( 1 ) = ( aComponentOfBodyUndergoingTorque - cComponentOfBodyUndergoingTorque ) * auxiliaryQuantities.fxzFunction +
            ixzComponentOfBodyUndergoingTorque * auxiliaryQuantities.gxzFunction -
            ixyComponentOfBodyUndergoingTorque * auxiliaryQuantities.fyzFunction +
            iyzComponentOfBodyUndergoingTorque * auxiliaryQuantities.fxyFunction;
    torqueFunctionVector( 2 ) = ( bComponentOfBodyUndergoingTorque - aComponentOfBodyUndergoingTorque ) * auxiliaryQuantities.fxyFunction +
            ixzComponentOfBodyUndergoingTorque * auxiliaryQuantities.fyzFunction +
            ixyComponentOfBodyUndergoingTorque * auxiliaryQuantities.gxyFunction -
            iyzComponentOfBodyUndergoingTorque * auxiliaryQuantities.fxzFunction;
    return torqueFunctionVector;
}

Eigen::Matrix< double, 3, 6 > computePartialOfTorqueFunctionVectorWrtAuxiliaryFunctions(
        const Eigen::Matrix< double, 6, 1 >& independentInertiaTensorComponentsOfBodyUndergoingTorque )
{
    const double aComponentOfBodyUndergoingTorque = independentInertiaTensorComponentsOfBodyUndergoingTorque( aComponentIndex );
    const double bComponentOfBodyUndergoingTorque = independentInertiaTensorComponentsOfBodyUndergoingTorque( bComponentIndex );
    const double cComponentOfBodyUndergoingTorque = independentInertiaTensorComponentsOfBodyUndergoingTorque( cComponentIndex );
    const double ixyComponentOfBodyUndergoingTorque = independentInertiaTensorComponentsOfBodyUndergoingTorque( ixyComponentIndex );
    const double ixzComponentOfBodyUndergoingTorque = independentInertiaTensorComponentsOfBodyUndergoingTorque( ixzComponentIndex );
    const double iyzComponentOfBodyUndergoingTorque = independentInertiaTensorComponentsOfBodyUndergoingTorque( iyzComponentIndex );

    Eigen::Matrix< double, 3, 6 > partialOfTorqueFunctionVectorWrtAuxiliaryFunctions = Eigen::Matrix< double, 3, 6 >::Zero( );

    partialOfTorqueFunctionVectorWrtAuxiliaryFunctions( 0, fyzFunctionIndex ) =
            cComponentOfBodyUndergoingTorque - bComponentOfBodyUndergoingTorque;
    partialOfTorqueFunctionVectorWrtAuxiliaryFunctions( 0, fxzFunctionIndex ) = ixyComponentOfBodyUndergoingTorque;
    partialOfTorqueFunctionVectorWrtAuxiliaryFunctions( 0, fxyFunctionIndex ) = -ixzComponentOfBodyUndergoingTorque;
    partialOfTorqueFunctionVectorWrtAuxiliaryFunctions( 0, gyzFunctionIndex ) = iyzComponentOfBodyUndergoingTorque;

    partialOfTorqueFunctionVectorWrtAuxiliaryFunctions( 1, fyzFunctionIndex ) = -ixyComponentOfBodyUndergoingTorque;
    partialOfTorqueFunctionVectorWrtAuxiliaryFunctions( 1, fxzFunctionIndex ) =
            aComponentOfBodyUndergoingTorque - cComponentOfBodyUndergoingTorque;
    partialOfTorqueFunctionVectorWrtAuxiliaryFunctions( 1, fxyFunctionIndex ) = iyzComponentOfBodyUndergoingTorque;
    partialOfTorqueFunctionVectorWrtAuxiliaryFunctions( 1, gxzFunctionIndex ) = ixzComponentOfBodyUndergoingTorque;

    partialOfTorqueFunctionVectorWrtAuxiliaryFunctions( 2, fyzFunctionIndex ) = ixzComponentOfBodyUndergoingTorque;
    partialOfTorqueFunctionVectorWrtAuxiliaryFunctions( 2, fxzFunctionIndex ) = -iyzComponentOfBodyUndergoingTorque;
    partialOfTorqueFunctionVectorWrtAuxiliaryFunctions( 2, fxyFunctionIndex ) =
            bComponentOfBodyUndergoingTorque - aComponentOfBodyUndergoingTorque;
    partialOfTorqueFunctionVectorWrtAuxiliaryFunctions( 2, gxyFunctionIndex ) = ixyComponentOfBodyUndergoingTorque;

    return partialOfTorqueFunctionVectorWrtAuxiliaryFunctions;
}

Eigen::Matrix< double, 3, 6 > computePartialOfTorqueWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque(
        const FourthDegreeTorqueAuxiliaryQuantities& auxiliaryQuantities )
{
    Eigen::Matrix< double, 3, 6 > partialOfTorqueWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque =
            Eigen::Matrix< double, 3, 6 >::Zero( );

    partialOfTorqueWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque( 0, bComponentIndex ) =
            -auxiliaryQuantities.torquePrefactor * auxiliaryQuantities.fyzFunction;
    partialOfTorqueWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque( 0, cComponentIndex ) =
            auxiliaryQuantities.torquePrefactor * auxiliaryQuantities.fyzFunction;
    partialOfTorqueWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque( 0, ixyComponentIndex ) =
            auxiliaryQuantities.torquePrefactor * auxiliaryQuantities.fxzFunction;
    partialOfTorqueWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque( 0, ixzComponentIndex ) =
            -auxiliaryQuantities.torquePrefactor * auxiliaryQuantities.fxyFunction;
    partialOfTorqueWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque( 0, iyzComponentIndex ) =
            auxiliaryQuantities.torquePrefactor * auxiliaryQuantities.gyzFunction;

    partialOfTorqueWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque( 1, aComponentIndex ) =
            auxiliaryQuantities.torquePrefactor * auxiliaryQuantities.fxzFunction;
    partialOfTorqueWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque( 1, cComponentIndex ) =
            -auxiliaryQuantities.torquePrefactor * auxiliaryQuantities.fxzFunction;
    partialOfTorqueWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque( 1, ixyComponentIndex ) =
            -auxiliaryQuantities.torquePrefactor * auxiliaryQuantities.fyzFunction;
    partialOfTorqueWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque( 1, ixzComponentIndex ) =
            auxiliaryQuantities.torquePrefactor * auxiliaryQuantities.gxzFunction;
    partialOfTorqueWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque( 1, iyzComponentIndex ) =
            auxiliaryQuantities.torquePrefactor * auxiliaryQuantities.fxyFunction;

    partialOfTorqueWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque( 2, aComponentIndex ) =
            -auxiliaryQuantities.torquePrefactor * auxiliaryQuantities.fxyFunction;
    partialOfTorqueWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque( 2, bComponentIndex ) =
            auxiliaryQuantities.torquePrefactor * auxiliaryQuantities.fxyFunction;
    partialOfTorqueWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque( 2, ixyComponentIndex ) =
            auxiliaryQuantities.torquePrefactor * auxiliaryQuantities.gxyFunction;
    partialOfTorqueWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque( 2, ixzComponentIndex ) =
            auxiliaryQuantities.torquePrefactor * auxiliaryQuantities.fyzFunction;
    partialOfTorqueWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque( 2, iyzComponentIndex ) =
            -auxiliaryQuantities.torquePrefactor * auxiliaryQuantities.fxzFunction;

    return partialOfTorqueWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque;
}

double computePartialOfContractedInertiaTensorOfBodyExertingTorqueWrtCoordinate(
        const FourthDegreeTorqueAuxiliaryQuantities& auxiliaryQuantities,
        const int coordinateIndex )
{
    switch( coordinateIndex )
    {
        case 0:
            return 2.0 * auxiliaryQuantities.aComponentOfBodyExertingTorque * auxiliaryQuantities.xCoordinate -
                    2.0 * auxiliaryQuantities.ixyComponentOfBodyExertingTorque * auxiliaryQuantities.yCoordinate -
                    2.0 * auxiliaryQuantities.ixzComponentOfBodyExertingTorque * auxiliaryQuantities.zCoordinate;
        case 1:
            return 2.0 * auxiliaryQuantities.bComponentOfBodyExertingTorque * auxiliaryQuantities.yCoordinate -
                    2.0 * auxiliaryQuantities.ixyComponentOfBodyExertingTorque * auxiliaryQuantities.xCoordinate -
                    2.0 * auxiliaryQuantities.iyzComponentOfBodyExertingTorque * auxiliaryQuantities.zCoordinate;
        case 2:
            return 2.0 * auxiliaryQuantities.cComponentOfBodyExertingTorque * auxiliaryQuantities.zCoordinate -
                    2.0 * auxiliaryQuantities.ixzComponentOfBodyExertingTorque * auxiliaryQuantities.xCoordinate -
                    2.0 * auxiliaryQuantities.iyzComponentOfBodyExertingTorque * auxiliaryQuantities.yCoordinate;
        default:
            throw std::runtime_error( "Error when computing position partial of fourth-degree torque: invalid coordinate index." );
    }
}

Eigen::Matrix< double, 6, 1 > computePartialOfAuxiliaryFunctionsWrtPositionCoordinate(
        const FourthDegreeTorqueAuxiliaryQuantities& auxiliaryQuantities,
        const int coordinateIndex )
{
    const double coordinateValue = ( coordinateIndex == 0 ) ? auxiliaryQuantities.xCoordinate
            : ( coordinateIndex == 1 )                      ? auxiliaryQuantities.yCoordinate
                                                            : auxiliaryQuantities.zCoordinate;
    const double derivativeOfInverseRelativeDistanceSquaredWrtCoordinate =
            -2.0 * coordinateValue * auxiliaryQuantities.inverseRelativeDistanceToFourthPower;
    const double derivativeOfContractedInertiaTensorOfBodyExertingTorqueWrtCoordinate =
            computePartialOfContractedInertiaTensorOfBodyExertingTorqueWrtCoordinate( auxiliaryQuantities, coordinateIndex );
    const double derivativeOfWPrimeQuantityWrtCoordinate =
            7.5 * auxiliaryQuantities.traceOfInertiaTensorOfBodyExertingTorque * derivativeOfInverseRelativeDistanceSquaredWrtCoordinate -
            17.5 *
                    ( derivativeOfContractedInertiaTensorOfBodyExertingTorqueWrtCoordinate *
                              auxiliaryQuantities.inverseRelativeDistanceToFourthPower +
                      2.0 * auxiliaryQuantities.contractedInertiaTensorOfBodyExertingTorque *
                              auxiliaryQuantities.inverseRelativeDistanceSquared *
                              derivativeOfInverseRelativeDistanceSquaredWrtCoordinate );

    const double derivativeOfxyTermWrtCoordinate = ( coordinateIndex == 0 ) ? auxiliaryQuantities.yCoordinate
            : ( coordinateIndex == 1 )                                      ? auxiliaryQuantities.xCoordinate
                                                                            : 0.0;
    const double derivativeOfxzTermWrtCoordinate = ( coordinateIndex == 0 ) ? auxiliaryQuantities.zCoordinate
            : ( coordinateIndex == 2 )                                      ? auxiliaryQuantities.xCoordinate
                                                                            : 0.0;
    const double derivativeOfyzTermWrtCoordinate = ( coordinateIndex == 1 ) ? auxiliaryQuantities.zCoordinate
            : ( coordinateIndex == 2 )                                      ? auxiliaryQuantities.yCoordinate
                                                                            : 0.0;

    const double derivativeOfxSquaredPluszSquaredWrtCoordinate = ( coordinateIndex == 0 ) ? 2.0 * auxiliaryQuantities.xCoordinate
            : ( coordinateIndex == 2 )                                                    ? 2.0 * auxiliaryQuantities.zCoordinate
                                                                                          : 0.0;
    const double derivativeOfxSquaredPlusySquaredWrtCoordinate = ( coordinateIndex == 0 ) ? 2.0 * auxiliaryQuantities.xCoordinate
            : ( coordinateIndex == 1 )                                                    ? 2.0 * auxiliaryQuantities.yCoordinate
                                                                                          : 0.0;
    const double derivativeOfySquaredPluszSquaredWrtCoordinate = ( coordinateIndex == 1 ) ? 2.0 * auxiliaryQuantities.yCoordinate
            : ( coordinateIndex == 2 )                                                    ? 2.0 * auxiliaryQuantities.zCoordinate
                                                                                          : 0.0;

    const double derivativeOfzSquaredMinusySquaredWrtCoordinate = ( coordinateIndex == 2 ) ? 2.0 * auxiliaryQuantities.zCoordinate
            : ( coordinateIndex == 1 )                                                     ? -2.0 * auxiliaryQuantities.yCoordinate
                                                                                           : 0.0;
    const double derivativeOfxSquaredMinuszSquaredWrtCoordinate = ( coordinateIndex == 0 ) ? 2.0 * auxiliaryQuantities.xCoordinate
            : ( coordinateIndex == 2 )                                                     ? -2.0 * auxiliaryQuantities.zCoordinate
                                                                                           : 0.0;
    const double derivativeOfySquaredMinusxSquaredWrtCoordinate = ( coordinateIndex == 1 ) ? 2.0 * auxiliaryQuantities.yCoordinate
            : ( coordinateIndex == 0 )                                                     ? -2.0 * auxiliaryQuantities.xCoordinate
                                                                                           : 0.0;

    const double aPlusBMinusC = auxiliaryQuantities.aComponentOfBodyExertingTorque + auxiliaryQuantities.bComponentOfBodyExertingTorque -
            auxiliaryQuantities.cComponentOfBodyExertingTorque;
    const double aMinusBPlusC = auxiliaryQuantities.aComponentOfBodyExertingTorque - auxiliaryQuantities.bComponentOfBodyExertingTorque +
            auxiliaryQuantities.cComponentOfBodyExertingTorque;
    const double minusAPlusBPlusC = -auxiliaryQuantities.aComponentOfBodyExertingTorque +
            auxiliaryQuantities.bComponentOfBodyExertingTorque + auxiliaryQuantities.cComponentOfBodyExertingTorque;

    Eigen::Matrix< double, 6, 1 > partialOfAuxiliaryFunctionsWrtPositionCoordinate;

    partialOfAuxiliaryFunctionsWrtPositionCoordinate( fyzFunctionIndex ) = derivativeOfyzTermWrtCoordinate *
                    ( auxiliaryQuantities.wPrimeQuantity -
                      5.0 * auxiliaryQuantities.aComponentOfBodyExertingTorque * auxiliaryQuantities.inverseRelativeDistanceSquared ) +
            auxiliaryQuantities.yzTerm *
                    ( derivativeOfWPrimeQuantityWrtCoordinate -
                      5.0 * auxiliaryQuantities.aComponentOfBodyExertingTorque * derivativeOfInverseRelativeDistanceSquaredWrtCoordinate ) -
            5.0 * auxiliaryQuantities.ixzComponentOfBodyExertingTorque *
                    ( derivativeOfxyTermWrtCoordinate * auxiliaryQuantities.inverseRelativeDistanceSquared +
                      auxiliaryQuantities.xyTerm * derivativeOfInverseRelativeDistanceSquaredWrtCoordinate ) -
            5.0 * auxiliaryQuantities.ixyComponentOfBodyExertingTorque *
                    ( derivativeOfxzTermWrtCoordinate * auxiliaryQuantities.inverseRelativeDistanceSquared +
                      auxiliaryQuantities.xzTerm * derivativeOfInverseRelativeDistanceSquaredWrtCoordinate ) +
            auxiliaryQuantities.iyzComponentOfBodyExertingTorque *
                    ( -5.0 * derivativeOfySquaredPluszSquaredWrtCoordinate * auxiliaryQuantities.inverseRelativeDistanceSquared -
                      5.0 * ( auxiliaryQuantities.yCoordinateSquared + auxiliaryQuantities.zCoordinateSquared ) *
                              derivativeOfInverseRelativeDistanceSquaredWrtCoordinate );

    partialOfAuxiliaryFunctionsWrtPositionCoordinate( fxzFunctionIndex ) = derivativeOfxzTermWrtCoordinate *
                    ( auxiliaryQuantities.wPrimeQuantity -
                      5.0 * auxiliaryQuantities.bComponentOfBodyExertingTorque * auxiliaryQuantities.inverseRelativeDistanceSquared ) +
            auxiliaryQuantities.xzTerm *
                    ( derivativeOfWPrimeQuantityWrtCoordinate -
                      5.0 * auxiliaryQuantities.bComponentOfBodyExertingTorque * derivativeOfInverseRelativeDistanceSquaredWrtCoordinate ) +
            auxiliaryQuantities.ixzComponentOfBodyExertingTorque *
                    ( -5.0 * derivativeOfxSquaredPluszSquaredWrtCoordinate * auxiliaryQuantities.inverseRelativeDistanceSquared -
                      5.0 * ( auxiliaryQuantities.xCoordinateSquared + auxiliaryQuantities.zCoordinateSquared ) *
                              derivativeOfInverseRelativeDistanceSquaredWrtCoordinate ) -
            5.0 * auxiliaryQuantities.iyzComponentOfBodyExertingTorque *
                    ( derivativeOfxyTermWrtCoordinate * auxiliaryQuantities.inverseRelativeDistanceSquared +
                      auxiliaryQuantities.xyTerm * derivativeOfInverseRelativeDistanceSquaredWrtCoordinate ) -
            5.0 * auxiliaryQuantities.ixyComponentOfBodyExertingTorque *
                    ( derivativeOfyzTermWrtCoordinate * auxiliaryQuantities.inverseRelativeDistanceSquared +
                      auxiliaryQuantities.yzTerm * derivativeOfInverseRelativeDistanceSquaredWrtCoordinate );

    partialOfAuxiliaryFunctionsWrtPositionCoordinate( fxyFunctionIndex ) = derivativeOfxyTermWrtCoordinate *
                    ( auxiliaryQuantities.wPrimeQuantity -
                      5.0 * auxiliaryQuantities.cComponentOfBodyExertingTorque * auxiliaryQuantities.inverseRelativeDistanceSquared ) +
            auxiliaryQuantities.xyTerm *
                    ( derivativeOfWPrimeQuantityWrtCoordinate -
                      5.0 * auxiliaryQuantities.cComponentOfBodyExertingTorque * derivativeOfInverseRelativeDistanceSquaredWrtCoordinate ) -
            5.0 * auxiliaryQuantities.iyzComponentOfBodyExertingTorque *
                    ( derivativeOfxzTermWrtCoordinate * auxiliaryQuantities.inverseRelativeDistanceSquared +
                      auxiliaryQuantities.xzTerm * derivativeOfInverseRelativeDistanceSquaredWrtCoordinate ) +
            auxiliaryQuantities.ixyComponentOfBodyExertingTorque *
                    ( -5.0 * derivativeOfxSquaredPlusySquaredWrtCoordinate * auxiliaryQuantities.inverseRelativeDistanceSquared -
                      5.0 * ( auxiliaryQuantities.xCoordinateSquared + auxiliaryQuantities.yCoordinateSquared ) *
                              derivativeOfInverseRelativeDistanceSquaredWrtCoordinate ) -
            5.0 * auxiliaryQuantities.ixzComponentOfBodyExertingTorque *
                    ( derivativeOfyzTermWrtCoordinate * auxiliaryQuantities.inverseRelativeDistanceSquared +
                      auxiliaryQuantities.yzTerm * derivativeOfInverseRelativeDistanceSquaredWrtCoordinate );

    partialOfAuxiliaryFunctionsWrtPositionCoordinate( gyzFunctionIndex ) =
            derivativeOfzSquaredMinusySquaredWrtCoordinate * auxiliaryQuantities.wPrimeQuantity +
            ( auxiliaryQuantities.zCoordinateSquared - auxiliaryQuantities.yCoordinateSquared ) * derivativeOfWPrimeQuantityWrtCoordinate -
            10.0 * auxiliaryQuantities.ixzComponentOfBodyExertingTorque *
                    ( derivativeOfxzTermWrtCoordinate * auxiliaryQuantities.inverseRelativeDistanceSquared +
                      auxiliaryQuantities.xzTerm * derivativeOfInverseRelativeDistanceSquaredWrtCoordinate ) -
            10.0 * auxiliaryQuantities.ixyComponentOfBodyExertingTorque *
                    ( derivativeOfxyTermWrtCoordinate * auxiliaryQuantities.inverseRelativeDistanceSquared +
                      auxiliaryQuantities.xyTerm * derivativeOfInverseRelativeDistanceSquaredWrtCoordinate ) -
            20.0 * auxiliaryQuantities.iyzComponentOfBodyExertingTorque *
                    ( derivativeOfyzTermWrtCoordinate * auxiliaryQuantities.inverseRelativeDistanceSquared +
                      auxiliaryQuantities.yzTerm * derivativeOfInverseRelativeDistanceSquaredWrtCoordinate ) -
            5.0 * derivativeOfInverseRelativeDistanceSquaredWrtCoordinate *
                    ( auxiliaryQuantities.zCoordinateSquared * aPlusBMinusC + auxiliaryQuantities.yCoordinateSquared * aMinusBPlusC ) -
            5.0 * auxiliaryQuantities.inverseRelativeDistanceSquared *
                    ( ( coordinateIndex == 2 ? 2.0 * auxiliaryQuantities.zCoordinate : 0.0 ) * aPlusBMinusC +
                      ( coordinateIndex == 1 ? 2.0 * auxiliaryQuantities.yCoordinate : 0.0 ) * aMinusBPlusC );

    partialOfAuxiliaryFunctionsWrtPositionCoordinate( gxzFunctionIndex ) =
            derivativeOfxSquaredMinuszSquaredWrtCoordinate * auxiliaryQuantities.wPrimeQuantity +
            ( auxiliaryQuantities.xCoordinateSquared - auxiliaryQuantities.zCoordinateSquared ) * derivativeOfWPrimeQuantityWrtCoordinate -
            20.0 * auxiliaryQuantities.ixzComponentOfBodyExertingTorque *
                    ( derivativeOfxzTermWrtCoordinate * auxiliaryQuantities.inverseRelativeDistanceSquared +
                      auxiliaryQuantities.xzTerm * derivativeOfInverseRelativeDistanceSquaredWrtCoordinate ) -
            10.0 * auxiliaryQuantities.ixyComponentOfBodyExertingTorque *
                    ( derivativeOfxyTermWrtCoordinate * auxiliaryQuantities.inverseRelativeDistanceSquared +
                      auxiliaryQuantities.xyTerm * derivativeOfInverseRelativeDistanceSquaredWrtCoordinate ) -
            10.0 * auxiliaryQuantities.iyzComponentOfBodyExertingTorque *
                    ( derivativeOfyzTermWrtCoordinate * auxiliaryQuantities.inverseRelativeDistanceSquared +
                      auxiliaryQuantities.yzTerm * derivativeOfInverseRelativeDistanceSquaredWrtCoordinate ) -
            5.0 * derivativeOfInverseRelativeDistanceSquaredWrtCoordinate *
                    ( auxiliaryQuantities.xCoordinateSquared * minusAPlusBPlusC + auxiliaryQuantities.zCoordinateSquared * aPlusBMinusC ) -
            5.0 * auxiliaryQuantities.inverseRelativeDistanceSquared *
                    ( ( coordinateIndex == 0 ? 2.0 * auxiliaryQuantities.xCoordinate : 0.0 ) * minusAPlusBPlusC +
                      ( coordinateIndex == 2 ? 2.0 * auxiliaryQuantities.zCoordinate : 0.0 ) * aPlusBMinusC );

    partialOfAuxiliaryFunctionsWrtPositionCoordinate( gxyFunctionIndex ) =
            derivativeOfySquaredMinusxSquaredWrtCoordinate * auxiliaryQuantities.wPrimeQuantity +
            ( auxiliaryQuantities.yCoordinateSquared - auxiliaryQuantities.xCoordinateSquared ) * derivativeOfWPrimeQuantityWrtCoordinate -
            10.0 * auxiliaryQuantities.ixzComponentOfBodyExertingTorque *
                    ( derivativeOfxzTermWrtCoordinate * auxiliaryQuantities.inverseRelativeDistanceSquared +
                      auxiliaryQuantities.xzTerm * derivativeOfInverseRelativeDistanceSquaredWrtCoordinate ) -
            20.0 * auxiliaryQuantities.ixyComponentOfBodyExertingTorque *
                    ( derivativeOfxyTermWrtCoordinate * auxiliaryQuantities.inverseRelativeDistanceSquared +
                      auxiliaryQuantities.xyTerm * derivativeOfInverseRelativeDistanceSquaredWrtCoordinate ) -
            10.0 * auxiliaryQuantities.iyzComponentOfBodyExertingTorque *
                    ( derivativeOfyzTermWrtCoordinate * auxiliaryQuantities.inverseRelativeDistanceSquared +
                      auxiliaryQuantities.yzTerm * derivativeOfInverseRelativeDistanceSquaredWrtCoordinate ) -
            5.0 * derivativeOfInverseRelativeDistanceSquaredWrtCoordinate *
                    ( auxiliaryQuantities.yCoordinateSquared * aMinusBPlusC + auxiliaryQuantities.xCoordinateSquared * minusAPlusBPlusC ) -
            5.0 * auxiliaryQuantities.inverseRelativeDistanceSquared *
                    ( ( coordinateIndex == 1 ? 2.0 * auxiliaryQuantities.yCoordinate : 0.0 ) * aMinusBPlusC +
                      ( coordinateIndex == 0 ? 2.0 * auxiliaryQuantities.xCoordinate : 0.0 ) * minusAPlusBPlusC );

    return partialOfAuxiliaryFunctionsWrtPositionCoordinate;
}

Eigen::Matrix3d computePartialOfTorqueWrtBodyFixedRelativePosition(
        const FourthDegreeTorqueAuxiliaryQuantities& auxiliaryQuantities,
        const Eigen::Matrix< double, 6, 1 >& independentInertiaTensorComponentsOfBodyUndergoingTorque )
{
    Eigen::Matrix3d partialOfTorqueWrtBodyFixedRelativePosition = Eigen::Matrix3d::Zero( );
    const Eigen::Vector3d torqueFunctionVector =
            computeFourthDegreeTorqueFunctionVector( auxiliaryQuantities, independentInertiaTensorComponentsOfBodyUndergoingTorque );
    const Eigen::Matrix< double, 3, 6 > partialOfTorqueFunctionVectorWrtAuxiliaryFunctions =
            computePartialOfTorqueFunctionVectorWrtAuxiliaryFunctions( independentInertiaTensorComponentsOfBodyUndergoingTorque );

    const std::array< double, 3 > coordinateValues = { auxiliaryQuantities.xCoordinate,
                                                       auxiliaryQuantities.yCoordinate,
                                                       auxiliaryQuantities.zCoordinate };

    for( int coordinateIndex = 0; coordinateIndex < 3; coordinateIndex++ )
    {
        const double derivativeOfTorquePrefactorWrtCoordinate = -5.0 * coordinateValues.at( coordinateIndex ) *
                auxiliaryQuantities.inverseRelativeDistanceSquared * auxiliaryQuantities.torquePrefactor;
        const Eigen::Matrix< double, 6, 1 > partialOfAuxiliaryFunctionsWrtPositionCoordinate =
                computePartialOfAuxiliaryFunctionsWrtPositionCoordinate( auxiliaryQuantities, coordinateIndex );
        partialOfTorqueWrtBodyFixedRelativePosition.col( coordinateIndex ) =
                derivativeOfTorquePrefactorWrtCoordinate * torqueFunctionVector +
                auxiliaryQuantities.torquePrefactor * partialOfTorqueFunctionVectorWrtAuxiliaryFunctions *
                        partialOfAuxiliaryFunctionsWrtPositionCoordinate;
    }

    return partialOfTorqueWrtBodyFixedRelativePosition;
}

Eigen::Matrix< double, 6, 6 > computePartialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque(
        const FourthDegreeTorqueAuxiliaryQuantities& auxiliaryQuantities )
{
    Eigen::Matrix< double, 6, 6 > partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque =
            Eigen::Matrix< double, 6, 6 >::Zero( );

    const std::array< double, 6 > partialOfWPrimeQuantityWrtIndependentInertiaTensorComponentsOfBodyExertingTorque = {
        7.5 * auxiliaryQuantities.inverseRelativeDistanceSquared -
                17.5 * auxiliaryQuantities.xCoordinateSquared * auxiliaryQuantities.inverseRelativeDistanceToFourthPower,
        7.5 * auxiliaryQuantities.inverseRelativeDistanceSquared -
                17.5 * auxiliaryQuantities.yCoordinateSquared * auxiliaryQuantities.inverseRelativeDistanceToFourthPower,
        7.5 * auxiliaryQuantities.inverseRelativeDistanceSquared -
                17.5 * auxiliaryQuantities.zCoordinateSquared * auxiliaryQuantities.inverseRelativeDistanceToFourthPower,
        35.0 * auxiliaryQuantities.xyTerm * auxiliaryQuantities.inverseRelativeDistanceToFourthPower,
        35.0 * auxiliaryQuantities.xzTerm * auxiliaryQuantities.inverseRelativeDistanceToFourthPower,
        35.0 * auxiliaryQuantities.yzTerm * auxiliaryQuantities.inverseRelativeDistanceToFourthPower
    };

    for( int inertiaComponentIndex = 0; inertiaComponentIndex < 6; inertiaComponentIndex++ )
    {
        partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( fyzFunctionIndex, inertiaComponentIndex ) =
                auxiliaryQuantities.yzTerm *
                partialOfWPrimeQuantityWrtIndependentInertiaTensorComponentsOfBodyExertingTorque.at( inertiaComponentIndex );
        partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( fxzFunctionIndex, inertiaComponentIndex ) =
                auxiliaryQuantities.xzTerm *
                partialOfWPrimeQuantityWrtIndependentInertiaTensorComponentsOfBodyExertingTorque.at( inertiaComponentIndex );
        partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( fxyFunctionIndex, inertiaComponentIndex ) =
                auxiliaryQuantities.xyTerm *
                partialOfWPrimeQuantityWrtIndependentInertiaTensorComponentsOfBodyExertingTorque.at( inertiaComponentIndex );
        partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( gyzFunctionIndex, inertiaComponentIndex ) =
                ( auxiliaryQuantities.zCoordinateSquared - auxiliaryQuantities.yCoordinateSquared ) *
                partialOfWPrimeQuantityWrtIndependentInertiaTensorComponentsOfBodyExertingTorque.at( inertiaComponentIndex );
        partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( gxzFunctionIndex, inertiaComponentIndex ) =
                ( auxiliaryQuantities.xCoordinateSquared - auxiliaryQuantities.zCoordinateSquared ) *
                partialOfWPrimeQuantityWrtIndependentInertiaTensorComponentsOfBodyExertingTorque.at( inertiaComponentIndex );
        partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( gxyFunctionIndex, inertiaComponentIndex ) =
                ( auxiliaryQuantities.yCoordinateSquared - auxiliaryQuantities.xCoordinateSquared ) *
                partialOfWPrimeQuantityWrtIndependentInertiaTensorComponentsOfBodyExertingTorque.at( inertiaComponentIndex );
    }

    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( fyzFunctionIndex, aComponentIndex ) +=
            -5.0 * auxiliaryQuantities.yzTerm * auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( fyzFunctionIndex, ixyComponentIndex ) +=
            -5.0 * auxiliaryQuantities.xzTerm * auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( fyzFunctionIndex, ixzComponentIndex ) +=
            -5.0 * auxiliaryQuantities.xyTerm * auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( fyzFunctionIndex, iyzComponentIndex ) += 1.0 -
            5.0 * ( auxiliaryQuantities.yCoordinateSquared + auxiliaryQuantities.zCoordinateSquared ) *
                    auxiliaryQuantities.inverseRelativeDistanceSquared;

    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( fxzFunctionIndex, bComponentIndex ) +=
            -5.0 * auxiliaryQuantities.xzTerm * auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( fxzFunctionIndex, ixyComponentIndex ) +=
            -5.0 * auxiliaryQuantities.yzTerm * auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( fxzFunctionIndex, ixzComponentIndex ) += 1.0 -
            5.0 * ( auxiliaryQuantities.xCoordinateSquared + auxiliaryQuantities.zCoordinateSquared ) *
                    auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( fxzFunctionIndex, iyzComponentIndex ) +=
            -5.0 * auxiliaryQuantities.xyTerm * auxiliaryQuantities.inverseRelativeDistanceSquared;

    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( fxyFunctionIndex, cComponentIndex ) +=
            -5.0 * auxiliaryQuantities.xyTerm * auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( fxyFunctionIndex, ixyComponentIndex ) += 1.0 -
            5.0 * ( auxiliaryQuantities.xCoordinateSquared + auxiliaryQuantities.yCoordinateSquared ) *
                    auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( fxyFunctionIndex, ixzComponentIndex ) +=
            -5.0 * auxiliaryQuantities.yzTerm * auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( fxyFunctionIndex, iyzComponentIndex ) +=
            -5.0 * auxiliaryQuantities.xzTerm * auxiliaryQuantities.inverseRelativeDistanceSquared;

    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( gyzFunctionIndex, aComponentIndex ) += -5.0 *
            ( auxiliaryQuantities.zCoordinateSquared + auxiliaryQuantities.yCoordinateSquared ) *
            auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( gyzFunctionIndex, bComponentIndex ) += 1.0 -
            5.0 * auxiliaryQuantities.zCoordinateSquared * auxiliaryQuantities.inverseRelativeDistanceSquared +
            5.0 * auxiliaryQuantities.yCoordinateSquared * auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( gyzFunctionIndex, cComponentIndex ) += -1.0 +
            5.0 * auxiliaryQuantities.zCoordinateSquared * auxiliaryQuantities.inverseRelativeDistanceSquared -
            5.0 * auxiliaryQuantities.yCoordinateSquared * auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( gyzFunctionIndex, ixyComponentIndex ) +=
            -10.0 * auxiliaryQuantities.xyTerm * auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( gyzFunctionIndex, ixzComponentIndex ) +=
            -10.0 * auxiliaryQuantities.xzTerm * auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( gyzFunctionIndex, iyzComponentIndex ) +=
            -20.0 * auxiliaryQuantities.yzTerm * auxiliaryQuantities.inverseRelativeDistanceSquared;

    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( gxzFunctionIndex, aComponentIndex ) += -1.0 +
            5.0 * auxiliaryQuantities.xCoordinateSquared * auxiliaryQuantities.inverseRelativeDistanceSquared -
            5.0 * auxiliaryQuantities.zCoordinateSquared * auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( gxzFunctionIndex, bComponentIndex ) += -5.0 *
            ( auxiliaryQuantities.xCoordinateSquared + auxiliaryQuantities.zCoordinateSquared ) *
            auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( gxzFunctionIndex, cComponentIndex ) += 1.0 -
            5.0 * auxiliaryQuantities.xCoordinateSquared * auxiliaryQuantities.inverseRelativeDistanceSquared +
            5.0 * auxiliaryQuantities.zCoordinateSquared * auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( gxzFunctionIndex, ixyComponentIndex ) +=
            -10.0 * auxiliaryQuantities.xyTerm * auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( gxzFunctionIndex, ixzComponentIndex ) +=
            -20.0 * auxiliaryQuantities.xzTerm * auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( gxzFunctionIndex, iyzComponentIndex ) +=
            -10.0 * auxiliaryQuantities.yzTerm * auxiliaryQuantities.inverseRelativeDistanceSquared;

    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( gxyFunctionIndex, aComponentIndex ) += 1.0 -
            5.0 * auxiliaryQuantities.yCoordinateSquared * auxiliaryQuantities.inverseRelativeDistanceSquared +
            5.0 * auxiliaryQuantities.xCoordinateSquared * auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( gxyFunctionIndex, bComponentIndex ) += -1.0 +
            5.0 * auxiliaryQuantities.yCoordinateSquared * auxiliaryQuantities.inverseRelativeDistanceSquared -
            5.0 * auxiliaryQuantities.xCoordinateSquared * auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( gxyFunctionIndex, cComponentIndex ) += -5.0 *
            ( auxiliaryQuantities.yCoordinateSquared + auxiliaryQuantities.xCoordinateSquared ) *
            auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( gxyFunctionIndex, ixyComponentIndex ) +=
            -20.0 * auxiliaryQuantities.xyTerm * auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( gxyFunctionIndex, ixzComponentIndex ) +=
            -10.0 * auxiliaryQuantities.xzTerm * auxiliaryQuantities.inverseRelativeDistanceSquared;
    partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( gxyFunctionIndex, iyzComponentIndex ) +=
            -10.0 * auxiliaryQuantities.yzTerm * auxiliaryQuantities.inverseRelativeDistanceSquared;

    return partialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque;
}

Eigen::Matrix< double, 3, 6 > computePartialOfTorqueWrtIndependentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque(
        const FourthDegreeTorqueAuxiliaryQuantities& auxiliaryQuantities,
        const Eigen::Matrix< double, 6, 1 >& independentInertiaTensorComponentsOfBodyUndergoingTorque )
{
    return auxiliaryQuantities.torquePrefactor *
            computePartialOfTorqueFunctionVectorWrtAuxiliaryFunctions( independentInertiaTensorComponentsOfBodyUndergoingTorque ) *
            computePartialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque( auxiliaryQuantities );
}

}  // namespace detail

FourthDegreeFullTwoBodyGravitationalTorquePartial::FourthDegreeFullTwoBodyGravitationalTorquePartial(
        const std::shared_ptr< gravitation::FourthDegreeFullTwoBodyGravitationalTorqueModel > torqueModel,
        const std::shared_ptr< gravitation::SphericalHarmonicsGravityField > gravityFieldOfBodyUndergoingTorque,
        const std::shared_ptr< gravitation::SphericalHarmonicsGravityField > gravityFieldOfBodyExertingTorque,
        const std::string& acceleratedBody,
        const std::string& acceleratingBody ):
    TorquePartial( acceleratedBody, acceleratingBody, basic_astrodynamics::fourth_degree_full_two_body_gravitational_torque ),
    torqueModel_( torqueModel ), gravityFieldOfBodyUndergoingTorque_( gravityFieldOfBodyUndergoingTorque ),
    gravityFieldOfBodyExertingTorque_( gravityFieldOfBodyExertingTorque ),
    currentPartialWrtQuaternionOfBodyUndergoingTorque_( Eigen::Matrix< double, 3, 4 >::Zero( ) ),
    currentPartialWrtQuaternionOfBodyExertingTorque_( Eigen::Matrix< double, 3, 4 >::Zero( ) ),
    currentPartialWrtPositionOfBodyUndergoingTorque_( Eigen::Matrix3d::Zero( ) ),
    currentPartialWrtPositionOfBodyExertingTorque_( Eigen::Matrix3d::Zero( ) ),
    currentPartialWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque_( Eigen::Matrix< double, 3, 6 >::Zero( ) ),
    currentPartialWrtIndependentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque_(
            Eigen::Matrix< double, 3, 6 >::Zero( ) ),
    currentRotationFromInertialToBodyFixedFrameOfBodyUndergoingTorque_( Eigen::Matrix3d::Identity( ) ),
    currentRotationFromInertialToBodyFixedFrameOfBodyExertingTorque_( Eigen::Matrix3d::Identity( ) ),
    currentRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorque_( Eigen::Matrix3d::Identity( ) ),
    currentRelativePositionOfBodyExertingTorqueInInertialFrame_( Eigen::Vector3d::Zero( ) ),
    currentRelativePositionOfBodyExertingTorqueInBodyFixedFrameOfBodyUndergoingTorque_( Eigen::Vector3d::Zero( ) ),
    currentInertiaTensorOfBodyExertingTorque_( Eigen::Matrix3d::Zero( ) )
{
    if( torqueModel_ == nullptr )
    {
        throw std::runtime_error( "Error when creating FourthDegreeFullTwoBodyGravitationalTorquePartial, torque model is nullptr." );
    }
}

std::pair< std::function< void( Eigen::MatrixXd& ) >, int > FourthDegreeFullTwoBodyGravitationalTorquePartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > )
{
    return std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );
}

std::pair< std::function< void( Eigen::MatrixXd& ) >, int > FourthDegreeFullTwoBodyGravitationalTorquePartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    using namespace estimatable_parameters;

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunction =
            std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );

    if( parameter->getParameterName( ).second.first == bodyUndergoingTorque_ )
    {
        switch( parameter->getParameterName( ).first )
        {
            case spherical_harmonics_cosine_coefficient_block: {
                std::shared_ptr< SphericalHarmonicsCosineCoefficients > coefficientsParameter =
                        std::dynamic_pointer_cast< SphericalHarmonicsCosineCoefficients >( parameter );
                int c20Index, c21Index, c22Index;
                coefficientsParameter->getDegreeTwoEntries( c20Index, c21Index, c22Index );
                if( c20Index >= 0 || c21Index >= 0 || c22Index >= 0 )
                {
                    partialFunction = std::make_pair( std::bind( &FourthDegreeFullTwoBodyGravitationalTorquePartial::
                                                                         wrtCosineSphericalHarmonicCoefficientsOfBodyUndergoingTorque,
                                                                 this,
                                                                 std::placeholders::_1,
                                                                 c20Index,
                                                                 c21Index,
                                                                 c22Index ),
                                                      coefficientsParameter->getParameterSize( ) );
                }
                break;
            }
            case spherical_harmonics_sine_coefficient_block: {
                std::shared_ptr< SphericalHarmonicsSineCoefficients > coefficientsParameter =
                        std::dynamic_pointer_cast< SphericalHarmonicsSineCoefficients >( parameter );
                int s21Index, s22Index;
                coefficientsParameter->getDegreeTwoEntries( s21Index, s22Index );
                if( s21Index >= 0 || s22Index >= 0 )
                {
                    partialFunction = std::make_pair( std::bind( &FourthDegreeFullTwoBodyGravitationalTorquePartial::
                                                                         wrtSineSphericalHarmonicCoefficientsOfBodyUndergoingTorque,
                                                                 this,
                                                                 std::placeholders::_1,
                                                                 s21Index,
                                                                 s22Index ),
                                                      coefficientsParameter->getParameterSize( ) );
                }
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
            case spherical_harmonics_cosine_coefficient_block: {
                std::shared_ptr< SphericalHarmonicsCosineCoefficients > coefficientsParameter =
                        std::dynamic_pointer_cast< SphericalHarmonicsCosineCoefficients >( parameter );
                int c20Index, c21Index, c22Index;
                coefficientsParameter->getDegreeTwoEntries( c20Index, c21Index, c22Index );
                if( c20Index >= 0 || c21Index >= 0 || c22Index >= 0 )
                {
                    partialFunction = std::make_pair( std::bind( &FourthDegreeFullTwoBodyGravitationalTorquePartial::
                                                                         wrtCosineSphericalHarmonicCoefficientsOfBodyExertingTorque,
                                                                 this,
                                                                 std::placeholders::_1,
                                                                 c20Index,
                                                                 c21Index,
                                                                 c22Index ),
                                                      coefficientsParameter->getParameterSize( ) );
                }
                break;
            }
            case spherical_harmonics_sine_coefficient_block: {
                std::shared_ptr< SphericalHarmonicsSineCoefficients > coefficientsParameter =
                        std::dynamic_pointer_cast< SphericalHarmonicsSineCoefficients >( parameter );
                int s21Index, s22Index;
                coefficientsParameter->getDegreeTwoEntries( s21Index, s22Index );
                if( s21Index >= 0 || s22Index >= 0 )
                {
                    partialFunction = std::make_pair( std::bind( &FourthDegreeFullTwoBodyGravitationalTorquePartial::
                                                                         wrtSineSphericalHarmonicCoefficientsOfBodyExertingTorque,
                                                                 this,
                                                                 std::placeholders::_1,
                                                                 s21Index,
                                                                 s22Index ),
                                                      coefficientsParameter->getParameterSize( ) );
                }
                break;
            }
            default:
                break;
        }
    }

    return partialFunction;
}

void FourthDegreeFullTwoBodyGravitationalTorquePartial::wrtOrientationOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
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

void FourthDegreeFullTwoBodyGravitationalTorquePartial::wrtOrientationOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
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

bool FourthDegreeFullTwoBodyGravitationalTorquePartial::isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
        const std::pair< std::string, std::string >& stateReferencePoint,
        const propagators::IntegratedStateType integratedStateType )
{
    return ( integratedStateType == propagators::translational_state &&
             ( stateReferencePoint.first == bodyUndergoingTorque_ || stateReferencePoint.first == bodyExertingTorque_ ) );
}

void FourthDegreeFullTwoBodyGravitationalTorquePartial::wrtNonRotationalStateOfAdditionalBody(
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

Eigen::Matrix< double, 6, 1 > FourthDegreeFullTwoBodyGravitationalTorquePartial::getIndependentInertiaTensorComponentsFromMatrixDerivative(
        const Eigen::Matrix3d& inertiaTensorDerivative ) const
{
    return detail::getIndependentInertiaTensorComponentsFromMatrix( inertiaTensorDerivative );
}

Eigen::Matrix3d FourthDegreeFullTwoBodyGravitationalTorquePartial::getInertiaTensorPartialWrtNormalizedCosineCoefficient(
        const unsigned int order,
        const std::shared_ptr< gravitation::SphericalHarmonicsGravityField >& gravityField ) const
{
    if( gravityField == nullptr )
    {
        throw std::runtime_error(
                "Error when computing fourth-degree torque partial w.r.t. cosine coefficient: gravity field is nullptr." );
    }

    const double scalingFactor = basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, order ) *
            gravityField->getInertiaTensorNormalizationFactor( );

    Eigen::Matrix3d unscaledInertiaTensorPartial = Eigen::Matrix3d::Zero( );
    switch( order )
    {
        case 0:
            unscaledInertiaTensorPartial = UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C20;
            break;
        case 1:
            unscaledInertiaTensorPartial = UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C21;
            break;
        case 2:
            unscaledInertiaTensorPartial = UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_C22;
            break;
        default:
            throw std::runtime_error(
                    "Error when computing fourth-degree torque partial w.r.t. cosine coefficient: invalid degree-2 order." );
    }

    return scalingFactor * unscaledInertiaTensorPartial;
}

Eigen::Matrix3d FourthDegreeFullTwoBodyGravitationalTorquePartial::getInertiaTensorPartialWrtNormalizedSineCoefficient(
        const unsigned int order,
        const std::shared_ptr< gravitation::SphericalHarmonicsGravityField >& gravityField ) const
{
    if( gravityField == nullptr )
    {
        throw std::runtime_error( "Error when computing fourth-degree torque partial w.r.t. sine coefficient: gravity field is nullptr." );
    }

    const double scalingFactor = basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, order ) *
            gravityField->getInertiaTensorNormalizationFactor( );

    Eigen::Matrix3d unscaledInertiaTensorPartial = Eigen::Matrix3d::Zero( );
    switch( order )
    {
        case 1:
            unscaledInertiaTensorPartial = UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_S21;
            break;
        case 2:
            unscaledInertiaTensorPartial = UNSCALED_INERTIAL_TENSOR_PARTIAL_WRT_S22;
            break;
        default:
            throw std::runtime_error(
                    "Error when computing fourth-degree torque partial w.r.t. sine coefficient: invalid degree-2 order." );
    }

    return scalingFactor * unscaledInertiaTensorPartial;
}

Eigen::Matrix< double, 6, 1 >
FourthDegreeFullTwoBodyGravitationalTorquePartial::getIndependentInertiaTensorComponentsPartialWrtNormalizedCosineCoefficient(
        const unsigned int order,
        const std::shared_ptr< gravitation::SphericalHarmonicsGravityField >& gravityField ) const
{
    return getIndependentInertiaTensorComponentsFromMatrixDerivative(
            getInertiaTensorPartialWrtNormalizedCosineCoefficient( order, gravityField ) );
}

Eigen::Matrix< double, 6, 1 >
FourthDegreeFullTwoBodyGravitationalTorquePartial::getIndependentInertiaTensorComponentsPartialWrtNormalizedSineCoefficient(
        const unsigned int order,
        const std::shared_ptr< gravitation::SphericalHarmonicsGravityField >& gravityField ) const
{
    return getIndependentInertiaTensorComponentsFromMatrixDerivative(
            getInertiaTensorPartialWrtNormalizedSineCoefficient( order, gravityField ) );
}

void FourthDegreeFullTwoBodyGravitationalTorquePartial::wrtCosineSphericalHarmonicCoefficientsOfBodyUndergoingTorque(
        Eigen::MatrixXd& partialMatrix,
        const int c20Index,
        const int c21Index,
        const int c22Index )
{
    partialMatrix.setZero( );

    if( c20Index >= 0 )
    {
        partialMatrix.col( c20Index ) = currentPartialWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque_ *
                getIndependentInertiaTensorComponentsPartialWrtNormalizedCosineCoefficient( 0, gravityFieldOfBodyUndergoingTorque_ );
    }
    if( c21Index >= 0 )
    {
        partialMatrix.col( c21Index ) = currentPartialWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque_ *
                getIndependentInertiaTensorComponentsPartialWrtNormalizedCosineCoefficient( 1, gravityFieldOfBodyUndergoingTorque_ );
    }
    if( c22Index >= 0 )
    {
        partialMatrix.col( c22Index ) = currentPartialWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque_ *
                getIndependentInertiaTensorComponentsPartialWrtNormalizedCosineCoefficient( 2, gravityFieldOfBodyUndergoingTorque_ );
    }
}

void FourthDegreeFullTwoBodyGravitationalTorquePartial::wrtSineSphericalHarmonicCoefficientsOfBodyUndergoingTorque(
        Eigen::MatrixXd& partialMatrix,
        const int s21Index,
        const int s22Index )
{
    partialMatrix.setZero( );

    if( s21Index >= 0 )
    {
        partialMatrix.col( s21Index ) = currentPartialWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque_ *
                getIndependentInertiaTensorComponentsPartialWrtNormalizedSineCoefficient( 1, gravityFieldOfBodyUndergoingTorque_ );
    }
    if( s22Index >= 0 )
    {
        partialMatrix.col( s22Index ) = currentPartialWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque_ *
                getIndependentInertiaTensorComponentsPartialWrtNormalizedSineCoefficient( 2, gravityFieldOfBodyUndergoingTorque_ );
    }
}

void FourthDegreeFullTwoBodyGravitationalTorquePartial::wrtCosineSphericalHarmonicCoefficientsOfBodyExertingTorque(
        Eigen::MatrixXd& partialMatrix,
        const int c20Index,
        const int c21Index,
        const int c22Index )
{
    partialMatrix.setZero( );

    if( c20Index >= 0 )
    {
        const Eigen::Matrix3d inertiaTensorPartialInBodyFixedFrameOfBodyUndergoingTorque =
                currentRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorque_ *
                getInertiaTensorPartialWrtNormalizedCosineCoefficient( 0, gravityFieldOfBodyExertingTorque_ ) *
                currentRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorque_.transpose( );
        partialMatrix.col( c20Index ) =
                currentPartialWrtIndependentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque_ *
                getIndependentInertiaTensorComponentsFromMatrixDerivative( inertiaTensorPartialInBodyFixedFrameOfBodyUndergoingTorque );
    }
    if( c21Index >= 0 )
    {
        const Eigen::Matrix3d inertiaTensorPartialInBodyFixedFrameOfBodyUndergoingTorque =
                currentRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorque_ *
                getInertiaTensorPartialWrtNormalizedCosineCoefficient( 1, gravityFieldOfBodyExertingTorque_ ) *
                currentRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorque_.transpose( );
        partialMatrix.col( c21Index ) =
                currentPartialWrtIndependentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque_ *
                getIndependentInertiaTensorComponentsFromMatrixDerivative( inertiaTensorPartialInBodyFixedFrameOfBodyUndergoingTorque );
    }
    if( c22Index >= 0 )
    {
        const Eigen::Matrix3d inertiaTensorPartialInBodyFixedFrameOfBodyUndergoingTorque =
                currentRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorque_ *
                getInertiaTensorPartialWrtNormalizedCosineCoefficient( 2, gravityFieldOfBodyExertingTorque_ ) *
                currentRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorque_.transpose( );
        partialMatrix.col( c22Index ) =
                currentPartialWrtIndependentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque_ *
                getIndependentInertiaTensorComponentsFromMatrixDerivative( inertiaTensorPartialInBodyFixedFrameOfBodyUndergoingTorque );
    }
}

void FourthDegreeFullTwoBodyGravitationalTorquePartial::wrtSineSphericalHarmonicCoefficientsOfBodyExertingTorque(
        Eigen::MatrixXd& partialMatrix,
        const int s21Index,
        const int s22Index )
{
    partialMatrix.setZero( );

    if( s21Index >= 0 )
    {
        const Eigen::Matrix3d inertiaTensorPartialInBodyFixedFrameOfBodyUndergoingTorque =
                currentRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorque_ *
                getInertiaTensorPartialWrtNormalizedSineCoefficient( 1, gravityFieldOfBodyExertingTorque_ ) *
                currentRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorque_.transpose( );
        partialMatrix.col( s21Index ) =
                currentPartialWrtIndependentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque_ *
                getIndependentInertiaTensorComponentsFromMatrixDerivative( inertiaTensorPartialInBodyFixedFrameOfBodyUndergoingTorque );
    }
    if( s22Index >= 0 )
    {
        const Eigen::Matrix3d inertiaTensorPartialInBodyFixedFrameOfBodyUndergoingTorque =
                currentRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorque_ *
                getInertiaTensorPartialWrtNormalizedSineCoefficient( 2, gravityFieldOfBodyExertingTorque_ ) *
                currentRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorque_.transpose( );
        partialMatrix.col( s22Index ) =
                currentPartialWrtIndependentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque_ *
                getIndependentInertiaTensorComponentsFromMatrixDerivative( inertiaTensorPartialInBodyFixedFrameOfBodyUndergoingTorque );
    }
}

void FourthDegreeFullTwoBodyGravitationalTorquePartial::update( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        torqueModel_->updateMembers( currentTime );

        currentRotationFromInertialToBodyFixedFrameOfBodyUndergoingTorque_ =
                torqueModel_->getCurrentRotationToBodyFixedFrameOfBodyUndergoingTorque( ).toRotationMatrix( );
        currentRotationFromInertialToBodyFixedFrameOfBodyExertingTorque_ =
                torqueModel_->getCurrentRotationToBodyFixedFrameOfBodyExertingTorque( ).toRotationMatrix( );
        currentRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorque_ =
                currentRotationFromInertialToBodyFixedFrameOfBodyUndergoingTorque_ *
                currentRotationFromInertialToBodyFixedFrameOfBodyExertingTorque_.transpose( );
        currentRelativePositionOfBodyExertingTorqueInInertialFrame_ = torqueModel_->getCurrentRelativePositionInInertialFrame( );
        currentRelativePositionOfBodyExertingTorqueInBodyFixedFrameOfBodyUndergoingTorque_ =
                torqueModel_->getCurrentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque( );
        currentInertiaTensorOfBodyExertingTorque_ = torqueModel_->getCurrentInertiaTensorOfBodyExertingTorque( );

        const Eigen::Matrix< double, 6, 1 > independentInertiaTensorComponentsOfBodyUndergoingTorque =
                detail::getIndependentInertiaTensorComponentsFromMatrix( torqueModel_->getCurrentInertiaTensorOfBodyUndergoingTorque( ) );
        const Eigen::Matrix< double, 6, 1 > independentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque =
                detail::getIndependentInertiaTensorComponentsFromMatrix(
                        torqueModel_->getCurrentInertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( ) );

        const detail::FourthDegreeTorqueAuxiliaryQuantities auxiliaryQuantities = detail::computeFourthDegreeTorqueAuxiliaryQuantities(
                currentRelativePositionOfBodyExertingTorqueInBodyFixedFrameOfBodyUndergoingTorque_,
                torqueModel_->getCurrentMassOfBodyExertingTorque( ),
                independentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque );

        currentPartialWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque_ =
                detail::computePartialOfTorqueWrtIndependentInertiaTensorComponentsOfBodyUndergoingTorque( auxiliaryQuantities );
        currentPartialWrtIndependentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque_ =
                detail::computePartialOfTorqueWrtIndependentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque(
                        auxiliaryQuantities, independentInertiaTensorComponentsOfBodyUndergoingTorque );

        const Eigen::Matrix3d partialOfTorqueWrtBodyFixedRelativePosition = detail::computePartialOfTorqueWrtBodyFixedRelativePosition(
                auxiliaryQuantities, independentInertiaTensorComponentsOfBodyUndergoingTorque );
        currentPartialWrtPositionOfBodyExertingTorque_ =
                partialOfTorqueWrtBodyFixedRelativePosition * currentRotationFromInertialToBodyFixedFrameOfBodyUndergoingTorque_;
        currentPartialWrtPositionOfBodyUndergoingTorque_ = -currentPartialWrtPositionOfBodyExertingTorque_;

        const std::array< Eigen::Matrix3d, 4 > derivativeOfBodyUndergoingTorqueRotationFromBodyFixedToInertialWrtQuaternion =
                detail::getDerivativeOfBodyFixedToInertialRotationMatrixWrtQuaternion(
                        torqueModel_->getCurrentRotationToBodyFixedFrameOfBodyUndergoingTorque( ) );
        const std::array< Eigen::Matrix3d, 4 > derivativeOfBodyExertingTorqueRotationFromBodyFixedToInertialWrtQuaternion =
                detail::getDerivativeOfBodyFixedToInertialRotationMatrixWrtQuaternion(
                        torqueModel_->getCurrentRotationToBodyFixedFrameOfBodyExertingTorque( ) );

        for( int i = 0; i < 4; i++ )
        {
            const Eigen::Vector3d partialOfBodyFixedRelativePositionWrtQuaternionOfBodyUndergoingTorque =
                    derivativeOfBodyUndergoingTorqueRotationFromBodyFixedToInertialWrtQuaternion.at( i ).transpose( ) *
                    currentRelativePositionOfBodyExertingTorqueInInertialFrame_;

            const Eigen::Matrix3d
                    partialOfRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorqueWrtQuaternionOfBodyUndergoingTorque =
                            derivativeOfBodyUndergoingTorqueRotationFromBodyFixedToInertialWrtQuaternion.at( i ).transpose( ) *
                    currentRotationFromInertialToBodyFixedFrameOfBodyExertingTorque_.transpose( );

            const Eigen::Matrix3d partialOfInertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorqueWrtQuaternionOfBodyUndergoingTorque =
                    partialOfRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorqueWrtQuaternionOfBodyUndergoingTorque *
                            currentInertiaTensorOfBodyExertingTorque_ *
                            currentRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorque_.transpose( ) +
                    currentRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorque_ *
                            currentInertiaTensorOfBodyExertingTorque_ *
                            partialOfRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorqueWrtQuaternionOfBodyUndergoingTorque
                                    .transpose( );

            currentPartialWrtQuaternionOfBodyUndergoingTorque_.col( i ) =
                    partialOfTorqueWrtBodyFixedRelativePosition * partialOfBodyFixedRelativePositionWrtQuaternionOfBodyUndergoingTorque +
                    currentPartialWrtIndependentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque_ *
                            getIndependentInertiaTensorComponentsFromMatrixDerivative(
                                    partialOfInertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorqueWrtQuaternionOfBodyUndergoingTorque );

            const Eigen::Matrix3d
                    partialOfRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorqueWrtQuaternionOfBodyExertingTorque =
                            currentRotationFromInertialToBodyFixedFrameOfBodyUndergoingTorque_ *
                    derivativeOfBodyExertingTorqueRotationFromBodyFixedToInertialWrtQuaternion.at( i );

            const Eigen::Matrix3d partialOfInertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorqueWrtQuaternionOfBodyExertingTorque =
                    partialOfRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorqueWrtQuaternionOfBodyExertingTorque *
                            currentInertiaTensorOfBodyExertingTorque_ *
                            currentRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorque_.transpose( ) +
                    currentRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorque_ *
                            currentInertiaTensorOfBodyExertingTorque_ *
                            partialOfRotationFromBodyFixedFrameOfBodyExertingTorqueToBodyFixedFrameOfBodyUndergoingTorqueWrtQuaternionOfBodyExertingTorque
                                    .transpose( );

            currentPartialWrtQuaternionOfBodyExertingTorque_.col( i ) =
                    currentPartialWrtIndependentInertiaTensorComponentsOfBodyExertingTorqueInFrameOfBodyUndergoingTorque_ *
                    getIndependentInertiaTensorComponentsFromMatrixDerivative(
                            partialOfInertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorqueWrtQuaternionOfBodyExertingTorque );
        }

        currentTime_ = currentTime;
    }
}

}  // namespace acceleration_partials

}  // namespace tudat
