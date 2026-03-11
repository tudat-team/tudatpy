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

#include <cmath>

#include "tudat/astro/basic_astro/physicalConstants.h"

namespace tudat
{

namespace gravitation
{

namespace
{

inline double getMomentOfInertiaTerm( const Eigen::Matrix3d& inertiaTensor, const int axis )
{
    return inertiaTensor( axis, axis );
}

inline double getProductOfInertiaTerm( const Eigen::Matrix3d& inertiaTensor, const int firstAxis, const int secondAxis )
{
    return -0.5 * ( inertiaTensor( firstAxis, secondAxis ) + inertiaTensor( secondAxis, firstAxis ) );
}

}  // namespace

Eigen::Vector3d calculateFourthDegreeFullTwoBodyGravitationalTorque(
        const Eigen::Vector3d& relativePositionOfBodyExertingTorqueInBodyFixedFrame,
        const double massOfBodyExertingTorque,
        const Eigen::Matrix3d& inertiaTensorOfBodyUndergoingTorque,
        const Eigen::Matrix3d& inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque )
{
    const double x = relativePositionOfBodyExertingTorqueInBodyFixedFrame.x( );
    const double y = relativePositionOfBodyExertingTorqueInBodyFixedFrame.y( );
    const double z = relativePositionOfBodyExertingTorqueInBodyFixedFrame.z( );

    const double rSquared = relativePositionOfBodyExertingTorqueInBodyFixedFrame.squaredNorm( );
    if( rSquared <= 0.0 )
    {
        throw std::runtime_error(
                "Error when computing fourth-degree full two-body gravitational torque: relative distance is zero." );
    }

    const double inverseRSquared = 1.0 / rSquared;
    const double inverseRFifth = inverseRSquared * inverseRSquared / std::sqrt( rSquared );

    const double xSquaredOverRSquared = x * x * inverseRSquared;
    const double ySquaredOverRSquared = y * y * inverseRSquared;
    const double zSquaredOverRSquared = z * z * inverseRSquared;
    const double xyOverRSquared = x * y * inverseRSquared;
    const double xzOverRSquared = x * z * inverseRSquared;
    const double yzOverRSquared = y * z * inverseRSquared;

    const double A = getMomentOfInertiaTerm( inertiaTensorOfBodyUndergoingTorque, 0 );
    const double B = getMomentOfInertiaTerm( inertiaTensorOfBodyUndergoingTorque, 1 );
    const double C = getMomentOfInertiaTerm( inertiaTensorOfBodyUndergoingTorque, 2 );
    const double Ixy = getProductOfInertiaTerm( inertiaTensorOfBodyUndergoingTorque, 0, 1 );
    const double Ixz = getProductOfInertiaTerm( inertiaTensorOfBodyUndergoingTorque, 0, 2 );
    const double Iyz = getProductOfInertiaTerm( inertiaTensorOfBodyUndergoingTorque, 1, 2 );

    const double APrime =
            getMomentOfInertiaTerm( inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque, 0 );
    const double BPrime =
            getMomentOfInertiaTerm( inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque, 1 );
    const double CPrime =
            getMomentOfInertiaTerm( inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque, 2 );
    const double IxyPrime =
            getProductOfInertiaTerm( inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque, 0, 1 );
    const double IxzPrime =
            getProductOfInertiaTerm( inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque, 0, 2 );
    const double IyzPrime =
            getProductOfInertiaTerm( inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque, 1, 2 );

    const double QPrime = APrime + BPrime + CPrime;
    const double IPrimeAlongRelativePosition = ( APrime * x * x + BPrime * y * y + CPrime * z * z -
                                                  2.0 * IxyPrime * x * y - 2.0 * IyzPrime * y * z -
                                                  2.0 * IxzPrime * x * z ) *
            inverseRSquared;
    const double WPrime =
            massOfBodyExertingTorque + 15.0 * QPrime * inverseRSquared / 2.0 -
            35.0 * IPrimeAlongRelativePosition * inverseRSquared / 2.0;

    const double fyz = y * z * ( WPrime - 5.0 * APrime * inverseRSquared ) - 5.0 * IxzPrime * xyOverRSquared -
            5.0 * IxyPrime * xzOverRSquared + IyzPrime * ( 1.0 - 5.0 * ( ySquaredOverRSquared + zSquaredOverRSquared ) );
    const double fxz = x * z * ( WPrime - 5.0 * BPrime * inverseRSquared ) +
            IxzPrime * ( 1.0 - 5.0 * ( xSquaredOverRSquared + zSquaredOverRSquared ) ) -
            5.0 * IxyPrime * yzOverRSquared - 5.0 * IyzPrime * xyOverRSquared;
    const double fxy = x * y * ( WPrime - 5.0 * CPrime * inverseRSquared ) - 5.0 * IyzPrime * xzOverRSquared +
            IxyPrime * ( 1.0 - 5.0 * ( xSquaredOverRSquared + ySquaredOverRSquared ) ) - 5.0 * IxzPrime * yzOverRSquared;

    const double gyz = ( z * z - y * y ) * WPrime + BPrime - CPrime - 10.0 * IxzPrime * xzOverRSquared -
            10.0 * IxyPrime * xyOverRSquared - 20.0 * IyzPrime * yzOverRSquared -
            5.0 * zSquaredOverRSquared * ( APrime + BPrime - CPrime ) -
            5.0 * ySquaredOverRSquared * ( APrime - BPrime + CPrime );
    const double gxz = ( x * x - z * z ) * WPrime + CPrime - APrime - 20.0 * IxzPrime * xzOverRSquared -
            10.0 * IxyPrime * xyOverRSquared - 10.0 * IyzPrime * yzOverRSquared -
            5.0 * xSquaredOverRSquared * ( -APrime + BPrime + CPrime ) -
            5.0 * zSquaredOverRSquared * ( APrime + BPrime - CPrime );
    const double gxy = ( y * y - x * x ) * WPrime + APrime - BPrime - 10.0 * IxzPrime * xzOverRSquared -
            20.0 * IxyPrime * xyOverRSquared - 10.0 * IyzPrime * yzOverRSquared -
            5.0 * ySquaredOverRSquared * ( APrime - BPrime + CPrime ) -
            5.0 * xSquaredOverRSquared * ( -APrime + BPrime + CPrime );

    const double preMultiplier = 3.0 * physical_constants::GRAVITATIONAL_CONSTANT * inverseRFifth;

    return Eigen::Vector3d(
            preMultiplier * ( ( C - B ) * fyz - Ixz * fxy + Ixy * fxz + Iyz * gyz ),
            preMultiplier * ( ( A - C ) * fxz + Ixz * gxz - Ixy * fyz + Iyz * fxy ),
            preMultiplier * ( ( B - A ) * fxy + Ixz * fyz + Ixy * gxy - Iyz * fxz ) );
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
{
}

void FourthDegreeFullTwoBodyGravitationalTorqueModel::updateMembers( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        currentRotationToBodyFixedFrameOfBodyUndergoingTorque_ = rotationToBodyFixedFrameOfBodyUndergoingTorqueFunction_( );
        currentRotationToBodyFixedFrameOfBodyExertingTorque_ = rotationToBodyFixedFrameOfBodyExertingTorqueFunction_( );
        currentRelativePositionInInertialFrame_ =
                positionOfBodyExertingTorqueFunction_( ) - positionOfBodyUndergoingTorqueFunction_( );
        currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_ =
                currentRotationToBodyFixedFrameOfBodyUndergoingTorque_ * currentRelativePositionInInertialFrame_;

        currentMassOfBodyExertingTorque_ = massOfBodyExertingTorqueFunction_( );
        currentInertiaTensorOfBodyUndergoingTorque_ = inertiaTensorOfBodyUndergoingTorqueFunction_( );
        currentInertiaTensorOfBodyExertingTorque_ = inertiaTensorOfBodyExertingTorqueFunction_( );

        const Eigen::Matrix3d rotationFromBody2ToBody1 =
                currentRotationToBodyFixedFrameOfBodyUndergoingTorque_.toRotationMatrix( ) *
                currentRotationToBodyFixedFrameOfBodyExertingTorque_.toRotationMatrix( ).transpose( );
        currentInertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque_ =
                rotationFromBody2ToBody1 * currentInertiaTensorOfBodyExertingTorque_ * rotationFromBody2ToBody1.transpose( );

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
