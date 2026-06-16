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

namespace tudat
{

namespace gravitation
{

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
    const double inverseR2 = 1.0 / r2;
    const double r5 = r2 * r2 * std::sqrt( r2 );

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

    // Eq. (11): torque components in body-1-fixed frame.
    const double prefactor = ( 3.0 * physical_constants::GRAVITATIONAL_CONSTANT ) / r5;
    Eigen::Vector3d torque;
    const double torqueXFromDiagonalTerm = ( C - B ) * fyz;
    const double torqueXFromIxzTerm = -Ixz * fxy;
    const double torqueXFromIxyTerm = Ixy * fxz;
    const double torqueXFromIyzTerm = Iyz * gyz;

    const double torqueYFromDiagonalTerm = ( A - C ) * fxz;
    const double torqueYFromIxzTerm = Ixz * gxz;
    const double torqueYFromIxyTerm = -Ixy * fyz;
    const double torqueYFromIyzTerm = Iyz * fxy;

    const double torqueZFromDiagonalTerm = ( B - A ) * fxy;
    const double torqueZFromIxzTerm = Ixz * fyz;
    const double torqueZFromIxyTerm = Ixy * gxy;
    const double torqueZFromIyzTerm = -Iyz * fxz;

    torque( 0 ) = prefactor * ( torqueXFromDiagonalTerm + torqueXFromIxzTerm + torqueXFromIxyTerm + torqueXFromIyzTerm );
    torque( 1 ) = prefactor * ( torqueYFromDiagonalTerm + torqueYFromIxzTerm + torqueYFromIxyTerm + torqueYFromIyzTerm );
    torque( 2 ) = prefactor * ( torqueZFromDiagonalTerm + torqueZFromIxzTerm + torqueZFromIxyTerm + torqueZFromIyzTerm );

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
