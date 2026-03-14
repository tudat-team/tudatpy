/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_FOURTHDEGREEFULLTWOBODYGRAVITATIONALTORQUE_H
#define TUDAT_FOURTHDEGREEFULLTWOBODYGRAVITATIONALTORQUE_H

#include <cmath>
#include <functional>

#include <Eigen/Geometry>

#include "tudat/astro/basic_astro/torqueModel.h"
#include "tudat/astro/basic_astro/physicalConstants.h"

namespace tudat
{

namespace gravitation
{

template< typename ScalarType >
Eigen::Matrix< ScalarType, 3, 1 > calculateFourthDegreeFullTwoBodyGravitationalTorqueFromTensorComponents(
        const Eigen::Matrix< ScalarType, 3, 1 >& relativePositionOfBodyExertingTorqueInBodyFixedFrame,
        const ScalarType massOfBodyExertingTorque,
        const Eigen::Matrix< ScalarType, 3, 3 >& inertiaTensorOfBodyUndergoingTorque,
        const Eigen::Matrix< ScalarType, 3, 3 >& inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque )
{
    const ScalarType x = relativePositionOfBodyExertingTorqueInBodyFixedFrame( 0 );
    const ScalarType y = relativePositionOfBodyExertingTorqueInBodyFixedFrame( 1 );
    const ScalarType z = relativePositionOfBodyExertingTorqueInBodyFixedFrame( 2 );
    const ScalarType x2 = x * x;
    const ScalarType y2 = y * y;
    const ScalarType z2 = z * z;
    const ScalarType xy = x * y;
    const ScalarType xz = x * z;
    const ScalarType yz = y * z;

    const ScalarType r2 = x2 + y2 + z2;
    const ScalarType inverseR2 = ScalarType( 1.0 ) / r2;
    using std::sqrt;
    const ScalarType r5 = r2 * r2 * sqrt( r2 );

    const ScalarType A = inertiaTensorOfBodyUndergoingTorque( 0, 0 );
    const ScalarType B = inertiaTensorOfBodyUndergoingTorque( 1, 1 );
    const ScalarType C = inertiaTensorOfBodyUndergoingTorque( 2, 2 );
    const ScalarType Ixy = -inertiaTensorOfBodyUndergoingTorque( 0, 1 );
    const ScalarType Ixz = -inertiaTensorOfBodyUndergoingTorque( 0, 2 );
    const ScalarType Iyz = -inertiaTensorOfBodyUndergoingTorque( 1, 2 );

    const ScalarType Aprime = inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( 0, 0 );
    const ScalarType Bprime = inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( 1, 1 );
    const ScalarType Cprime = inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( 2, 2 );
    const ScalarType IxyPrime = -inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( 0, 1 );
    const ScalarType IxzPrime = -inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( 0, 2 );
    const ScalarType IyzPrime = -inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( 1, 2 );

    const ScalarType Qprime = Aprime + Bprime + Cprime;
    const ScalarType IellPrime =
            ( Aprime * x2 + Bprime * y2 + Cprime * z2 - ScalarType( 2.0 ) * IxyPrime * xy -
              ScalarType( 2.0 ) * IxzPrime * xz - ScalarType( 2.0 ) * IyzPrime * yz ) *
            inverseR2;
    const ScalarType Wprime =
            massOfBodyExertingTorque + ScalarType( 7.5 ) * Qprime * inverseR2 - ScalarType( 17.5 ) * IellPrime * inverseR2;

    const ScalarType fyz = yz * ( Wprime - ScalarType( 5.0 ) * Aprime * inverseR2 ) -
            ScalarType( 5.0 ) * IxzPrime * xy * inverseR2 - ScalarType( 5.0 ) * IxyPrime * xz * inverseR2 +
            IyzPrime * ( ScalarType( 1.0 ) - ScalarType( 5.0 ) * ( y2 + z2 ) * inverseR2 );
    const ScalarType fxz = xz * ( Wprime - ScalarType( 5.0 ) * Bprime * inverseR2 ) +
            IxzPrime * ( ScalarType( 1.0 ) - ScalarType( 5.0 ) * ( x2 + z2 ) * inverseR2 ) -
            ScalarType( 5.0 ) * IyzPrime * xy * inverseR2 - ScalarType( 5.0 ) * IxyPrime * yz * inverseR2;
    const ScalarType fxy = xy * ( Wprime - ScalarType( 5.0 ) * Cprime * inverseR2 ) -
            ScalarType( 5.0 ) * IyzPrime * xz * inverseR2 +
            IxyPrime * ( ScalarType( 1.0 ) - ScalarType( 5.0 ) * ( x2 + y2 ) * inverseR2 ) -
            ScalarType( 5.0 ) * IxzPrime * yz * inverseR2;

    const ScalarType gyz = ( z2 - y2 ) * Wprime + Bprime - Cprime -
            ScalarType( 10.0 ) * IxzPrime * xz * inverseR2 - ScalarType( 10.0 ) * IxyPrime * xy * inverseR2 -
            ScalarType( 20.0 ) * IyzPrime * yz * inverseR2 -
            ScalarType( 5.0 ) * z2 * ( Aprime + Bprime - Cprime ) * inverseR2 -
            ScalarType( 5.0 ) * y2 * ( Aprime - Bprime + Cprime ) * inverseR2;
    const ScalarType gxz = ( x2 - z2 ) * Wprime + Cprime - Aprime -
            ScalarType( 20.0 ) * IxzPrime * xz * inverseR2 - ScalarType( 10.0 ) * IxyPrime * xy * inverseR2 -
            ScalarType( 10.0 ) * IyzPrime * yz * inverseR2 -
            ScalarType( 5.0 ) * x2 * ( -Aprime + Bprime + Cprime ) * inverseR2 -
            ScalarType( 5.0 ) * z2 * ( Aprime + Bprime - Cprime ) * inverseR2;
    const ScalarType gxy = ( y2 - x2 ) * Wprime + Aprime - Bprime -
            ScalarType( 10.0 ) * IxzPrime * xz * inverseR2 - ScalarType( 20.0 ) * IxyPrime * xy * inverseR2 -
            ScalarType( 10.0 ) * IyzPrime * yz * inverseR2 -
            ScalarType( 5.0 ) * y2 * ( Aprime - Bprime + Cprime ) * inverseR2 -
            ScalarType( 5.0 ) * x2 * ( -Aprime + Bprime + Cprime ) * inverseR2;

    const ScalarType prefactor = ScalarType( 3.0 * physical_constants::GRAVITATIONAL_CONSTANT ) / r5;
    Eigen::Matrix< ScalarType, 3, 1 > torque;
    torque( 0 ) = prefactor * ( ( C - B ) * fyz - Ixz * fxy + Ixy * fxz + Iyz * gyz );
    torque( 1 ) = prefactor * ( ( A - C ) * fxz + Ixz * gxz - Ixy * fyz + Iyz * fxy );
    torque( 2 ) = prefactor * ( ( B - A ) * fxy + Ixz * fyz + Ixy * gxy - Iyz * fxz );
    return torque;
}

//! Compute the fourth-degree full two-body gravitational torque from Schutz (1981), Eq. (11).
/*!
 * Computes the torque on body 1 (the body undergoing torque), in body-1-fixed coordinates.
 * The input inertia tensor of body 2 must be expressed in body-1-fixed coordinates.
 * The torque is evaluated directly from tensor components as given by the fourth-degree two-body formulation.
 * \param relativePositionOfBodyExertingTorqueInBodyFixedFrame Position of body 2 w.r.t. body 1, in body-1-fixed frame.
 * \param massOfBodyExertingTorque Mass of body 2.
 * \param inertiaTensorOfBodyUndergoingTorque Inertia tensor of body 1, in body-1-fixed frame.
 * \param inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque Inertia tensor of body 2 in body-1-fixed frame.
 * \return Fourth-degree two-body gravitational torque on body 1.
 */
Eigen::Vector3d calculateFourthDegreeFullTwoBodyGravitationalTorque(
        const Eigen::Vector3d& relativePositionOfBodyExertingTorqueInBodyFixedFrame,
        const double massOfBodyExertingTorque,
        const Eigen::Matrix3d& inertiaTensorOfBodyUndergoingTorque,
        const Eigen::Matrix3d& inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque );

//! Torque model implementing the direct fourth-degree two-body tensor formulation for two finite-sized rigid bodies.
class FourthDegreeFullTwoBodyGravitationalTorqueModel : public basic_astrodynamics::TorqueModel
{
public:
    //! Constructor.
    FourthDegreeFullTwoBodyGravitationalTorqueModel(
            const std::function< Eigen::Vector3d( ) > positionOfBodyUndergoingTorqueFunction,
            const std::function< Eigen::Vector3d( ) > positionOfBodyExertingTorqueFunction,
            const std::function< double( ) > massOfBodyExertingTorqueFunction,
            const std::function< Eigen::Matrix3d( ) > inertiaTensorOfBodyUndergoingTorqueFunction,
            const std::function< Eigen::Matrix3d( ) > inertiaTensorOfBodyExertingTorqueFunction,
            const std::function< Eigen::Quaterniond( ) > rotationToBodyFixedFrameOfBodyUndergoingTorqueFunction,
            const std::function< Eigen::Quaterniond( ) > rotationToBodyFixedFrameOfBodyExertingTorqueFunction );

    //! Get currently computed torque.
    Eigen::Vector3d getTorque( )
    {
        return currentTorque_;
    }

    //! Update member variables used by the torque model.
    void updateMembers( const double currentTime = TUDAT_NAN );

    //! Get current body-2 inertia tensor expressed in body-1-fixed frame.
    Eigen::Matrix3d getCurrentInertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( ) const
    {
        return currentInertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque_;
    }

    Eigen::Quaterniond getCurrentRotationToBodyFixedFrameOfBodyUndergoingTorque( ) const
    {
        return currentRotationToBodyFixedFrameOfBodyUndergoingTorque_;
    }

    Eigen::Quaterniond getCurrentRotationToBodyFixedFrameOfBodyExertingTorque( ) const
    {
        return currentRotationToBodyFixedFrameOfBodyExertingTorque_;
    }

    Eigen::Vector3d getCurrentRelativePositionInInertialFrame( ) const
    {
        return currentRelativePositionInInertialFrame_;
    }

    Eigen::Vector3d getCurrentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque( ) const
    {
        return currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_;
    }

    double getCurrentMassOfBodyExertingTorque( ) const
    {
        return currentMassOfBodyExertingTorque_;
    }

    Eigen::Matrix3d getCurrentInertiaTensorOfBodyUndergoingTorque( ) const
    {
        return currentInertiaTensorOfBodyUndergoingTorque_;
    }

    Eigen::Matrix3d getCurrentInertiaTensorOfBodyExertingTorque( ) const
    {
        return currentInertiaTensorOfBodyExertingTorque_;
    }

private:
    std::function< Eigen::Vector3d( ) > positionOfBodyUndergoingTorqueFunction_;
    std::function< Eigen::Vector3d( ) > positionOfBodyExertingTorqueFunction_;
    std::function< double( ) > massOfBodyExertingTorqueFunction_;
    std::function< Eigen::Matrix3d( ) > inertiaTensorOfBodyUndergoingTorqueFunction_;
    std::function< Eigen::Matrix3d( ) > inertiaTensorOfBodyExertingTorqueFunction_;
    std::function< Eigen::Quaterniond( ) > rotationToBodyFixedFrameOfBodyUndergoingTorqueFunction_;
    std::function< Eigen::Quaterniond( ) > rotationToBodyFixedFrameOfBodyExertingTorqueFunction_;

    Eigen::Quaterniond currentRotationToBodyFixedFrameOfBodyUndergoingTorque_;
    Eigen::Quaterniond currentRotationToBodyFixedFrameOfBodyExertingTorque_;
    Eigen::Vector3d currentRelativePositionInInertialFrame_;
    Eigen::Vector3d currentRelativePositionInBodyFixedFrameOfBodyUndergoingTorque_;
    double currentMassOfBodyExertingTorque_ = TUDAT_NAN;
    Eigen::Matrix3d currentInertiaTensorOfBodyUndergoingTorque_ = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d currentInertiaTensorOfBodyExertingTorque_ = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d currentInertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque_ = Eigen::Matrix3d::Zero( );
    Eigen::Vector3d currentTorque_ = Eigen::Vector3d::Zero( );
};

}  // namespace gravitation

}  // namespace tudat

#endif  // TUDAT_FOURTHDEGREEFULLTWOBODYGRAVITATIONALTORQUE_H
