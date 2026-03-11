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

#include <functional>

#include <Eigen/Geometry>

#include "tudat/astro/basic_astro/torqueModel.h"

namespace tudat
{

namespace gravitation
{

//! Compute the fourth-degree full two-body gravitational torque from Schutz (1981), Eq. (11).
/*!
 * Computes the torque on body 1 (the body undergoing torque), in body-1-fixed coordinates.
 * The input inertia tensor of body 2 must be expressed in body-1-fixed coordinates.
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

//! Torque model implementing Schutz (1981), Eq. (11), for two finite-sized rigid bodies.
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
