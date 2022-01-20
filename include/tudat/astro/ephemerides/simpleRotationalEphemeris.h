/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_SIMPLE_ROTATIONAL_EPHEMERIS_H
#define TUDAT_SIMPLE_ROTATIONAL_EPHEMERIS_H

#include "tudat/astro/ephemerides/rotationalEphemeris.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"

namespace tudat
{
namespace ephemerides
{

//! Class to model the rotation of bodies with constant rotation axis and rate.
/*!
 * Class to model the rotation of bodies with constant rotation axis and rate, i.e., those for
 * which the angular momentum vector remains constant in time (according to an inertial observer).
 */
class SimpleRotationalEphemeris : public RotationalEphemeris
{
public:

    //! Constructor from initial rotational state.
    /*!
     * Constructor from initial rotational state that is empty because all member variables are
     * initialized through the constructor arguments.
     * \param initialRotationToTargetFrame Rotation from base to target frame at initial time
     *          (specified by 2nd and 3rd parameter).
     * \param rotationRate Constant rotation rate of body.
     * \param initialSecondsSinceEpoch Seconds since epoch at which primeMeridianOfDate  is valid.
     * \param baseFrameOrientation Base frame identifier.
     * \param targetFrameOrientation Target frame identifier.
     */
    SimpleRotationalEphemeris( const Eigen::Quaterniond& initialRotationToTargetFrame,
                               const double rotationRate,
                               const double initialSecondsSinceEpoch,
                               const std::string& baseFrameOrientation = "",
                               const std::string& targetFrameOrientation = "" )
        : RotationalEphemeris( baseFrameOrientation, targetFrameOrientation ),
          rotationRate_( rotationRate ),
          initialRotationToTargetFrame_( initialRotationToTargetFrame ),
          initialSecondsSinceEpoch_( initialSecondsSinceEpoch )
    {
        initialEulerAngles_ = reference_frames::calculateInertialToPlanetFixedRotationAnglesFromMatrix(
                    Eigen::Matrix3d( initialRotationToTargetFrame_ ) );
    }

    //! Constructor from rotation state angles.
    /*!
     * Constructor from rotation state angles, i.e., right ascension and declination of pole
     * (both assumed constant) and the position of the prime meridian at a given epoch
     * (reference Julian day and seconds since that Julian day). Calculates member variables from
     * constructor arguments.
     * \param poleRightAscension Right ascension of body's pole in base frame.
     * \param poleDeclination Declination of body's pole in base frame.
     * \param primeMeridianOfDate Position of prime meridian at given time.
     * \param rotationRate Constant rotation rate of body.
     * \param initialSecondsSinceEpoch Seconds since epoch at which primeMeridianOfDate  is valid.
     * \param baseFrameOrientation Base frame identifier.
     * \param targetFrameOrientation Target frame identifier.
     */
    SimpleRotationalEphemeris( const double poleRightAscension,
                               const double poleDeclination,
                               const double primeMeridianOfDate,
                               const double rotationRate,
                               const double initialSecondsSinceEpoch,
                               const std::string& baseFrameOrientation = "",
                               const std::string& targetFrameOrientation = ""  )
        : RotationalEphemeris( baseFrameOrientation, targetFrameOrientation ),
          rotationRate_( rotationRate ),
          initialRotationToTargetFrame_(
              reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion(
                  poleDeclination, poleRightAscension, primeMeridianOfDate ) ),
          initialSecondsSinceEpoch_( initialSecondsSinceEpoch )
    {
        // Set vector of euler angles
        initialEulerAngles_.x( ) = poleRightAscension;
        initialEulerAngles_.y( ) = poleDeclination;
        initialEulerAngles_.z( ) = primeMeridianOfDate;
    }

    //! Destructor
    ~SimpleRotationalEphemeris( ){ }

    //! Calculate rotation quaternion from target frame to base frame.
    /*!
     * Pure virtual function that calculates the rotation quaternion from target frame to base
     * frame at specified time.
     * \param secondsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated.
     * \return Rotation quaternion computed.
     */
    Eigen::Quaterniond getRotationToBaseFrame(
            const double secondsSinceEpoch )
    {
        return getRotationToTargetFrame( secondsSinceEpoch ).inverse( );
    }

    //! Get rotation quaternion to target frame from base frame.
    /*!
     * Returns the rotation quaternion to target frame from base frame at specified time.
     * \param secondsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated.
     * \return Rotation quaternion computed.
     */
    Eigen::Quaterniond getRotationToTargetFrame(
            const double secondsSinceEpoch );

    //! Function to calculate the derivative of the rotation matrix from target frame to original
    //! frame.
    /*!
     *  Function to calculate the derivative of the rotation matrix from target frame to original
     *  frame at specified time.
     *  \param secondsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated.
     *  \return Derivative of rotation from target (typically local) to original (typically global)
     *          frame at specified time.
     */
    Eigen::Matrix3d getDerivativeOfRotationToBaseFrame(
            const double secondsSinceEpoch )
    {
        return getDerivativeOfRotationToTargetFrame( secondsSinceEpoch ).
                transpose( );
    }

    //! Function to calculate the derivative of the rotation matrix from original frame to target
    //! frame.
    /*!
     *  Function to calculate the derivative of the rotation matrix from original frame to target
     *  frame at specified time.
     * \param secondsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated.
     *  \return Derivative of rotation from original (typically global) to target (typically local)
     *          frame at specified time.
     */
    Eigen::Matrix3d getDerivativeOfRotationToTargetFrame(
            const double secondsSinceEpoch );

    //! Get rotation from target to base frame at initial time.
    /*!
     * Returns rotation from target to base frame at initial time.
     * \return Rotation quaternion at initial time.
     */
    Eigen::Quaterniond getInitialRotationToTargetFrame( ) { return initialRotationToTargetFrame_; }

    //! Get seconds since epoch at which initialRotationToTargetFrame_ is valid.
    /*!
     * Returns seconds since epoch at which initialRotationToTargetFrame_ is valid.
     * \return Seconds since Julian day epoch [s].
     */
    double getInitialSecondsSinceEpoch( ) { return initialSecondsSinceEpoch_; }

    //! Get rotation rate of body.
    /*!
     * Returns rotation rate of body.
     * \return Rotation rate [rad/s].
     */
    double getRotationRate( ) { return rotationRate_; }

    //! Function to reset the rotation rate of the body.
    /*!
     * Function to reset the rotation rate of the body.
     * \param rotationRate New rotation rate [rad/s].
     */
    void resetRotationRate( const double rotationRate ) { rotationRate_ = rotationRate; }

    //! Function to get vector of euler angles at initialSecondsSinceEpoch_
    /*!
     *  Function to get vector of euler angles at initialSecondsSinceEpoch_, in order right ascension, declination,
     *  prime meridian.
     *  \return Vector of euler angles at initialSecondsSinceEpoch_, in order right ascension, declination, prime meridian.
     */
    Eigen::Vector3d getInitialEulerAngles( )
    {
        return initialEulerAngles_;
    }

    //! Function to reset the right ascension and declination of body's north pole.
    /*!
     *  Function to reset the right ascension and declination of body's north pole,
     *  recalculates the initialRotationToOriginalFrame_ member.
     *  \param rightAscension New right ascension of north pole.
     *  \param declination New declination of north pole.
     */
    void resetInitialPoleRightAscensionAndDeclination( const double rightAscension,
                                                       const double declination );


private:

    //! Rotation rate of body (about local z-axis).
    /*!
     * Rotation rate of body (about local z-axis).
     */
    double rotationRate_;

    //! Rotation from target to base frame at initial time.
    /*!
     * Rotation from target to base frame at initial time.
     */
    Eigen::Quaterniond initialRotationToTargetFrame_;

    //! Seconds since epoch at which initialRotationToTargetFrame_ is valid.
    /*!
     * Seconds since epoch at which initialRotationToTargetFrame is valid.
     */
    double initialSecondsSinceEpoch_;


    //! Initial Euler angles describing the rotational state of the body at initialSecondsSinceEpoch_
    /*!
     *  Initial Euler angles describing the rotational state of the body at initialSecondsSinceEpoch_. Order of the vector
     *  is: right ascension (alpha), declination (delta), prime meridian of date (W)
     */
    Eigen::Vector3d initialEulerAngles_;
};

} // namespace ephemerides
} // namespace tudat

#endif // TUDAT_SIMPLE_ROTATIONAL_EPHEMERIS_H
