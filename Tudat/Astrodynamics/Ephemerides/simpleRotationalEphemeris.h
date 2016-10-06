/*   Copyright (c) 2010-2015, Delft University of Technology
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or without modification, are
 *   permitted provided that the following conditions are met:
 *     - Redistributions of source code must retain the above copyright notice, this list of
 *       conditions and the following disclaimer.
 *     - Redistributions in binary form must reproduce the above copyright notice, this list of
 *       conditions and the following disclaimer in the documentation and/or other materials
 *       provided with the distribution.
 *     - Neither the name of the Delft University of Technology nor the names of its contributors
 *       may be used to endorse or promote products derived from this software without specific
 *       prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *   OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *   COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *   GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *   OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *   Changelog
 *     YYMMDD    Author            Comment
 *     130219    D. Dirkx          Migrated from personal code.
 *     130227    R.C.A. Boon       Changed include guard, changed header indentation, minor
 *                                 textual changes, improved commenting.
 *
 *   References
 *
 *   Notes
 *
 */

#ifndef TUDAT_SIMPLE_ROTATIONAL_EPHEMERIS_H
#define TUDAT_SIMPLE_ROTATIONAL_EPHEMERIS_H

#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

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
        auxiliaryMatrix_<< 0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0;

        initialEulerAngles_ = reference_frames::calculateInertialToPlanetFixedRotationAnglesFromMatrix(
                    Eigen::Matrix3d( initialRotationToTargetFrame_ ) );
    }

    //! Constructor from rotation state angles.
    /*!
     * Constructor from rotation state angles, i.e., right ascension and declination of pole
     * (both assumuned constant) and the position of the prime meridian at a given epoch
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
        auxiliaryMatrix_<< 0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0;

        // Set vector of euler angles
        initialEulerAngles_.x( ) = poleRightAscension;
        initialEulerAngles_.y( ) = poleDeclination;
        initialEulerAngles_.z( ) = primeMeridianOfDate;
    }

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
\
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

    //! Auxiliary matrix used to calculate the time derivative of a rotation matrix.
    Eigen::Matrix3d auxiliaryMatrix_;
};

} // namespace tudat
} // namespace ephemerides

#endif // TUDAT_SIMPLE_ROTATIONAL_EPHEMERIS_H
