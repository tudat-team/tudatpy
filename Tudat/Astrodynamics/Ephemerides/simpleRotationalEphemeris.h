/*   Copyright (c) 2010-2013, Delft University of Technology
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
     * \param initialSecondsSinceEpoch Seconds since epoch at which initialRotationToTargetFrame
     *          is valid. Epoch is given by next parameter.
     * \param inputReferenceJulianDay Julian day of epoch since which initialSecondsSinceEpoch is
     *          counted.
     * \param baseFrameOrientation Base frame identifier.
     * \param targetFrameOrientation Target frame identifier.
     */
    SimpleRotationalEphemeris( const Eigen::Quaterniond& initialRotationToTargetFrame,
                               const double rotationRate,
                               const double initialSecondsSinceEpoch,
                               const double inputReferenceJulianDay,
                               const std::string& baseFrameOrientation = "",
                               const std::string& targetFrameOrientation = "" )
        : RotationalEphemeris( baseFrameOrientation, targetFrameOrientation ),
          rotationRate_( rotationRate ),
          initialRotationToTargetFrame_( initialRotationToTargetFrame ),
          initialSecondsSinceEpoch_( initialSecondsSinceEpoch ),
          inputReferenceJulianDay_( inputReferenceJulianDay )
    { }

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
     * \param initialSecondsSinceEpoch Seconds since epoch at which primeMeridianOfDate
     *          is valid. Epoch is gievn by next parameter.
     * \param inputReferenceJulianDay Julian day of epoch since which previous parameter is
     *          counted.
     * \param baseFrameOrientation Base frame identifier.
     * \param targetFrameOrientation Target frame identifier.
     */
    SimpleRotationalEphemeris( const double poleRightAscension,
                               const double poleDeclination,
                               const double primeMeridianOfDate,
                               const double rotationRate,
                               const double initialSecondsSinceEpoch,
                               const double inputReferenceJulianDay,
                               const std::string& baseFrameOrientation = "",
                               const std::string& targetFrameOrientation = ""  )
        : RotationalEphemeris( baseFrameOrientation, targetFrameOrientation ),
          rotationRate_( rotationRate ),
          initialRotationToTargetFrame_(
              reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion(
                  poleDeclination, poleRightAscension, primeMeridianOfDate ) ),
          initialSecondsSinceEpoch_( initialSecondsSinceEpoch ),
          inputReferenceJulianDay_( inputReferenceJulianDay )
    { }

    //! Calculate rotation quaternion from target frame to base frame.
    /*!
     * Pure virtual function that calculates the rotation quaternion from target frame to base
     * frame at specified time.
     * \param secondsSinceEpoch Seconds since Julian day epoch specified by 2nd argument.
     * \param julianDayAtEpoch Reference epoch in Julian days from which number of seconds are
     *          counted.
     * \return Rotation quaternion computed.
     */
    Eigen::Quaterniond getRotationToBaseFrame(
            const double secondsSinceEpoch, const double julianDayAtEpoch )
    {
        return getRotationToTargetFrame( secondsSinceEpoch, julianDayAtEpoch ).inverse( );
    }

    //! Get rotation quaternion to target frame from base frame.
    /*!
     * Returns the rotation quaternion to target frame from base frame at specified time.
     * \param secondsSinceEpoch Seconds since Julian day epoch specified by 2nd argument
     * \param julianDayAtEpoch Reference epoch in Julian days from which number of seconds are
     *          counted.
     * \return Rotation quaternion computed.
     */
    virtual Eigen::Quaterniond getRotationToTargetFrame(
            const double secondsSinceEpoch, const double julianDayAtEpoch );

    //! Get rotation from target to base frame at initial time.
    /*!
     * Returns rotation from target to base frame at initial time.
     * \return Rotation quaternion at initial time.
     */
    Eigen::Quaterniond getInitialRotationalState( ) { return initialRotationToTargetFrame_; }

    //! Get seconds since epoch at which initialRotationToTargetFrame_ is valid.
    /*!
     * Returns seconds since epoch at which initialRotationToTargetFrame_ is valid.
     * \return Seconds since Julian day epoch [s].
     */
    double getInitialSecondsSinceEpoch( ) { return initialSecondsSinceEpoch_; }

    //! Get Julian day of reference epoch.
    /*!
     * Returns Julian day of reference epoch.
     * \return Julian day of reference epoch.
     */
    double getInputReferenceJulianDay( ) { return inputReferenceJulianDay_; }

    //! Get rotation rate of body.
    /*!
     * Returns rotation rate of body.
     * \return Rotation rate [rad/s].
     */
    double getRotationRate( ) { return rotationRate_; }

private:

    //! Rotation rate of body (about local z-axis).
    /*!
     * Rotation rate of body (about local z-axis).
     */
     const double rotationRate_;

    //! Rotation from target to base frame at initial time.
    /*!
     * Rotation from target to base frame at initial time.
     */
    const Eigen::Quaterniond initialRotationToTargetFrame_;

    //! Seconds since epoch at which initialRotationToTargetFrame_ is valid.
    /*!
     * Seconds since epoch at which initialRotationToTargetFrame is valid.
     */
    const double initialSecondsSinceEpoch_;

    //! Julian day of reference epoch.
    /*!
     * Julian day of reference epoch.
     */
   const  double inputReferenceJulianDay_;
};

} // namespace tudat
} // namespace ephemerides

#endif // TUDAT_SIMPLE_ROTATIONAL_EPHEMERIS_H
