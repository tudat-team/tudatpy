/*    Copyright (c) 2010-2015, Delft University of Technology
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
   *     130227    R.C.A. Boon       Changed include guard, improved commenting.
   *
   *   References
   *
   *   Notes
   *
   */

#ifndef TUDAT_ROTATIONAL_EPHEMERIS_H
#define TUDAT_ROTATIONAL_EPHEMERIS_H

#include <string>

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace tudat
{
namespace ephemerides
{

//! Base class for rotational ephemerides of bodies
/*!
 * Base class for rotational ephemerides of bodies. The rotation (quaternion) between two frames
 * specified by member variable ids can be calculated as a function of time in the manner
 * determined by the derived class.
 */
class RotationalEphemeris
{
public:

    //! Constructor.
    /*!
     * Constructor, sets frames between which rotation is determined.
     * \param baseFrameOrientation Base frame identifier.
     * \param targetFrameOrientation Target frame identifier.
     */
    RotationalEphemeris( const std::string& baseFrameOrientation = "",
                         const std::string& targetFrameOrientation = "" )
        : baseFrameOrientation_( baseFrameOrientation ),
          targetFrameOrientation_( targetFrameOrientation )
    { }

    //! Virtual destructor.
    /*!
     * Virtual destructor.
     */
    virtual ~RotationalEphemeris( ) { }

    //! Get rotation quaternion from target frame to base frame.
    /*!
     * Pure virtual function to calculate and return the rotation quaternion from target frame to
     * base frame at specified time.
     * \param secondsSinceEpoch Seconds since Julian day epoch specified by 2nd argument
     * \param julianDayAtEpoch Reference epoch in Julian days from which number of seconds
     *          are counted.
     * \return Rotation quaternion computed.
     */
    virtual Eigen::Quaterniond getRotationToBaseFrame(
            const double secondsSinceEpoch, const double julianDayAtEpoch ) = 0;

    //! Get rotation quaternion to target frame from base frame.
    /*!
     * Pure virtual function to calculate and return the rotation quaternion to target frame from
     * base frame at specified time.
     * \param secondsSinceEpoch Seconds since Julian day epoch specified by 2nd argument
     * \param julianDayAtEpoch Reference epoch in Julian days from which number of seconds are
     *          counted.
     * \return Rotation quaternion computed.
     */
    virtual Eigen::Quaterniond getRotationToTargetFrame(
            const double secondsSinceEpoch, const double julianDayAtEpoch ) = 0;

    //! Get base reference frame orientation.
    /*!
     * Function to retrieve the base reference frame orientation.
     * \return Base reference frame orientation.
     */
    std::string getBaseFrameOrientation( ) { return baseFrameOrientation_; }

    //! Get target reference frame orientation.
    /*!
     * Function to retrieve the target reference frame orientation.
     * \return Target reference frame orientation.
     */
    std::string getTargetFrameOrientation( ) { return targetFrameOrientation_; }

protected:

    //! Base reference frame orientation.
    /*!
     * Base reference frame orientation.
     */
    const std::string baseFrameOrientation_;

    //! Target reference frame orientation.
    /*!
     * Target reference frame orientation.
     */
    const std::string targetFrameOrientation_;

};

} // namespace tudat
} // namespace ephemerides

#endif // TUDAT_ROTATIONAL_EPHEMERIS_H
