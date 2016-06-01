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

#include <boost/function.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

namespace tudat
{
namespace ephemerides
{


//! Function to calculate the rotational velocity vector of frame B w.r.t frame A.
/*!
 *  Function to calculate the rotational velocity vector of frame B (local) w.r.t frame A (global)
 *  from the rotation matrix between the frames, as well as its time derivative.
 *  \param rotationToTargetFrame Rotation matrix from frame A to frame B.
 *  \param rotationMatrixToGlobalFrameDerivative Time derivative if rotation matrix from
 *  frame B to frame A.
 *  \return Angular velocity vector of frame B, expressed in frame A.
 */
Eigen::Vector3d getRotationalVelocityVectorInBaseFrameFromMatrices(
        const Eigen::Matrix3d& rotationToTargetFrame,
        const Eigen::Matrix3d& rotationMatrixToGlobalFrameDerivative );

//! Function to calculate the time derivative of rotation matrix from frame A to frame B.
/*!
 *  Function to calculate the time derivative of rotation matrix from frame A (global) to frame B
 *  (local) from the rotation matrix between the frames, as well as the angular velocity
 *  vector of frame B w.r.t. frame A.
 *  \param rotationToTargetFrame Rotation matrix from frame A to frame B.
 *  \param rotationalVelocityVectorInTargetFrame Angular velocity vector of frame B,
 *  expressed in frame A.
 *  \return Time derivative if rotation matrix from frame A to frame B.
 */
Eigen::Matrix3d getDerivativeOfRotationMatrixToFrame(
        const Eigen::Matrix3d& rotationToTargetFrame,
        const Eigen::Vector3d& rotationalVelocityVectorInTargetFrame );

//! Transform a state (Cartesian position and velocity) from one frame to another.
/*!
 *  Transform a state (Cartesian position and velocity) from one frame to another, taking into
 *  account both the instantaneous rotational state of the two frames, and the rotational
 *  rate of one frame w.r.t. the other.
 *  \param stateInBaseFrame State that is to be transformed from base to target frame.
 *  \param rotationToFrame Rotation from base to target frame.
 *  \param rotationMatrixToFrameDerivative Time derivative of rotation matrix from base to target
 *  frame.
 *  \return State (Cartesian position and velocity) in target frame.
 */
template< typename StateScalarType >
Eigen::Matrix< StateScalarType, 6, 1 > transformStateToFrame(
        const Eigen::Matrix< StateScalarType, 6, 1 >& stateInBaseFrame,
        const Eigen::Quaterniond& rotationToFrame,
        const Eigen::Matrix3d& rotationMatrixToFrameDerivative )
{
    Eigen::Matrix< StateScalarType, 6, 1 >stateInTargetFrame;
    stateInTargetFrame.segment( 0, 3 ) =
            rotationToFrame.template cast< StateScalarType >( ) * stateInBaseFrame.segment( 0, 3 );
    stateInTargetFrame.segment( 3, 3 ) =
            rotationMatrixToFrameDerivative.template cast< StateScalarType >( ) * stateInBaseFrame.segment( 0, 3 ) +
            rotationToFrame.template cast< StateScalarType >( ) * stateInBaseFrame.segment( 3, 3 );
    return stateInTargetFrame;

}

//! Transform a state (Cartesian position and velocity) from one frame to another.
/*!
 *  Transform a state (Cartesian position and velocity) from one frame to another, taking into
 *  account both the instantaneous rotational state of the two frames, and the rotational
 *  rate of one frame w.r.t. the other.
 *  \param stateInBaseFrame State that is to be transformed from base to target frame.
 *  \param rotationToFrameFunction Function returning rotation from base to target frame.
 *  \param rotationMatrixToFrameDerivativeFunction Function returning time derivative of rotation
 *   matrix from base to target frame.
 *  \return State (Cartesian position and velocity) in target frame.
 */
template< typename StateScalarType >
Eigen::Matrix< StateScalarType, 6, 1 > transformStateToFrame(
        const Eigen::Matrix< StateScalarType, 6, 1 >& stateInBaseFrame,
        const boost::function< Eigen::Quaterniond( ) > rotationToFrameFunction,
        const boost::function< Eigen::Matrix3d( ) > rotationMatrixToFrameDerivativeFunction )
{
    return transformStateToFrame(
                stateInBaseFrame, rotationToFrameFunction( ),
                rotationMatrixToFrameDerivativeFunction( ) );
}

//! Transform a state (Cartesian position and velocity) from one frame to another.
/*!
 *  Transform a state (Cartesian position and velocity) from one frame to another, taking into
 *  account both the instantaneous rotational state of the two frames, and the rotational
 *  rate of one frame w.r.t. the other.
 *  \param stateInBaseFrame State that is to be transformed from base to target frame.
 *  \param currentTime Input time for rotationToFrameFunction and rotationMatrixToFrameDerivativeFunction
 *  \param rotationToFrameFunction Function returning rotation from base to target frame.
 *  \param rotationMatrixToFrameDerivativeFunction Function returning time derivative of rotation
 *   matrix from base to target frame.
 *  \return State (Cartesian position and velocity) in target frame.
 */
template< typename StateScalarType >
Eigen::Matrix< StateScalarType, 6, 1 > transformStateToFrameFromStateFunctions(
        const Eigen::Matrix< StateScalarType, 6, 1 >& stateInBaseFrame,
        const double currentTime,
        const boost::function< Eigen::Quaterniond( const double ) > rotationToFrameFunction,
        const boost::function< Eigen::Matrix3d( const double ) > rotationMatrixToFrameDerivativeFunction )
{
    return transformStateToFrame(
                stateInBaseFrame, rotationToFrameFunction( currentTime ),
                rotationMatrixToFrameDerivativeFunction( currentTime ) );
}

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
     *          are counted (ddefault J2000).
     * \return Rotation quaternion computed.
     */
    virtual Eigen::Quaterniond getRotationToBaseFrame(
            const double secondsSinceEpoch,
            const double julianDayAtEpoch = basic_astrodynamics::JULIAN_DAY_ON_J2000 ) = 0;

    //! Get rotation quaternion to target frame from base frame.
    /*!
     * Pure virtual function to calculate and return the rotation quaternion to target frame from
     * base frame at specified time.
     * \param secondsSinceEpoch Seconds since Julian day epoch specified by 2nd argument
     * \param julianDayAtEpoch Reference epoch in Julian days from which number of seconds are
     *          counted (ddefault J2000).
     * \return Rotation quaternion computed.
     */
    virtual Eigen::Quaterniond getRotationToTargetFrame(
            const double secondsSinceEpoch,
            const double julianDayAtEpoch = basic_astrodynamics::JULIAN_DAY_ON_J2000 ) = 0;

    //! Function to calculate the derivative of the rotation matrix from target frame to original
    //! frame.
    /*!
     *  Function to calculate the derivative of the rotation matrix from target frame to original
     *  frame at specified time, to be implemented by derived class.
     *  \param secondsSinceEpoch Seconds since Julian day epoch specified by 2nd argument
     *  \param julianDayAtEpoch Reference epoch in Julian days from which number of seconds are
     *          counted.
     *  \return Derivative of rotation from target (typically local) to original (typically global)
     *          frame at specified time.
     */
    virtual Eigen::Matrix3d getDerivativeOfRotationToBaseFrame(
            const double secondsSinceEpoch, const double julianDayAtEpoch =
            basic_astrodynamics::JULIAN_DAY_ON_J2000 ) = 0;

    //! Function to calculate the derivative of the rotation matrix from original frame to target
    //! frame.
    /*!
     *  Function to calculate the derivative of the rotation matrix from original frame to target
     *  frame at specified time, to be implemented by derived class.
     * \param secondsSinceEpoch Seconds since Julian day epoch specified by 2nd argument
     * \param julianDayAtEpoch Reference epoch in Julian days from which number of seconds are
     *          counted (default J2000).
     *  \return Derivative of rotation from original (typically global) to target (typically local)
     *          frame at specified time.
     */
    virtual Eigen::Matrix3d getDerivativeOfRotationToTargetFrame(
            const double secondsSinceEpoch, const double julianDayAtEpoch =
            basic_astrodynamics::JULIAN_DAY_ON_J2000 ) = 0;

    //! Function to retrieve the angular velocity vector, expressed in base frame.
    /*!
     * Function to retrieve the angular velocity vector, expressed in base frame. Thsi function uses the functions
     * that calculate the rotation matrix and its time derivative. It may be redefined in a derived class, to
     * calculate the angular velocity vector directly
     * \param secondsSinceEpoch Seconds since Julian day epoch specified by 2nd argument
     * \param julianDayAtEpoch Reference epoch in Julian days from which number of seconds are
     * counted (default J2000).
     * \return Angular velocity vector, expressed in base frame.
     */
    virtual Eigen::Vector3d getRotationalVelocityVectorInBaseFrame(
            const double secondsSinceEpoch, const double julianDayAtEpoch =
            basic_astrodynamics::JULIAN_DAY_ON_J2000 )
    {
        return getRotationalVelocityVectorInBaseFrameFromMatrices(
                    Eigen::Matrix3d( getRotationToTargetFrame( secondsSinceEpoch, julianDayAtEpoch ) ),
                    getDerivativeOfRotationToTargetFrame( secondsSinceEpoch, julianDayAtEpoch ) );
    }

    //! Function to retrieve the angular velocity vector, expressed in target frame.
    /*!
     * Function to retrieve the angular velocity vector, expressed in target frame. This function calls the function
     * that calculates the angular velocioty vector in the base frame, and rotates it to teh target frame.
     * It may be redefined in a derived class, to calculate the angular velocity vector directly
     * \param secondsSinceEpoch Seconds since Julian day epoch specified by 2nd argument
     * \param julianDayAtEpoch Reference epoch in Julian days from which number of seconds are
     * counted (default J2000).
     * \return Angular velocity vector, expressed in target frame.
     */
    virtual Eigen::Vector3d getRotationalVelocityVectorInTargetFrame(
            const double secondsSinceEpoch, const double julianDayAtEpoch =
            basic_astrodynamics::JULIAN_DAY_ON_J2000  )
    {
        return getRotationToTargetFrame( secondsSinceEpoch, julianDayAtEpoch ) *
                getRotationalVelocityVectorInBaseFrame( secondsSinceEpoch, julianDayAtEpoch );
    }

    //! Function to calculate the full rotational state at given time
    /*!
     * Function to calculate the full rotational state at given time (rotation matrix, derivative of rotation matrix
     * and angular velocity vector).
     * \param currentRotationToLocalFrame Current rotation to local frame (returned by reference)
     * \param currentRotationToLocalFrameDerivative Current derivative of rotation matrix to local frame
     * (returned by reference)
     * \param currentAngularVelocityVectorInGlobalFrame Current angular velocity vector, expressed in global frame
     * (returned by reference)
     * \param secondsSinceEpoch Seconds since Julian day epoch specified by 2nd argument
     * \param julianDayAtEpoch Reference epoch in Julian days from which number of seconds are
     * counted (default J2000).
     */
    virtual void getFullRotationalQuantitiesToTargetFrame(
            Eigen::Quaterniond& currentRotationToLocalFrame,
            Eigen::Matrix3d& currentRotationToLocalFrameDerivative,
            Eigen::Vector3d& currentAngularVelocityVectorInGlobalFrame,
            const double secondsSinceEpoch, const double julianDayAtEpoch =
            basic_astrodynamics::JULIAN_DAY_ON_J2000 )
    {
        currentRotationToLocalFrame = getRotationToTargetFrame( secondsSinceEpoch, julianDayAtEpoch );
        currentRotationToLocalFrameDerivative = getDerivativeOfRotationToTargetFrame( secondsSinceEpoch, julianDayAtEpoch );
        currentAngularVelocityVectorInGlobalFrame = getRotationalVelocityVectorInBaseFrameFromMatrices(
                    Eigen::Matrix3d( currentRotationToLocalFrame ), currentRotationToLocalFrameDerivative.transpose( ) );
    }

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
