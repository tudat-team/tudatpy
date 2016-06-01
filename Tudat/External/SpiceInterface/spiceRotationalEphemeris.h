#ifndef TUDAT_SPICEROTATIONALEPHEMERIS_H
#define TUDAT_SPICEROTATIONALEPHEMERIS_H

#include <string>

#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"

namespace tudat
{

namespace ephemerides
{

//! Class to directly calculate rotational state of body from spice kernel(s)
/*!
 *  Class to directly calculate rotational state of body from spice kernel(s).
 *  Note that if an object of this class is called often, it may become a performance bottleneck,
 *  using the an interpolated rotation will generally be more efficient, but at the expense of
 *  interpolation errors.
 */
class SpiceRotationalEphemeris : public RotationalEphemeris
{
public:
    //! Constructor for spice ephemeris
    /*!
     *  Constructor for spice ephemeris, sets frames between which rotation is determined.
     * \param baseFrameOrientation Base frame identifier.
     * \param targetFrameOrientation Target frame identifier.
     */
    SpiceRotationalEphemeris( const std::string& baseFrameOrientation = "ECLIPJ2000",
                              const std::string& targetFrameOrientation = ""  ):
        RotationalEphemeris( baseFrameOrientation, targetFrameOrientation )
    { }

    //! Destructor
    /*!
     *  Destructor.
     */
    ~SpiceRotationalEphemeris( ){ }

    //! Function to calculate the rotation quaternion from target frame to original frame.
    /*!
     *  Function to calculate the rotation quaternion from target frame to original frame at
     *  specified time, calculated by call to pxform_c spice function
     *  (through computeRotationQuaternionBetweenFrames function)
     *  \param secondsSinceEpoch Seconds since Julian day epoch specified by 2nd argument
     *  at which rotation is to be calculated.
     *  \param julianDayAtEpoch Reference epoch in Julian days from which number of seconds are
     *  in secondsSinceEpoch are counted.
     *  \return Rotation from target (typically local) to original (typically global) frame at
     *  specified time.
     */
    Eigen::Quaterniond getRotationToBaseFrame(
            const double secondsSinceEpoch,
            const double julianDayAtEpoch = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    //! Function to calculate the rotation quaternion from original frame to target frame.
    /*!
     *  Function to calculate the rotation quaternion from original frame to target frame at
     *  specified time, calculated by call to pxform_c spice function
     *  (through computeRotationQuaternionBetweenFrames function)
     *  \param secondsSinceEpoch Seconds since Julian day epoch specified by 2nd argument
     *  at which rotation is to be calculated.
     *  \param julianDayAtEpoch Reference epoch in Julian days from which number of seconds are
     *  in secondsSinceEpoch are counted.
     *  \return Rotation from original (typically global) to target (typically local)
     *  frame at specified time.
     */
    Eigen::Quaterniond getRotationToTargetFrame(
            const double secondsSinceEpoch,
            const double julianDayAtEpoch = basic_astrodynamics::JULIAN_DAY_ON_J2000 )
    {
        return getRotationToBaseFrame( secondsSinceEpoch, julianDayAtEpoch ).inverse( );
    }

    //! Function to calculate the derivative of the rotation matrix from target frame to original
    //! frame.
    /*!
     *  Function to calculate the derivative of the rotation matrix from target frame to original
     *  frame at specified time, calculated by SPICE sx_form function.
     *  \param secondsSinceEpoch Seconds since Julian day epoch specified by 2nd argument
     *  \param julianDayAtEpoch Reference epoch in Julian days from which number of seconds are
     *          counted.
     *  \return Derivative of rotation from target (typically local) to original (typically global)
     *          frame at specified time.
     */
    Eigen::Matrix3d getDerivativeOfRotationToBaseFrame(
            const double secondsSinceEpoch, const double julianDayAtEpoch =
            basic_astrodynamics::JULIAN_DAY_ON_J2000 );


    //! Function to calculate the derivative of the rotation matrix from original frame to target
    //! frame.
    /*!
     *  Function to calculate the derivative of the rotation matrix from target frame to original
     *  frame at specified time, calculated by SPICE sx_form function.
     *  \param secondsSinceEpoch Seconds since Julian day epoch specified by 2nd argument
     *  \param julianDayAtEpoch Reference epoch in Julian days from which number of seconds are
     *          counted.
     *  \return Derivative of rotation from original (typically global) to target (typically local)
     *          frame at specified time.
     */
    Eigen::Matrix3d getDerivativeOfRotationToTargetFrame(
            const double secondsSinceEpoch, const double julianDayAtEpoch =
            basic_astrodynamics::JULIAN_DAY_ON_J2000 )
    {
        return getDerivativeOfRotationToBaseFrame( secondsSinceEpoch, julianDayAtEpoch ).
                transpose( );
    }

    //! Function to calculate the full rotational state at given time
    /*!
     * Function to calculate the full rotational state at given time (rotation matrix, derivative of
     * rotation matrix and angular velocity vector).
     * \param currentRotationToLocalFrame Current rotation to local frame (returned by reference)
     * \param currentRotationToLocalFrameDerivative Current derivative of rotation matrix to local
     * frame (returned by reference)
     * \param currentAngularVelocityVectorInGlobalFrame Current angular velocity vector, expressed
     * in global frame (returned by reference)
     * \param secondsSinceEpoch Seconds since Julian day epoch specified by 2nd argument
     * \param julianDayAtEpoch Reference epoch in Julian days from which number of seconds are
     * counted (default J2000).
     */
    void getFullRotationalQuantitiesToTargetFrame(
            Eigen::Quaterniond& currentRotationToLocalFrame,
            Eigen::Matrix3d& currentRotationToLocalFrameDerivative,
            Eigen::Vector3d& currentAngularVelocityVectorInGlobalFrame,
            const double secondsSinceEpoch, const double julianDayAtEpoch =
            basic_astrodynamics::JULIAN_DAY_ON_J2000 );

};

} // namespace ephemerides

} // namespace tudat
#endif // TUDAT_SPICEROTATIONALEPHEMERIS_H
