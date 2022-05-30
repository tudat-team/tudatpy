/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ROTATIONAL_EPHEMERIS_H
#define TUDAT_ROTATIONAL_EPHEMERIS_H

#include <string>

//https://stackoverflow.com/questions/28914711/getting-error-shared-ptr-in-namespace-std-does-not-name-a-type
#include <memory>

#include <functional>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/basics/timeType.h"

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
Eigen::Matrix< StateScalarType, 6, 1 > transformStateToFrameFromRotations(
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

template< typename StateScalarType >
Eigen::Matrix< StateScalarType, 6, 1 > transformStateToFrameFromRotations(
        const Eigen::Matrix< StateScalarType, 6, 1 >& stateInBaseFrame,
        const Eigen::Matrix3d& rotationToFrame,
        const Eigen::Matrix3d& rotationMatrixToFrameDerivative )
{
    return transformStateToFrameFromRotations< StateScalarType >(
                stateInBaseFrame, Eigen::Quaterniond( rotationToFrame ), rotationMatrixToFrameDerivative );
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
Eigen::Matrix< StateScalarType, 6, 1 > transformStateToFrameFromRotationFunctions(
        const Eigen::Matrix< StateScalarType, 6, 1 >& stateInBaseFrame,
        const std::function< Eigen::Quaterniond( ) > rotationToFrameFunction,
        const std::function< Eigen::Matrix3d( ) > rotationMatrixToFrameDerivativeFunction )
{
    return transformStateToFrameFromRotations< StateScalarType >(
                stateInBaseFrame, rotationToFrameFunction( ),
                rotationMatrixToFrameDerivativeFunction( ) );
}

//! Transform a relative state (Cartesian position and velocity) from one frame to another.
/*!
 *  Transform a relative state (Cartesian position and velocity) from one frame to another, taking into
 *  account both the instantaneous rotational state of the two frames, and the rotational
 *  rate of one frame w.r.t. the other.
 *  \param stateInBaseFrame State that is to be transformed from base to target frame.
 *  \param centralBodyStateInBaseFrame State of central body w.r.t. which returned state is to be computed.
 *  State returned by this function must be in frame with same orientation as that returned by stateInBaseFrame.
 *  \param rotationToFrameFunction Function returning rotation from base to target frame.
 *  \param rotationMatrixToFrameDerivativeFunction Function returning time derivative of rotation
 *   matrix from base to target frame.
 *  \return State (Cartesian position and velocity) in target frame.
 */
template< typename StateScalarType >
Eigen::Matrix< StateScalarType, 6, 1 > transformRelativeStateToFrame(
        const std::function< Eigen::Matrix< StateScalarType, 6, 1 >( ) > stateInBaseFrame,
        const std::function< Eigen::Matrix< StateScalarType, 6, 1 >( ) > centralBodyStateInBaseFrame,
        const std::function< Eigen::Quaterniond( ) > rotationToFrameFunction,
        const std::function< Eigen::Matrix3d( ) > rotationMatrixToFrameDerivativeFunction )
{
    return transformStateToFrameFromRotations< StateScalarType >(
                std::move( stateInBaseFrame( ) ) -
                std::move( centralBodyStateInBaseFrame( ) ),
                std::move( rotationToFrameFunction( ) ),
                std::move( rotationMatrixToFrameDerivativeFunction( ) ) );
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
template< typename StateScalarType, typename TimeType >
Eigen::Matrix< StateScalarType, 6, 1 > transformStateToFrameFromRotationTimeFunctions(
        const Eigen::Matrix< StateScalarType, 6, 1 >& stateInBaseFrame,
        const double currentTime,
        const std::function< Eigen::Quaterniond( const TimeType ) > rotationToFrameFunction,
        const std::function< Eigen::Matrix3d( const TimeType ) > rotationMatrixToFrameDerivativeFunction )
{
    return transformStateToFrameFromRotations(
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
     * \param secondsSinceEpoch Seconds since Julian day reference epoch.
     * \return Rotation quaternion computed.
     */
    virtual Eigen::Quaterniond getRotationToBaseFrame(
            const double secondsSinceEpoch ) = 0;


    Eigen::Matrix3d getRotationMatrixToBaseFrame( const double secondsSinceEpoch )
    {
        return Eigen::Matrix3d( getRotationToBaseFrame( secondsSinceEpoch ) );
    }

    //! Get rotation quaternion from target frame to base frame in Time precision.
    /*!
     * Pure virtual function to calculate and return the rotation quaternion from target frame to
     * base frame at specified time.
     * \param timeSinceEpoch Seconds since Julian day reference epoch.
     * \return Rotation quaternion computed.
     */
    virtual Eigen::Quaterniond getRotationToBaseFrameFromExtendedTime(
            const Time timeSinceEpoch )
    {
        return getRotationToBaseFrame( timeSinceEpoch.getSeconds< double >( ) );
    }

    //! Get rotation quaternion from target frame to base frame in templated precision.
    /*!
     * Pure virtual function to calculate and return the rotation quaternion from target frame to
     * base frame at specified time.
     * \param timeSinceEpoch Seconds since Julian day reference epoch.
     * \return Rotation quaternion computed.
     */
    template< typename TimeType >
    Eigen::Quaterniond getRotationToBaseFrameTemplated(
            const TimeType timeSinceEpoch );

    //! Get rotation quaternion to target frame from base frame.
    /*!
     * Pure virtual function to calculate and return the rotation quaternion to target frame from
     * base frame at specified time.
     * \param secondsSinceEpoch Seconds since Julian day reference epoch.
     * \return Rotation quaternion computed.
     */
    virtual Eigen::Quaterniond getRotationToTargetFrame(
            const double secondsSinceEpoch ) = 0;

    Eigen::Matrix3d getRotationMatrixToTargetFrame( const double secondsSinceEpoch )
    {
        return Eigen::Matrix3d( getRotationToTargetFrame( secondsSinceEpoch ) );
    }


    //! Get rotation quaternion to target frame from base frame in Time precision.
    /*!
     * Pure virtual function to calculate and return the rotation quaternion to target frame from
     * base frame at specified time.
     * \param timeSinceEpoch Seconds since Julian day reference epoch.
     * \return Rotation quaternion computed.
     */
    virtual Eigen::Quaterniond getRotationToTargetFrameFromExtendedTime(
            const Time timeSinceEpoch )
    {
        return getRotationToTargetFrame( timeSinceEpoch.getSeconds< double >( ) );
    }

    //! Get rotation quaternion to target frame from base frame in templated precision
    /*!
     * Pure virtual function to calculate and return the rotation quaternion to target frame from
     * base frame at specified time.
     * \param secondsSinceEpoch Seconds since Julian day reference epoch.
     * \return Rotation quaternion computed.
     */
    template< typename TimeType >
    Eigen::Quaterniond getRotationToTargetFrameTemplated(
            const TimeType secondsSinceEpoch );

    //! Function to calculate the derivative of the rotation matrix from target frame to original frame.
    /*!
     *  Function to calculate the derivative of the rotation matrix from target frame to original
     *  frame at specified time, to be implemented by derived class.
     * \param secondsSinceEpoch Seconds since Julian day reference epoch.
     *  \return Derivative of rotation from target (typically local) to original (typically global)
     *          frame at specified time.
     */
    virtual Eigen::Matrix3d getDerivativeOfRotationToBaseFrame(
            const double secondsSinceEpoch ) = 0;

    //! Function to calculate the derivative of the rotation matrix from target frame to original frame in Time precision.
    /*!
     *  Function to calculate the derivative of the rotation matrix from target frame to original
     *  frame at specified time, to be implemented by derived class.
     * \param timeSinceEpoch Seconds since Julian day reference epoch.
     *  \return Derivative of rotation from target (typically local) to original (typically global)
     *          frame at specified time.
     */
    virtual Eigen::Matrix3d getDerivativeOfRotationToBaseFrameFromExtendedTime(
            const Time timeSinceEpoch )
    {
        return getDerivativeOfRotationToBaseFrame( timeSinceEpoch.getSeconds< double >( ) );
    }

    //! Function to calculate the derivative of the rotation matrix from target frame to original frame in templated precision
    /*!
     *  Function to calculate the derivative of the rotation matrix from target frame to original
     *  frame at specified time, to be implemented by derived class.
     * \param timeSinceEpoch Seconds since Julian day reference epoch.
     *  \return Derivative of rotation from target (typically local) to original (typically global)
     *          frame at specified time.
     */
    template< typename TimeType >
    Eigen::Matrix3d getDerivativeOfRotationToBaseFrameTemplated(
            const TimeType timeSinceEpoch );

    //! Function to calculate the derivative of the rotation matrix from original frame to target frame.
    /*!
     *  Function to calculate the derivative of the rotation matrix from original frame to target
     *  frame at specified time, to be implemented by derived class.
     * \param secondsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated.
     *  \return Derivative of rotation from original (typically global) to target (typically local)
     *          frame at specified time.
     */
    virtual Eigen::Matrix3d getDerivativeOfRotationToTargetFrame(
            const double secondsSinceEpoch ) = 0;

    //! Function to calculate the derivative of the rotation matrix from original frame to target frame in Time precision.
    /*!
     *  Function to calculate the derivative of the rotation matrix from original frame to target
     *  frame at specified time, to be implemented by derived class.
     * \param timeSinceEpoch Seconds since epoch at which ephemeris is to be evaluated.
     *  \return Derivative of rotation from original (typically global) to target (typically local)
     *          frame at specified time.
     */
    virtual Eigen::Matrix3d getDerivativeOfRotationToTargetFrameFromExtendedTime(
            const Time timeSinceEpoch )
    {
        return getDerivativeOfRotationToTargetFrame( timeSinceEpoch.getSeconds< double >( ) );
    }

    //! Function to calculate the derivative of the rotation matrix from original frame to target frame in templated precision
    /*!
     *  Function to calculate the derivative of the rotation matrix from original frame to target
     *  frame at specified time, to be implemented by derived class.
     * \param secondsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated.
     *  \return Derivative of rotation from original (typically global) to target (typically local)
     *          frame at specified time.
     */
    template< typename TimeType >
    Eigen::Matrix3d getDerivativeOfRotationToTargetFrameTemplated(
            const TimeType secondsSinceEpoch );

    //! Function to retrieve the angular velocity vector, expressed in base frame.
    /*!
     * Function to retrieve the angular velocity vector, expressed in base frame. Thsi function uses the functions
     * that calculate the rotation matrix and its time derivative. It may be redefined in a derived class, to
     * calculate the angular velocity vector directly
     * \param secondsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated.
     * \return Angular velocity vector, expressed in base frame.
     */
    virtual Eigen::Vector3d getRotationalVelocityVectorInBaseFrame(
            const double secondsSinceEpoch )
    {
        return getRotationalVelocityVectorInBaseFrameFromMatrices(
                    Eigen::Matrix3d( getRotationToTargetFrame( secondsSinceEpoch ) ),
                    getDerivativeOfRotationToBaseFrame( secondsSinceEpoch ) );
    }

    //! Function to retrieve the angular velocity vector, expressed in target frame.
    /*!
     * Function to retrieve the angular velocity vector, expressed in target frame. This function calls the function
     * that calculates the angular velocioty vector in the base frame, and rotates it to the target frame.
     * It may be redefined in a derived class, to calculate the angular velocity vector directly
     * \param secondsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated.
     * \return Angular velocity vector, expressed in target frame.
     */
    virtual Eigen::Vector3d getRotationalVelocityVectorInTargetFrame(
            const double secondsSinceEpoch )
    {
        return getRotationToTargetFrame( secondsSinceEpoch ) *
                getRotationalVelocityVectorInBaseFrame( secondsSinceEpoch );
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
     * \param secondsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated.
     */
    virtual void getFullRotationalQuantitiesToTargetFrame(
            Eigen::Quaterniond& currentRotationToLocalFrame,
            Eigen::Matrix3d& currentRotationToLocalFrameDerivative,
            Eigen::Vector3d& currentAngularVelocityVectorInGlobalFrame,
            const double secondsSinceEpoch )
    {
        currentRotationToLocalFrame = getRotationToTargetFrame( secondsSinceEpoch );
        currentRotationToLocalFrameDerivative = getDerivativeOfRotationToTargetFrame( secondsSinceEpoch );
        currentAngularVelocityVectorInGlobalFrame = getRotationalVelocityVectorInBaseFrameFromMatrices(
                    Eigen::Matrix3d( currentRotationToLocalFrame ), currentRotationToLocalFrameDerivative.transpose( ) );
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
     * \param timeSinceEpoch Seconds since epoch at which ephemeris is to be evaluated.
     */
    virtual void getFullRotationalQuantitiesToTargetFrameFromExtendedTime(
            Eigen::Quaterniond& currentRotationToLocalFrame,
            Eigen::Matrix3d& currentRotationToLocalFrameDerivative,
            Eigen::Vector3d& currentAngularVelocityVectorInGlobalFrame,
            const Time timeSinceEpoch )
    {
        currentRotationToLocalFrame = getRotationToTargetFrameFromExtendedTime( timeSinceEpoch );
        currentRotationToLocalFrameDerivative = getDerivativeOfRotationToTargetFrameFromExtendedTime( timeSinceEpoch );
        currentAngularVelocityVectorInGlobalFrame = getRotationalVelocityVectorInBaseFrameFromMatrices(
                    Eigen::Matrix3d( currentRotationToLocalFrame ), currentRotationToLocalFrameDerivative.transpose( ) );
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
     * \param timeSinceEpoch Seconds since epoch at which ephemeris is to be evaluated.
     */
    template< typename TimeType >
    void getFullRotationalQuantitiesToTargetFrameTemplated(
            Eigen::Quaterniond& currentRotationToLocalFrame,
            Eigen::Matrix3d& currentRotationToLocalFrameDerivative,
            Eigen::Vector3d& currentAngularVelocityVectorInGlobalFrame,
            const TimeType timeSinceEpoch );

    //! Function to retrieve the current rotation as a state vector
    /*!
     * Function to retrieve the current rotation as a state vector (quaternion entries of body-fixed to base frame, and
     * angular velocity vector in body-fixed frame.
     * \param time Seconds since epoch at which ephemeris is to be evaluated.
     * \return Current rotation as a state vector
     */
    Eigen::Vector7d getRotationStateVector( const double time )
    {
        Eigen::Vector7d rotationalState;
        rotationalState.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat(
                    getRotationToBaseFrame( time ) );
        rotationalState.segment( 4, 3 ) = getRotationalVelocityVectorInTargetFrame( time );
        return rotationalState;
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

    virtual void resetCurrentTime( ){ }

    virtual void setIsBodyInPropagation( const bool isBodyInPropagation )
    { }
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

//! Function to transform a state from the target to base frame of a rotational ephemeris
/*!
 *  Function to transform a state from the target (body-fixed) to base (inertial) frame of a rotational ephemeris
 *  \param stateInLocalFrame State in body-fixed frame (target frame of rotational ephemeris)
 *  \param currentTime Time at which rotational ephemeris is to be evaluated
 *  \param rotationalEphemeris Rotational ephemeris object to compute the rotation.
 *  \return stateInLocalFrame State in inertial frame (base frame of rotational ephemeris)
 */
template< typename StateScalarType, typename TimeType >
Eigen::Matrix< StateScalarType, 6, 1 > transformStateToInertialOrientation(
        const Eigen::Matrix< StateScalarType, 6, 1 >& stateInLocalFrame,
        const TimeType currentTime,
        const std::shared_ptr< RotationalEphemeris > rotationalEphemeris )
{
    return transformStateToFrameFromRotations< StateScalarType >(
                stateInLocalFrame, rotationalEphemeris->getRotationToBaseFrameTemplated< TimeType >( currentTime ),
                rotationalEphemeris->getDerivativeOfRotationToBaseFrameTemplated< TimeType >( currentTime ) );

}

//! Function to transform a state from the base to target frame of a rotational ephemeris
/*!
 *  Function to transform a state from the base (inertial) to target (body-fixed) frame of a rotational ephemeris
 *  \param stateInGlobalFrame State in inertial frame (base frame of rotational ephemeris)
 *  \param currentTime Time at which rotational ephemeris is to be evaluated
 *  \param rotationalEphemeris Rotational ephemeris object to compute the rotation.
 *  \return State in body-fixed frame (target frame of rotational ephemeris)
 */
template< typename StateScalarType, typename TimeType >
Eigen::Matrix< StateScalarType, 6, 1 > transformStateToTargetFrame(
        const Eigen::Matrix< StateScalarType, 6, 1 >& stateInGlobalFrame,
        const TimeType currentTime,
        const std::shared_ptr< RotationalEphemeris > rotationalEphemeris )
{
    return transformStateToFrameFromRotations< StateScalarType >(
                stateInGlobalFrame, rotationalEphemeris->getRotationToTargetFrameTemplated< TimeType >( currentTime ),
                rotationalEphemeris->getDerivativeOfRotationToTargetFrameTemplated< TimeType >( currentTime ) );

}

} // namespace ephemerides

} // namespace tudat

#endif // TUDAT_ROTATIONAL_EPHEMERIS_H
