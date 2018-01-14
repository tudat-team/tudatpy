/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      J.S. Torok, Analytical Mechanics, Wiley-Interscience, 2000.
 *
 */

#ifndef TUDAT_APPARENT_ACCELERATION_MODEL_H
#define TUDAT_APPARENT_ACCELERATION_MODEL_H

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include <Eigen/Dense>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

namespace tudat
{
namespace reference_frames
{

//! Compute apparent acceleration due to non-inertiality of reference frame.
/*!
 * Computes apparent acceleration experienced by an object, of which the position and velocity
 * are provided in a non-inertial (i.e. rotating and/or accelerating) reference frame.
 * The apparent acceleration of a particle due to the non-inertiality of the reference frame in
 * which its state is given, is given by following equation:
 * \f[
 *      \boldsymbol{a}_{apparent} =
 *      -\boldsymbol{a}_B - \dot{\boldsymbol{\omega}} \times \boldsymbol{r}_{ni} -
 *      \boldsymbol{\omega} \times (\boldsymbol{\omega} \times \boldsymbol{r}_{ni}) -
 *      2(\boldsymbol{\omega} \times \boldsymbol{v}_{ni})
 * \f]
 * where \f$\boldsymbol{a}_B\f$ is the acceleration of the non-inertial frame with respect to
 * an inertial reference frame, \f$\boldsymbol{\omega}\f$ is the respective rotation rate
 * and \f$\boldsymbol{r}_{ni}\f$ and \f$\boldsymbol{v}_{ni}\f$ are the object's position and
 * velocity vector in the non-inertial frame in which the apparent acceleration is computed.
 * \param accelerationOfNonInertialReferenceFrame Acceleration vector of the non-inertial
 *          frame with respect to an inertial reference frame [m s^-2].
 * \param angularVelocityOfNonInertialReferenceFrame Angular velocity vector of the non-inertial
 *          frame with respect to an inertial reference frame [rad s^-1].
 * \param angularAccelerationOfNonInertialReferenceFrame Angular acceleration vector of the
 *          non-inertial frame with respect to an inertial reference frame [rad s^-1].
 * \param positionOfBodyInNonInertialReferenceFrame Position vector of body in the non-inertial
 *          frame of reference in which the apparent acceleration is computed [m].
 * \param velocityOfBodyInNonInertialReferenceFrame Velocity vector of body in the non-inertial
 *          frame of reference in which the apparent acceleration is computed [m s^-1].
 * \return Apparent acceleration as seen by an observer in the rotating and accelerating frame.
 */
Eigen::Vector3d computeApparentAcceleration(
        const Eigen::Vector3d& accelerationOfNonInertialReferenceFrame,
        const Eigen::Vector3d& angularVelocityOfNonInertialReferenceFrame,
        const Eigen::Vector3d& angularAccelerationOfNonInertialReferenceFrame,
        const Eigen::Vector3d& positionOfBodyInNonInertialReferenceFrame,
        const Eigen::Vector3d& velocityOfBodyInNonInertialReferenceFrame );

//! Compute centripetal acceleration due to non-inertiality of reference frame.
/*!
  * Computes centripetal acceleration experienced by an object, of which the position is
  * provided in a rotating reference frame.
  * The centripetal acceleration of a particle due to the non-inertiality of the reference frame in
  * which its state is given, is given by following equation:
  * \f[
  *      \boldsymbol{a}_{centripetal} =
  *      -\boldsymbol{\omega} \times (\boldsymbol{\omega} \times \boldsymbol{r}_{ni})
  * \f]
  * where \f$\boldsymbol{\omega}\f$ is the rotation rate of the non-inertial frame with respect to
  * an inertial reference frame  and \f$\boldsymbol{r}_{ni}\f$ is the object's position in the
  * non-inertial frame in which the apparent acceleration is computed.
  * \param angularVelocityOfNonInertialReferenceFrame Angular velocity vector of the non-inertial
  *          frame with respect to an inertial reference frame [rad s^-1].
  * \param positionOfBodyInNonInertialReferenceFrame Position vector of body in the non-inertial
  *          frame of reference in which the apparent acceleration is computed [m].
  * \return Centripetal acceleration as seen by an observer in the rotating frame.
  */
Eigen::Vector3d computeCentripetalAcceleration(
        const Eigen::Vector3d& angularVelocityOfNonInertialReferenceFrame,
        const Eigen::Vector3d& positionOfBodyInNonInertialReferenceFrame );

//! Compute Coriolis acceleration due to non-inertiality of reference frame.
/*!
  * Computes Coriolis acceleration experienced by an object, of which the position is
  * provided in a rotating reference frame.
  * The Coriolis acceleration of a particle due to the non-inertiality of the reference frame in
  * which its state is given, is given by following equation:
  * \f[
  *      \boldsymbol{a}_{Coriolis} =
  *      -2(\boldsymbol{\omega} \times \boldsymbol{v}_{ni})
  * \f]
  * where \f$\boldsymbol{\omega}\f$ is the rotation rate of the non-inertial frame with respect to
  * an inertial reference frame and \f$\boldsymbol{v}_{ni}\f$ is the object's velocity in the
  * non-inertial frame in which the apparent acceleration is computed.
  * \param angularVelocityOfNonInertialReferenceFrame Angular velocity vector of the non-inertial
  *          frame with respect to an inertial reference frame [rad s^-1].
  * \param velocityOfBodyInNonInertialReferenceFrame Velocity vector of body in the non-inertial
  *          frame of reference in which the apparent acceleration is computed [m s^-1].
  * \return Coriolis acceleration as seen by an observer in the rotating frame.
  */
Eigen::Vector3d computeCoriolisAcceleration(
        const Eigen::Vector3d& angularVelocityOfNonInertialReferenceFrame,
        const Eigen::Vector3d& velocityOfBodyInNonInertialReferenceFrame );

//! Compute Euler acceleration due to non-inertiality of reference frame.
/*!
  * Computes Euler acceleration experienced by an object, of which the position is
  * provided in a rotationally accelerating reference frame.
  * The Euler acceleration of a particle due to the non-inertiality of the reference frame in
  * which its state is given, is given by following equation:
  * \f[
  *      \boldsymbol{a}_{Euler} =
  *      - \dot{\boldsymbol{\omega}} \times \boldsymbol{r}_{ni}
  * \f]
  * where \f$\dot{\boldsymbol{\omega}}\f$ is the rate of change of the rotation rate of the
  * non-inertial frame with respect to an inertial reference frame and \f$\boldsymbol{r}_{ni}\f$
  * is the object's position in the non-inertial frame in which the apparent acceleration is
  * computed.
  * \param angularAccelerationOfNonInertialReferenceFrame Angular acceleration vector of the
  *          non-inertial frame with respect to an inertial reference frame [rad s^-1].
  * \param positionOfBodyInNonInertialReferenceFrame Position vector of body in the non-inertial
  *          frame of reference in which the apparent acceleration is computed [m].
  * \return Euler acceleration as seen by an observer in the rotationally accelerating frame.
  */
Eigen::Vector3d computeEulerAcceleration(
        const Eigen::Vector3d& angularAccelerationOfNonInertialReferenceFrame,
        const Eigen::Vector3d& positionOfBodyInNonInertialReferenceFrame );

//! Apparent acceleration model class.
/*!
 * Implementation of apparent acceleration due to non-inertiality of a reference system in which
 * the equations of motion are evaluated. It evaluates the total apparent acceleration due to the
 * acceleration of the reference frame, as well as the Coriolis, centripetal and Euler
 * accelerations.
 */
class ApparentAccelerationModel : public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{
private:

    //! Typedef for Eigen::Vector3d returning function.
    typedef boost::function< Eigen::Vector3d( ) > Vector3dReturningFunction;

public:

    //! Class constructor.
    /*!
     * Constructor for apparent acceleration model.
     * \param accelerationOfNonInertialReferenceFrameFunction Pointer to a function returning
     *         the acceleration vector of the non-intertial reference frame w.r.t. an
     *         inertial frame of reference.
     * \param angularVelocityOfNonInertialReferenceFrameFunction Pointer to a function returning
     *         the angular velocity vector of the non-intertial reference frame w.r.t. an
     *         inertial frame of reference.
     * \param angularAccelerationOfNonInertialReferenceFrameFunction Pointer to a function 
     *         returning the angular acceleration vector of the non-inertial frame with respect to 
     *         an inertial reference frame.
     * \param positionOfBodyInNonInertialReferenceFrameFunction Pointer to a function returning
     *         the position vector in the non-inertial frame of reference in which the apparent
     *         acceleration is computed.
     * \param velocityOfBodyInNonInertialReferenceFrameFunction Pointer to a function returning
     *          the velocity vector in the non-inertial frame of reference in which the apparent
     *          acceleration is computed.
     */
    ApparentAccelerationModel(
            Vector3dReturningFunction accelerationOfNonInertialReferenceFrameFunction,
            Vector3dReturningFunction angularVelocityOfNonInertialReferenceFrameFunction,
            Vector3dReturningFunction angularAccelerationOfNonInertialReferenceFrameFunction,
            Vector3dReturningFunction positionOfBodyInNonInertialReferenceFrameFunction,
            Vector3dReturningFunction velocityOfBodyInNonInertialReferenceFrameFunction )
        : accelerationOfNonInertialReferenceFrameFunction_(
              accelerationOfNonInertialReferenceFrameFunction ),
          angularVelocityOfNonInertialReferenceFrameFunction_(
              angularVelocityOfNonInertialReferenceFrameFunction ),
          angularAccelerationOfNonInertialReferenceFrameFunction_(
              angularAccelerationOfNonInertialReferenceFrameFunction ),
          positionOfBodyInNonInertialReferenceFrameFunction_(
              positionOfBodyInNonInertialReferenceFrameFunction ),
          velocityOfBodyInNonInertialReferenceFrameFunction_(
              velocityOfBodyInNonInertialReferenceFrameFunction )
    {
        updateMembers( );
    }

    //! Get apparent acceleration.
    /*!
     * Computes and returns the acceleration.
     * \return Vector of the apparent acceleration in the non-inertial Cartesian frame in which the
     *          object's state is defined.
     */
    Eigen::Vector3d getAcceleration( )
    {
        return computeApparentAcceleration( currentAccelerationOfNonInertialReferenceFrame_,
                                            currentAngularVelocityOfNonInertialReferenceFrame_,
                                            currentAngularAccelerationOfNonInertialReferenceFrame_,
                                            currentPositionOfBodyInNonInertialReferenceFrame_,
                                            currentVelocityOfBodyInNonInertialReferenceFrame_ );
    }

    //! Update member variables used by apparent acceleration model.
    /*!
     * Function to update member variables used by this acceleration model. The variables
     * that are required as input for the evaluation of the accelerations are retrieved from
     * boost::functions that may or may not give constant return values, depending on the
     * user input. This function sets the member variables using these functions.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN );

protected:

private:

    //! Function returning translational acceleration vector of non-inertial reference frame.
    /*!
     * Function returning translational acceleration vector of non-inertial reference frame.
     */
    const Vector3dReturningFunction accelerationOfNonInertialReferenceFrameFunction_;

    //! Function returning angular velocity vector of non-inertial reference frame.
    /*!
     * Function returning angular velocity vector of non-inertial reference frame.
     */
    const Vector3dReturningFunction angularVelocityOfNonInertialReferenceFrameFunction_;

    //! Function returning angular acceleration vector of non-inertial reference frame.
    /*!
     * Function returning angular acceleration vector of non-inertial reference frame.
     */
    const Vector3dReturningFunction angularAccelerationOfNonInertialReferenceFrameFunction_;

    //! Function returning position vector in non-inertial reference frame.
    /*!
     * Function returning position vector in non-inertial reference frame.
     */
    const Vector3dReturningFunction positionOfBodyInNonInertialReferenceFrameFunction_;

    //! Function returning velocity vector in non-inertial reference frame.
    /*!
     * Function returning velocity vector in non-inertial reference frame.
     */
    const Vector3dReturningFunction velocityOfBodyInNonInertialReferenceFrameFunction_;

    //! Current translational acceleration vector of non-inertial reference frame.
    /*!
     * Current translational acceleration vector of non-inertial reference frame. This vector
     * is updated by the updateMembers( ) function.
     */
    Eigen::Vector3d currentAccelerationOfNonInertialReferenceFrame_;

    //! Current angular velocity vector of non-inertial reference frame.
    /*!
     * Current angular velocity vector of non-inertial reference frame. This vector
     * is updated by the updateMembers( ) function.
     */
    Eigen::Vector3d currentAngularVelocityOfNonInertialReferenceFrame_;

    //! Current angular acceleration vector of non-inertial reference frame.
    /*!
     * Current angular acceleration vector of non-inertial reference frame. This vector
     * is updated by the updateMembers( ) function.
     */
    Eigen::Vector3d currentAngularAccelerationOfNonInertialReferenceFrame_;

    //! Current position vector in non-inertial reference frame.
    /*!
     * Current position vector in non-inertial reference frame. This vector
     * is updated by the updateMembers( ) function.
     */
    Eigen::Vector3d currentPositionOfBodyInNonInertialReferenceFrame_;

    //! Current velocity vector in non-inertial reference frame.
    /*!
     * Current velocity vector in non-inertial reference frame. This vector
     * is updated by the updateMembers( ) function.
     */
    Eigen::Vector3d currentVelocityOfBodyInNonInertialReferenceFrame_;
};

//! Typedef for shared-pointer to ApparentAccelerationModel object.
typedef boost::shared_ptr< ApparentAccelerationModel > ApparentAccelerationModelPointer;

}   // namespace reference_frames
}   // namespace tudat

#endif // TUDAT_APPARENT_ACCELERATION_MODEL_H
