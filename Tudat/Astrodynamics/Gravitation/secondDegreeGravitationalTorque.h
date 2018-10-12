/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SECONDDEGREEGRAVITATIONALTORQUE_H
#define TUDAT_SECONDDEGREEGRAVITATIONALTORQUE_H

#include <iostream>

#include <functional>

#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModel.h"

namespace tudat
{

namespace gravitation
{

//! Function to calculate the gravitational torque exerted by a point mass on a body with degree two gravity field
/*!
 * Function to calculate the gravitational torque exerted by a point mass on a body with degree two gravity field, which is
 * provided here as an inertia tensor. Higher order terms of the torque are omitted.
 * \param relativePositionOfBodySubjectToTorque Position of body exerting torque, w.r.t. body undergoing torque (typically
 * expressed in frame fixed to body undergoing torque).
 * \param premultipler Torque pre-multipler (3*mu/r^5).
 * \param inertiaTensorTimesRelativePositionOfBody The inertia tensor of the body undergoing the torque, multiplied by
 * relativePositionOfBodySubjectToTorque
 * \return Gravitational torque of point mass on second-degree body.
 */
Eigen::Vector3d calculateSecondDegreeGravitationalTorque(
        const Eigen::Vector3d& relativePositionOfBodySubjectToTorque,
        const double premultipler,
        const Eigen::Vector3d& inertiaTensorTimesRelativePositionOfBody );

//! Function to calculate the gravitational torque exerted by a point mass on a body with degree two gravity field
/*!
 * Function to calculate the gravitational torque exerted by a point mass on a body with degree two gravity field, which is
 * provided here as an inertia tensor. Higher order terms of the torque are omitted.
 * \param relativePositionOfBodySubjectToTorque Position of body exerting torque, w.r.t. body undergoing torque (typically
 * expressed in frame fixed to body undergoing torque).
 * \param gravitationalParameterOfAttractingBody Tha gravitational parameter of the body that exerts the torque
 * \param inertiaTensorOfRotatingBody The inertia tensor of the body undergoing the torqie, in the same frame as
 * relativePositionOfBodySubjectToTorque (typically frame fixed to body undergoing torque)
 * \return Gravitational torque of point mass on second-degree body.
 */
Eigen::Vector3d calculateSecondDegreeGravitationalTorque(
        const Eigen::Vector3d& relativePositionOfBodySubjectToTorque,
        const double gravitationalParameterOfAttractingBody,
        const Eigen::Matrix3d& inertiaTensorOfRotatingBody );

//! Class to compute the second degree gravitational torque
/*!
 *  Class to compute the second degree gravitational torque: the gravitational torque exerted by a point mass on a body with
 *  degree two gravity field (parameterized by inertia tensor).
 */
class SecondDegreeGravitationalTorqueModel: public basic_astrodynamics::TorqueModel
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param positionOfBodySubjectToTorqueFunction Function returning the position of the body that is subject to the torque,
     * in inertial frame.
     * \param gravitationalParameterOfAttractingBodyFunction Function returning the gravitational parameter of the body that is
     * exerting the torque.
     * \param inertiaTensorOfRotatingBodyFunction Function returning the inertia tensor of the body that is subject to the
     * torque, in the frame fixed to that body.
     * \param positionOfBodyExertingTorqueFunction Function returning the position of the body that is exerting to the torque,
     * in inertial frame.
     * \param rotationToBodyFixedFrameFunction Function returning the rotation from inertial frame to frame fixed to body
     *  undergoing torque.
     */
    SecondDegreeGravitationalTorqueModel(
            const std::function< Eigen::Vector3d( ) > positionOfBodySubjectToTorqueFunction,
            const std::function< double( ) > gravitationalParameterOfAttractingBodyFunction,
            const std::function< Eigen::Matrix3d( ) > inertiaTensorOfRotatingBodyFunction,
            const std::function< Eigen::Vector3d( ) > positionOfBodyExertingTorqueFunction =
            [ ]( ){ return Eigen::Vector3d::Zero( ); },
            const std::function< Eigen::Quaterniond( ) > rotationToBodyFixedFrameFunction =
            [ ]( ){ return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ); } ):
        positionOfBodySubjectToTorqueFunction_( positionOfBodySubjectToTorqueFunction ),
        gravitationalParameterOfAttractingBodyFunction_( gravitationalParameterOfAttractingBodyFunction ),
        inertiaTensorOfRotatingBodyFunction_( inertiaTensorOfRotatingBodyFunction ),
        positionOfBodyExertingTorqueFunction_( positionOfBodyExertingTorqueFunction ),
        rotationToBodyFixedFrameFunction_( rotationToBodyFixedFrameFunction ){  }

    //! Get gravitational torque.
    /*!
     * Returns the gravitational torque. All data required for the computation is taken
     * from member variables, which are set to their latest values by the last call of the
     * updateMembers function.
     * \return Gravitational torque.
     * \sa updateMembers().
     */
    Eigen::Vector3d getTorque( )
    {
        return currentTorque_;
    }

    //! Update member variables used by the gravitational torque model.
    /*!
     * Updates member variables used by the gravitational accfeleration model.
     * Function pointers to retrieve the current values of quantities from which the
     * torque is to be calculated are set by constructor. This function calls
     * them to update the associated variables to their current state.
     * \param currentTime Time at which torque model is to be updated.
     */
    void updateMembers( const double currentTime )
    {
        currentRotationToBodyFixedFrame_ = rotationToBodyFixedFrameFunction_( );
        currentRelativePositionOfBodySubjectToTorque_ = positionOfBodyExertingTorqueFunction_( ) - positionOfBodySubjectToTorqueFunction_( );
        currentRotationToBodyFixedFrameFunction_ = currentRotationToBodyFixedFrame_;
        currentRelativeBodyFixedPositionOfBodySubjectToTorque_ =
                currentRotationToBodyFixedFrameFunction_ * currentRelativePositionOfBodySubjectToTorque_;

        currentGravitationalParameterOfAttractingBody_ = gravitationalParameterOfAttractingBodyFunction_( );
        currentInertiaTensorOfRotatingBody_ = inertiaTensorOfRotatingBodyFunction_( );

        currentInertiaTensorTimesRelativePositionOfBody_ = currentInertiaTensorOfRotatingBody_ *
                currentRelativeBodyFixedPositionOfBodySubjectToTorque_;
        currentTorqueMagnitudePremultiplier_ =
                3.0 * currentGravitationalParameterOfAttractingBody_ /
                std::pow(
                    currentRelativePositionOfBodySubjectToTorque_.norm( ), 5.0 );

        currentTorque_ = calculateSecondDegreeGravitationalTorque(
                    currentRelativeBodyFixedPositionOfBodySubjectToTorque_,
                    currentTorqueMagnitudePremultiplier_,
                    currentInertiaTensorTimesRelativePositionOfBody_ );
    }

    //! Function to retrieve current rotation from inertial to body-fixed frame
    /*!
     * Function to retrieve current rotation from inertial to body-fixed frame
     * \return Current rotation from inertial to body-fixed frame
     */
    Eigen::Quaterniond& getCurrentRotationToBodyFixedFrame( )
    {
        return currentRotationToBodyFixedFrame_;
    }

    //! Function to retrieve current position of body exerting torque, w.r.t. body undergoing torque in inertial frame
    /*!
     * Function to retrieve current position of body exerting torque, w.r.t. body undergoing torque in inertial frame
     * \return Current position of body exerting torque, w.r.t. body undergoing torque in inertial frame
     */
    Eigen::Vector3d& getCurrentRelativePositionOfBodySubjectToTorque( )
    {
        return currentRelativePositionOfBodySubjectToTorque_;
    }

    //! Function to retrieve current position of body exerting torque, w.r.t. body undergoing torque in frame fixed to body
    /*!
     * Function to retrieve current position of body exerting torque, w.r.t. body undergoing torque in frame fixed to body
     * undergoing torque,
     * \return Current position of body exerting torque, w.r.t. body undergoing torque in frame fixed to body undergoing torque,
     */
    Eigen::Vector3d& getCurrentRelativeBodyFixedPositionOfBodySubjectToTorque( )
    {
        return currentRelativeBodyFixedPositionOfBodySubjectToTorque_;
    }

    //! Function to retrieve current inertia tensor of the body that is subject to the torque, in the frame fixed to that body
    /*!
     * Function to retrieve current inertia tensor of the body that is subject to the torque, in the frame fixed to that body
     * \return Current inertia tensor of the body that is subject to the torque, in the frame fixed to that body
     */
    Eigen::Matrix3d& getCurrentInertiaTensorOfRotatingBody( )
    {
        return currentInertiaTensorOfRotatingBody_;
    }

    //! Function to retrieve current torque pre-multipler
    /*!
     * Function to retrieve current torque pre-multipler
     * \return Current torque pre-multipler
     */
    double getCurrentTorqueMagnitudePremultiplier( )
    {
        return currentTorqueMagnitudePremultiplier_;
    }

    //! Function to retrieve current gravitational parameter of the body that is exerting the torque,
    /*!
     * Function to retrieve current gravitational parameter of the body that is exerting the torque,
     * \return Current gravitational parameter of the body that is exerting the torque,
     */
    double getCurrentGravitationalParameterOfAttractingBody( )
    {
        return currentGravitationalParameterOfAttractingBody_;
    }

    //! Function to retrieve The current inertia tensor of the body undergoing the torque, premultiplied by state
    /*!
     *  Function to retrieve The current inertia tensor of the body undergoing the torque, multiplied by relative
     *  currentRelativeBodyFixedPositionOfBodySubjectToTorque_
     *  \return The current inertia tensor of the body undergoing the torque, multiplied by relative
     *  currentRelativeBodyFixedPositionOfBodySubjectToTorque_
     */
    Eigen::Vector3d& getCurrentInertiaTensorTimesRelativePositionOfBody( )
    {
        return currentInertiaTensorTimesRelativePositionOfBody_;
    }

protected:

private:

    //! Function returning the position of the body that is subject to the torque, in inertial frame.
    std::function< Eigen::Vector3d( ) > positionOfBodySubjectToTorqueFunction_;

    //! Function returning the gravitational parameter of the body that is exerting the torque.
    std::function< double( ) > gravitationalParameterOfAttractingBodyFunction_;

    //! Function returning the inertia tensor of the body that is subject to the torque, in the frame fixed to that body.
    std::function< Eigen::Matrix3d( ) > inertiaTensorOfRotatingBodyFunction_;

    //! Function returning the position of the body that is exerting to the torque, in inertial frame.
    std::function< Eigen::Vector3d( ) > positionOfBodyExertingTorqueFunction_;

    //! Function returning the rotation from inertial frame to frame fixed to body undergoing torque.
    std::function< Eigen::Quaterniond( ) > rotationToBodyFixedFrameFunction_;

    //! Current rotation from inertial to body-fixed frame
    Eigen::Quaterniond currentRotationToBodyFixedFrame_;

    //! Current position of body exerting torque, w.r.t. body undergoing torque in inertial frame, as set
    //! by updateMembers function.
    Eigen::Vector3d currentRelativePositionOfBodySubjectToTorque_;

    //! Current position of body exerting torque, w.r.t. body undergoing torque in frame fixed to body undergoing torque, as set
    //! by updateMembers function.
    Eigen::Vector3d currentRelativeBodyFixedPositionOfBodySubjectToTorque_;

    //! Current gravitational parameter of the body that is exerting the torque, as set by updateMembers function.
    double currentGravitationalParameterOfAttractingBody_;

    //! Current inertia tensor of the body that is subject to the torque, in the frame fixed to that body, as set by updateMembers
    //! function.
    Eigen::Matrix3d currentInertiaTensorOfRotatingBody_;

    //! Rotation from inertial frame to frame fixed to body undergoing torque, as set by updateMembers function.
    Eigen::Quaterniond currentRotationToBodyFixedFrameFunction_;

    //! The current inertia tensor of the body undergoing the torque, multiplied by relative
    //! currentRelativeBodyFixedPositionOfBodySubjectToTorque_
    Eigen::Vector3d currentInertiaTensorTimesRelativePositionOfBody_;

    //! Current torque pre-multipler (3*mu/r^5).
    double currentTorqueMagnitudePremultiplier_;

    //! Current gravitational torque, as set by updateMembers function.
    Eigen::Vector3d currentTorque_;
};

}

}


#endif // TUDAT_SECONDDEGREEGRAVITATIONALTORQUE_H
