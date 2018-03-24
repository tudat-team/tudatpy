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


#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

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
 * \param gravitationalParameterOfAttractingBody Tha gravitational parameter of teh body that exerts the torque
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
            const boost::function< Eigen::Vector3d( ) > positionOfBodySubjectToTorqueFunction,
            const boost::function< double( ) > gravitationalParameterOfAttractingBodyFunction,
            const boost::function< Eigen::Matrix3d( ) > inertiaTensorOfRotatingBodyFunction,
            const boost::function< Eigen::Vector3d( ) > positionOfBodyExertingTorqueFunction =
            boost::lambda::constant( Eigen::Vector3d::Zero( ) ),
            const boost::function< Eigen::Quaterniond( ) > rotationToBodyFixedFrameFunction =
            boost::lambda::constant( Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) ) ):
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
        currentRelativePositionOfBodySubjectToTorque_ = positionOfBodyExertingTorqueFunction_( ) - positionOfBodySubjectToTorqueFunction_( );
        currentRotationToBodyFixedFrameFunction_ = rotationToBodyFixedFrameFunction_( );

        currentGravitationalParameterOfAttractingBody_ = gravitationalParameterOfAttractingBodyFunction_( );
        currentInertiaTensorOfRotatingBody_ = inertiaTensorOfRotatingBodyFunction_( );

        currentTorque_ = calculateSecondDegreeGravitationalTorque(
                    currentRotationToBodyFixedFrameFunction_ * currentRelativePositionOfBodySubjectToTorque_,
                    currentGravitationalParameterOfAttractingBody_,
                    currentInertiaTensorOfRotatingBody_ );
    }

protected:

private:

    //! Function returning the position of the body that is subject to the torque, in inertial frame.
    boost::function< Eigen::Vector3d( ) > positionOfBodySubjectToTorqueFunction_;

    //! Function returning the gravitational parameter of the body that is exerting the torque.
    boost::function< double( ) > gravitationalParameterOfAttractingBodyFunction_;

    //! Function returning the inertia tensor of the body that is subject to the torque, in the frame fixed to that body.
    boost::function< Eigen::Matrix3d( ) > inertiaTensorOfRotatingBodyFunction_;

    //! Function returning the position of the body that is exerting to the torque, in inertial frame.
    boost::function< Eigen::Vector3d( ) > positionOfBodyExertingTorqueFunction_;

    //! Function returning the rotation from inertial frame to frame fixed to body undergoing torque.
    const boost::function< Eigen::Quaterniond( ) > rotationToBodyFixedFrameFunction_;


    //! Current [osition of body exerting torque, w.r.t. body undergoing torque in frame fixed to body undergoing torque, as set
    //! by updateMembers function.
    Eigen::Vector3d currentRelativePositionOfBodySubjectToTorque_;

    //! Current gravitational parameter of the body that is exerting the torque, as set by updateMembers function.
    double currentGravitationalParameterOfAttractingBody_;

    //! Current inertia tensor of the body that is subject to the torque, in the frame fixed to that body, as set by updateMembers
    //! function.
    Eigen::Matrix3d currentInertiaTensorOfRotatingBody_;

    //! Rotation from inertial frame to frame fixed to body undergoing torque, as set by updateMembers function.
    Eigen::Quaterniond currentRotationToBodyFixedFrameFunction_;

    //! Current gravitational torque, as set by updateMembers function.
    Eigen::Vector3d currentTorque_;
};

}

}

#endif // TUDAT_SECONDDEGREEGRAVITATIONALTORQUE_H
