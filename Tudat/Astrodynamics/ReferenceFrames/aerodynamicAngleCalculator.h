/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#ifndef TUDAT_AERODYNAMICANGLECALCULATOR_H
#define TUDAT_AERODYNAMICANGLECALCULATOR_H

#include <vector>
#include <map>

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

namespace tudat
{

namespace reference_frames
{

//! Enum to define ids for various reference frames for calculating between inertial and body-fixed
//! frame, using transformation chain via aerodynamic frame.
enum AerodynamicsReferenceFrames
{
    inertial_frame = -1,
    corotating_frame = 0,
    vertical_frame = 1,
    trajectory_frame = 2,
    aerodynamic_frame = 3,
    body_frame = 4
};

//! Enum to define ids for various angles needed for converting between inertial and body-fixed
//! frame, using transformation chain via aerodynamic frame.
enum AerodynamicsReferenceFrameAngles
{
    latitude_angle = 0,
    longitude_angle = 1,
    heading_angle = 2,
    flight_path_angle = 3,
    angle_of_attack = 4,
    angle_of_sideslip = 5,
    bank_angle = 6
};

//! Object to calculate aerodynamic orientation angles from current vehicle state.
/*!
 *  Object to calculate aerodynamic orientation angles from current vehicle state.
 */
class AerodynamicAngleCalculator
{
public:

    //! Constructor
    /*!
     *  Constructor.
     *  \param bodyFixedStateFunction Vehicle state in a frame fixed w.r.t. the central body.  Note
     *  that this state is w.r.t. the body itself, not w.r.t. the local atmosphere
     *  \param angleOfAttackFunction Function to determine the angle of attack of the vehicle.
     *  \param angleOfSideslipFunction Function to determine the angle of sideslip of the vehicle.
     *  \param bankAngleFunction Function to determine the bank angle of the vehicle.
     *  \param calculateVerticalToAerodynamicFrame Boolean to determine whether to determine
     *  vertical <-> aerodynamic frame conversion when calling update function.
     */
    AerodynamicAngleCalculator(
            const boost::function< basic_mathematics::Vector6d( ) > bodyFixedStateFunction,
            const boost::function< double( ) > angleOfAttackFunction =
            boost::lambda::constant ( 0.0 ),
            const boost::function< double( ) > angleOfSideslipFunction =
            boost::lambda::constant ( 0.0 ),
            const boost::function< double( ) > bankAngleFunction =
            boost::lambda::constant ( 0.0 ),
            const bool calculateVerticalToAerodynamicFrame = 0 ):
        bodyFixedStateFunction_( bodyFixedStateFunction ),
        angleOfAttackFunction_( angleOfAttackFunction ),
        angleOfSideslipFunction_( angleOfSideslipFunction ),
        bankAngleFunction_( bankAngleFunction ),
        calculateVerticalToAerodynamicFrame_( calculateVerticalToAerodynamicFrame ){ }

    //! Function to update the orientation angles to the current state.
    /*!
     *  This function updates all requires orientation angles to the current state of the vehicle.
     *  The current state is retrieved from the bodyFixedStateFunction_ member variable
     *  function pointer.
     */
    void update( );

    //! Function to get the rotation quaternion between two frames
    /*!
     * Function to get the rotation quaternion between two frames. This function uses the values
     * calculated by the previous call of the update( ) function.
     * \param originalFrame Id for 'current' frame
     * \param targetFrame Id for frame to which transformation object should transfrom a vector
     * (from originalFrame).
     * \return Rotation quaternion from originalFrame to targetFrame.
     */
    Eigen::Quaterniond getRotationQuaternionBetweenFrames(
            const AerodynamicsReferenceFrames originalFrame,
            const AerodynamicsReferenceFrames targetFrame );

    //! Function to get a single orientation angle.
    /*!
     * Function to get a single orientation angle, as calculated by previous call to update( )
     * function.
     * \param angleId Id of requested angle.
     * \return Value of requested angle.
     */
    double getAerodynamicAngle( const AerodynamicsReferenceFrameAngles angleId );

    void setOrientationAngleFunctions(
            const boost::function< double( ) > angleOfAttackFunction =
            boost::lambda::constant ( 0.0 ),
            const boost::function< double( ) > angleOfSideslipFunction =
            boost::lambda::constant ( 0.0 ),
            const boost::function< double( ) > bankAngleFunction =
            boost::lambda::constant ( 0.0 ) )
    {
        angleOfAttackFunction_ = angleOfAttackFunction;
        angleOfSideslipFunction_ = angleOfSideslipFunction;
        bankAngleFunction_ = bankAngleFunction;
    }

private:

    //! Map of current angles, as calculated by previous call to update( ) function.
    std::map< AerodynamicsReferenceFrameAngles, double > currentAerodynamicAngles_;

    //! Map of current transformation quaternions, as calculated since previous call to update( )
    //! function.
    std::map< std::pair< AerodynamicsReferenceFrames, AerodynamicsReferenceFrames >,
    Eigen::Quaterniond > currentRotationMatrices_;

    //! Current body-fixed state of vehicle, as set by previous call to update( ).
    basic_mathematics::Vector6d currentBodyFixedState_;

    //! Vehicle state in a frame fixed w.r.t. the central body.
    /*!
     *  Vehicle state in a frame fixed w.r.t. the central body.
     *  Note that this state is w.r.t. the body itself, not w.r.t. the local atmosphere
     */
    boost::function< basic_mathematics::Vector6d( ) > bodyFixedStateFunction_;

    //! Function to determine the angle of attack of the vehicle.
    boost::function< double( ) > angleOfAttackFunction_;

    //! Function to determine the angle of sideslip of the vehicle.
    boost::function< double( ) > angleOfSideslipFunction_;

    //! Function to determine the bank angle of the vehicle.
    boost::function< double( ) > bankAngleFunction_;

    //! Boolean to determine whether to determine vertical <-> aerodynamic frame conversion
    //! when calling update function.
    bool calculateVerticalToAerodynamicFrame_;

};

//! Get a function to transform aerodynamic force from local to propagation frame.
/*!
 *  Get a function to transform aerodynamic force from local to propagation frame. The returned
 *  function takes the force in the accelerationFrame, and returns this force in the
 *  propagationFrame.
 *  \param aerodynamicAngleCalculator Object to calculate aerodynamic angles and rotation
 *  quaternions.
 *  \param accelerationFrame Id of frame in which aerodynamic force is calculated
 *  \param bodyFixedToInertialFrameFunction Function to get the central-body corotating frame to
 *  inertial frame quaternion.
 *  \param propagationFrame Id of frame in which orbit propagation is done.
 *  \return Aerodynamic force conversion function.
 */
boost::function< Eigen::Vector3d( const Eigen::Vector3d& ) >
getAerodynamicForceTransformationFunction(
        const boost::shared_ptr< AerodynamicAngleCalculator > aerodynamicAngleCalculator,
        const AerodynamicsReferenceFrames accelerationFrame,
        const boost::function< Eigen::Quaterniond( ) > bodyFixedToInertialFrameFunction =
        boost::lambda::constant( Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) ),
        const AerodynamicsReferenceFrames propagationFrame = inertial_frame);
} // namespace reference_frames

} // namespace tudat
#endif // TUDAT_AERODYNAMICANGLECALCULATOR_H
