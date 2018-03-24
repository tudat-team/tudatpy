/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#ifndef TUDAT_THRUSTGUIDANCE_H
#define TUDAT_THRUSTGUIDANCE_H
#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"
#include "Tudat/Astrodynamics/ReferenceFrames/dependentOrientationCalculator.h"
#include "Tudat/Basics/basicTypedefs.h"
namespace tudat
{

namespace propulsion
{

//! Base class for computing the direction of a given force/acceleration
/*!
 *  Base class for computing the direction of a given force/acceleration. The computation is done directly from some
 *  a priori imposed rule, possibly as a function of dependent, independent or state variables.
 *  A derived class for using the current orientation of the vehicle (as computed by some other method and
 *  retrieved from the body class) may be set using the OrientationBasedForceGuidance class.
 */
class BodyFixedForceDirectionGuidance : public reference_frames::DependentOrientationCalculator
{
public:

    //! Constructor
    /*!
     * Constructor, sets the direction of the force in the body-fixed frame.
     * \param bodyFixedForceDirection Function returning the direction of the force in the body-fixed frame.
     */
    BodyFixedForceDirectionGuidance (
            const boost::function< Eigen::Vector3d( ) > bodyFixedForceDirection ):
    DependentOrientationCalculator( ), bodyFixedForceDirection_( bodyFixedForceDirection ){ }

    //! Destructor
    virtual ~BodyFixedForceDirectionGuidance ( ){ }

    //! Function to get the force/acceleration direction in the propagation frame
    /*!
     *  Function to get the force/acceleration direction in the propagation frame. Function is pure virtual and to be
     *  implemented in derived class.
     *  \return Direction of force (as unit vector) expressed in propagation frame.
     */
    virtual Eigen::Vector3d getCurrentForceDirectionInPropagationFrame( ) = 0;

    //! Function to compute the rotation from the propagation frame to the body-fixed frame
    /*!
     *  Function to compute the rotation from the propagation frame to the body-fixed frame using the algorithm implemented
     *  in the derived class. The derived class must implement the function for the inverse rotation
     *  (getRotationToGlobalFrame).
     *  \return Quaternion that provides the rotation from the propagation frame to the body-fixed frame.
     */
    Eigen::Quaterniond getRotationToLocalFrame( )
    {
       return getRotationToGlobalFrame( ).inverse( );
    }

    //! Function to update the object to the current time.
    /*!
     *  Function to update the object to the current time. This function updates only this base class, derived class
     *  must implement the updateForceDirection function that obdates the full object.
     *  \param time Time to which object is to be updated.
     */
    void updateCalculator( const double time )
    {
        currentBodyFixedForceDirection_ = ( bodyFixedForceDirection_( ) ).normalized( );
        updateForceDirection( time );
    }

protected:

    //! Function to update the force/acceleration direction to the current time.
    /*!
     *  Function to update the force/acceleration direction to the current time. This function is to be implemented in the
     *  derived class.
     *  \param time Time to which object is to be updated.
     */
    virtual void updateForceDirection( const double time ) = 0;

    //! Function returning the direction of the force in the body-fixed frame.
    boost::function< Eigen::Vector3d( ) > bodyFixedForceDirection_;

    //! Current direction of the force in the body-fixed frame as set by last call to updateCalculator function.
    Eigen::Vector3d currentBodyFixedForceDirection_;
};


//! Function to get the unit vector colinear with velocity segment of a translational state.
/*!
 * Function to get the unit vector colinear with velocity segment of a translational state.
 * \param currentStateFunction Function returning (by reference) translational Cartesian state from which the unit velocity
 * vector is to be retrieved.
 * \param currentTime Time at which computation is to be done (not used here; included for interface compatibility).
 * \param putForceInOppositeDirection Boolean denoting whether the output vector should be in opposite (if true) or same
 * direction (if false) as velocity segment of currentState
 * \return Unit vector colinear with velocity segment of currentState.
 */
Eigen::Vector3d getForceDirectionColinearWithVelocity(
        const boost::function< void( Eigen::Vector6d& ) > currentStateFunction,
        const double currentTime, const bool putForceInOppositeDirection );

//! Function to get the unit vector colinear with position segment of a translational state.
/*!
 * Function to get the unit vector colinear with position segment of a translational state.
 * \param currentStateFunction Function returning (by reference) translational Cartesian state from which the unit position
 * vector is to be retrieved.
 * \param currentTime Time at which computation is to be done (not used here; included for interface compatibility).
 * \param putForceInOppositeDirection Boolean denoting whether the output vector should be in opposite (if true) or same
 * direction (if false) as position segment of current state
 * \return Unit vector colinear with position segment of current state.
 */
Eigen::Vector3d getForceDirectionColinearWithPosition(
        const boost::function< void( Eigen::Vector6d& ) > currentStateFunction,
        const double currentTime, const bool putForceInOppositeDirection );

//! Function to get the force direction from a time-only function.
/*!
 * Function to get the force direction from a time-only function.
 * \param currentTime Current time.
 * \param timeOnlyFunction Function returning unit vector (thrust direction) as a funtion of time.
 * \return Thrust direction.
 */
Eigen::Vector3d getForceDirectionFromTimeOnlyFunction(
        const double currentTime,
        const boost::function< Eigen::Vector3d( const double ) > timeOnlyFunction );

//! Class for computing the force direction directly from a function returning the associated unit vector of direction.
class DirectionBasedForceGuidance: public BodyFixedForceDirectionGuidance
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param forceDirectionFunction Function returning thrust-direction (represented in the relevant propagation frame)
     * as a function of time.
     * \param centralBody Name of central body
     * \param bodyFixedForceDirection Function returning the unit-vector of the force direction in a body-fixed frame (e.g.
     * engine thrust pointing in body-fixed frame).
     */
    DirectionBasedForceGuidance(
            const boost::function< Eigen::Vector3d( const double ) > forceDirectionFunction,
            const std::string& centralBody,
            const boost::function< Eigen::Vector3d( ) > bodyFixedForceDirection =
            boost::lambda::constant( Eigen::Vector3d::UnitX( ) ) ):
        BodyFixedForceDirectionGuidance ( bodyFixedForceDirection ),
        forceDirectionFunction_( forceDirectionFunction ),
        centralBody_( centralBody ){ }

    //! Function returning the current force direction, as computed by last call to updateCalculator/updateForceDirection.
    /*!
     *  Function returning the current force direction, as computed by last call to updateCalculator/updateForceDirection.
     *  \return Current force direction, expressed in propagation frame.
     */
    Eigen::Vector3d getCurrentForceDirectionInPropagationFrame( )
    {
        return currentForceDirection_;
    }

    //! Function to get the rotation from body-fixed to inertial frame.
    /*!
     *  Function to get the rotation from body-fixed to inertial frame. NOT YET IMPLEMENTED IN THIS DERIVED CLASS.
     *  \return NOT YET IMPLEMENTED IN THIS DERIVED CLASS.
     */
    Eigen::Quaterniond getRotationToGlobalFrame( )
    {
        throw std::runtime_error( "Error, body-fixed frame to propagation frame not yet implemented for DirectionBasedForceGuidance." );
    }

    //! Function to return the name of the central body.
    /*!
     * Function to return the name of the central body.
     * \return Name of the central body.
     */
    std::string getCentralBody( )
    {
        return centralBody_;
    }

protected:

    //! Function to update the force direction to the current time.
    /*!
     *  Function to update the force direction to the current time.
     *  \param time Time to which object is to be updated.
     */
    void updateForceDirection( const double time )
    {
        currentForceDirection_ = forceDirectionFunction_( time ).normalized( );
    }

    //! Function returning thrust-direction (represented in the relevant propagation frame) as a function of time.
    boost::function< Eigen::Vector3d( const double ) > forceDirectionFunction_;

    //! Current force direction, as computed by last call to updateCalculator/updateForceDirection.
    Eigen::Vector3d currentForceDirection_;

    //! Name of central body
    std::string centralBody_;

};

//! Class for computing the force direction using the rotation matrix from the propagation to the body-fixed frame.
class OrientationBasedForceGuidance: public BodyFixedForceDirectionGuidance
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param bodyFixedFrameToBaseFrameFunction Function returning rotation from the body-fixed frame to the relevant
     * propagation frame as a function of time.
     * \param bodyFixedForceDirection Function returning the unit-vector of the force direction in a body-fixed frame (e.g.
     * engine thrust pointing in body-fixed frame).
     */
    OrientationBasedForceGuidance(
            const boost::function< Eigen::Quaterniond( const double ) > bodyFixedFrameToBaseFrameFunction,
            const boost::function< Eigen::Vector3d( ) > bodyFixedForceDirection =
            boost::lambda::constant( Eigen::Vector3d::UnitX( ) ) ):
        BodyFixedForceDirectionGuidance ( bodyFixedForceDirection ),
        bodyFixedFrameToBaseFrameFunction_( bodyFixedFrameToBaseFrameFunction ){  }

    //! Function returning the current force direction, as computed by last call to updateCalculator/updateForceDirection.
    /*!
     *  Function returning the current force direction, as computed by last call to updateCalculator/updateForceDirection.
     *  \return Current force direction, expressed in propagation frame.
     */
    Eigen::Vector3d getCurrentForceDirectionInPropagationFrame( )
    {
       return ( getRotationToGlobalFrame( ) * currentBodyFixedForceDirection_ ).normalized( );
    }

    //! Function to get the rotation from body-fixed to inertial frame.
    /*!
     *  Function to get the rotation from body-fixed to inertial frame. For thsi derived class, this simply means returning
     *  the currentBodyFixedFrameToBaseFrame_ rotation.
     *  \return The rotation from body-fixed to inertial frame.
     */
    Eigen::Quaterniond getRotationToGlobalFrame( )
    {
        return currentBodyFixedFrameToBaseFrame_ ;
    }

protected:

    //! Function to update the force direction to the current time.
    /*!
     *  Function to update the force direction to the current time.
     *  \param time Time to which object is to be updated.
     */
    void updateForceDirection( const double time )
    {
        currentBodyFixedFrameToBaseFrame_ = bodyFixedFrameToBaseFrameFunction_( time );
    }

    //! Function returning rotation from the body-fixed frame to the relevant propagation frame as a function of time.
    boost::function< Eigen::Quaterniond( const double ) > bodyFixedFrameToBaseFrameFunction_;

    //! The rotation from body-fixed to inertial frame.
    Eigen::Quaterniond currentBodyFixedFrameToBaseFrame_;

};



} // namespace propulsion

} // namespace tudat


#endif // TUDAT_THRUSTGUIDANCE_H
