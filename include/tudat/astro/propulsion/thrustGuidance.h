/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <iostream>
#include <functional>

#include "tudat/astro/ephemerides/directionBasedRotationalEphemeris.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
//#include "tudat/astro/reference_frames/dependentOrientationCalculator.h"

#include "tudat/basics/basicTypedefs.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/astro/ephemerides/rotationalEphemeris.h"

namespace tudat
{

namespace propulsion
{



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
Eigen::Vector3d getDirectionColinearWithVelocity(
        const std::function< void( Eigen::Vector6d& ) > currentStateFunction,
        const double currentTime, const bool putVectorInOppositeDirection );

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
Eigen::Vector3d getDirectionColinearWithPosition(
        const std::function< void( Eigen::Vector6d& ) > currentStateFunction,
        const double currentTime, const bool putVectorInOppositeDirection );

//! Function to get the force direction from a time-only function.
/*!
 * Function to get the force direction from a time-only function.
 * \param currentTime Current time.
 * \param timeOnlyFunction Function returning unit vector (thrust direction) as a funtion of time.
 * \return Thrust direction.
 */
Eigen::Vector3d getForceDirectionFromTimeOnlyFunction(
        const double currentTime,
        const std::function< Eigen::Vector3d( const double ) > timeOnlyFunction );

class ThrustDirectionWrapper
{
public:
    ThrustDirectionWrapper( ){ }

    virtual ~ThrustDirectionWrapper( ){ }

    void resetCurrentTime( const double currentTime = TUDAT_NAN )
    {
        currentTime_ = currentTime;
        resetDerivedClassCurrentTime( currentTime );
    }

    virtual void resetDerivedClassCurrentTime( const double currentTime = TUDAT_NAN )
    { }

    virtual void update( const double time ) = 0;

    virtual Eigen::Vector3d getCurrentThrustDirection( )
    {
        return currentThrustDirection_;
    }

protected:

    double currentTime_;

    Eigen::Vector3d currentThrustDirection_;

};

class DirectThrustDirectionWrapper: public ThrustDirectionWrapper
{
public:
    DirectThrustDirectionWrapper(
            const std::shared_ptr< ephemerides::DirectionBasedRotationalEphemeris > directionBasedRotationModel,
            const bool inPositiveDirection ):
    directionBasedRotationModel_( directionBasedRotationModel ),
    thrustMultiplier_( inPositiveDirection ? 1.0 : -1.0 ){ }

    virtual ~DirectThrustDirectionWrapper( ){ }


    virtual void resetDerivedClassCurrentTime( const double currentTime = TUDAT_NAN )
    {
        directionBasedRotationModel_->update( currentTime );
    }

    virtual void update( const double time )
    {
        if( time != currentTime_  && time == time )
        {
            currentThrustDirection_ = (
                        thrustMultiplier_ * directionBasedRotationModel_->getCurrentInertialDirection( time ) ).normalized( );
        }
        currentTime_ = time;
    }


protected:

    std::shared_ptr< ephemerides::DirectionBasedRotationalEphemeris > directionBasedRotationModel_;

    double thrustMultiplier_;
};

class OrientationBasedThrustDirectionWrapper: public ThrustDirectionWrapper
{
public:
    OrientationBasedThrustDirectionWrapper(
            const std::function< Eigen::Quaterniond( ) > rotationFunction,
            const std::function< Eigen::Vector3d( ) > bodyFixedThrustDirection ):
        ThrustDirectionWrapper( ),
    rotationFunction_( rotationFunction ), bodyFixedThrustDirection_( bodyFixedThrustDirection ){ }

    virtual ~OrientationBasedThrustDirectionWrapper( ){ }

    virtual void update( const double time )
    {
        if( time != currentTime_ && time == time )
        {
            currentThrustDirection_ = ( rotationFunction_( ) * bodyFixedThrustDirection_( ) ).normalized( );
        }
        currentTime_ = time;

    }

protected:

    const std::function< Eigen::Quaterniond( ) > rotationFunction_;

    const std::function< Eigen::Vector3d( ) > bodyFixedThrustDirection_;
};

////! Class for computing the force direction directly from a function returning the associated unit vector of direction.
//class DirectionBasedForceGuidance: public BodyFixedForceDirectionGuidance
//{
//public:

//    //! Constructor.
//    /*!
//     *  Constructor.
//     *  \param forceDirectionFunction Function returning thrust-direction (represented in the relevant propagation frame)
//     *      as a function of time.
//     *  \param centralBody Name of central body
//     *  \param bodyFixedForceDirection Function returning the unit-vector of the force direction in a body-fixed frame (e.g.
//     *      engine thrust pointing in body-fixed frame).
//     */
//    DirectionBasedForceGuidance(
//            const std::function< Eigen::Vector3d( const double ) > forceDirectionFunction,
//            const std::string& centralBody,
//            const std::function< Eigen::Vector3d( ) > bodyFixedForceDirection =
//            [ ]( ){ return Eigen::Vector3d::UnitX( ); } ):
//        BodyFixedForceDirectionGuidance ( bodyFixedForceDirection ),
//        forceDirectionFunction_( forceDirectionFunction ),
//        centralBody_( centralBody ),
//        hasWarningBeenGiven_( false ){ }

//    //! Function returning the current force direction, as computed by last call to updateCalculator/updateForceDirection.
//    /*!
//     *  Function returning the current force direction, as computed by last call to updateCalculator/updateForceDirection.
//     *  \return Current force direction, expressed in propagation frame.
//     */
//    Eigen::Vector3d getCurrentForceDirectionInPropagationFrame( )
//    {
//        return currentForceDirection_;
//    }

//    //! Function to get the rotation from body-fixed to inertial frame.
//    /*!
//     *  Function to get the rotation from body-fixed to inertial frame.
//     *  \return Current quaternion representing rotation from body-fixed to inertial frame.
//     */
//    Eigen::Quaterniond getRotationToGlobalFrame( )
//    {
//        if ( !hasWarningBeenGiven_ )
//        {
//            hasWarningBeenGiven_ = true;
//            std::cerr << "Warning, body-fixed frame to propagation frame not yet implemented for "
//                         "DirectionBasedForceGuidance. The function will return a unit rotation to represent the "
//                         "rotation from body-fixed to inertial frame." << std::endl;
//        }
//        return Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) );
//    }

//    //! Function to return the name of the central body.
//    /*!
//     * Function to return the name of the central body.
//     * \return Name of the central body.
//     */
//    std::string getCentralBody( )
//    {
//        return centralBody_;
//    }

//protected:

//    //! Function to update the force direction to the current time.
//    /*!
//     *  Function to update the force direction to the current time.
//     *  \param time Time to which object is to be updated.
//     */
//    void updateForceDirection( const double time )
//    {
//        currentForceDirection_ = forceDirectionFunction_( time ).normalized( );
//    }

//    //! Function returning thrust-direction (represented in the relevant propagation frame) as a function of time.
//    std::function< Eigen::Vector3d( const double ) > forceDirectionFunction_;

//    //! Current force direction, as computed by last call to updateCalculator/updateForceDirection.
//    Eigen::Vector3d currentForceDirection_;

//    //! Name of central body
//    std::string centralBody_;

//private:

//    //! Boolean denoting whether user has been warned of incomplete function.
//    bool hasWarningBeenGiven_;

//};

////! Class for computing the force direction using the rotation matrix from the propagation to the body-fixed frame.
//class OrientationBasedForceGuidance: public BodyFixedForceDirectionGuidance
//{
//public:

//    //! Constructor
//    /*!
//     * Constructor
//     * \param bodyFixedFrameToBaseFrameFunction Function returning rotation from the body-fixed frame to the relevant
//     * propagation frame as a function of time.
//     * \param bodyFixedForceDirection Function returning the unit-vector of the force direction in a body-fixed frame (e.g.
//     * engine thrust pointing in body-fixed frame).
//     */
//    OrientationBasedForceGuidance(
//            const std::function< Eigen::Quaterniond( const double ) > bodyFixedFrameToBaseFrameFunction,
//            const std::function< Eigen::Vector3d( ) > bodyFixedForceDirection =
//            [ ]( ){ return  Eigen::Vector3d::UnitX( ); } ):
//        BodyFixedForceDirectionGuidance ( bodyFixedForceDirection ),
//        bodyFixedFrameToBaseFrameFunction_( bodyFixedFrameToBaseFrameFunction ){  }

//    //! Function returning the current force direction, as computed by last call to updateCalculator/updateForceDirection.
//    /*!
//     *  Function returning the current force direction, as computed by last call to updateCalculator/updateForceDirection.
//     *  \return Current force direction, expressed in propagation frame.
//     */
//    Eigen::Vector3d getCurrentForceDirectionInPropagationFrame( )
//    {
//       return ( getRotationToGlobalFrame( ) * currentBodyFixedForceDirection_ ).normalized( );
//    }

//    //! Function to get the rotation from body-fixed to inertial frame.
//    /*!
//     *  Function to get the rotation from body-fixed to inertial frame. For thsi derived class, this simply means returning
//     *  the currentBodyFixedFrameToBaseFrame_ rotation.
//     *  \return The rotation from body-fixed to inertial frame.
//     */
//    Eigen::Quaterniond getRotationToGlobalFrame( )
//    {
//        return currentBodyFixedFrameToBaseFrame_ ;
//    }

//protected:

//    //! Function to update the force direction to the current time.
//    /*!
//     *  Function to update the force direction to the current time.
//     *  \param time Time to which object is to be updated.
//     */
//    void updateForceDirection( const double time )
//    {
//        currentBodyFixedFrameToBaseFrame_ = bodyFixedFrameToBaseFrameFunction_( time );
//    }

//    //! Function returning rotation from the body-fixed frame to the relevant propagation frame as a function of time.
//    std::function< Eigen::Quaterniond( const double ) > bodyFixedFrameToBaseFrameFunction_;

//    //! The rotation from body-fixed to inertial frame.
//    Eigen::Quaterniond currentBodyFixedFrameToBaseFrame_;

//};



} // namespace propulsion

} // namespace tudat


#endif // TUDAT_THRUSTGUIDANCE_H
