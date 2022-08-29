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
#include "tudat/astro/system_models/engineModel.h"

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

class ThrustDirectionCalculator
{
public:
    ThrustDirectionCalculator( ):
    currentTime_( TUDAT_NAN ){ }

    virtual ~ThrustDirectionCalculator( ){ }

    void resetCurrentTime( )
    {
        currentTime_ = TUDAT_NAN;
        resetDerivedClassCurrentTime( );
    }

    virtual void resetDerivedClassCurrentTime( ){ }

    virtual void update( const double time ) = 0;


    virtual Eigen::Vector3d getInertialThrustDirection(
            const std::shared_ptr< system_models::EngineModel > engineModel ) = 0;

protected:

    double currentTime_;

};

class DirectThrustDirectionCalculator: public ThrustDirectionCalculator
{
public:
    DirectThrustDirectionCalculator(
            const std::shared_ptr< ephemerides::DirectionBasedRotationalEphemeris > directionBasedRotationModel ):
        directionBasedRotationModel_( directionBasedRotationModel ),
    currentQuaterionTime_( TUDAT_NAN ){ }

    virtual ~DirectThrustDirectionCalculator( ){ }

    void resetDerivedClassCurrentTime( )
    {
        directionBasedRotationModel_->resetCurrentTime( );
        currentQuaterionTime_ = TUDAT_NAN;
    }

    void update( const double time )
    {
        if( time != currentTime_  && time == time )
        {
            currentInertialDirection_ = directionBasedRotationModel_->getCurrentInertialDirection( time );
        }
        currentTime_ = time;
    }


    void updateQuaternion( const double time );


    Eigen::Vector3d getInertialThrustDirection(
            const std::shared_ptr< system_models::EngineModel > engineModel );


protected:

    std::shared_ptr< ephemerides::DirectionBasedRotationalEphemeris > directionBasedRotationModel_;

    Eigen::Vector3d currentInertialDirection_;

    Eigen::Quaterniond currentRotationToBaseFrame_;

    double currentQuaterionTime_;
};

class OrientationBasedThrustDirectionCalculator: public ThrustDirectionCalculator
{
public:
    OrientationBasedThrustDirectionCalculator(
            const std::function< Eigen::Quaterniond( ) > rotationFunction ):
        ThrustDirectionCalculator( ),
        rotationFunction_( rotationFunction ){ }

    virtual ~OrientationBasedThrustDirectionCalculator( ){ }

    virtual void update( const double time )
    {
        if( time != currentTime_ && time == time )
        {
            currentRotation_ = rotationFunction_( );
        }
        currentTime_ = time;
    }


    Eigen::Vector3d getInertialThrustDirection(
            const std::shared_ptr< system_models::EngineModel > engineModel )
    {
        return ( currentRotation_ * engineModel->getBodyFixedThrustDirection( ) ).normalized( );
    }

    Eigen::Quaterniond getCurrentRotation( )
    {
        return currentRotation_;
    }

protected:

    const std::function< Eigen::Quaterniond( ) > rotationFunction_;

    Eigen::Quaterniond currentRotation_;
};



} // namespace propulsion

} // namespace tudat


#endif // TUDAT_THRUSTGUIDANCE_H
