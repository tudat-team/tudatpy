/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#ifndef THRUSTGUIDANCE_H
#define THRUSTGUIDANCE_H

#include <iostream>

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include <Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h>

namespace tudat
{

namespace propulsion
{

//! Base class for computing the direction of a thrust acceleration
/*!
 *  Base class for computing the direction of a thrust acceleration. The computation is done directly from some
 *  a priori imposed rule, possibly as a function of dependent, independent or state variables.
 *  A derived class for using the current orientation of the vehicle (as computed by some other method and
 *  retrieved from the body class) may be set using the OrientationBasedThrustGuidance class.
 */
class ThrustDirectionGuidance: public reference_frames::DependentOrientationCalculator
{
public:

    //! Constructor
    /*!
     * Constructor, sets the direction of the force in the body-fixed frame.
     * \param bodyFixedThrustDirection Function returning the direction of the force in the body-fixed frame.
     */
    ThrustDirectionGuidance(
            const boost::function< Eigen::Vector3d( ) > bodyFixedThrustDirection ):
    DependentOrientationCalculator( ){ }

    //! Destructor
    virtual ~ThrustDirectionGuidance( ){ }

    //! Function to get the thrust force direction in the propagation frame
    /*!
     *  Function to get the thrust force direction in the propagation frame. Function is pure virtual and to be
     *  implemented in derived class.
     *  \return Direction of force (as unit vector) expressed in propagation frame.
     */
    virtual Eigen::Vector3d getCurrentThrustDirectionInPropagationFrame( ) = 0;

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
     *  must implement the updateThrustDirection function that obdates the full object.
     *  \param time Time to which object is to be updated.
     */
    void updateCalculator( const double time )
    {
        currentBodyFixedThrustDirection_ = bodyFixedThrustDirection_( );
        updateThrustDirection( time );
    }

protected:

    //! Function to update the thrust direction to the current time.
    /*!
     *  Function to update the thrust direction to the current time. This function is to be implemented in the derived
     *  class.
     *  \param time Time to which object is to be updated.
     */
    virtual void updateThrustDirection( const double time ) = 0;

    //! Function returning the direction of the force in the body-fixed frame.
    boost::function< Eigen::Vector3d( ) > bodyFixedThrustDirection_;

    //! Current direction of the force in the body-fixed frame as set by last call to updateCalculator function.
    Eigen::Vector3d currentBodyFixedThrustDirection_;
};

Eigen::Vector3d getThrustDirectionColinearWithVelocity(
        const basic_mathematics::Vector6d& currentState, const double currentTime, const bool putThrustInOppositeDirection );

Eigen::Vector3d getThrustDirectionColinearWithPosition(
        const basic_mathematics::Vector6d& currentState, const double currentTime, const bool putThrustInOppositeDirection );


Eigen::Vector3d getThrustDirectionFromTimeOnlyFunction(
        const basic_mathematics::Vector6d& currentState, const double currentTime,
        const boost::function< Eigen::Vector3d( const double ) > timeOnlyFunction );

class StateBasedThrustGuidance: public ThrustDirectionGuidance
{
public:
    StateBasedThrustGuidance(
            const boost::function< Eigen::Vector3d( const basic_mathematics::Vector6d&, const double ) > thrustDirectionFunction,
            const boost::function< basic_mathematics::Vector6d( ) > bodyStateFunction,
            const std::string& centralBody,
            const boost::function< Eigen::Vector3d( ) > bodyFixedThrustDirection =
            boost::lambda::constant( Eigen::Vector3d::UnitX( ) ) ):
        ThrustDirectionGuidance( bodyFixedThrustDirection ),
        thrustDirectionFunction_( thrustDirectionFunction ),
        bodyStateFunction_( bodyStateFunction ),
        centralBody_( centralBody ){ }

    Eigen::Vector3d getCurrentThrustDirectionInPropagationFrame( )
    {
        return currentThrustDirection_;
    }

    Eigen::Quaterniond getRotationToGlobalFrame( )
    {
        throw std::runtime_error( "Error, body-fixed frame to propagation frame not yet implemented for StateBasedThrustGuidance." );
    }
     std::string getCentralBody( )
     {
         return centralBody_;
     }

protected:

     void updateThrustDirection( const double time )
     {
         currentThrustDirection_ = thrustDirectionFunction_( bodyStateFunction_( ), time );
     }

    boost::function< Eigen::Vector3d( const basic_mathematics::Vector6d&, const double ) > thrustDirectionFunction_;

    boost::function< basic_mathematics::Vector6d( ) > bodyStateFunction_;

    Eigen::Vector3d currentThrustDirection_;


    std::string centralBody_;

};

class OrientationBasedThrustGuidance: public ThrustDirectionGuidance
{
public:

    OrientationBasedThrustGuidance(
            const boost::function< Eigen::Quaterniond( const double ) > thrustFrameToBaseFrameFunction,
            const boost::function< Eigen::Vector3d( ) > bodyFixedThrustDirection_ =
            boost::lambda::constant( Eigen::Vector3d::UnitX( ) ) ):
        ThrustDirectionGuidance( bodyFixedThrustDirection_ ),
        thrustFrameToBaseFrameFunction( thrustFrameToBaseFrameFunction ){  }

    Eigen::Vector3d getCurrentThrustDirectionInPropagationFrame( )
    {
       return ( getRotationToGlobalFrame( ) * currentBodyFixedThrustDirection_ ).normalized( );
    }

    Eigen::Quaterniond getRotationToGlobalFrame( )
    {
        return currentThrustFrameToBaseFrame_ ;
    }

protected:

    void updateThrustDirection( const double time )
    {
        currentThrustFrameToBaseFrame_ = thrustFrameToBaseFrameFunction( time );
    }

    boost::function< Eigen::Quaterniond( const double ) > thrustFrameToBaseFrameFunction;

    boost::function< Eigen::Vector3d( ) > thrustFrameThrustDirectionFunction_;

    Eigen::Quaterniond currentThrustFrameToBaseFrame_;

};



} // namespace propulsion

} // namespace tudat


#endif // THRUSTGUIDANCE_H
