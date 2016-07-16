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

class ThrustDirectionGuidance: public reference_frames::DependentOrientationCalculator
{
public:

    ThrustDirectionGuidance( ){ }

    virtual ~ThrustDirectionGuidance( ){ }

    virtual Eigen::Vector3d getCurrentThrustDirectionInPropagationFrame( ) = 0;

    Eigen::Quaterniond getRotationToLocalFrame( )
    {
       return getRotationToGlobalFrame( ).inverse( );
    }

protected:

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
            const boost::function< Eigen::Vector3d( const double ) > bodyFixedThrustDirection =
            boost::lambda::constant( Eigen::Vector3d::UnitX( ) ) ):
        ThrustDirectionGuidance( ),
        thrustDirectionFunction_( thrustDirectionFunction ),
        bodyStateFunction_( bodyStateFunction ),
        centralBody_( centralBody ),
        bodyFixedThrustDirection_( bodyFixedThrustDirection ){ }

    void updateCalculator( const double time )
    {
        currentThrustDirection_ = thrustDirectionFunction_( bodyStateFunction_( ), time );
        currentBodyFixedThrustDirection_ = bodyFixedThrustDirection_( time );

    }

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

    boost::function< Eigen::Vector3d( const basic_mathematics::Vector6d&, const double ) > thrustDirectionFunction_;

    boost::function< basic_mathematics::Vector6d( ) > bodyStateFunction_;

    Eigen::Vector3d currentThrustDirection_;

    Eigen::Vector3d currentBodyFixedThrustDirection_;


    std::string centralBody_;

    boost::function< Eigen::Vector3d( const double ) > bodyFixedThrustDirection_;
};

class OrientationBasedThrustGuidance: public ThrustDirectionGuidance
{
public:

    OrientationBasedThrustGuidance(
            const boost::function< Eigen::Quaterniond( const double ) > thrustFrameToBaseFrameFunction,
            const Eigen::Vector3d thrustFrameThrustDirection = Eigen::Vector3d::UnitX( ) ):
        ThrustDirectionGuidance( ), thrustFrameToBaseFrameFunction( thrustFrameToBaseFrameFunction ),
        thrustFrameThrustDirection_(  thrustFrameThrustDirection ){  }

    Eigen::Vector3d getCurrentThrustDirectionInPropagationFrame( )
    {
       return ( getRotationToGlobalFrame( ) * thrustFrameThrustDirection_ ).normalized( );
    }

    Eigen::Quaterniond getRotationToGlobalFrame( )
    {
        return currentThrustFrameToBaseFrame_ * currentThrustFrameToBodyFixedFrame.inverse( );
    }

    virtual void updateThrustAngles( const double time ) = 0;

    void updateCalculator( const double time )
    {
        currentThrustFrameToBaseFrame_ = thrustFrameToBaseFrameFunction( time );

        updateThrustAngles( time );
    }

protected:

    boost::function< Eigen::Quaterniond( const double ) > thrustFrameToBaseFrameFunction;

    Eigen::Quaterniond currentThrustFrameToBodyFixedFrame;

    Eigen::Quaterniond currentThrustFrameToBaseFrame_;

    Eigen::Vector3d thrustFrameThrustDirection_;

};

class DirectOrientationBasedThrustGuidance: public OrientationBasedThrustGuidance
{
public:

    DirectOrientationBasedThrustGuidance(
            const boost::function< Eigen::Quaterniond( const double ) > thrustFrameToBaseFrameFunction,
            const boost::function< Eigen::Quaterniond( const double ) > thrustFrameToBodyFixedFrameFunction =
            boost::lambda::constant( Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) ),
            const Eigen::Vector3d thrustFrameThrustDirection = Eigen::Vector3d::UnitX( ) ):
        OrientationBasedThrustGuidance( thrustFrameToBaseFrameFunction, thrustFrameThrustDirection ),
    thrustFrameToBodyFixedFrameFunction_( thrustFrameToBodyFixedFrameFunction ){  }

    virtual void updateThrustAngles( const double time )
    {
        currentThrustFrameToBodyFixedFrame = thrustFrameToBodyFixedFrameFunction_( time );
    }

private:

    boost::function< Eigen::Quaterniond( const double ) > thrustFrameToBodyFixedFrameFunction_;

};


} // namespace propulsion

} // namespace tudat


#endif // THRUSTGUIDANCE_H
