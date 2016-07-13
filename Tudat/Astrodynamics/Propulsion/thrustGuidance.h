#ifndef THRUSTGUIDANCE_H
#define THRUSTGUIDANCE_H

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
            boost::function< Eigen::Vector3d( const basic_mathematics::Vector6d&, const double ) > thrustDirectionFunction,
            boost::function< basic_mathematics::Vector6d( ) > bodyStateFunction ):
        ThrustDirectionGuidance( ), thrustDirectionFunction_( thrustDirectionFunction ), bodyStateFunction_( bodyStateFunction ){ }

    void updateCalculator( const double time )
    {
        currentThrustDirection_ = thrustDirectionFunction_( bodyStateFunction_( ), time );
    }

    Eigen::Vector3d getCurrentThrustDirectionInPropagationFrame( )
    {
        return currentThrustDirection_;
    }

    Eigen::Quaterniond getRotationToGlobalFrame( )
    {
        throw std::runtime_error( "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB" );
    }


protected:

    boost::function< Eigen::Vector3d( const basic_mathematics::Vector6d&, const double ) > thrustDirectionFunction_;

    boost::function< basic_mathematics::Vector6d( ) > bodyStateFunction_;

    Eigen::Vector3d currentThrustDirection_;
};

class OrientationBasedThrustGuidance: public ThrustDirectionGuidance
{
public:

    OrientationBasedThrustGuidance(
            const boost::function< Eigen::Quaterniond( ) > guidanceBaseFrameToPropagationFrameFunction,
            const Eigen::Vector3d thrustFrameThrustDirection = -Eigen::Vector3d::UnitX( ) ):
        ThrustDirectionGuidance( ), guidanceBaseFrameToPropagationFrameFunction_( guidanceBaseFrameToPropagationFrameFunction ),
        thrustFrameThrustDirection_(  thrustFrameThrustDirection ){  }

    Eigen::Vector3d getCurrentThrustDirectionInPropagationFrame( )
    {
        return getRotationToGlobalFrame( ) * thrustFrameThrustDirection_;
    }

    Eigen::Quaterniond getRotationToGlobalFrame( )
    {
        return currentBaseFrameToPropagationFrame_ * currentThrustFrameToBaseFrameRotation_;
    }

    virtual void updateThrustAngles( const double time ) = 0;

    void updateCalculator( const double time )
    {
        currentBaseFrameToPropagationFrame_ = guidanceBaseFrameToPropagationFrameFunction_( );

        updateThrustAngles( time );
    }

protected:

    boost::function< Eigen::Quaterniond( ) > guidanceBaseFrameToPropagationFrameFunction_;

    Eigen::Quaterniond currentThrustFrameToBaseFrameRotation_;

    Eigen::Quaterniond currentBaseFrameToPropagationFrame_;

    Eigen::Vector3d thrustFrameThrustDirection_;

};

class DirectOrientationBasedThrustGuidance: public OrientationBasedThrustGuidance
{
public:

    DirectOrientationBasedThrustGuidance(
            const boost::function< Eigen::Quaterniond( const double ) > guidanceThrustFrameToBaseFrameFunction,
            const boost::function< Eigen::Quaterniond( ) > guidanceBaseFrameToPropagationFrameFunction =
            boost::lambda::constant( Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) ),
            const Eigen::Vector3d thrustFrameThrustDirection = -Eigen::Vector3d::UnitX( ) ):
        OrientationBasedThrustGuidance( guidanceBaseFrameToPropagationFrameFunction, thrustFrameThrustDirection ),
    guidanceThrustFrameToBaseFrameFunction_( guidanceThrustFrameToBaseFrameFunction ){  }

    virtual void updateThrustAngles( const double time )
    {
        currentThrustFrameToBaseFrameRotation_ = guidanceThrustFrameToBaseFrameFunction_( time );
    }

private:

    boost::function< Eigen::Quaterniond(  const double ) > guidanceThrustFrameToBaseFrameFunction_;


};


}

}


#endif // THRUSTGUIDANCE_H
