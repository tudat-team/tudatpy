#ifndef THRUSTGUIDANCE_H
#define THRUSTGUIDANCE_H

#include <boost/function.hpp>

#include <Tudat/Astrodynamics/ReferenceFrames/aerodynamicAngleCalculator.h>



namespace tudat
{

namespace basic_astrodynamics
{

class ThrustGuidance
{
public:

    ThrustGuidance( ){ }

    virtual ~ThrustGuidance( ){ }

    virtual void updateThrustDirection( const double time ) = 0;

    virtual Eigen::Vector3d getCurrentThrustDirectionInPropagationFrame( ) = 0;

protected:

};

Eigen::Vector3d getThrustDirectionColinearWithVelocity(
        const basic_mathematics::Vector6d& currentState, const double currentTime, const bool putThrustInOppositeDirection )
{
    return ( ( putThrustInOppositeDirection == 1 ) ? -1.0 : 1.0 ) * ( currentState.segment( 3, 3 ) ).normalize( );
}

Eigen::Vector3d getThrustDirectionColinearWithPosition(
        const basic_mathematics::Vector6d& currentState, const double currentTime, const bool putThrustInOppositeDirection )
{
    return ( ( putThrustInOppositeDirection == 1 ) ? -1.0 : 1.0 ) * ( currentState.segment( 0, 3 ) ).normalize( );
}

class StateBasedThrustGuidance: public ThrustGuidance
{

    StateBasedThrustGuidance(
            boost::function< Eigen::Vector3d( const basic_mathematics::Vector6d&, const double ) > thrustDirectionFunction,
            boost::function< basic_mathematics::Vector6d( ) > bodyStateFunction ):
        ThrustGuidance( ), thrustDirectionFunction_( thrustDirectionFunction ), bodyStateFunction_( bodyStateFunction ){ }

    void updateThrustDirection( const double time )
    {
        currentThrustDirection_ = thrustDirectionFunction_( bodyStateFunction_( ), time );
    }

    Eigen::Vector3d getCurrentThrustDirectionInPropagationFrame( )
    {
        return currentThrustDirection_;
    }

protected:

    boost::function< Eigen::Vector3d( const basic_mathematics::Vector6d&, const double ) > thrustDirectionFunction_;

    boost::function< basic_mathematics::Vector6d( ) > bodyStateFunction_;

    Eigen::Vector3d currentThrustDirection_;
};

class OrientationBasedThrustGuidance: public ThrustGuidance
{
public:

    OrientationBasedThrustGuidance(
            const boost::function< Eigen::Quaterniond( ) > guidanceBaseFrameToPropagationFrameFunction,
            const Eigen::Vector3d thrustFrameThrustDirection = -Eigen::Vector3d::UnitX( ) ):
        ThrustGuidance( ), guidanceBaseFrameToPropagationFrameFunction_( guidanceBaseFrameToPropagationFrameFunction ),
        thrustFrameThrustDirection_(  thrustFrameThrustDirection ){  }

    Eigen::Vector3d getCurrentThrustDirectionInPropagationFrame( )
    {
        return getCurrentThrustFrameToPropagationFrameRotation( ) * thrustFrameThrustDirection_;
    }

    Eigen::Quaterniond getCurrentThrustFrameToPropagationFrameRotation( )
    {
        return currentBaseFrameToPropagationFrame_ * currentThrustFrameToBaseFrameRotation_;
    }

    virtual void updateThrustAngles( const double time ) = 0;

    void updateThrustDirection( const double time )
    {
        currentBaseFrameToPropagationFrame_ = guidanceBaseFrameToPropagationFrameFunction_( );

        updateThrustAngles( );
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
            const boost::function< Eigen::Quaterniond( ) > guidanceThrustFrameToBaseFrameFunction,
            const boost::function< Eigen::Quaterniond( ) > guidanceBaseFrameToPropagationFrameFunction =
            boost::lambda::constant( Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) ),
            const Eigen::Vector3d thrustFrameThrustDirection = -Eigen::Vector3d::UnitX( ) ):
        OrientationBasedThrustGuidance( guidanceBaseFrameToPropagationFrameFunction, thrustFrameThrustDirection ),
    guidanceThrustFrameToBaseFrameFunction_( guidanceThrustFrameToBaseFrameFunction ){  }

    virtual void updateThrustAngles( const double time )
    {
        currentThrustFrameToBaseFrameRotation_ = guidanceThrustFrameToBaseFrameFunction_( );
    }

private:

    boost::function< Eigen::Quaterniond( ) > guidanceThrustFrameToBaseFrameFunction_;


};


}

}


#endif // THRUSTGUIDANCE_H
