#include "Tudat/Astrodynamics/BasicAstrodynamics/thrustGuidance.h"

namespace tudat
{

namespace basic_astrodynamics
{

Eigen::Vector3d getThrustDirectionColinearWithVelocity(
        const basic_mathematics::Vector6d& currentState, const double currentTime, const bool putThrustInOppositeDirection )
{
    return ( ( putThrustInOppositeDirection == 1 ) ? -1.0 : 1.0 ) * ( currentState.segment( 3, 3 ) ).normalized( );
}

Eigen::Vector3d getThrustDirectionColinearWithPosition(
        const basic_mathematics::Vector6d& currentState, const double currentTime, const bool putThrustInOppositeDirection )
{
    return ( ( putThrustInOppositeDirection == 1 ) ? -1.0 : 1.0 ) * ( currentState.segment( 0, 3 ) ).normalized( );
}

Eigen::Vector3d getThrustDirectionFromTimeOnlyFunction(
        const basic_mathematics::Vector6d& currentState, const double currentTime,
                 const boost::function< Eigen::Vector3d( const double ) > timeOnlyFunction )
{
    return timeOnlyFunction( currentTime ).normalized( );
}

}

}
