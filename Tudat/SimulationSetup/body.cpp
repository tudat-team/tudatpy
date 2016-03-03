#include "Tudat/SimulationSetup/body.h"

namespace tudat
{

namespace simulation_setup
{

template< >
void Body::setTemplatedState( const Eigen::Matrix< double, 6, 1 >& state )
{
    setState( state );
}

template< >
void Body::setTemplatedState( const Eigen::Matrix< long double, 6, 1 >& state )
{
    setLongState( state );
}

template<  >
Eigen::Matrix< double, 6, 1 > Body::getTemplatedStateInBaseFrameFromEphemeris( const double& time )
{
    return getStateInBaseFrameFromEphemeris( time );
}

template<  >
Eigen::Matrix< long double, 6, 1 > Body::getTemplatedStateInBaseFrameFromEphemeris( const double& time )
{
    return getLongStateInBaseFrameFromEphemeris( time );
}

template< >
void Body::setTemplatedStateFromEphemeris< double, double >( const double& time )
{
    currentState = getTemplatedStateInBaseFrameFromEphemeris< double, double >( time );
}

template< >
void Body::setTemplatedStateFromEphemeris< long double, double >( const double& time )
{
    currentLongState = getTemplatedStateInBaseFrameFromEphemeris< long double, double >( time );
    currentState = currentLongState.cast< double >( );
}

}

}

