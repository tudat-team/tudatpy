#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"

namespace tudat
{
namespace ephemerides
{

basic_mathematics::Vector6d getDifferenceBetweenStates(
        const boost::function< basic_mathematics::Vector6d( const double ) > stateFunction,
        const boost::function< basic_mathematics::Vector6d( const double ) > centralBodyStateFunction,
        const double time )
{
    return stateFunction( time ) - centralBodyStateFunction( time );
}


Eigen::Vector3d calculateAccelerationFromEphemeris(
        const boost::function< basic_mathematics::Vector6d( const double ) > stateFunction, const double time, const double timePerturbation )
{
    return ( stateFunction( time + timePerturbation ) - stateFunction( time - timePerturbation ) ).segment( 3, 3 ) / ( 2.0 * timePerturbation );
}


template<  >
Eigen::Matrix< double, 6, 1 > Ephemeris::getTemplatedStateFromEphemeris( const double& time )
{
    return getCartesianStateFromEphemeris( time, basic_astrodynamics::JULIAN_DAY_ON_J2000 );
}

template<  >
Eigen::Matrix< long double, 6, 1 > Ephemeris::getTemplatedStateFromEphemeris( const double& time )
{
    return getCartesianLongStateFromEphemeris( time, basic_astrodynamics::JULIAN_DAY_ON_J2000 );
}


}

}

