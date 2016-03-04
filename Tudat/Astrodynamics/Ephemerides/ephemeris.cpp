#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"

namespace tudat
{
namespace ephemerides
{

//! Get state from ephemeris, with state scalar as template type (double specialization).
template<  >
Eigen::Matrix< double, 6, 1 > Ephemeris::getTemplatedStateFromEphemeris( const double& time )
{
    return getCartesianStateFromEphemeris( time, basic_astrodynamics::JULIAN_DAY_ON_J2000 );
}

//! Get state from ephemeris, with state scalar as template type (long double specialization).
template<  >
Eigen::Matrix< long double, 6, 1 > Ephemeris::getTemplatedStateFromEphemeris( const double& time )
{
    return getCartesianLongStateFromEphemeris( time, basic_astrodynamics::JULIAN_DAY_ON_J2000 );
}


}

}

