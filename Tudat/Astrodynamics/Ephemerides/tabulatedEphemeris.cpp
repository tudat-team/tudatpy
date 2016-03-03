#include "Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h"

namespace tudat
{

namespace ephemerides
{

template< >
basic_mathematics::Vector6d TabulatedCartesianEphemeris< double, double >::getCartesianStateFromEphemeris(
        const double ephemerisTime, const double julianDayAtEpoch )
{
    if( julianDayAtEpoch != julianDayAtEpoch_ )
    {
        throw std::runtime_error(
                    "Error in Tabulated Ephemeris, reference epochs are inconsistent" );
    }

    return interpolator_->interpolate( ephemerisTime );
}

template< >
Eigen::Matrix< long double, 6, 1 > TabulatedCartesianEphemeris< double, double >::getCartesianLongStateFromEphemeris(
        const double secondsSinceEpoch, const double julianDayAtEpoch )
{
    if( julianDayAtEpoch != julianDayAtEpoch_ )
    {
        throw std::runtime_error(
                    "Error in Tabulated Ephemeris, reference epochs are inconsistent" );
    }

    return interpolator_->interpolate( secondsSinceEpoch ).cast< long double >( );
}



template< >
basic_mathematics::Vector6d TabulatedCartesianEphemeris< long double, double >::getCartesianStateFromEphemeris(
        const double ephemerisTime, const double julianDayAtEpoch )
{
    if( julianDayAtEpoch != julianDayAtEpoch_ )
    {
        throw std::runtime_error(
                    "Error in Tabulated Ephemeris, reference epochs are inconsistent" );
    }

    return interpolator_->interpolate( ephemerisTime ).cast< double >( );
}

template< >
Eigen::Matrix< long double, 6, 1 > TabulatedCartesianEphemeris< long double, double >::getCartesianLongStateFromEphemeris(
        const double secondsSinceEpoch, const double julianDayAtEpoch )
{
    if( julianDayAtEpoch != julianDayAtEpoch_ )
    {
        throw std::runtime_error(
                    "Error in Tabulated Ephemeris, reference epochs are inconsistent" );
    }

    return interpolator_->interpolate( secondsSinceEpoch );
}


}

}

