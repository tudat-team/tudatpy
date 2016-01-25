#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

namespace tudat
{

namespace physical_constants
{

template< >
double getJulianDay< double >( )
{
    return JULIAN_DAY;
}

template< >
long double getJulianDay< long double >( )
{
    return JULIAN_DAY_LONG;
}

template< >
double getSpeedOfLight< double >( )
{
    return SPEED_OF_LIGHT;
}

template< >
long double getSpeedOfLight< long double >( )
{
    return LONG_SPEED_OF_LIGHT;
}

template< >
double getLgTimeRateTerm< double >( )
{
    return LG_TIME_RATE_TERM;
}

template< >
long double getLgTimeRateTerm< long double >( )
{
    return LG_TIME_RATE_TERM_LONG;
}

template< >
double getLbTimeRateTerm< double >( )
{
    return LB_TIME_RATE_TERM;
}

template< >
long double getLbTimeRateTerm< long double >( )
{
    return LB_TIME_RATE_TERM_LONG;
}

template< >
double getLcTimeRateTerm< double >( )
{
    return LC_TIME_RATE_TERM;
}

template< >
long double getLcTimeRateTerm< long double >( )
{
    return LC_TIME_RATE_TERM_LONG;
}



}

}
