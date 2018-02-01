/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

namespace tudat
{

namespace physical_constants
{

//! Function to get the length of a Julian day in seconds, with double precision.
template< >
double getJulianDay< double >( )
{
    return JULIAN_DAY;
}

//! Function to get the length of a Julian day in seconds, with long double precision.
template< >
long double getJulianDay< long double >( )
{
    return JULIAN_DAY_LONG;
}

//! Function to get the length of a Julian year in days, with double precision.
template< >
double getJulianYearInDays< double >( )
{
    return JULIAN_YEAR_IN_DAYS;
}

//! Function to get the length of a Julian year in days, with long double precision.
template< >
long double getJulianYearInDays< long double >( )
{
    return JULIAN_YEAR_IN_DAYS_LONG;
}

//! Function to get the speed of light in meters per second with double precision.
template< >
double getSpeedOfLight< double >( )
{
    return SPEED_OF_LIGHT;
}

//! Function to get the speed of light in meters per second with long double precision.
template< >
long double getSpeedOfLight< long double >( )
{
    return SPEED_OF_LIGHT_LONG;
}

//! Function to get the TCG and TT relative rate difference with double precision.
template< >
double getLgTimeRateTerm< double >( )
{
    return LG_TIME_RATE_TERM;
}

//! Function to get the TCG and TT relative rate difference with long double precision.
template< >
long double getLgTimeRateTerm< long double >( )
{
    return LG_TIME_RATE_TERM_LONG;
}

//! Function to get the TCB and TDB relative rate difference with double precision.
template< >
double getLbTimeRateTerm< double >( )
{
    return LB_TIME_RATE_TERM;
}

//! Function to get the TCB and TDB relative rate difference with long double precision.
template< >
long double getLbTimeRateTerm< long double >( )
{
    return LB_TIME_RATE_TERM_LONG;
}

} // namespace physical_constants

} // namespace tudat
