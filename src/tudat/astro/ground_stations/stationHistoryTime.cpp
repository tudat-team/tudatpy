/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/astro/ground_stations/stationCoordinateHistory.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"
#include <stdexcept>

namespace tudat
{
namespace ground_stations
{
namespace
{
//! Modified Julian day of 1 January of a year.
int januaryFirstMjd( int year )
{
    double origin, mjd;
    iauCal2jd( year, 1, 1, &origin, &mjd );
    return static_cast< int >( mjd );
}

//! Midnight (start) of a modified Julian day, with the day taken as a TAI day.
Time getMidnight( int mjd )
{
    return Time( ( mjd - 51544 ) * 24 - 12, 0.0L );
}
}  // namespace

Time convertSinexCivilUtcToTai( const input_output::SinexCivilUtc& civil )
{
    if( civil.year < 1972 )
    {
        throw std::runtime_error( "Station history: UTC epochs before 1972 are not supported" );
    }
    const int mjd = januaryFirstMjd( civil.year ) + civil.dayOfYear - 1;
    int year, month, day;
    double fraction, taiMinusUtc;
    iauJd2cal( 2400000.5, mjd, &year, &month, &day, &fraction );
    iauDat( year, month, day, 0.0, &taiMinusUtc );
    return getMidnight( mjd ) + ( static_cast< double >( civil.secondsOfDay ) + taiMinusUtc );
}

long double getCalendarDecimalYear( const Time& tai )
{
    int year, month, day;
    double fraction;
    iauJd2cal( 2451545.0 + tai.fullDaysSinceEpoch( ),
               static_cast< double >( tai.secondsSinceNoon( ) / 86400.0L ),
               &year,
               &month,
               &day,
               &fraction );
    const int start = januaryFirstMjd( year );
    const long double daysInYear = januaryFirstMjd( year + 1 ) - start;
    return year + static_cast< long double >( tai - getMidnight( start ) ) / ( daysInYear * 86400.0L );
}
}  // namespace ground_stations
}  // namespace tudat
