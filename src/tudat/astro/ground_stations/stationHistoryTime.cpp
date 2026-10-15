/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/astro/ground_stations/stationHistoryTrajectory.h"
#include "tudat/astro/ground_stations/earthDeformationTime.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace tudat
{
namespace ground_stations
{
namespace
{
void checkSofa( int status )
{
    if( status != 0 )
    {
        throw std::runtime_error( "IGN nominal UTC: invalid or dubious date under the installed SOFA leap-second table" );
    }
}
bool leapYear( int year )
{
    return year % 4 == 0 && ( year % 100 != 0 || year % 400 == 0 );
}

//! A UTC calendar day, with its leap-second information.
struct CivilDay {
    //! Calendar date and modified Julian day.
    int year;
    int month;
    int day;
    int mjd;

    //! TAI - UTC (seconds) at the start of the day, and the start of the day as TAI.
    double offset;
    Time midnightTai;
};
CivilDay civilDay( int mjd )
{
    CivilDay result;
    result.mjd = mjd;
    double fraction;
    checkSofa( iauJd2cal( 2400000.5, mjd, &result.year, &result.month, &result.day, &fraction ) );
    checkSofa( iauDat( result.year, result.month, result.day, 0.0, &result.offset ) );
    result.midnightTai = Time( ( mjd - 51544 ) * 24 - 12, static_cast< long double >( result.offset ) );
    return result;
}
int januaryFirstMjd( int year )
{
    double origin, mjd;
    checkSofa( iauCal2jd( year, 1, 1, &origin, &mjd ) );
    return static_cast< int >( mjd );
}
}  // namespace

Time convertSinexCivilUtcToTai( const input_output::SinexCivilUtc& civil )
{
    if( civil.year < 1972 || civil.dayOfYear < 1 || civil.dayOfYear > ( leapYear( civil.year ) ? 366 : 365 ) || civil.secondsOfDay < 0 ||
        civil.secondsOfDay > 86400 )
    {
        throw std::runtime_error( "IGN nominal UTC: invalid civil date or unsupported pre-1972 UTC frequency drift" );
    }
    const auto day = civilDay( januaryFirstMjd( civil.year ) + civil.dayOfYear - 1 );
    const auto next = civilDay( day.mjd + 1 );
    const Time time = day.midnightTai + static_cast< double >( civil.secondsOfDay );
    if( time >= next.midnightTai )
    {
        throw std::runtime_error( "IGN nominal UTC: civil inserted-second label does not exist" );
    }
    return time;
}

StationHistoryTime getIgnNominalUtcStationTimeFromTai( const Time& tai, std::optional< double > taiRateWithRespectToTdb )
{
    if( !std::isfinite( static_cast< long double >( tai ) ) )
    {
        throw std::runtime_error( "IGN nominal UTC: nonfinite TAI epoch" );
    }
    double utc1, utc2, fraction, origin, mjd;
    int year, month, date;
    checkSofa(
            iauTaiutc( 2451545.0 + tai.fullDaysSinceEpoch( ), static_cast< double >( tai.secondsSinceNoon( ) / 86400.0L ), &utc1, &utc2 ) );
    checkSofa( iauJd2cal( utc1, utc2, &year, &month, &date, &fraction ) );
    checkSofa( iauCal2jd( year, month, date, &origin, &mjd ) );
    auto day = civilDay( static_cast< int >( mjd ) );

    // SOFA supplies a civil-day candidate only. Precise TAI comparisons repair
    // rounding of a sub-nanosecond sample across midnight and detect leap seconds.
    if( tai < day.midnightTai )
    {
        day = civilDay( day.mjd - 1 );
    }
    else if( tai >= civilDay( day.mjd + 1 ).midnightTai )
    {
        day = civilDay( day.mjd + 1 );
    }
    if( day.year < 1972 )
    {
        throw std::runtime_error( "IGN nominal UTC: pre-1972 UTC frequency drift is unsupported" );
    }
    const long double seconds = static_cast< long double >( tai - day.midnightTai );
    if( seconds < 0.0L || seconds >= 86400.0L )
    {
        throw std::runtime_error( "IGN nominal UTC source coordinates cannot represent an inserted second" );
    }
    StationHistoryTime result;
    result.tai = tai;
    result.mjdDay = day.mjd;
    result.mjdFraction = seconds / 86400.0L;
    result.calendarYear = day.year;
    const int daysInYear = leapYear( day.year ) ? 366 : 365;
    result.yearFraction = ( day.mjd - januaryFirstMjd( day.year ) + result.mjdFraction ) / daysInYear;

    // A year fraction cannot retain the full sub-femtosecond Time resolution;
    // keep a rounded pre-boundary fraction inside its genuine calendar year.
    if( result.yearFraction == 1.0L )
    {
        result.yearFraction = std::nextafter( 1.0L, 0.0L );
    }
    result.convention = "IGN nominal UTC";
    if( taiRateWithRespectToTdb )
    {
        if( !std::isfinite( *taiRateWithRespectToTdb ) || *taiRateWithRespectToTdb <= 0.0 )
        {
            throw std::runtime_error( "IGN nominal UTC: invalid TAI derivative with respect to TDB" );
        }
        if( tai == day.midnightTai && civilDay( day.mjd - 1 ).offset != day.offset )
        {
            throw std::runtime_error( "IGN nominal UTC source rate is undefined at a UTC offset step" );
        }
        result.mjdRate = *taiRateWithRespectToTdb / 86400.0;
        result.decimalYearRate = *taiRateWithRespectToTdb / ( daysInYear * 86400.0 );
    }
    return result;
}

StationHistoryTime getIgnNominalUtcStationTime( const input_output::SinexCivilUtc& civil )
{
    if( civil.secondsOfDay == 86400 )
    {
        throw std::runtime_error( "IGN nominal UTC reference/event epoch cannot be an inserted-second label" );
    }
    return getIgnNominalUtcStationTimeFromTai( convertSinexCivilUtcToTai( civil ) );
}
StationHistoryTime getIgnNominalUtcStationTimeFromTdb( const Time& tdb,
                                                       std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter,
                                                       bool rates )
{
    if( !converter )
    {
        converter = basic_astrodynamics::getGeocentricTimeConverter( );
    }
    return getIgnNominalUtcStationTimeFromTai(
            converter->getTaiTime( basic_astrodynamics::tdb_scale, tdb ),
            rates ? std::optional< double >( converter->getGeocentricTtRateWithRespectToTdb( tdb ) ) : std::nullopt );
}
}  // namespace ground_stations
}  // namespace tudat
