/*    Copyright (c) 2010-2014, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      130218    D. Dirkx          File created from personal application.
 *
 *    References
 *
 *    Notes
 *
 */

#include <boost/date_time/gregorian/gregorian.hpp>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/physicalConstants.h>

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"

namespace tudat
{
namespace basic_astrodynamics
{

//! Compute number of seconds since a reference Julian day.
double convertJulianDayToSecondsSinceEpoch( const double julianDay,
                                            const double epochSinceJulianDayZero )
{
    return ( julianDay - epochSinceJulianDayZero ) * physical_constants::JULIAN_DAY;
}

//! Compute Julian day from seconds since reference Julian day epoch.
double convertSecondsSinceEpochToJulianDay( const double secondsSinceEpoch,
                                            const double epochSinceJulianDayZero )
{
    return ( secondsSinceEpoch / physical_constants::JULIAN_DAY + epochSinceJulianDayZero );
}

//! Compute the Julian day from the calendar date and time.
double convertCalendarDateToJulianDay( const int calendarYear,
                                       const int calendarMonth,
                                       const int calendarDay,
                                       const int calendarHour,
                                       const int calendarMinutes,
                                       const double calendarSeconds )
{
    // Calculate julian day of calendar date.
    double julianDay = boost::gregorian::date( calendarYear, calendarMonth, calendarDay ).julian_day( );

    //Compute day fraction
    const double dayFraction = static_cast< double >( calendarHour ) / 24.0 +
                               static_cast< double >( calendarMinutes ) / ( 24.0 * 60.0 ) +
                               calendarSeconds / ( 24.0 * 3600.0 );

    // Compute Julian day by adding day fraction and subtracting 0.5 to reference to midnight
    // instead of noon..
    return julianDay + dayFraction - 0.5;
}

} // namespace basic_astrodynamics
} // namespace tudat
