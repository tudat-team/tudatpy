/*    Copyright (c) 2010-2015, Delft University of Technology
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

#ifndef TUDAT_TIME_CONVERSIONS_H
#define TUDAT_TIME_CONVERSIONS_H

namespace tudat
{
namespace basic_astrodynamics
{

//! Julian day at J2000.
/*!
 * Julian day at J2000, i.e. 01-01-2000, at 12:00 (in TT).
 */
const static double JULIAN_DAY_ON_J2000 = 2451545.0;

//! Julian day at Modified Julain Date 0.
/*!
 * Julian day at Modified Julain Date 0, i.e. Nov 17, 1858, 00:00.
 */
const static double JULIAN_DAY_AT_0_MJD = 2400000.5;

//! Compute number of seconds since a reference Julian day.
/*!
 * Computes the number of seconds since a reference Julian day from
 * Julian day.
 * \param julianDay Julian day at which number of seconds since epoch is to be determined.
 * \param epochSinceJulianDayZero Epoch in Julian day.
 * \return Number of seconds since epoch.
 */
double convertJulianDayToSecondsSinceEpoch( const double julianDay,
                                            const double epochSinceJulianDayZero );

//! Compute Julian day from seconds since reference Julian day epoch.
/*!
 * Computes the Julian day bsaed on seconds since reference Julian day epoch provided.
 * \param secondsSinceEpoch Seconds since epoch.
 * \param epochSinceJulianDayZero Epoch in Julian day.
 * \return Number of Julian days since epoch.
 */
double convertSecondsSinceEpochToJulianDay(
                                      const double secondsSinceEpoch,
                                      const double epochSinceJulianDayZero = JULIAN_DAY_ON_J2000 );

//! Compute Julian day from given date and time
/*!
  * Computes the Julian day from given year, month, day, hour, minutes, seconds as used in everyday
  * life. The function uses the internal calcualtions of the boost::date_time::gregorian class.
  *
  * \param calendarYear Year of the standard calendar in years.
  * \param calendarMonth Month of the standard calendar in months.
  * \param calendarDay Day of the standard calendar in days.
  * \param calendarHour Hour of the time of this day in hours.
  * \param calendarMinutes Minutes of the time of this day in minutes.
  * \param calendarSeconds Seconds of the time of this day in seconds.
  */
double convertCalendarDateToJulianDay( const int calendarYear,
                                       const int calendarMonth,
                                       const int calendarDay,
                                       const int calendarHour,
                                       const int calendarMinutes,
                                       const double calendarSeconds );

//! Dummy time conversion function (time remains constant)
/*!
 * Dummy time conversion function (time remains constant)
 * \param inputTime Time given as input.
 * \return Same as inputTime.
 */
template< typename TimeType >
TimeType doDummyTimeConversion( const TimeType inputTime )
{
    return inputTime;
}


} // namespace basic_astrodynamics
} // tudat

#endif // TUDAT_TIME_CONVERSIONS_H
