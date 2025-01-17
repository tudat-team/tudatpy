/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SOFATIMECONVERSIONS_H
#define TUDAT_SOFATIMECONVERSIONS_H

#include <vector>

#include <Eigen/Core>

extern "C"
{
    #include <sofa/sofa.h>
    #include <sofa/sofam.h>
}

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/timeConversions.h"

namespace tudat
{

namespace sofa_interface
{


//! Function to calculate number of leap seconds from UTC input
/*!
 *  Function to calculate number of leap seconds from UTC (Universal Coordinated Time) input.
 *  \param utcInJulianDays Time in UTC; in julianDays since J2000.
 *  \return Number of leap seconds at requested time.
 */
double getDeltaAtFromUtc( const double utcInJulianDays );

//! Function to calculate number of leap seconds from TAI input
/*!
 *  Function to calculate number of leap seconds from UTC (International Atomic Time) input. Note: this function will not
 *  work properly near the introduction of a leap second (specifically if TAI is less than the number of leap seconds from
 *  the leap second introduction) . In case of incorrect conversion, an error message is produced.
 *  \param taiInJulianDays Time in TAI; in seconds since J2000.
 *  \return Number of leap seconds at requested time.
 */
double getDeltaAtFromTai( const double taiInJulianDays );


//! Function to convert TAI to UTC
/*!
 *  Function to convert TAI (International Atomic Time) to UTC (Universal Coordinated Time) by subtracting bias (
 *  i.e. number of leap seconds since reference epoch) as determined by Sofa. Note: this function will not
 *  work properly near the introduction of a leap second,
 *  \sa getDeltaAtFromTai
 *  \param taiSeconds Time in TAI; in seconds since J2000.
 *  \return Time in UTC; in seconds since J2000
 */
template< typename TimeType >
TimeType convertTAItoUTC( const TimeType taiSeconds )
{
    // Retrieve number of leap seconds from Sofa, assuming TAI=UTC
    double deltaAt = getDeltaAtFromUtc( static_cast< double >( taiSeconds ) / physical_constants::JULIAN_DAY );

    // Update correction in case conversion is close to leap second introduction.
    TimeType utc = taiSeconds - static_cast< TimeType >( deltaAt );
    deltaAt = getDeltaAtFromUtc( static_cast< double >( utc ) / physical_constants::JULIAN_DAY );

    // Return converted time
    return taiSeconds - static_cast< TimeType >( deltaAt );
}

//! Function to convert UTC to TAI
/*!
 *  Function to convert UTC (Universal Coordinated Time) to TAI (International Atomic Time) by adding bias
 *  (i.e. number of leap seconds since reference epoch) as determined by Sofa.
 *  \param utcSeconds Time in UTC; in seconds since J2000
 *  \return Time in TAI; in seconds since J2000
 */
template< typename TimeType >
TimeType convertUTCtoTAI( const TimeType utcSeconds )
{
    double deltaAt = getDeltaAtFromUtc( static_cast< double >( utcSeconds ) / physical_constants::JULIAN_DAY );
    return utcSeconds + static_cast< TimeType >( deltaAt );
}


//! Function to convert TT to UTC
/*!
 *  Function to convert TT (Terrestrial Time) to UTC (Universal Coordinated Time)  by determining bias
 *  (i.e. number of leap seconds  since reference epoch and constant bias of TAI wrt TT) as determined by Sofa.
 *  Note: this function will not work properly near the introduction of a leap second,
 *  \sa getDeltaAtFromTai
 *  \param tt Time in TT; in seconds since J2000
 *  \return Time in UTC; in seconds since J2000
 */
template< typename TimeType >
TimeType convertTTtoUTC( const TimeType tt )
{
    return convertTAItoUTC( basic_astrodynamics::convertTTtoTAI( tt ) );
}

//! Function to convert UTC to TT
/*!
 *  Function to convert UTC (Universal Coordinated Time) to TT (Terrestrial Time) by determining bias
 *  (i.e. number of leap seconds  since reference epoch and constant bias of TAI wrt TT) as determined by Sofa.
 *  \param utc Time in UTC in seconds since J2000
 *  \return Time in TT; in seconds since J2000
 */
template< typename TimeType >
TimeType convertUTCtoTT( const TimeType utc )
{
    return basic_astrodynamics::convertTAItoTT( convertUTCtoTAI( utc ) );
}


////! Function to convert UTC to UT1
///*!
// * Function to convert UTC (Universal Coordinated Time) to UT1 (Universal Time version 1) by looking up the difference
// * in tables from SOFA. UTC closely follows (within 0.9 s) UT1, but for some applications an accurate conversion between
// * the two is necessary.
// * @param utc Time in UTC in seconds since J2000
// * @return Time in UT1 in seconds since J2000
// */
//double convertUTCtoUT1( const double utc );

//! Function to calculate difference between TDB and TT.
/*!
 *  Function to calculate difference between TDB (Dynamical Barycentric Time) and TT (Terrestrial Time),
 *  from Sofa function, which is based  on decomposition of time ephemeris generated by Fairhead an Bretagnon (1990).
 *  Note that the time ephemeris should be generated by the direct solution of differential equations of proper time
 *  using custom ephemerides that, if very high consistency is required.
 *  \param tdbTime TDB in seconds since J2000.
 *  \param universalTimeFractionOfDay UT1 in fraction of current day.
 *  \param stationLongitude Longitude of point on Earth where difference is to be calculated
 *  \param distanceFromSpinAxis Distance from Earth spin axis where difference is to be calculated
 *  \param distanceFromEquatorialPlane Distance from Earth equatorial plane where difference is to be calculated
 *  \return Difference between TDB and TT at requested position and TDB
 */
double getTDBminusTT( const double tdbTime, const double universalTimeFractionOfDay, const double stationLongitude,
                      const double distanceFromSpinAxis, const double distanceFromEquatorialPlane );

//! Function to calculate difference between TDB and TT.
/*!
 *  Function to calculate difference between TDB (Dynamical Barycentric Time) and TT (Terrestrial Time),
 *  from Sofa function, which is based  on decomposition of time ephemeris generated by Fairhead an Bretagnon (1990).
 *  Note that the time ephemeris should be generated by the direct solution of differential equations of proper time
 *  using custom ephemerides that, if very high consistency is required.
 *  \param tdbTime TDB in seconds since J2000.
 *  \param universalTimeFractionOfDay UT1 in fraction of current day.
 *  \param stationCartesianPosition Earth-fixed Cartesian position of evaluation point
 *  \return Difference between TDB and TT at requested position and TDB
 */
double getTDBminusTT( const double tdbTime, const double universalTimeFractionOfDay,
                      const Eigen::Vector3d& stationCartesianPosition );

//! Function to calculate difference between TDB and TT
/*!
 *  Function to calculate difference between TDB (Dynamical Barycentric Time) and TT (Terrestrial Time), from Sofa function,
 *  which is based on decomposition of time ephemeris generated by Fairhead an Bretagnon (1990).
 *  This function uses a simplification for the input to the algorithm: it is assumed that UTC and UT1 coincide for the
 *  purposes of the time difference calculation (as needed for the topocentric term of the time
 *  conversion). Also, the first input argument may be TDB or TT and the difference between the two is neglected in the
 *  internal calculations (TDB is formally required). These simplifications make this function unsuitable for moderately
 *  high precision applications.
 *  \param ttOrTdbSinceJ2000 TDB or TT in seconds since J2000.
 *  \param stationLongitude Longitude of point on Earth where difference is to be calculated
 *  \param distanceFromSpinAxis Distance from Earth spin axis where difference is to be calculated
 *  \param distanceFromEquatorialPlane Distance from Earth equatorial plane where difference is to be calculated
 *  \return Difference between TDB and TT at requested position and TDB
 */
double getTDBminusTT( const double ttOrTdbSinceJ2000, const double stationLongitude, const double distanceFromSpinAxis,
                      const double distanceFromEquatorialPlane );

//! Function to calculate difference between TDB and TT
/*!
 *  Function to calculate difference between TDB (Dynamical Barycentric Time) and TT (Terrestrial Time), from Sofa function,
 *  which is based on decomposition of time ephemeris generated by Fairhead an Bretagnon (1990).
 *  This function uses a simplification for the input to the algorithm: it is assumed that UTC and UT1 coincide for the
 *  purposes of the time difference calculation (as needed for the topocentric term of the time
 *  conversion). Also, the first input argument may be TDB or TT and the difference between the two is neglected in the
 *  internal calculations (TDB is formally required). These simplifications make this function unsuitable for moderately
 *  high precision applications.
 *  \param ttOrTdbSinceJ2000 TDB or TT in seconds since J2000.
 *  \param stationCartesianPosition Earth-fixed Cartesian position of evaluation point
 *  \return Difference between TDB and TT at requested position and TDB
 */
double getTDBminusTT( const double ttOrTdbSinceJ2000, const Eigen::Vector3d& stationCartesianPosition );

//! Determine the number of seconds that have passed in the current year.
/*!
 * Determine the number of seconds that have passed in the current year.
 * @param secondsSinceEpoch Seconds since reference epoch.
 * @param epochSinceJulianDayZero Reference epoch in julian days.
 * @return Seconds since beginning of calendar year
 */
template< typename TimeScalarType = double >
TimeScalarType convertSecondsSinceEpochToSecondsOfYear(
        const TimeScalarType secondsSinceEpoch,
        const TimeScalarType epochSinceJulianDayZero = basic_astrodynamics::getJulianDayOnJ2000< TimeScalarType >( ) )
{
    // Get calendar date
    int year, month, day;
    double fractionOfDay;
    iauJd2cal( static_cast< double >( epochSinceJulianDayZero ),
               static_cast< double >( secondsSinceEpoch ) / physical_constants::JULIAN_DAY, &year, &month, &day, &fractionOfDay );

    // Get time since start of the calendar year
    double timeStartOfYear = basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
            year, 1, 1, 0, 0, 0.0, epochSinceJulianDayZero ) * physical_constants::JULIAN_DAY;

    return secondsSinceEpoch - static_cast< TimeScalarType >( timeStartOfYear );
}



template< typename TimeType >
TimeType convertTTtoTDB( const TimeType ttTime, const Eigen::Vector3d& stationCartesianPosition )
{
    return ttTime + static_cast< TimeType >( sofa_interface::getTDBminusTT( ttTime, stationCartesianPosition ) );
}


template< typename TimeType >
TimeType convertTDBtoTT( const TimeType tdbTime, const Eigen::Vector3d& stationCartesianPosition )
{
    return tdbTime - static_cast< TimeType >( sofa_interface::getTDBminusTT( tdbTime, stationCartesianPosition ) );
}




} // namespace sofa_interfaces

} // namespace tudat

#endif // TUDAT_SOFATIMECONVERSIONS_H
