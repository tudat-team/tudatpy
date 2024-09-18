/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/aerodynamics/nrlmsise00Atmosphere.h"
#include "tudat/astro/aerodynamics/nrlmsise00InputFunctions.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/basic_astro/timeConversions.h"


// Tudat library namespace.
namespace tudat
{
namespace aerodynamics
{

//! Function to convert Eigen::VectorXd to std::vector<double>
std::vector< double >  eigenToStlVector(
        const Eigen::VectorXd& vector)
{
    std::vector< double > stdVector( vector.rows( ) );
    for( int i = 0; i < vector.rows( ); i++)
    {
        stdVector[ i ] = vector( i );
    }
    return stdVector;
}

//! NRLMSISE00Input function
NRLMSISE00Input nrlmsiseInputFunction( const double altitude, const double longitude,
                                       const double latitude, const double time,
                                       const tudat::input_output::solar_activity::SolarActivityDataMap& solarActivityMap,
                                       const bool adjustSolarTime,
                                       const double localSolarTime ) {
    using namespace tudat::input_output::solar_activity;

    // Declare input data class member
    NRLMSISE00Input nrlmsiseInputData;

    // Julian dates
    double julianDate = tudat::basic_astrodynamics::convertSecondsSinceEpochToJulianDay(
                time, basic_astrodynamics::JULIAN_DAY_ON_J2000 );
    double julianDay = std::floor( julianDate - 0.5 ) + 0.5;

    // Check if solar activity is found for current day.
    SolarActivityDataPtr solarActivity;

    if( solarActivityMap.count( julianDay ) == 0 )
    {
        auto activityIterator = solarActivityMap.lower_bound( julianDay );
        if( ( julianDay - activityIterator->first > 30 ) )
        {
            throw std::runtime_error( "Error when retrieving solar activity data at JD" + std::to_string( julianDate ) +
                                      ", most recent data is more than 30 days old" );
        }
        solarActivity = activityIterator->second;
    }
    else
    {
        solarActivity = solarActivityMap.at( julianDay );
    }


    // Compute julian date at the first of januari
    double julianDate1Jan = tudat::basic_astrodynamics::convertCalendarDateToJulianDay(
                solarActivity->year, 1, 1, 0, 0, 0.0 );

    nrlmsiseInputData.year = solarActivity->year; // int
    nrlmsiseInputData.dayOfTheYear = julianDay - julianDate1Jan + 1;
    nrlmsiseInputData.secondOfTheDay = time -
            tudat::basic_astrodynamics::convertJulianDayToSecondsSinceEpoch( julianDay,
                                                            tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    if( solarActivity->fluxQualifier == 1 )
    { // requires adjustment
        nrlmsiseInputData.f107 = solarActivity->solarRadioFlux107Adjusted;
        nrlmsiseInputData.f107a = solarActivity->centered81DaySolarRadioFlux107Adjusted;
    }
    else
    { // no adjustment required
        nrlmsiseInputData.f107 = solarActivity->solarRadioFlux107Observed;
        nrlmsiseInputData.f107a = solarActivity->centered81DaySolarRadioFlux107Observed;
    }
    nrlmsiseInputData.apDaily = solarActivity->planetaryEquivalentAmplitudeAverage;
    nrlmsiseInputData.apVector = eigenToStlVector( solarActivity->planetaryEquivalentAmplitudeVector );

    // Compute local solar time
    // Hrs since begin of the day at longitude 0 (GMT) + Hrs passed at current longitude
    if( adjustSolarTime )
    {
        nrlmsiseInputData.localSolarTime = localSolarTime;
    }
    else
    {
        nrlmsiseInputData.localSolarTime = nrlmsiseInputData.secondOfTheDay / 3600.0
                + longitude / ( tudat::mathematical_constants::PI / 12.0 );
    }

    return nrlmsiseInputData;
}

}  // namespace aerodynamics
}  // namespace tudat
