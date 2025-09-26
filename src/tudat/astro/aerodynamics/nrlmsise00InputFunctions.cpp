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
std::vector< double > eigenToStlVector( const Eigen::VectorXd& vector )
{
    std::vector< double > stdVector( vector.rows( ) );
    for( int i = 0; i < vector.rows( ); i++ )
    {
        stdVector[ i ] = vector( i );
    }
    return stdVector;
}

//! NRLMSISE00Input function
NRLMSISE00Input nrlmsiseInputFunctionFromMap( const double altitude,
                                              const double longitude,
                                              const double latitude,
                                              const double time,
                                              const tudat::input_output::solar_activity::SolarActivityDataMap& solarActivityMap,
                                              const bool adjustSolarTime,
                                              const double localSolarTime )
{
    return nrlmsiseInputFunction( altitude,
                                  longitude,
                                  latitude,
                                  time,
                                  tudat::input_output::solar_activity::SolarActivityContainer( solarActivityMap ),
                                  adjustSolarTime,
                                  localSolarTime );
}

NRLMSISE00Input nrlmsiseInputFunction( const double altitude,
                                       const double longitude,
                                       const double latitude,
                                       const double time,
                                       const tudat::input_output::solar_activity::SolarActivityContainer& solarActivityContainer,
                                       const bool adjustSolarTime,
                                       const double localSolarTime,
                                       const int geomagneticActivity )
{
    using namespace tudat::input_output::solar_activity;

    // Declare input data class member
    NRLMSISE00Input nrlmsiseInputData;

    // Julian dates
    double julianDate = tudat::basic_astrodynamics::convertSecondsSinceEpochToJulianDay( time, basic_astrodynamics::JULIAN_DAY_ON_J2000 );
    double julianDay = std::floor( julianDate - 0.5 ) + 0.5;

    // Check if solar activity is found for current day.
    SolarActivityDataPtr solarActivity = solarActivityContainer.getSolarActivityData( time );
    SolarActivityDataPtr solarActivityDayOld = solarActivityContainer.getDelayedSolarActivityData( time, 1.0 );

    // Compute julian date at the first of januari
    double julianDate1Jan = tudat::basic_astrodynamics::convertCalendarDateToJulianDay( solarActivity->year, 1, 1, 0, 0, 0.0 );

    nrlmsiseInputData.year = solarActivity->year;  // int
    nrlmsiseInputData.dayOfTheYear = julianDay - julianDate1Jan + 1;
    nrlmsiseInputData.secondOfTheDay = time -
            tudat::basic_astrodynamics::convertJulianDayToSecondsSinceEpoch( julianDay, tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    if( solarActivity->fluxQualifier == 1 )
    {
        // requires adjustment
        nrlmsiseInputData.f107 = solarActivityDayOld->solarRadioFlux107Adjusted;
        nrlmsiseInputData.f107a = solarActivity->centered81DaySolarRadioFlux107Adjusted;
    }
    else
    {
        // no adjustment required
        nrlmsiseInputData.f107 = solarActivityDayOld->solarRadioFlux107Observed;
        nrlmsiseInputData.f107a = solarActivity->centered81DaySolarRadioFlux107Observed;
    }
    nrlmsiseInputData.apDaily = solarActivity->planetaryEquivalentAmplitudeAverage;

    // Create switches vector
    nrlmsiseInputData.switches = std::vector< int >( 24, 1 );
    nrlmsiseInputData.switches[ 0 ] = 0;

    if( geomagneticActivity == -1 )
    {
        // Custom behavior: maybe just use apDaily in vector, or do not fill the vector
        solarActivityContainer.getDelayedApValues( time, nrlmsiseInputData.apVector );
        nrlmsiseInputData.switches[ 9 ] = -1;
    }
    else if( geomagneticActivity == 1 )
    {
        // Example: zero ap vector
        nrlmsiseInputData.apVector = std::vector< double >( 6, 0.0 );
    }

    // Compute local solar time
    // Hrs since begin of the day at longitude 0 (GMT) + Hrs passed at current longitude
    if( adjustSolarTime )
    {
        nrlmsiseInputData.localSolarTime = localSolarTime;
    }
    else
    {
        nrlmsiseInputData.localSolarTime =
                nrlmsiseInputData.secondOfTheDay / 3600.0 + longitude / ( tudat::mathematical_constants::PI / 12.0 );
    }

    return nrlmsiseInputData;
}

}  // namespace aerodynamics
}  // namespace tudat
