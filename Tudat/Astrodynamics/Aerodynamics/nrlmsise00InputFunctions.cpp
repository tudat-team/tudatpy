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
 *
 *    References
 *
 *    Notes
 *
 */
#include <iostream>

#include "Tudat/Astrodynamics/Aerodynamics/nrlmsise00Atmosphere.h"
#include "Tudat/Astrodynamics/Aerodynamics/nrlmsise00InputFunctions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"


//! Tudat library namespace.
namespace tudat
{
namespace aerodynamics
{

//! NRLMSISE00Input function
NRLMSISE00Input nrlmsiseInputFunction(double altitude, double longitude,
                                      double latitude, double time,
                                      tudat::input_output::solar_activity::SolarActivityDataMap& solarActivityMap) {
    // Declare input data class member
    NRLMSISE00Input nrlmsiseInputData ;

    // Julian dates
    double julianDate = tudat::basic_astrodynamics::convertSecondsSinceEpochToJulianDay( time );
    double julianDay = std::floor(julianDate-0.5)+0.5 ;

    // Move not always find, only if day changes
    // Only use datamap with this day + some extra days..
    // Find solar activity data for current date
    SolarActivityDataPtr solarActivity = solarActivityMap[ julianDay ] ;
    if(!solarActivity){
        std::cerr << "Solar activity data could not be found for this date.." << std::endl;
    }

    // Compute julian date at the first of januari
    double julianDate1Jan = tudat::basic_astrodynamics::convertCalendarDateToJulianDay(
                solarActivity->year, 1, 1, 0, 0, 0.0) ;

    nrlmsiseInputData.year = solarActivity->year ; // int
    nrlmsiseInputData.dayOfTheYear = julianDay - julianDate1Jan + 1 ;
    nrlmsiseInputData.secondOfTheDay = time -
            tudat::basic_astrodynamics::convertJulianDayToSecondsSinceEpoch( julianDay,
                                                            tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000 ) ;

    if( solarActivity->fluxQualifier == 1 ){ // requires adjustment
        nrlmsiseInputData.f107 = solarActivity->solarRadioFlux107Adjusted ;
        nrlmsiseInputData.f107a = solarActivity->centered81DaySolarRadioFlux107Adjusted ;
    }
    else{ // no adjustment required
        nrlmsiseInputData.f107 = solarActivity->solarRadioFlux107Observed ;
        nrlmsiseInputData.f107a = solarActivity->centered81DaySolarRadioFlux107Observed ;
    }
    nrlmsiseInputData.apDaily = solarActivity->planetaryEquivalentAmplitudeAverage ;
    nrlmsiseInputData.apVector = EigenToStdVector( solarActivity->planetaryEquivalentAmplitudeVector ) ; // std vector

    // Compute local solar time
    // Hrs since begin of the day at longitude 0 (GMT) + Hrs passed at current longitude
//    nrlmsiseInputData.localSolarTime = 16 ; // use this for unit test
    nrlmsiseInputData.localSolarTime = nrlmsiseInputData.secondOfTheDay/3600.0
            + longitude/(tudat::mathematical_constants::PI/12.0);

    return nrlmsiseInputData;
}

}  // namespace aerodynamics
}  // namespace tudat
