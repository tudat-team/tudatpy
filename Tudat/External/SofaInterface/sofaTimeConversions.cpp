
#include <cmath>
#include <iostream>
#include <iomanip>

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/External/SofaInterface/sofaTimeConversions.h"

namespace tudat
{

namespace sofa_interface
{

template< >
double getTTMinusTai< double >( )
{
    return 32.184;
}

template< >
long double getTTMinusTai< long double >( )
{
    return 32.184L;
}

double approximateConvertTTtoTDB( const double tt, const double earthMeanAnomaly )
{
    return tt + 0.001657  * sin( earthMeanAnomaly );
}


//! Function to calculate number of leap seconds from UTC input
double getDeltaAtFromUtc( const double utcInJulianDays )
{
    // Declare return variables (by reference) from sofa calendar function.
    int year, month, day;
    double fractionOfDay;

    // Get calendar date and check feasibility of input.
    if( iauJd2cal( basic_astrodynamics::JULIAN_DAY_ON_J2000, utcInJulianDays, &year, &month, &day, &fractionOfDay ) != 0 )
    {
        std::cerr<<"Provided julian date too small to convert to calendar date 1"<<year<<" "<<month<<" "<<day<<std::endl;
    }

    // Get number of leap seconds and check feasibility of calculation
    double deltaAt;
    int deltaAtReturn = iauDat( year, month, day, fractionOfDay, &deltaAt);
    if( deltaAtReturn != 0 )
    {
        std::cerr<<"Provided caledar date cannot properly give Delta AT 1 2"<<year<<" "<<month<<" "<<day<<std::endl;
    }

    return deltaAt;
}

//! Function to calculate number of leap seconds from TAI input
double getDeltaAtFromTai( const double taiInJulianDays )
{
    // Declare return variables (by reference) from sofa calendar function.
    int year, month, day;
    double fractionOfDay;

    // Get calendar date and check feasibility of input.
    if( iauJd2cal( basic_astrodynamics::JULIAN_DAY_ON_J2000, taiInJulianDays, &year, &month, &day, &fractionOfDay ) != 0 )
    {
        std::cerr<<"Provided julian date too small to convert to calendar date"<<std::endl;
    }

    // Estimate number of leap seconds (by assuming TAI = UTC for Sofa input) and check feasibility of calculation
    double deltaAt;
    int deltaAtReturn = iauDat( year, month, day, fractionOfDay, &deltaAt);
    if( deltaAtReturn != 0 )
    {
        std::cerr<<"Provided caledar date cannot properly give Delta AT 2"<<std::endl;
    }

    // Reperform calculation with converted utc time and check consistency with previous calculation.
    double utcInJulianDays = taiInJulianDays - deltaAt / physical_constants::JULIAN_DAY;
    double deltaAtCheck = getDeltaAtFromUtc( utcInJulianDays );
    if( deltaAt != deltaAtCheck )
    {
        deltaAt--;
        std::cout<<"Warning, Delta TAI calculation encountered error, iteration not yet implemented"<<std::endl;
    }


    return deltaAt;
}

//! Function to calculate difference between TDB and TT
template< >
double getTDBminusTT< double >( const double ephemerisTime, const double universalTimeFractionOfDay, const double stationLongitude,
                                const double distanceFromSpinAxis, const double distanceFromEquatorialPlane )
{
    return iauDtdb( basic_astrodynamics::JULIAN_DAY_ON_J2000, ephemerisTime / physical_constants::JULIAN_DAY,
                    universalTimeFractionOfDay, stationLongitude,
                    distanceFromSpinAxis / 1000.0, distanceFromEquatorialPlane / 1000.0 );
}

} // namespace sofa_interfaces

} // namespace tudat
