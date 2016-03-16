
#include <cmath>
#include <iostream>
#include <iomanip>

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/geodeticCoordinateConversions.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/External/SofaInterface/sofaTimeConversions.h"

namespace tudat
{

namespace sofa_interface
{

//! Function to calculate number of leap seconds from UTC input
double getDeltaAtFromUtc( const double utcInJulianDays )
{
    // Declare return variables (by reference) from sofa calendar function.
    int year, month, day;
    double fractionOfDay;

    // Get calendar date and check feasibility of input.
    if( iauJd2cal( basic_astrodynamics::JULIAN_DAY_ON_J2000, utcInJulianDays, &year, &month, &day, &fractionOfDay ) != 0 )
    {
        throw std::runtime_error(  "Provided julian date too small to convert to calendar date" +
                                   boost::lexical_cast< std::string >( year ) + " " +
                                   boost::lexical_cast< std::string >( month ) + " " +
                                   boost::lexical_cast< std::string >( day ) );
    }


    // Get number of leap seconds and check feasibility of calculation
    double deltaAt;
    int deltaAtReturn = iauDat( year, month, day, fractionOfDay, &deltaAt);
    if( deltaAtReturn != 0 )
    {
        throw std::runtime_error(  "Provided caledar date cannot properly give Delta AT" +
                                   boost::lexical_cast< std::string >( year ) + " " +
                                   boost::lexical_cast< std::string >( month ) + " " +
                                   boost::lexical_cast< std::string >( day ) );
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
        throw std::runtime_error(  "Provided julian date too small to convert to calendar date" );
    }

    // Estimate number of leap seconds (by assuming TAI = UTC for Sofa input) and check feasibility of calculation
    double deltaAt;
    int deltaAtReturn = iauDat( year, month, day, fractionOfDay, &deltaAt);
    if( deltaAtReturn != 0 )
    {
        throw std::runtime_error( "Provided caledar date cannot properly give Delta AT" );
    }

    // Reperform calculation with converted utc time and check consistency with previous calculation.
    double utcInJulianDays = taiInJulianDays - deltaAt / physical_constants::JULIAN_DAY;
    double deltaAtCheck = getDeltaAtFromUtc( utcInJulianDays );
    if( deltaAt != deltaAtCheck )
    {
        deltaAt--;
        throw std::runtime_error( "Warning, Delta TAI calculation encountered error, iteration not yet implemented" );
    }


    return deltaAt;
}

//! Function to calculate difference between TDB and TT
double getTDBminusTT( const double tdbTime, const double universalTimeFractionOfDay, const double stationLongitude,
                      const double distanceFromSpinAxis, const double distanceFromEquatorialPlane )
{
    return iauDtdb( basic_astrodynamics::JULIAN_DAY_ON_J2000, tdbTime / physical_constants::JULIAN_DAY,
                    universalTimeFractionOfDay, stationLongitude,
                    distanceFromSpinAxis / 1000.0, distanceFromEquatorialPlane / 1000.0 );
}

//! Function to calculate difference between TDB and TT.
double getTDBminusTT( const double tdbTime, const double universalTimeFractionOfDay,
                      const Eigen::Vector3d& stationCartesianPosition )
{
    double stationLongitude = std::atan2( stationCartesianPosition( 1 ),
                                          stationCartesianPosition( 0 ) );
    return getTDBminusTT( tdbTime, universalTimeFractionOfDay, stationLongitude,
                          std::sqrt( stationCartesianPosition.x( ) * stationCartesianPosition.x( ) +
                                     stationCartesianPosition.y( ) * stationCartesianPosition.y( ) ),
                          stationCartesianPosition.z( ) );
}

//! Function to calculate difference between TDB and TT.
double getTDBminusTT( const double ttOrTdbSinceJ2000, const double stationLongitude, const double distanceFromSpinAxis,
                      const double distanceFromEquatorialPlane )
{
    // Calculate current TAI (approximately if input is in TDB)
    double tai = basic_astrodynamics::convertTTtoTAI< double >( ttOrTdbSinceJ2000 );

    // Calculate current UT1 (by assuming it equal to UTC)
    double ut1 = static_cast< double >( convertTAItoUTC< double >( tai ) );
    double ut1FractionOfDay = ( ut1 / physical_constants::JULIAN_DAY + 0.5 ) -
            static_cast< double >( std::floor( ut1 / physical_constants::JULIAN_DAY  + 0.5 ) );

    // Calculate and return difference (introducing addition approximation if input is in TT, by assuming TDB is equal to TT)
    return getTDBminusTT( ttOrTdbSinceJ2000, ut1FractionOfDay, stationLongitude, distanceFromSpinAxis,
                          distanceFromEquatorialPlane );
}

//! Function to calculate difference between TDB and TT.
double getTDBminusTT( const double ttOrTdbSinceJ2000, const Eigen::Vector3d& stationCartesianPosition )
{
    double stationLongitude = std::atan2( stationCartesianPosition( 1 ),
                                          stationCartesianPosition( 0 ) );
    return getTDBminusTT( ttOrTdbSinceJ2000, stationLongitude,
                          std::sqrt( stationCartesianPosition.x( ) * stationCartesianPosition.x( ) +
                                     stationCartesianPosition.y( ) * stationCartesianPosition.y( ) ),
                          stationCartesianPosition.z( ) );
}


} // namespace sofa_interfaces

} // namespace tudat
