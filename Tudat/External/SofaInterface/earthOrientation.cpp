#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/External/SofaInterface/earthOrientation.h"


namespace tudat
{

namespace sofa_interface
{

//! Function to calculate CIP and CIO locator according to requested IAU conventions
std::pair< Eigen::Vector2d, double > getPositionOfCipInGcrs(
        const double terrestrialTime, const double julianDaysEpochShift, const IAUConventions precessionNutationTheory )
{
    // Declare Sofa function return arguments (by reference)
    double xAngle, yAngle;
    double originLocator;

    // Check for IAU convention and retrieve requested values.
    switch( precessionNutationTheory )
    {
    case iau_2000_a:
        iauXys00a( julianDaysEpochShift, terrestrialTime / physical_constants::JULIAN_DAY,
                   &xAngle, &yAngle, &originLocator );
        break;

    case iau_2000_b:
        iauXys00b( julianDaysEpochShift, terrestrialTime / physical_constants::JULIAN_DAY,
                   &xAngle, &yAngle, &originLocator );
        break;

    case iau_2006:
        iauXys06a( julianDaysEpochShift, terrestrialTime / physical_constants::JULIAN_DAY,
                   &xAngle, &yAngle, &originLocator );
        break;
    default:
        throw std::runtime_error( "Warning, precession nutation theory selection not recongnized" );

    }

    // Set and return requested values.
    Eigen::Vector2d cioPosition;
    cioPosition << xAngle, yAngle;
    return std::pair< Eigen::Vector2d, double >( cioPosition, originLocator );
}

//! Function to calculate GMST according to requested IAU conventions
double calculateGreenwichMeanSiderealTime(
        const double terrestrialTime, const double universalTime1JulianDaysSinceJ2000,
        const double referenceJulianDay, const IAUConventions iauConvention )
{
    // Declare GMST variable
    double gmst = TUDAT_NAN;

    // Check for IAU convention and retrieve requested GMST
    switch( iauConvention )
    {
    case iau_2000_a:
        gmst = iauGmst00( referenceJulianDay, universalTime1JulianDaysSinceJ2000 / physical_constants::JULIAN_DAY,
                          referenceJulianDay, terrestrialTime / physical_constants::JULIAN_DAY );
        break;

    case iau_2000_b:
        gmst = iauGmst00( referenceJulianDay, universalTime1JulianDaysSinceJ2000 / physical_constants::JULIAN_DAY,
                          referenceJulianDay, terrestrialTime / physical_constants::JULIAN_DAY );
        break;

    case iau_2006:
        gmst = iauGmst06( referenceJulianDay, universalTime1JulianDaysSinceJ2000 / physical_constants::JULIAN_DAY,
                          referenceJulianDay, terrestrialTime / physical_constants::JULIAN_DAY );
        break;
    default:
       throw std::runtime_error( "Warning, iau convention for GMST calculation not recongnized" );

    }

    return gmst;

}

//! Function to calculate ERA (earth rotation angle)
double calculateEarthRotationAngle( const double ut1, const double julianDaysEpochShift )
{
    return iauEra00( julianDaysEpochShift, ut1 / physical_constants::JULIAN_DAY );
}

}

}
