#ifndef EARTHORIENTATION_H
#define EARTHORIENTATION_H

extern "C"
{
#include "sofa/src/sofa.h"
#include "sofa/src/sofam.h"
}

#include <map>
#include <iostream>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"


namespace tudat
{

enum IAUConventions
{
    iau_2000_a,
    iau_2000_b,
    iau_2006
};

namespace sofa_interface
{

//! Function to calculate CIP and CIO locator according to requested IAU conventions
/*!
 *  Function to calculate CIP (Celestial Intermediate Pole, typically denoted as X and Y) and CIO (Celestial Intermediate Origin, typically denoted s)
 *  locator according to requested IAU conventions using Sofa implementation. See Petit et al. (2010), chap. 5 for more details on requested variables.
 *  Note that this function does not include empirical corrections, as published as daily values by IERS.
 *  \param terrestrialTime Time in TT in seconds since julianDaysEpochShift julianDay.
 *  \param julianDaysEpochShift Julian day wrt which tt is referenced (default = J2000)
 *  \param precessionNutationTheory IAU conventions that are to be used for calculation (i.e. determining which Sofa function to call)
 *  \return Pair of first: Vector of entries X, Y (in that order) CIP values and second: CIO locator.
 */
std::pair< Eigen::Vector2d, double > getPositionOfCipInGcrs(
        const double terrestrialTime, const double julianDaysEpochShift = basic_astrodynamics::JULIAN_DAY_ON_J2000,
        const IAUConventions precessionNutationTheory = iau_2000_b );

//! Function to calculate GMST according to requested IAU conventions
/*!
 *  Function to calculate GMST (Greenwich Mean Sidereal Time) according to requested IAU conventions using Sofa implementation.
 *  Note that this function does not include empirical corrections, as published as daily values by IERS.
 *  \param terrestrialTimeJulianDaysSinceJ2000 Time in TT in seconds since J2000
 *  \param universalTime1JulianDaysSinceJ2000 Time in UT1 in seconds since J2000
 *  \param iauConvention IAU conventions that are to be used for calculation (i.e. determining which Sofa function to call)
 *  \return Current GMST
 */
double calculateGreenwichMeanSiderealTime(
        const double terrestrialTimeJulianDaysSinceJ2000, const double universalTime1JulianDaysSinceJ2000,
        const IAUConventions iauConvention = iau_2000_b );

//! Function to calculate ERA (earth rotation angle)
/*!
 *  Function to calculate GMST ERA (earth rotation angle) from current UT1. ERA represents one of the Euler angles (between CIO and TIO)
 *  for transforming from ITRS to GCRS, see Petit et al. chap. 5
 *  \param ut1 Time in UT1 in seconds since julianDaysEpochShift julianDay.
 *  \param julianDaysEpochShift Julian day wrt which ut1 is referenced (default = J2000)
 *  \return Current Earth rotation angle
 */
double calculateEarthRotationAngle(
        const double ut1, const double julianDaysEpochShift = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

}

}

#endif // EARTHORIENTATION_H
