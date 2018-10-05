/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SOFAEARTHORIENTATION_H
#define TUDAT_SOFAEARTHORIENTATION_H

extern "C"
{
    #include <sofa/src/sofa.h>
    #include <sofa/src/sofam.h>
}

#include <map>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"

namespace tudat
{

namespace basic_astrodynamics
{

//! Enum of IAU conventions for Earth rotation.
enum IAUConventions
{
    iau_2000_a,
    iau_2000_b,
    iau_2006
};

}
namespace sofa_interface
{

//! Function to calculate CIP and CIO locator according to requested IAU conventions
/*!
 *  Function to calculate CIP (Celestial Intermediate Pole, typically denoted as X and Y) and CIO locator
 *  (Celestial Intermediate Origin, typically denoted s) according to requested IAU conventions using Sofa implementation.
 *  See Petit et al. (2010), chap. 5 for more details on requested variables.
 *  Note that this function does not include empirical corrections, published as daily values by the IERS.
 *  \param terrestrialTime Time in TT in seconds since referenceJulianDay julianDay.
 *  \param referenceJulianDay Julian day wrt which tt is referenced (default = J2000)
 *  \param precessionNutationTheory IAU conventions that are to be used for calculation
 *  (i.e. determining which Sofa function to call)
 *  \return Pair of first: Vector of entries X, Y (in that order) CIP values and second: CIO locator.
 */
std::pair< Eigen::Vector2d, double > getPositionOfCipInGcrs(
        const double terrestrialTime, const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000,
        const basic_astrodynamics::IAUConventions precessionNutationTheory = basic_astrodynamics::iau_2000_b );

//! Function to calculate GMST according to requested IAU conventions
/*!
 *  Function to calculate GMST (Greenwich Mean Sidereal Time) according to requested IAU conventions
 *  using Sofa implementation.  Note that this function does not include empirical corrections, as published as daily
 *  values by IERS.
 *  \param terrestrialTime Time in TT in seconds since reference julian day.
 *  \param universalTime1 Time in UT1 in seconds since reference julian day.
 *  \param referenceJulianDay Julian day wrt which input times are referenced (default = J2000)
 *  \param iauConvention IAU conventions that are to be used for calculation (i.e. determining which Sofa function to call)
 *  \return Current GMST (normalized to [0,2 pi]).
 */
double calculateGreenwichMeanSiderealTime(
        const double terrestrialTime, const double universalTime1,
        const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000,
        const basic_astrodynamics::IAUConventions iauConvention = basic_astrodynamics::iau_2000_b );

//! Function to calculate ERA (earth rotation angle)
/*!
 *  Function to calculate ERA (earth rotation angle) from current UT1.
 *  The ERA represents one of the Euler angles (between CIO and TIO) for transforming from ITRS to GCRS,
 *  see Petit et al. chap. 5
 *  \param ut1 Time in UT1 in seconds since referenceJulianDay julianDay.
 *  \param referenceJulianDay Julian day wrt which ut1 is referenced (default = J2000)
 *  \return Current Earth rotation angle (normalized to [0,2 pi]).
 */
double calculateEarthRotationAngle(
        const double ut1, const double referenceJulianDay );

//! Function to calculate ERA (earth rotation angle) in templated precision
/*!
 *  Function to calculate ERA (earth rotation angle) from current UT1 in templated precision.
 *  The ERA represents one of the Euler angles (between CIO and TIO) for transforming from ITRS to GCRS,
 *  see Petit et al. chap. 5
 *  \param currentUt1 Time in UT1
 *  \return Current Earth rotation angle (normalized to [0,2 pi]).
 */
template< typename TimeType >
double calculateEarthRotationAngleTemplated( const TimeType currentUt1 );

//! Function to compute the rotation matrix from GCRS to J2000 at epoch
/*!
 * Function to compute the rotation matrix from GCRS to J2000 at epoch, e.g. frame bias
 * \param julianDaysSinceReference Days since reference epoch at which the frame bias is to be determined
 * \param precessionNutationTheory IAU precession theory to use
 * \param referenceJulianDay Reference Julian day for julianDaysSinceReference
 * \return
 */
Eigen::Matrix3d getFrameBias(
        const double julianDaysSinceReference,
        const basic_astrodynamics::IAUConventions precessionNutationTheory = basic_astrodynamics::iau_2006,
        const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000 );
}

}

#endif // TUDAT_SOFAEARTHORIENTATION_H
