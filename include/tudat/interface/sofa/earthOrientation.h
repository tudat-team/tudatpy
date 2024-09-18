/*    Copyright (c) 2010-2019, Delft University of Technology
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
    #include <sofa/sofa.h>
    #include <sofa/sofam.h>
}

#include <map>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/timeConversions.h"

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
Eigen::Vector3d getPositionOfCipInGcrs(
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


//! Function to calculate the equation of equinoxes at the given epoch according to requested IAU conventions
/*!
 * Function to calculate the equation of equinoxes at the given epoch according to requested IAU conventions
 * using a SOFA implementation. This equation gives the angular difference between the mean and true equinox
 * of date. This can in turn be used to e.g. convert GMST to GAST (and vice versa) or convert a state in the
 * TEME frame to the TOD frame.
 * @param terrestrialTime Time in TT in seconds since reference julian day.
 * @param referenceJulianDay Julian day wrt which input times are referenced (default = J2000).
 * @param iauConvention IAU conventions that are to be used for calculation (i.e. determining which Sofa function to call).
 * @return Equation of equinoxes (in radians).
 */
double calculateEquationOfEquinoxes(
		const double terrestrialTime, const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000,
		const basic_astrodynamics::IAUConventions iauConvention = basic_astrodynamics::iau_2000_b );

//! Function to calculate the combined precession/nutation rotation matrix to go from the GCRF to the True Of Date (TOD) frame.
/*!
 * Function to calculate the rotation matrix incorporating the combined effect of nutation and precession, based on the 1976 and 1980
 * IAU models. It is a function of the epoch that is the argument of this function. Since Tudat usually works with times measured from J2000,
 * it is normally sufficient to provide the first argument only, containing the number of seconds terrestrial time since the J2000 epoch.
 * If a different reference epoch is desired, it is mandatory to provide this epoch as a Julian day as the second argument.
 * @param terrestrialTime Epoch in seconds since the reference epoch at which the precession/nutation rotation matrix is to be retrieved from SOFA.
 * @param referenceJulianDay Julian Day (JD) of the reference epoch from which the terrestrial time argument is measured.
 * @return Eigen rotation matrix that can be used to rotate 3D vectors (position or velocity) from the J2000 to the True Of Date frame.
 */
Eigen::Matrix3d getPrecessionNutationMatrix(
		const double terrestrialTime, const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

//! Function to retrieve the IAU-1976 precession model angles zeta, z, and theta.
/*!
 * Function to retrieve the IAU-1976 precession model angles zeta, z, and theta. These angles embody the transformation from the Mean-Of-Date
 * frame to the GCRF frame. The definition of these angles can be found on pages 226 to 228 in Vallado (2013).
 * @param zeta Precession angle zeta (third rotation; about z-axis).
 * @param z Precession angle z (first rotation, about z-axis).
 * @param theta Precession angle Theta (second rotation, about y-axis).
 * @param terrestrialTime Epoch in seconds since reference day at which the angles are to be retrieved from SOFA.
 * @param referenceJulianDay Reference epoch for the former parameter, given in Julian Days (default = J2000).
 */
void getPrecessionAngles(double &zeta, double &z, double &theta, const double terrestrialTime,
						 const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000);

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
