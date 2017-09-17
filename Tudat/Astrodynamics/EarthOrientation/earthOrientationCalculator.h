/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_EARTHORIENTATIONCALCULATOR_H
#define TUDAT_EARTHORIENTATIONCALCULATOR_H

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/EarthOrientation/terrestrialTimeScaleConverter.h"
#include "Tudat/Astrodynamics/EarthOrientation/polarMotionCalculator.h"
#include "Tudat/Astrodynamics/EarthOrientation/precessionNutationCalculator.h"
#include "Tudat/Astrodynamics/EarthOrientation/eopReader.h"
#include "Tudat/Mathematics/Interpolators/createInterpolator.h"

namespace tudat
{

namespace earth_orientation
{

//! Function to get the value of the TIO locator, approximated as described in IERS 2010 conventions
/*!
 * Function to get the value of the TIO locator, approximated as described in IERS 2010 conventions, Eq. (5.13), evaluated using
 * iauSp00 SOFA function.
 * \param secondsSinceJ2000 Time at which TIO locator is to be computed. Formally, input should be in TT, but may be provided
 * in any common time scale, as variations in TIO are slow.
 * \return Value of TIO locator
 */
double getApproximateTioLocator( const double secondsSinceJ2000 );

//! Function to calculate rotation from CIRS to GCRS, i.e. applying CIO-based rotations due to nutation and precession
/*!
 * Function to calculate rotation from CIRS to GCRS, i.e. applying CIO-based rotations due to nutation and precession, as
 * described by IERS Conventions, Eq. (5.10)
 * \param celestialPoleXPosition Parameter X in  IERS Conventions 2010, Section 5.4.4
 * \param celestialPoleYPosition Parameter X in  IERS Conventions 2010, Section 5.4.4
 * \param cioLocator Celestial intermediate origin locator; parameter s in  IERS Conventions 2010, Section 5.4.4
 * \return Rotation from CIRS to GCRS
 */
Eigen::Quaterniond calculateRotationFromCirsToGcrs(
        const double celestialPoleXPosition, const double celestialPoleYPosition, const double cioLocator );

//! Calculate rotation from TIRS to CIRS.
/*!
 * Calculate rotation from TIRS to CIRS, i.e. rotating over earth rotation angle. Implementes Eq. (5.5) of IERS 2010
 * Conventions
 * \param earthRotationAngle Current Earth Rotation angle
 * \return Rotation from TIRS to CIRS
 */
Eigen::Quaterniond calculateRotationFromTirsToCirs( const double earthRotationAngle );

//! Calculate rotation from ITRS to TIRS,
/*!
 * Calculate rotation from ITRS to TIRS, according to Eq. (5.3) of IERS 2010 Conventions. Applies the effect of polar motion.
 * \param xPolePosition Polar motion parameter in x-direction (typically denoted x_{p})
 * \param yPolePosition Polar motion parameter in y-direction (typically denoted x_{p})
 * \param tioLocator TIO locator.
 * \return Rotation from ITRS to TIRS
 */
Eigen::Quaterniond calculateRotationFromItrsToTirs(
        const double xPolePosition, const double yPolePosition, const double tioLocator );

//! Calculate time-derivative of rotation matrix from ITRS to GCRS
/*!
 * Calculate time-derivative of rotation matrix from ITRS to GCRS. Function approximates derivative by only including derivative
 * of Earth rotation sub-matrix
 * \param celestialPoleXPosition Parameter X in  IERS Conventions 2010, Section 5.4.4
 * \param celestialPoleYPosition Parameter X in  IERS Conventions 2010, Section 5.4.4
 * \param cioLocator Celestial intermediate origin locator; parameter s in  IERS Conventions 2010, Section 5.4.4
 * \param earthRotationAngle Current Earth Rotation angle
 * \param xPolePosition Polar motion parameter in x-direction (typically denoted x_{p})
 * \param yPolePosition Polar motion parameter in y-direction (typically denoted x_{p})
 * \param tioLocator TIO locator.
 * \return Time-derivative of rotation matrix from ITRS to GCRS
 */
Eigen::Matrix3d calculateRotationRateFromItrsToGcrs(
        const double celestialPoleXPosition, const double celestialPoleYPosition, const double cioLocator,
        const double earthRotationAngle, const double xPolePosition, const double yPolePosition, const double tioLocator );

//! Calculate time-derivative of rotation matrix from ITRS to GCRS
/*!
 * Calculate time-derivative of rotation matrix from ITRS to GCRS. Function approximates derivative by only including derivative
 * of Earth rotation sub-matrix
 * \param rotationAngles Vector containing quantities (in IERS Conventions 2010 notation): X, Y, x, ERA, xp, yp.
 * \param secondsSinceJ2000 Current time in seconds since J2000, used for computing TIO locator.
 * \return Time-derivative of rotation matrix from ITRS to GCRS
 */
Eigen::Matrix3d calculateRotationRateFromItrsToGcrs(
        const Eigen::Vector6d& rotationAngles, const double secondsSinceJ2000 );

//! Calculate rotation from ITRS to GCRS
/*!
 * Calculate rotation from ITRS to GCRS.
 * \param celestialPoleXPosition Parameter X in  IERS Conventions 2010, Section 5.4.4
 * \param celestialPoleYPosition Parameter X in  IERS Conventions 2010, Section 5.4.4
 * \param cioLocator Celestial intermediate origin locator; parameter s in  IERS Conventions 2010, Section 5.4.4
 * \param earthRotationAngle Current Earth Rotation angle
 * \param xPolePosition Polar motion parameter in x-direction (typically denoted x_{p})
 * \param yPolePosition Polar motion parameter in y-direction (typically denoted x_{p})
 * \param tioLocator TIO locator.
 * \return Rotation from ITRS to GCRS
 */
Eigen::Quaterniond calculateRotationFromItrsToGcrs(
        const double celestialPoleXPosition, const double celestialPoleYPosition, const double cioLocator,
        const double earthRotationAngle, const double xPolePosition, const double yPolePosition, const double tioLocator );

//! Calculate rotation from ITRS to GCRS
/*!
 * Calculate rotation from ITRS to GCRS. Function approximates derivative by only including derivative
 * of Earth rotation sub-matrix
 * \param rotationAngles Vector containing quantities (in IERS Conventions 2010 notation): X, Y, x, ERA, xp, yp.
 * \param secondsSinceJ2000 Current time in seconds since J2000, used for computing TIO locator.
 * \return Rotation from ITRS to GCRS
 */
Eigen::Quaterniond calculateRotationFromItrsToGcrs(
        const Eigen::Vector6d& rotationAngles, const double secondsSinceJ2000 );


//! Class to calculate earth orientation angles, i.e. those used for transforming from ITRS to GCRS
class EarthOrientationAnglesCalculator
{
public:

    //! Constructor from objects calculating three sub-parts of earth orientation.
    /*!
     *  Constructor from objects calculating sub-parts of earth orientation and time conversion,
     *  i.e. polar motion, precession/nutation and conversion between TT,TDB,UTC and UT1.
     *  \param polarMotionCalculator Pointer to object for calcu2012lating position of pole in ITRS (polarm motion).
     *  \param precessionNutationCalculator Pointer to object for calculating position of pole in GCRS (precession/nutation).
     *  \param terrestrialTimeScaleConverter Pointer to object to convert between different time scales.
     */
    EarthOrientationAnglesCalculator(
            const boost::shared_ptr< PolarMotionCalculator > polarMotionCalculator,
            const boost::shared_ptr< PrecessionNutationCalculator > precessionNutationCalculator,
            const boost::shared_ptr< TerrestrialTimeScaleConverter > terrestrialTimeScaleConverter ):
        polarMotionCalculator_( polarMotionCalculator ),
        precessionNutationCalculator_( precessionNutationCalculator ),
        terrestrialTimeScaleConverter_( terrestrialTimeScaleConverter ) { }

    //! Calculate rotation angles from ITRS to GCRS at given time value.
    /*!
     *  Calculate rotation angles from ITRS to GCRS at given time value. Any time scale combined with any time value
     *  can be used as input. TIO locator is not included in output as its value is minute and cal be easily evaluated, with
     *  no regard for time conversions in input value.
     *  \param timeValue Number of seconds since J2000 at which orientation is to be evaluated.
     *  \param timeScale Time scale in which the timeValue is given. To be taken from TimeScales enum.
     *  \return Rotation angles for ITRS<->GCRS transformation at given epoch. Order is: X, Y, s, ERA, x_p, y_p
     */
    Eigen::Vector6d getRotationAnglesFromItrsToGcrs(
            const double timeValue,
            basic_astrodynamics::TimeScales timeScale = basic_astrodynamics::tt_scale );

    //! Function to get object that calculates polar motion.
    /*!
     *  Function to get object that calculates polar motion.
     */
    boost::shared_ptr< PolarMotionCalculator > getPolarMotionCalculator( )
    { return polarMotionCalculator_; }

    //! Function to get object that calculates precession/nutation.
    /*!
     *  Function to get object that calculates precession/nutation.
     */
    boost::shared_ptr< PrecessionNutationCalculator > getPrecessionNutationCalculator( )
    { return precessionNutationCalculator_; }

    //! Function to get object that converts between time scales.
    /*!
     *  Function to get object that converts between time scales.
     */
    boost::shared_ptr< TerrestrialTimeScaleConverter > getTerrestrialTimeScaleConverter( )
    { return terrestrialTimeScaleConverter_; }


private:
    //! Pointer to object for calculating position of pole in ITRS (polarm motion).
    /*!
     *  Pointer to object for calculating position of pole in ITRS (polarm motion).
     */
    boost::shared_ptr< PolarMotionCalculator > polarMotionCalculator_;

    //! Pointer to object for calculating position of pole in GCRS (precession/nutation).
    /*!
     *  Pointer to object for calculating position of pole in GCRS (precession/nutation).
     */
    boost::shared_ptr< PrecessionNutationCalculator > precessionNutationCalculator_;

    //! Pointer to object to convert between different time scales.
    /*!
     *  Pointer to object to convert between different time scales.
     */
    boost::shared_ptr< TerrestrialTimeScaleConverter > terrestrialTimeScaleConverter_;
};

//! Function to create an EarthOrientationAnglesCalculator object, with default settings
/*!
 * Function to create an EarthOrientationAnglesCalculator object, with default settings:
 * IAU 2006 theory for precession/nutation, all (sub-)diurnal corrections to UTC-UT1 and polar motion according to IERS 2010,
 * polar motion/nutation/UT1 daily corrections published by IERS (linearly interpolated in time)
 * \return Default Earth rotation parameter object.
 */
boost::shared_ptr< EarthOrientationAnglesCalculator > createStandardEarthOrientationCalculator( );

//! Function to compute the Earth rotation angle, without normalization to [0, 2pi) range.
/*!
 * Function to compute the Earth rotation angle, without normalization to [0, 2pi) range, according to IERS Convetions 2019,
 * Eq. (5.14).
 * \param ut1SinceEpoch Time since reference Julian day in UT1.
 * \param referenceJulianDay Reference time for UT1 input
 * \return Unnormalized Earth rotation at input time.
 */
double calculateUnnormalizedEarthRotationAngle( const double ut1SinceEpoch,
                                                const double referenceJulianDay );

//! Function to create an interpolator for the Earth orientation angles
/*!
 * Function to create an interpolator for the Earth orientation angles, to reduce computation time of Earth rotation during
 * orbit propagation/estimation
 * \param intervalStart Start of time interval where interpolation data is to be generated
 * \param intervalEnd End of time interval where interpolation data is to be generated
 * \param timeStep Time step between evaluations of rotation data
 * \param timeScale Time scale for evaluation data
 * \param earthOrientationCalculator Object from which Earth orientation data is to be retrieved
 * \param interpolatorSettings Settings for the interpolation proces (default Lagrange 6 point)
 * \return Interpolator for the Earth orientation angles. Interpolated vector contains quantities (in IERS Conventions 2010
 * notation): X, Y, x, ERA, xp, yp.
 */
boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 6,1 > > >
createInterpolatorForItrsToGcrsAngles(
        const double intervalStart, const double intervalEnd, const double timeStep,
        const basic_astrodynamics::TimeScales timeScale = basic_astrodynamics::tdb_scale,
        const boost::shared_ptr< EarthOrientationAnglesCalculator > earthOrientationCalculator =
        createStandardEarthOrientationCalculator( ),
        const boost::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
        boost::make_shared< interpolators::LagrangeInterpolatorSettings >( 6 ) );

}

}

#endif // TUDAT_EARTHORIENTATIONCALCULATOR_H
