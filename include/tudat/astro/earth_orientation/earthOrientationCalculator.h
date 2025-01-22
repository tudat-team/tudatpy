/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/astro/earth_orientation/polarMotionCalculator.h"
#include "tudat/astro/earth_orientation/precessionNutationCalculator.h"
#include "tudat/astro/earth_orientation/eopReader.h"
#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"

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
 * \param celestialPoleYPosition Parameter Y in  IERS Conventions 2010, Section 5.4.4
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
 * \param ut1 Current UT1 time, used to compute Earth Rotation angle
 * \param xPolePosition Polar motion parameter in x-direction (typically denoted x_{p})
 * \param yPolePosition Polar motion parameter in y-direction (typically denoted x_{p})
 * \param tioLocator TIO locator.
 * \return Time-derivative of rotation matrix from ITRS to GCRS
 */
template< typename TimeType >
Eigen::Matrix3d calculateRotationRateFromItrsToGcrs(
        const double celestialPoleXPosition, const double celestialPoleYPosition, const double cioLocator,
        const TimeType ut1, const double xPolePosition, const double yPolePosition, const double tioLocator )
{
    Eigen::Matrix3d auxiliaryMatrix = reference_frames::Z_AXIS_ROTATION_MATRIX_DERIVATIVE_PREMULTIPLIER *
            ( -2.0 * mathematical_constants::PI / 86400.0 * 1.00273781191135448 );

    return  ( calculateRotationFromCirsToGcrs( celestialPoleXPosition, celestialPoleYPosition, cioLocator )  *
              calculateRotationFromTirsToCirs( sofa_interface::calculateEarthRotationAngleTemplated< TimeType >( ut1 - 1.42E-3 ) )
              ).toRotationMatrix( ) * auxiliaryMatrix *
            calculateRotationFromItrsToTirs( xPolePosition, yPolePosition, tioLocator ).toRotationMatrix( );

}

//! Calculate time-derivative of rotation matrix from ITRS to GCRS
/*!
 * Calculate time-derivative of rotation matrix from ITRS to GCRS. Function approximates derivative by only including derivative
 * of Earth rotation sub-matrix
 * \param rotationAngles Vector containing quantities (in IERS Conventions 2010 notation): X, Y, s, xp, yp.
 * \param ut1 Current UT1 time, used to compute Earth Rotation angle
 * \param secondsSinceJ2000 Current time in seconds since J2000, used for computing TIO locator.
 * \return Time-derivative of rotation matrix from ITRS to GCRS
 */
template< typename TimeType >
Eigen::Matrix3d calculateRotationRateFromItrsToGcrs(
        const Eigen::Vector5d& rotationAngles, const TimeType ut1, const double secondsSinceJ2000 )
{
    return calculateRotationRateFromItrsToGcrs(
                rotationAngles[ 0 ], rotationAngles[ 1 ], rotationAngles[ 2 ],
            ut1, rotationAngles[ 3 ], rotationAngles[ 4 ],
            getApproximateTioLocator( secondsSinceJ2000 ) );
}

template< typename TimeType >
Eigen::Matrix3d calculateRotationRateFromItrsToGcrs(
        const std::pair< Eigen::Vector5d, TimeType > rotationAnglesAndUt1, const double secondsSinceJ2000 )
{
    return calculateRotationRateFromItrsToGcrs(
                rotationAnglesAndUt1.first[ 0 ], rotationAnglesAndUt1.first[ 1 ], rotationAnglesAndUt1.first[ 2 ],
            rotationAnglesAndUt1.second, rotationAnglesAndUt1.first[ 3 ], rotationAnglesAndUt1.first[ 4 ],
            getApproximateTioLocator( secondsSinceJ2000 ) );
}

//! Calculate rotation from ITRS to GCRS
/*!
 * Calculate rotation from ITRS to GCRS.
 * \param celestialPoleXPosition Parameter X in  IERS Conventions 2010, Section 5.4.4
 * \param celestialPoleYPosition Parameter X in  IERS Conventions 2010, Section 5.4.4
 * \param cioLocator Celestial intermediate origin locator; parameter s in  IERS Conventions 2010, Section 5.4.4
 * \param ut1 Current UT1 time, used to compute Earth Rotation angle
 * \param xPolePosition Polar motion parameter in x-direction (typically denoted x_{p})
 * \param yPolePosition Polar motion parameter in y-direction (typically denoted x_{p})
 * \param tioLocator TIO locator.
 * \return Rotation from ITRS to GCRS
 */
template< typename TimeType >
Eigen::Quaterniond calculateRotationFromItrsToGcrs(
        const double celestialPoleXPosition, const double celestialPoleYPosition, const double cioLocator,
        const TimeType ut1, const double xPolePosition, const double yPolePosition, const double tioLocator )
{
    double currentEra = sofa_interface::calculateEarthRotationAngleTemplated< TimeType >( ut1 );

    return  calculateRotationFromCirsToGcrs( celestialPoleXPosition, celestialPoleYPosition, cioLocator ) *
            calculateRotationFromTirsToCirs( currentEra ) *
            calculateRotationFromItrsToTirs( xPolePosition, yPolePosition, tioLocator );

}

//! Calculate rotation from ITRS to GCRS
/*!
 * Calculate rotation from ITRS to GCRS. Function approximates derivative by only including derivative
 * of Earth rotation sub-matrix
 * \param rotationAngles Vector containing quantities (in IERS Conventions 2010 notation): X, Y, s, xp, yp.
 * \param ut1 Current UT1 time, used to compute Earth Rotation angle
 * \param secondsSinceJ2000 Current time in seconds since J2000, used for computing TIO locator.
 * \return Rotation from ITRS to GCRS
 */
template< typename TimeType >
Eigen::Quaterniond calculateRotationFromItrsToGcrs(
        const Eigen::Vector5d& rotationAngles, const TimeType ut1, const double secondsSinceJ2000 )
{
    return calculateRotationFromItrsToGcrs(
                rotationAngles[ 0 ], rotationAngles[ 1 ], rotationAngles[ 2 ],
            ut1, rotationAngles[ 3 ], rotationAngles[ 4 ],
            getApproximateTioLocator( secondsSinceJ2000 ) );
}

template< typename TimeType >
Eigen::Quaterniond calculateRotationFromItrsToGcrs(
        const std::pair< Eigen::Vector5d, TimeType > rotationAnglesAndUt1, const double secondsSinceJ2000 )
{
    return calculateRotationFromItrsToGcrs(
                rotationAnglesAndUt1.first[ 0 ], rotationAnglesAndUt1.first[ 1 ], rotationAnglesAndUt1.first[ 2 ],
            rotationAnglesAndUt1.second, rotationAnglesAndUt1.first[ 3 ], rotationAnglesAndUt1.first[ 4 ],
            getApproximateTioLocator( secondsSinceJ2000 ) );
}



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
            const std::shared_ptr< PolarMotionCalculator > polarMotionCalculator,
            const std::shared_ptr< PrecessionNutationCalculator > precessionNutationCalculator,
            const std::shared_ptr< TerrestrialTimeScaleConverter > terrestrialTimeScaleConverter ):
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
     *  \return Rotation angles for ITRS<->GCRS transformation at given epoch. First pair entry is: X, Y, s, x_p, y_p. Second
     *  defines UT1.
     */
    template< typename TimeType >
    std::pair< Eigen::Vector5d, TimeType > getRotationAnglesFromItrsToGcrs(
            const TimeType timeValue,
            basic_astrodynamics::TimeScales timeScale = basic_astrodynamics::tt_scale )
    {
        // Compute required time values
        TimeType terrestrialTime = terrestrialTimeScaleConverter_->getCurrentTime< TimeType >(
                    timeScale, basic_astrodynamics::tt_scale, timeValue, Eigen::Vector3d::Zero( ) );
        TimeType utc = terrestrialTimeScaleConverter_->getCurrentTime< TimeType >(
                    timeScale, basic_astrodynamics::utc_scale, timeValue, Eigen::Vector3d::Zero( ) );
        TimeType ut1 = terrestrialTimeScaleConverter_->getCurrentTime< TimeType >(
                    timeScale, basic_astrodynamics::ut1_scale, timeValue, Eigen::Vector3d::Zero( ) );

        // Compute nutation/precession parameters
        Eigen::Vector3d positionOfCipInGcrs =
                precessionNutationCalculator_->getPositionOfCipInGcrs(
                    terrestrialTime, utc );

        // Compute polar motion values
        Eigen::Vector2d positionOfCipInItrs = polarMotionCalculator_->getPositionOfCipInItrs(
                    terrestrialTime, utc );

        // Return vector of angles.
        Eigen::Vector5d rotationAngles;
        rotationAngles[ 0 ] = positionOfCipInGcrs( 0 );
        rotationAngles[ 1 ] = positionOfCipInGcrs( 1 );
        rotationAngles[ 2 ] = positionOfCipInGcrs( 2 );
        rotationAngles[ 3 ] = positionOfCipInItrs.x( );
        rotationAngles[ 4 ] = positionOfCipInItrs.y( );
        return std::make_pair( rotationAngles, ut1 );
    }


    //! Function to get object that calculates polar motion.
    /*!
     *  Function to get object that calculates polar motion.
     */
    std::shared_ptr< PolarMotionCalculator > getPolarMotionCalculator( )
    { return polarMotionCalculator_; }

    //! Function to get object that calculates precession/nutation.
    /*!
     *  Function to get object that calculates precession/nutation.
     */
    std::shared_ptr< PrecessionNutationCalculator > getPrecessionNutationCalculator( )
    { return precessionNutationCalculator_; }

    //! Function to get object that converts between time scales.
    /*!
     *  Function to get object that converts between time scales.
     */
    std::shared_ptr< TerrestrialTimeScaleConverter > getTerrestrialTimeScaleConverter( )
    { return terrestrialTimeScaleConverter_; }


private:
    //! Pointer to object for calculating position of pole in ITRS (polarm motion).
    /*!
     *  Pointer to object for calculating position of pole in ITRS (polarm motion).
     */
    std::shared_ptr< PolarMotionCalculator > polarMotionCalculator_;

    //! Pointer to object for calculating position of pole in GCRS (precession/nutation).
    /*!
     *  Pointer to object for calculating position of pole in GCRS (precession/nutation).
     */
    std::shared_ptr< PrecessionNutationCalculator > precessionNutationCalculator_;

    //! Pointer to object to convert between different time scales.
    /*!
     *  Pointer to object to convert between different time scales.
     */
    std::shared_ptr< TerrestrialTimeScaleConverter > terrestrialTimeScaleConverter_;
};

//! Function to create an EarthOrientationAnglesCalculator object, with default settings
/*!
 * Function to create an EarthOrientationAnglesCalculator object, with default settings:
 * IAU 2006 theory for precession/nutation, all (sub-)diurnal corrections to UTC-UT1 and polar motion according to IERS 2010,
 * polar motion/nutation/UT1 daily corrections published by IERS (linearly interpolated in time)
 * \return Default Earth rotation parameter object.
 */
std::shared_ptr< EarthOrientationAnglesCalculator > createStandardEarthOrientationCalculator(
        const std::shared_ptr< EOPReader > eopReader = std::make_shared< EOPReader >( ) );

//! Function to create an interpolator for the Earth orientation angles and UT1
/*!
 * Function to create an interpolator for the Earth orientation angles and UT1, to reduce computation time of Earth rotation
 * during  orbit propagation/estimation
 * \param intervalStart Start of time interval where interpolation data is to be generated
 * \param intervalEnd End of time interval where interpolation data is to be generated
 * \param timeStep Time step between evaluations of rotation data
 * \param timeScale Time scale for evaluation data
 * \param earthOrientationCalculator Object from which Earth orientation data is to be retrieved
 * \param interpolatorSettings Settings for the interpolation proces (default Lagrange 6 point)
 * \return interpolators for the Earth orientation angles (first) and for UT1 (second). Interpolated angle vector contains
 * quantities (in IERS Conventions 2010 notation): X, Y, s, xp, yp.
 */
template< typename UT1ScalarType >
std::pair< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 5, 1 > > >,
std::shared_ptr< interpolators::OneDimensionalInterpolator< double, UT1ScalarType > > >
createInterpolatorsForItrsToGcrsAngles(
        const double intervalStart, const double intervalEnd, const double timeStep,
        const basic_astrodynamics::TimeScales timeScale = basic_astrodynamics::tdb_scale,
        const std::shared_ptr< EarthOrientationAnglesCalculator > earthOrientationCalculator =
        createStandardEarthOrientationCalculator( ),
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
        std::make_shared< interpolators::LagrangeInterpolatorSettings >( 6 ) )
{
    // Define interpolators
    std::map< double, Eigen::Matrix< double, 5, 1 > > anglesMap;
    std::map< double, UT1ScalarType > ut1Map;

    // Iterate over all times and compute rotation parameters
    std::pair< Eigen::Vector5d, UT1ScalarType > currentRotationValues;
    double currentTime = intervalStart;
    while( currentTime < intervalEnd )
    {
        currentRotationValues = earthOrientationCalculator->getRotationAnglesFromItrsToGcrs< UT1ScalarType >(
                    currentTime, timeScale );
        anglesMap[ currentTime ] = currentRotationValues.first;
        ut1Map[ currentTime ] = currentRotationValues.second;
        currentTime += timeStep;
    }

    // Create interpolator for angles
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 5, 1 > > > anglesInterpolator =
            interpolators::createOneDimensionalInterpolator( anglesMap, interpolatorSettings );

    // Create interpolator for UT1
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, UT1ScalarType > > ut1Interpolator =
            interpolators::createOneDimensionalInterpolator( ut1Map, interpolatorSettings );

    return std::make_pair( anglesInterpolator, ut1Interpolator );

}

}

}

#endif // TUDAT_EARTHORIENTATIONCALCULATOR_H
