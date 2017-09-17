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

#include "Tudat/Mathematics/Interpolators/createInterpolator.h"
#include "Tudat/Astrodynamics/EarthOrientation/earthOrientationCalculator.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

namespace tudat
{

namespace earth_orientation
{

//! Function to get the value of the TIO locator, approximated as described in IERS 2010 conventions
double getApproximateTioLocator( const double secondsSinceJ2000 )
{
    return iauSp00( basic_astrodynamics::JULIAN_DAY_ON_J2000, secondsSinceJ2000 / physical_constants::JULIAN_DAY );
}

//! Function to calculate rotation from CIRS to GCRS, i.e. applying CIO-based rotations due to nutation and precession
Eigen::Quaterniond calculateRotationFromCirsToGcrs(
        const double celestialPoleXPosition, const double celestialPoleYPosition, const double cioLocator )
{
    // Pre-compute reused parameters
    double xParameterSquared =celestialPoleXPosition * celestialPoleXPosition;
    double yParameterSquared =celestialPoleYPosition * celestialPoleYPosition;
    double xyCrossTerm =celestialPoleXPosition * celestialPoleYPosition;
    double parameterA = 0.5 + 0.125 * ( xParameterSquared + yParameterSquared );

    // Set up rotation.
    Eigen::Matrix3d rotationMatrix;
    rotationMatrix << 1.0 - parameterA * xParameterSquared, -parameterA * xyCrossTerm, celestialPoleXPosition,
            - parameterA * xyCrossTerm, 1.0 - parameterA * yParameterSquared, celestialPoleYPosition,
            -celestialPoleXPosition, -celestialPoleYPosition, 1.0 - parameterA * ( xParameterSquared + yParameterSquared );
    return Eigen::Quaterniond( rotationMatrix ) *
            Eigen::Quaterniond( Eigen::AngleAxisd( -cioLocator, Eigen::Vector3d::UnitZ( ) ) );
}

//! Calculate rotation from TIRS to CIRS.
Eigen::Quaterniond calculateRotationFromTirsToCirs( const double earthRotationAngle )
{
    return Eigen::Quaterniond( Eigen::AngleAxisd( earthRotationAngle, Eigen::Vector3d::UnitZ( ) ) );
}

//! Calculate rotation from ITRS to TIRS.
Eigen::Quaterniond calculateRotationFromItrsToTirs(
        const double xPolePosition, const double yPolePosition, const double tioLocator )
{
    return Eigen::Quaterniond( Eigen::AngleAxisd( tioLocator, Eigen::Vector3d::UnitZ( ) ) *
                               Eigen::AngleAxisd( -xPolePosition, Eigen::Vector3d::UnitY( ) ) *
                               Eigen::AngleAxisd( -yPolePosition, Eigen::Vector3d::UnitX( ) ) );
}

//! Calculate time-derivative of rotation matrix from ITRS to GCRS
Eigen::Matrix3d calculateRotationRateFromItrsToGcrs(
        const double celestialPoleXPosition, const double celestialPoleYPosition, const double cioLocator,
        const double earthRotationAngle, const double xPolePosition, const double yPolePosition, const double tioLocator )
{
    Eigen::Matrix3d auxiliaryMatrix = reference_frames::Z_AXIS_ROTATION_MATRIX_DERIVATIVE_PREMULTIPLIER *
            ( -2.0 * mathematical_constants::PI / 86400.0 * 1.002737909350795 );

    return  ( calculateRotationFromCirsToGcrs( celestialPoleXPosition, celestialPoleYPosition, cioLocator )  *
              calculateRotationFromTirsToCirs( earthRotationAngle ) ).toRotationMatrix( ) * auxiliaryMatrix *
            calculateRotationFromItrsToTirs( xPolePosition, yPolePosition, tioLocator ).toRotationMatrix( );

}

//! Calculate time-derivative of rotation matrix from ITRS to GCRS
Eigen::Matrix3d calculateRotationRateFromItrsToGcrs(
        const Eigen::Vector6d& rotationAngles, const double secondsSinceJ2000 )
{
    return calculateRotationRateFromItrsToGcrs(
                rotationAngles[ 0 ], rotationAngles[ 1 ], rotationAngles[ 2 ],
            rotationAngles[ 3 ], rotationAngles[ 4 ], rotationAngles[ 5 ],
            getApproximateTioLocator( secondsSinceJ2000 ) );
}

//! Calculate rotation from ITRS to GCRS
Eigen::Quaterniond calculateRotationFromItrsToGcrs(
        const double celestialPoleXPosition, const double celestialPoleYPosition, const double cioLocator,
        const double earthRotationAngle, const double xPolePosition, const double yPolePosition, const double tioLocator )
{
    return  calculateRotationFromCirsToGcrs( celestialPoleXPosition, celestialPoleYPosition, cioLocator ) *
            calculateRotationFromTirsToCirs( earthRotationAngle ) *
            calculateRotationFromItrsToTirs( xPolePosition, yPolePosition, tioLocator );

}

//! Calculate rotation from ITRS to GCRS
Eigen::Quaterniond calculateRotationFromItrsToGcrs(
        const Eigen::Vector6d& rotationAngles, const double secondsSinceJ2000 )
{
    return calculateRotationFromItrsToGcrs(
                rotationAngles[ 0 ], rotationAngles[ 1 ], rotationAngles[ 2 ],
            rotationAngles[ 3 ], rotationAngles[ 4 ], rotationAngles[ 5 ],
            getApproximateTioLocator( secondsSinceJ2000 ) );
}

//! Calculate rotation angles from ITRS to GCRS at given time value.
Eigen::Vector6d EarthOrientationAnglesCalculator::getRotationAnglesFromItrsToGcrs(
        const double timeValue, basic_astrodynamics::TimeScales timeScale )
{
    // Compute required time values
    double terrestrialTime = terrestrialTimeScaleConverter_->getCurrentTime(
                timeScale, basic_astrodynamics::tt_scale, timeValue, Eigen::Vector3d::Zero( ) );
    double utc = terrestrialTimeScaleConverter_->getCurrentTime(
                timeScale, basic_astrodynamics::utc_scale, timeValue, Eigen::Vector3d::Zero( ) );
    double ut1 = terrestrialTimeScaleConverter_->getCurrentTime(
                timeScale, basic_astrodynamics::ut1_scale, timeValue, Eigen::Vector3d::Zero( ) );

    // Compute nutation/precession parameters
    std::pair< Eigen::Vector2d, double > positionOfCipInGcrs =
            precessionNutationCalculator_->getPositionOfCipInGcrs(
                terrestrialTime, utc );

    // Compute polar motion values
    Eigen::Vector2d positionOfCipInItrs = polarMotionCalculator_->getPositionOfCipInItrs(
                terrestrialTime, utc );

    // Compute Earth Rotation Angles
    double earthRotationAngle = calculateUnnormalizedEarthRotationAngle( ut1, basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    // Return vector of angles.
    Eigen::Vector6d rotationAngles;
    rotationAngles[ 0 ] = positionOfCipInGcrs.first.x( );
    rotationAngles[ 1 ] = positionOfCipInGcrs.first.y( );
    rotationAngles[ 2 ] = positionOfCipInGcrs.second;
    rotationAngles[ 3 ] = earthRotationAngle;
    rotationAngles[ 4 ] = positionOfCipInItrs.x( );
    rotationAngles[ 5 ] = positionOfCipInItrs.y( );
    return rotationAngles;
}


//! Function to create an EarthOrientationAnglesCalculator object, with default settings
boost::shared_ptr< EarthOrientationAnglesCalculator > createStandardEarthOrientationCalculator( )
{
    // Load EOP file
    boost::shared_ptr< EOPReader > eopReader = boost::make_shared< EOPReader >( );

    // Load polar motion corrections
    boost::shared_ptr< interpolators::LinearInterpolator< double, Eigen::Vector2d > > cipInItrsInterpolator =
            boost::make_shared< interpolators::LinearInterpolator< double, Eigen::Vector2d > >(
                eopReader->getCipInItrsMapInSecondsSinceJ2000( ) );

    // Load nutation corrections
    boost::shared_ptr< interpolators::LinearInterpolator< double, Eigen::Vector2d > > cipInGcrsCorrectionInterpolator =
            boost::make_shared< interpolators::LinearInterpolator< double, Eigen::Vector2d > >(
                eopReader->getCipInGcrsCorrectionMapInSecondsSinceJ2000( ) ); // dX, dY

    // Load default polar motion correction (sub-diural frequencies) object
    boost::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > > shortPeriodPolarMotionCalculator =
            getDefaultPolarMotionCorrectionCalculator( );

    // Create full polar motion calculator
    boost::shared_ptr< PolarMotionCalculator > polarMotionCalculator = boost::make_shared< PolarMotionCalculator >
            ( cipInItrsInterpolator, shortPeriodPolarMotionCalculator );

    // Create IAU 2006 precession/nutation calculator
    boost::shared_ptr< PrecessionNutationCalculator > precessionNutationCalculator =
            boost::make_shared< PrecessionNutationCalculator >( iau_2006, cipInGcrsCorrectionInterpolator );

    // Create default time scale converter
    boost::shared_ptr< TerrestrialTimeScaleConverter > terrestrialTimeScaleConverter =
            createDefaultTimeConverter( eopReader );

    // Create EarthOrientationAnglesCalculator object
    return boost::make_shared< EarthOrientationAnglesCalculator >(
                polarMotionCalculator, precessionNutationCalculator, terrestrialTimeScaleConverter );
}

//! Function to compute the Earth rotation angle, without normalization to [0, 2pi) range.
double calculateUnnormalizedEarthRotationAngle( const double ut1SinceEpoch,
                                                const double referenceJulianDay )
{
    return 2.0 * mathematical_constants::PI *
            ( 0.7790572732640 + 1.00273781191135448 * ( ut1SinceEpoch / physical_constants::JULIAN_DAY +
                                                        ( referenceJulianDay - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) ) );
}

//! Function to create an interpolator for the Earth orientation angles
boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 6,1 > > >
createInterpolatorForItrsToGcrsAngles(
        const double intervalStart, const double intervalEnd, const double timeStep,
        const basic_astrodynamics::TimeScales timeScale,
        const boost::shared_ptr< EarthOrientationAnglesCalculator > earthOrientationCalculator,
        const boost::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
{
    // Interpolate Earth orientation angles
    std::map< double, Eigen::Matrix< double, 6,1 > > orientationMap;
    double currentTime = intervalStart;
    while( currentTime < intervalEnd )
    {
        orientationMap[ currentTime ] = earthOrientationCalculator->getRotationAnglesFromItrsToGcrs(
                    currentTime, timeScale );
        currentTime += timeStep;
    }

    return interpolators::createOneDimensionalInterpolator(
                orientationMap, interpolatorSettings );
}


}

}
