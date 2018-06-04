/*    Copyright (c) 2010-2018, Delft University of Technology
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

//! Function to create an EarthOrientationAnglesCalculator object, with default settings
std::shared_ptr< EarthOrientationAnglesCalculator > createStandardEarthOrientationCalculator(
        const std::shared_ptr< EOPReader > eopReader )
{
    // Load polar motion corrections
    std::shared_ptr< interpolators::LinearInterpolator< double, Eigen::Vector2d > > cipInItrsInterpolator =
            std::make_shared< interpolators::LinearInterpolator< double, Eigen::Vector2d > >(
                eopReader->getCipInItrsMapInSecondsSinceJ2000( ) );

    // Load nutation corrections
    std::shared_ptr< interpolators::LinearInterpolator< double, Eigen::Vector2d > > cipInGcrsCorrectionInterpolator =
            std::make_shared< interpolators::LinearInterpolator< double, Eigen::Vector2d > >(
                eopReader->getCipInGcrsCorrectionMapInSecondsSinceJ2000( ) );

    // Load default polar motion correction (sub-diural frequencies) object
    std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > > shortPeriodPolarMotionCalculator =
            getDefaultPolarMotionCorrectionCalculator( );

    // Create full polar motion calculator
    std::shared_ptr< PolarMotionCalculator > polarMotionCalculator = std::make_shared< PolarMotionCalculator >
            ( cipInItrsInterpolator, shortPeriodPolarMotionCalculator );

    // Create IAU 2006 precession/nutation calculator
    std::shared_ptr< PrecessionNutationCalculator > precessionNutationCalculator =
            std::make_shared< PrecessionNutationCalculator >( basic_astrodynamics::iau_2006, cipInGcrsCorrectionInterpolator );

    // Create default time scale converter
    std::shared_ptr< TerrestrialTimeScaleConverter > terrestrialTimeScaleConverter =
            createDefaultTimeConverter( eopReader );

    // Create EarthOrientationAnglesCalculator object
    return std::make_shared< EarthOrientationAnglesCalculator >(
                polarMotionCalculator, precessionNutationCalculator, terrestrialTimeScaleConverter );
}

////! Function to create an interpolator for the Earth orientation angles
//std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 6,1 > > >
//createInterpolatorForItrsToGcrsAngles(
//        const double intervalStart, const double intervalEnd, const double timeStep,
//        const basic_astrodynamics::TimeScales timeScale,
//        const std::shared_ptr< EarthOrientationAnglesCalculator > earthOrientationCalculator,
//        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
//{
//    // Interpolate Earth orientation angles
//    std::map< double, Eigen::Matrix< double, 6,1 > > orientationMap;
//    double currentTime = intervalStart;
//    while( currentTime < intervalEnd )
//    {
//        orientationMap[ currentTime ] = earthOrientationCalculator->getRotationAnglesFromItrsToGcrs(
//                    currentTime, timeScale );
//        currentTime += timeStep;
//    }

//    return interpolators::createOneDimensionalInterpolator(
//                orientationMap, interpolatorSettings );
//}


}

}
