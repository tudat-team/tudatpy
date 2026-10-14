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

#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/astro/earth_orientation/earthOrientationCalculator.h"

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
Eigen::Quaterniond calculateRotationFromCirsToGcrs( const double celestialPoleXPosition,
                                                    const double celestialPoleYPosition,
                                                    const double cioLocator )
{
    // Pre-compute reused parameters
    double xParameterSquared = celestialPoleXPosition * celestialPoleXPosition;
    double yParameterSquared = celestialPoleYPosition * celestialPoleYPosition;
    double xyCrossTerm = celestialPoleXPosition * celestialPoleYPosition;
    double parameterA = 0.5 + 0.125 * ( xParameterSquared + yParameterSquared );

    // Set up rotation.
    Eigen::Matrix3d rotationMatrix;
    rotationMatrix << 1.0 - parameterA * xParameterSquared, -parameterA * xyCrossTerm, celestialPoleXPosition, -parameterA * xyCrossTerm,
            1.0 - parameterA * yParameterSquared, celestialPoleYPosition, -celestialPoleXPosition, -celestialPoleYPosition,
            1.0 - parameterA * ( xParameterSquared + yParameterSquared );
    return Eigen::Quaterniond( rotationMatrix ) * Eigen::Quaterniond( Eigen::AngleAxisd( -cioLocator, Eigen::Vector3d::UnitZ( ) ) );
}

//! Calculate rotation from TIRS to CIRS.
Eigen::Quaterniond calculateRotationFromTirsToCirs( const double earthRotationAngle )
{
    return Eigen::Quaterniond( Eigen::AngleAxisd( earthRotationAngle, Eigen::Vector3d::UnitZ( ) ) );
}

//! Calculate rotation from ITRS to TIRS.
Eigen::Quaterniond calculateRotationFromItrsToTirs( const double xPolePosition, const double yPolePosition, const double tioLocator )
{
    return Eigen::Quaterniond( Eigen::AngleAxisd( tioLocator, Eigen::Vector3d::UnitZ( ) ) *
                               Eigen::AngleAxisd( -xPolePosition, Eigen::Vector3d::UnitY( ) ) *
                               Eigen::AngleAxisd( -yPolePosition, Eigen::Vector3d::UnitX( ) ) );
}

//! Function to create an EarthOrientationAnglesCalculator object, with default settings
std::shared_ptr< EarthOrientationAnglesCalculator > createStandardEarthOrientationCalculator(
        const std::shared_ptr< EOPReader > eopReader,
        const EarthOrientationEvaluationSettings& settings )
{
    const auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( eopReader->getRawEarthOrientationData( ), settings );
    const auto ut1 = std::make_shared< EarthOrientationDataInterpolator< double > >( evaluator, EarthOrientationDataField::ut1_minus_utc );
    const auto timeConverter =
            std::make_shared< TerrestrialTimeScaleConverter >( ut1, getDefaultUT1CorrectionCalculator( ), nullptr, evaluator );
    const auto polar =
            std::make_shared< EarthOrientationDataInterpolator< Eigen::Vector2d > >( evaluator, EarthOrientationDataField::polar_motion );
    const auto cpo =
            std::make_shared< EarthOrientationDataInterpolator< Eigen::Vector2d > >( evaluator, EarthOrientationDataField::celestial_pole );
    const auto polarModel = std::make_shared< PolarMotionCalculator >(
            polar, getDefaultPolarMotionCorrectionCalculator( ), evaluator, [ timeConverter ]( const Time& tai ) {
                return timeConverter->getUt1FromTai( tai );
            } );
    const auto precession = std::make_shared< PrecessionNutationCalculator >( basic_astrodynamics::iau_2006, cpo, nullptr, evaluator );
    return std::make_shared< EarthOrientationAnglesCalculator >( polarModel, precession, timeConverter );
}

////! Function to create an interpolator for the Earth orientation angles
// std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 6,1 > > >
// createInterpolatorForItrsToGcrsAngles(
//         const double intervalStart, const double intervalEnd, const double timeStep,
//         const basic_astrodynamics::TimeScales timeScale,
//         const std::shared_ptr< EarthOrientationAnglesCalculator > earthOrientationCalculator,
//         const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
//{
//     // Interpolate Earth orientation angles
//     std::map< double, Eigen::Matrix< double, 6,1 > > orientationMap;
//     double currentTime = intervalStart;
//     while( currentTime < intervalEnd )
//     {
//         orientationMap[ currentTime ] = earthOrientationCalculator->getRotationAnglesFromItrsToGcrs(
//                     currentTime, timeScale );
//         currentTime += timeStep;
//     }

//    return interpolators::createOneDimensionalInterpolator(
//                orientationMap, interpolatorSettings );
//}

}  // namespace earth_orientation

}  // namespace tudat
