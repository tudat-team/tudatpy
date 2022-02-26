/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/observation_partials/rotationMatrixPartial.h"

namespace tudat
{

namespace observation_partials
{

//! Function to calculate a partial of rotation matrix from a body-fixed to inertial frame w.r.t. a constant rotation rate.
Eigen::Matrix3d calculatePartialOfRotationMatrixFromLocalFrameWrtConstantRotationRate(
        const Eigen::Quaterniond initialBodyFixedToIntegrationFrame,
        const double rotationRate, const double timeSinceEpoch )
{
    double currentRotationAngle = rotationRate * timeSinceEpoch;

    // Compute partial of rotation term containing rotation rate
    Eigen::Matrix3d rotationMatrixDerivative;
    double sineOfAngle = sin( currentRotationAngle );
    double cosineOfAngle = cos( currentRotationAngle );
    rotationMatrixDerivative << -sineOfAngle, -cosineOfAngle, 0.0, cosineOfAngle, - sineOfAngle, 0.0, 0.0, 0.0, 0.0;

    return timeSinceEpoch * ( initialBodyFixedToIntegrationFrame.toRotationMatrix( ) ) * rotationMatrixDerivative;

}

//! Function to calculate partial of a rotation matrix derivative (body-fixed to inertial) w.r.t. a constant rotation rate.
Eigen::Matrix3d calculatePartialOfRotationMatrixFromLocalFrameDerivativeWrtConstantRotationRate(
        const Eigen::Matrix3d currentRotationFromLocalToGlobalFrame,
        const double rotationRate, const double timeSinceEpoch )
{
    return currentRotationFromLocalToGlobalFrame * reference_frames::Z_AXIS_ROTATION_MATRIX_DERIVATIVE_PREMULTIPLIER *
            ( rotationRate * timeSinceEpoch * reference_frames::Z_AXIS_ROTATION_MATRIX_DERIVATIVE_PREMULTIPLIER -
              Eigen::Matrix3d::Identity( ) );

}
//! Function to calculate a partial of rotation matrix from a body-fixed to inertial frame w.r.t. a constant pole right
//! ascension and declination
std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixFromLocalFrameWrtPoleOrientation(
        const Eigen::Vector3d initialOrientationAngles,
        const double rotationRate, const double timeSinceEpoch )
{
    Eigen::Matrix3d commonTerm = Eigen::AngleAxisd(
                -1.0 * ( -initialOrientationAngles.z( ) - rotationRate * timeSinceEpoch ),
                Eigen::Vector3d::UnitZ( ) ).toRotationMatrix( );


    double rightAscension = initialOrientationAngles.x( );
    double declination = initialOrientationAngles.y( );

    // Compute partial of rotation term containing right ascension.
    Eigen::Matrix3d rightAscensionPartial = -reference_frames::getDerivativeOfZAxisRotationWrtAngle(
                -( rightAscension + mathematical_constants::PI / 2.0 ) );


    // Compute partial of rotation term containing declination.
    Eigen::Matrix3d declinationPartial = reference_frames::getDerivativeOfXAxisRotationWrtAngle(
                - ( mathematical_constants::PI / 2.0 - declination ) );

    // Compute partials.
    std::vector< Eigen::Matrix3d > rotationMatrixPartials;
    rotationMatrixPartials.push_back(
                rightAscensionPartial * Eigen::Matrix3d(
                    Eigen::AngleAxisd( -( declination - mathematical_constants::PI / 2.0 ),
                                       Eigen::Vector3d::UnitX( ) ) ) * commonTerm );
    rotationMatrixPartials.push_back(
                Eigen::AngleAxisd( ( rightAscension + mathematical_constants::PI / 2.0 ),
                                   Eigen::Vector3d::UnitZ( ) ).toRotationMatrix( ) *
                declinationPartial * commonTerm );

    return rotationMatrixPartials;
}

//! Function to calculate a partial of rotation matrix derivative from a body-fixed to inertial frame w.r.t. a constant
//! pole right  ascension and declination.
std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixFromLocalFrameDerivativeWrtPoleOrientation(
        const Eigen::Vector3d initialOrientationAngles,
        const double rotationRate, const double timeSinceEpoch )
{
    std::vector< Eigen::Matrix3d > partialsOfRotationMatrix =
            calculatePartialOfRotationMatrixFromLocalFrameWrtPoleOrientation(
                initialOrientationAngles, rotationRate, timeSinceEpoch );
    partialsOfRotationMatrix[ 0 ] = -rotationRate *
            partialsOfRotationMatrix.at( 0 ) * reference_frames::Z_AXIS_ROTATION_MATRIX_DERIVATIVE_PREMULTIPLIER ;
    partialsOfRotationMatrix[ 1 ] = -rotationRate *
            partialsOfRotationMatrix.at( 1 )* reference_frames::Z_AXIS_ROTATION_MATRIX_DERIVATIVE_PREMULTIPLIER ;

    return partialsOfRotationMatrix;
}

//! Function to calculate a partial of rotation matrix from a body-fixed to inertial frame w.r.t.
//! the periodic spin variations.
std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixFromLocalFrameWrtPeriodicSpinVariations(
        const std::shared_ptr< ephemerides::PlanetaryOrientationAngleCalculator > planetaryOrientationCalculator,
        const Eigen::Quaterniond &rotationFromMeanOrbitToIcrf,
        const Eigen::Quaterniond& polarMotionRotation,
        const double ephemerisTime )
{
    double currentMeanAnomaly = planetaryOrientationCalculator->getBodyMeanAnomalyAtEpoch( ) +
            planetaryOrientationCalculator->getBodyMeanMotion( ) * ephemerisTime;

    Eigen::Vector3d currentAngleCorrections = planetaryOrientationCalculator->updateAndGetRotationAngles( ephemerisTime );

    double currentPhiAngle = currentAngleCorrections.z( );

    std::vector< Eigen::Matrix3d > rotationMatrixPartials;

    std::map< double, std::pair< double, double > > rotationrateCorrections = planetaryOrientationCalculator->getRotationRateCorrections();
    for( std::map< double, std::pair< double, double > >::iterator correctionIterator = rotationrateCorrections.begin( );
         correctionIterator != rotationrateCorrections.end( ); correctionIterator++ )
    {
        rotationMatrixPartials.push_back(
                    std::cos( correctionIterator->first * currentMeanAnomaly) *
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) *
                      Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) ).toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << -std::sin( currentPhiAngle ), -std::cos( currentPhiAngle ), 0.0,
                      std::cos( currentPhiAngle ), -std::sin( currentPhiAngle ), 0.0, 0.0, 0.0, 0.0 ).finished( ) *
                    polarMotionRotation.toRotationMatrix() );

        rotationMatrixPartials.push_back(
                    std::sin( correctionIterator->first * currentMeanAnomaly) *
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) *
                      Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) ).toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << -std::sin( currentPhiAngle ), -std::cos( currentPhiAngle ), 0.0,
                      std::cos( currentPhiAngle ), -std::sin( currentPhiAngle ), 0.0, 0.0, 0.0, 0.0 ).finished( ) *
                    polarMotionRotation.toRotationMatrix() );
    }

    return rotationMatrixPartials;
}

//! Function to calculate a partial of rotation matrix derivative from a body-fixed to inertial frame w.r.t.
//! the periodic spin variations.
std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixFromLocalFrameDerivativeWrtPeriodicSpinVariations(
        const std::shared_ptr< ephemerides::PlanetaryOrientationAngleCalculator > planetaryOrientationCalculator,
        const Eigen::Quaterniond &rotationFromMeanOrbitToIcrf,
        const Eigen::Quaterniond& polarMotionRotation,
        const double ephemerisTime )
{
    double meanMotion = planetaryOrientationCalculator->getBodyMeanMotion( );

    double currentMeanAnomaly = planetaryOrientationCalculator->getBodyMeanAnomalyAtEpoch( )
            + meanMotion * ephemerisTime;

    Eigen::Vector3d currentAngleCorrections = planetaryOrientationCalculator->updateAndGetRotationAngles( ephemerisTime );

    double currentPhiAngle = currentAngleCorrections.z( );

    double meanPhiAngleDerivative = planetaryOrientationCalculator->getcurrentMeanPhiAngleDerivative( ephemerisTime );

    std::vector< Eigen::Matrix3d > partialsOfRotationMatrix;

    std::map< double, std::pair< double, double > > rotationrateCorrections = planetaryOrientationCalculator->getRotationRateCorrections();
    for( std::map< double, std::pair< double, double > >::iterator correctionIterator = rotationrateCorrections.begin( );
         correctionIterator != rotationrateCorrections.end( ); correctionIterator++ )
    {
        partialsOfRotationMatrix.push_back(
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) *
                      Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) ).toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) <<
                      correctionIterator->first * meanMotion * std::sin( correctionIterator->first * currentMeanAnomaly ) *
                      std::sin( currentPhiAngle ) - std::cos( correctionIterator->first * currentMeanAnomaly ) *
                      std::cos( currentPhiAngle ) * meanPhiAngleDerivative,
                      correctionIterator->first * meanMotion * std::sin( correctionIterator->first * currentMeanAnomaly ) *
                      std::cos( currentPhiAngle ) + std::cos( correctionIterator->first * currentMeanAnomaly ) *
                      std::sin( currentPhiAngle ) * meanPhiAngleDerivative,
                      0.0,
                      -correctionIterator->first * meanMotion * std::sin( correctionIterator->first * currentMeanAnomaly ) *
                      std::cos( currentPhiAngle ) - std::cos( correctionIterator->first * currentMeanAnomaly ) *
                      std::sin( currentPhiAngle ) * meanPhiAngleDerivative,
                      correctionIterator->first * meanMotion * std::sin( correctionIterator->first * currentMeanAnomaly ) *
                      std::sin( currentPhiAngle ) - std::cos( correctionIterator->first * currentMeanAnomaly ) *
                      std::cos( currentPhiAngle ) * meanPhiAngleDerivative,
                      0.0,
                      0.0, 0.0, 0.0 ).finished( ) *
                    polarMotionRotation.toRotationMatrix() );


        partialsOfRotationMatrix.push_back(
                    rotationFromMeanOrbitToIcrf .toRotationMatrix( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) *
                      Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) ).toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) <<
                      -correctionIterator->first * meanMotion * std::cos( correctionIterator->first * currentMeanAnomaly ) *
                      std::sin( currentPhiAngle ) - std::cos( currentPhiAngle ) *
                      std::sin( correctionIterator->first * currentMeanAnomaly ) * meanPhiAngleDerivative,
                      -correctionIterator->first * meanMotion * std::cos( correctionIterator->first * currentMeanAnomaly ) *
                      std::cos( currentPhiAngle ) + std::sin( currentPhiAngle ) *
                      std::sin( correctionIterator->first * currentMeanAnomaly ) * meanPhiAngleDerivative,
                      0.0,
                      correctionIterator->first * meanMotion * std::cos( correctionIterator->first * currentMeanAnomaly ) *
                      std::cos( currentPhiAngle ) - std::sin( currentPhiAngle ) *
                      std::sin( correctionIterator->first * currentMeanAnomaly ) * meanPhiAngleDerivative,
                      -correctionIterator->first * meanMotion * std::cos( correctionIterator->first * currentMeanAnomaly ) *
                      std::sin( currentPhiAngle ) - std::cos( currentPhiAngle ) *
                      std::sin( correctionIterator->first * currentMeanAnomaly ) * meanPhiAngleDerivative,
                      0.0,
                      0.0, 0.0, 0.0 ).finished( ) *
                    polarMotionRotation.toRotationMatrix() );
    }

    return partialsOfRotationMatrix;

}

//! Function to calculate a partial of rotation matrix from a body-fixed to inertial frame w.r.t.
//! polar motion amplitude.
std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixFromLocalFrameWrtPolarMotionAmplitude(
        const std::shared_ptr< ephemerides::PlanetaryOrientationAngleCalculator > planetaryOrientationCalculator,
        const Eigen::Quaterniond& rotationFromMeanOrbitToIcrf,
        const Eigen::Quaterniond& rotationFromBodyFixedToIntermediateInertialFrame,
        const double ephemerisTime )
{
    // Get current meanm anomaly and polar motion
    double currentMeanAnomaly = planetaryOrientationCalculator->getBodyMeanAnomalyAtEpoch( )
            + planetaryOrientationCalculator->getBodyMeanMotion( ) * ephemerisTime;
    Eigen::Vector2d polarMotion = planetaryOrientationCalculator ->getPolarMotion( ephemerisTime );


    std::map< double, std::pair< double, double > > xPolarMotionCoefficients = planetaryOrientationCalculator->getXpolarMotionCoefficients();
    std::map< double, std::pair< double, double > > yPolarMotionCoefficients = planetaryOrientationCalculator->getYpolarMotionCoefficients();


//    for( std::map< double, std::pair< double, double > >::iterator correctionIterator = rotationrateCorrections.begin( );
//         correctionIterator != rotationrateCorrections.end( ); correctionIterator++ )

    if ( xPolarMotionCoefficients.size() != yPolarMotionCoefficients.size() )
    {
        throw std::runtime_error( "Error, unconsistent sizes when comparing x and y polar motion"
                                  "amplitude coefficients." );
    }

    std::vector< Eigen::Matrix3d > rotationMatrixPartials;
    for( std::map< double, std::pair< double, double > >::iterator xPolarMotionCoefficientIterator = xPolarMotionCoefficients.begin( );
            xPolarMotionCoefficientIterator != xPolarMotionCoefficients.end( ); xPolarMotionCoefficientIterator++ )
    {
        rotationMatrixPartials.push_back(
                    -std::cos( xPolarMotionCoefficientIterator->first * currentMeanAnomaly) *
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    rotationFromBodyFixedToIntermediateInertialFrame.toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << -std::sin( -polarMotion.x( ) ), 0.0, std::cos( -polarMotion.x( ) ),
                      0.0, 0.0, 0.0,
                      -std::cos( -polarMotion.x( ) ), 0.0, -std::sin( -polarMotion.x( ) ) ).finished( ) *
                    Eigen::AngleAxisd( -polarMotion.y( ), Eigen::Vector3d::UnitX( ) ) );

        rotationMatrixPartials.push_back(
                    -std::sin( xPolarMotionCoefficientIterator->first * currentMeanAnomaly) *
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    rotationFromBodyFixedToIntermediateInertialFrame.toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << -std::sin( -polarMotion.x( ) ), 0.0, std::cos( -polarMotion.x( ) ),
                      0.0, 0.0, 0.0,
                      -std::cos( -polarMotion.x( ) ), 0.0, -std::sin( -polarMotion.x( ) ) ).finished( ) *
                    Eigen::AngleAxisd( -polarMotion.y( ), Eigen::Vector3d::UnitX( ) ) );

        rotationMatrixPartials.push_back(
                    -std::cos( xPolarMotionCoefficientIterator->first * currentMeanAnomaly) *
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    rotationFromBodyFixedToIntermediateInertialFrame.toRotationMatrix( ) *
                    Eigen::AngleAxisd( -polarMotion.x( ), Eigen::Vector3d::UnitY( ) ) *
                    ( Eigen::Matrix3d( ) << 0.0, 0.0, 0.0,
                      0.0, -std::sin( -polarMotion.y( ) ), -std::cos( -polarMotion.y( ) ),
                      0.0, std::cos( -polarMotion.y( ) ), -std::sin( -polarMotion.y( ) ) ).finished( ) );

        rotationMatrixPartials.push_back(
                    -std::sin( xPolarMotionCoefficientIterator->first * currentMeanAnomaly) *
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    rotationFromBodyFixedToIntermediateInertialFrame.toRotationMatrix( ) *
                    Eigen::AngleAxisd( -polarMotion.x( ), Eigen::Vector3d::UnitY( ) ) *
                    ( Eigen::Matrix3d( ) << 0.0, 0.0, 0.0,
                      0.0, -std::sin( -polarMotion.y( ) ), -std::cos( -polarMotion.y( ) ),
                      0.0, std::cos( -polarMotion.y( ) ), -std::sin( -polarMotion.y( ) ) ).finished( ) );
    }

    return rotationMatrixPartials;
}

//! Function to calculate a partial of rotation matrix derivative from a body-fixed to inertial frame w.r.t.
//! polar motion amplitude.
std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixFromLocalFrameDerivativeWrtPolarMotionAmplitude(
        const std::shared_ptr< ephemerides::PlanetaryOrientationAngleCalculator > planetaryOrientationCalculator,
        const Eigen::Quaterniond &rotationFromMeanOrbitToIcrf,
        const double ephemerisTime )
{
    double currentMeanAnomaly = planetaryOrientationCalculator->getBodyMeanAnomalyAtEpoch( )
            + planetaryOrientationCalculator->getBodyMeanMotion( ) * ephemerisTime;
    Eigen::Vector3d currentAngleCorrections = planetaryOrientationCalculator->updateAndGetRotationAngles( ephemerisTime );
    double currentPhiAngle = currentAngleCorrections.z( );
    double meanPhiAngleDerivative = planetaryOrientationCalculator->getcurrentMeanPhiAngleDerivative( ephemerisTime );
    Eigen::Vector2d polarMotion = planetaryOrientationCalculator ->getPolarMotion( ephemerisTime );

    std::vector< Eigen::Matrix3d > partialsOfRotationMatrix;

    std::map< double, std::pair< double, double > > xPolarMotionCoefficients = planetaryOrientationCalculator->getXpolarMotionCoefficients();
    std::map< double, std::pair< double, double > > yPolarMotionCoefficients = planetaryOrientationCalculator->getYpolarMotionCoefficients();

    if ( xPolarMotionCoefficients.size() != yPolarMotionCoefficients.size() )
    {
        throw std::runtime_error( "Error, unconsistent sizes when comparing x and y polar motion"
                                  "amplitude coefficients." );
    }

    for( std::map< double, std::pair< double, double > >::iterator xPolarMotionCoefficientIterator = xPolarMotionCoefficients.begin( );
            xPolarMotionCoefficientIterator != xPolarMotionCoefficients.end( ); xPolarMotionCoefficientIterator++ )
    {
        partialsOfRotationMatrix.push_back(
                    -meanPhiAngleDerivative *
                    std::cos( xPolarMotionCoefficientIterator->first * currentMeanAnomaly) *
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) *
                      Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) ).toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << -std::sin( currentPhiAngle ), -std::cos( currentPhiAngle ), 0.0,
                      std::cos( currentPhiAngle ), -std::sin( currentPhiAngle ), 0.0, 0.0, 0.0, 0.0 ).finished( ) *
                    ( Eigen::Matrix3d( ) << -std::sin( -polarMotion.x( ) ), 0.0, std::cos( -polarMotion.x( ) ),
                      0.0, 0.0, 0.0,
                      -std::cos( -polarMotion.x( ) ), 0.0, -std::sin( -polarMotion.x( ) ) ).finished( ) *
                    Eigen::AngleAxisd( -polarMotion.y( ), Eigen::Vector3d::UnitX( ) ) );

        partialsOfRotationMatrix.push_back(
                    -meanPhiAngleDerivative *
                    std::sin( xPolarMotionCoefficientIterator->first * currentMeanAnomaly) *
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) *
                      Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) ).toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << -std::sin( currentPhiAngle ), -std::cos( currentPhiAngle ), 0.0,
                      std::cos( currentPhiAngle ), -std::sin( currentPhiAngle ), 0.0, 0.0, 0.0, 0.0 ).finished( ) *
                    ( Eigen::Matrix3d( ) << -std::sin( -polarMotion.x( ) ), 0.0, std::cos( -polarMotion.x( ) ),
                      0.0, 0.0, 0.0,
                      -std::cos( -polarMotion.x( ) ), 0.0, -std::sin( -polarMotion.x( ) ) ).finished( ) *
                    Eigen::AngleAxisd( -polarMotion.y( ), Eigen::Vector3d::UnitX( ) ) );

        partialsOfRotationMatrix.push_back(
                    -meanPhiAngleDerivative *
                    std::cos( xPolarMotionCoefficientIterator->first * currentMeanAnomaly) *
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) *
                      Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) ).toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << -std::sin( currentPhiAngle ), -std::cos( currentPhiAngle ), 0.0,
                      std::cos( currentPhiAngle ), -std::sin( currentPhiAngle ), 0.0, 0.0, 0.0, 0.0 ).finished( ) *
                    Eigen::AngleAxisd( -polarMotion.x( ), Eigen::Vector3d::UnitY( ) ) *
                    ( Eigen::Matrix3d( ) << 0.0, 0.0, 0.0,
                      0.0, -std::sin( -polarMotion.y( ) ), -std::cos( -polarMotion.y( ) ),
                      0.0, std::cos( -polarMotion.y( ) ), -std::sin( -polarMotion.y( ) ) ).finished( ) );

        partialsOfRotationMatrix.push_back(
                    -meanPhiAngleDerivative *
                    std::sin( xPolarMotionCoefficientIterator->first * currentMeanAnomaly) *
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) *
                      Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) ).toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << -std::sin( currentPhiAngle ), -std::cos( currentPhiAngle ), 0.0,
                      std::cos( currentPhiAngle ), -std::sin( currentPhiAngle ), 0.0, 0.0, 0.0, 0.0 ).finished( ) *
                    Eigen::AngleAxisd( -polarMotion.x( ), Eigen::Vector3d::UnitY( ) ) *
                    ( Eigen::Matrix3d( ) << 0.0, 0.0, 0.0,
                      0.0, -std::sin( -polarMotion.y( ) ), -std::cos( -polarMotion.y( ) ),
                      0.0, std::cos( -polarMotion.y( ) ), -std::sin( -polarMotion.y( ) ) ).finished( ) );
    }

    return partialsOfRotationMatrix;

}

//! Function to calculate a partial of rotation matrix from a body-fixed to inertial frame w.r.t.
//! core factor.
Eigen::Matrix3d calculatePartialOfRotationMatrixFromLocalFrameWrtCoreFactor(
        const std::shared_ptr< ephemerides::PlanetaryOrientationAngleCalculator > planetaryOrientationCalculator,
        const Eigen::Quaterniond &rotationFromMeanOrbitToIcrf,
        const Eigen::Quaterniond& polarMotionRotation,
        const double ephemerisTime )
{
    double meanMotion = planetaryOrientationCalculator->getBodyMeanMotion( );

    double bodyMeanAnomalyAtEpoch = planetaryOrientationCalculator->getBodyMeanAnomalyAtEpoch( );

    std::map< double, std::pair< double, double > > meanMotionDirectNutationCorrections =
            planetaryOrientationCalculator->getMeanMotionDirectNutationCorrections( );

    std::vector< std::map< double, std::pair< double, double > > > meanMotionTimeDependentPhaseNutationCorrections =
            planetaryOrientationCalculator->getMeanMotionTimeDependentPhaseNutationCorrections( );

    std::vector< std::function< double( const double ) > > phaseAngleCorrectionFunctions =
            planetaryOrientationCalculator->getphaseAngleCorrectionFunctions( );

    double coreFactor = planetaryOrientationCalculator->getCorefactor( );

    double freeCoreNutationRate = planetaryOrientationCalculator->getFreeCoreNutationRate( );

    double angleIAtEpoch = planetaryOrientationCalculator->getAngleIAtEpoch( );

    Eigen::Vector3d currentAngleCorrections = planetaryOrientationCalculator->updateAndGetRotationAngles( ephemerisTime );

    Eigen::Matrix3d partialsOfRotationMatrix = Eigen::Matrix3d::Zero( );

    for( std::map< double, std::pair< double, double > >::iterator correctionIterator = meanMotionDirectNutationCorrections.begin( );
        correctionIterator != meanMotionDirectNutationCorrections.end( ); correctionIterator++ )
    {
        double a_m = correctionIterator->first * meanMotion;

        double Psi_m = correctionIterator->second.second + coreFactor * a_m /
                ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                ( a_m * correctionIterator->second.second + freeCoreNutationRate *
                  correctionIterator->second.first / std::sin( angleIAtEpoch ) );

        partialsOfRotationMatrix +=
                std::sin(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) ) *
                a_m / ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                ( a_m * correctionIterator->second.second + freeCoreNutationRate *
                  correctionIterator->second.first / std::sin( angleIAtEpoch ) ) *
                rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                ( Eigen::Matrix3d( ) << -std::sin( currentAngleCorrections.x( ) ), -std::cos( currentAngleCorrections.x( ) ), 0.0,
                  std::cos( currentAngleCorrections.x( ) ), -std::sin( currentAngleCorrections.x( ) ), 0.0,
                  0.0, 0.0, 0.0 ).finished( ) *
                ( Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) *
                  Eigen::AngleAxisd( currentAngleCorrections.z( ), Eigen::Vector3d::UnitZ( ) ) ).toRotationMatrix( ) *
                polarMotionRotation.toRotationMatrix( );

        partialsOfRotationMatrix +=
                std::cos(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) ) *
                a_m / ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                ( a_m * correctionIterator->second.first + freeCoreNutationRate *
                  correctionIterator->second.second * std::sin( angleIAtEpoch ) ) *
                rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) ).toRotationMatrix( ) *
                ( Eigen::Matrix3d( ) << 0.0, 0.0, 0.0,
                  0.0, -std::sin( currentAngleCorrections.y( ) ), -std::cos( currentAngleCorrections.y( ) ),
                  0.0, std::cos( currentAngleCorrections.y( ) ), -std::sin( currentAngleCorrections.y( ) ) ).finished( ) *
                ( Eigen::AngleAxisd( currentAngleCorrections.z( ), Eigen::Vector3d::UnitZ( ) ) ).toRotationMatrix( ) *
                polarMotionRotation.toRotationMatrix( );

        partialsOfRotationMatrix +=
                ( Psi_m * std::sin( correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) ) *
                  std::sin( currentAngleCorrections.y( ) ) *
                  std::cos(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) ) *
                  a_m / ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                  ( a_m * correctionIterator->second.first + freeCoreNutationRate *
                    correctionIterator->second.second * std::sin( angleIAtEpoch ) ) -
                  std::sin(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) )*
                  a_m / ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                  ( a_m * correctionIterator->second.second + freeCoreNutationRate *
                    correctionIterator->second.first / std::sin( angleIAtEpoch ) ) *
                  std::cos( currentAngleCorrections.y( ) ) ) *
                rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) *
                  Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) ).toRotationMatrix( ) *
                ( Eigen::Matrix3d( ) << -std::sin( currentAngleCorrections.z( ) ), -std::cos( currentAngleCorrections.z( ) ), 0.0,
                  std::cos( currentAngleCorrections.z( ) ), -std::sin( currentAngleCorrections.z( ) ), 0.0,
                  0.0, 0.0, 0.0 ).finished( ) *
                polarMotionRotation.toRotationMatrix( );
    }

    double currentPhaseAngleCorrection;

    for( unsigned int i = 0; i < meanMotionTimeDependentPhaseNutationCorrections.size( ); i++ )
    {
        currentPhaseAngleCorrection = phaseAngleCorrectionFunctions[ i ]( ephemerisTime );
        for( std::map< double, std::pair< double, double > >::iterator correctionIterator =
            meanMotionTimeDependentPhaseNutationCorrections[ i ].begin( );
            correctionIterator != meanMotionTimeDependentPhaseNutationCorrections[ i ].end( ); correctionIterator++ )
        {
            double a_m = correctionIterator->first * meanMotion;

            double Psi_m = correctionIterator->second.second + coreFactor * a_m /
                    ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                    ( a_m * correctionIterator->second.second + freeCoreNutationRate *
                      correctionIterator->second.first / std::sin( angleIAtEpoch ) );

            partialsOfRotationMatrix +=
                    std::sin( correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) +
                              currentPhaseAngleCorrection ) *
                    a_m / ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                    ( a_m * correctionIterator->second.second + freeCoreNutationRate *
                      correctionIterator->second.first / std::sin( angleIAtEpoch ) ) *
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << -std::sin( currentAngleCorrections.x( ) ), -std::cos( currentAngleCorrections.x( ) ), 0.0,
                      std::cos( currentAngleCorrections.x( ) ), -std::sin( currentAngleCorrections.x( ) ), 0.0,
                      0.0, 0.0, 0.0 ).finished( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) *
                      Eigen::AngleAxisd( currentAngleCorrections.z( ), Eigen::Vector3d::UnitZ( ) ) ).toRotationMatrix( ) *
                    polarMotionRotation.toRotationMatrix( );

            partialsOfRotationMatrix +=
                    std::cos( correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) +
                              currentPhaseAngleCorrection ) *
                    a_m / ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                    ( a_m * correctionIterator->second.first + freeCoreNutationRate *
                      correctionIterator->second.second * std::sin( angleIAtEpoch ) ) *
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) ).toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << 0.0, 0.0, 0.0,
                      0.0, -std::sin( currentAngleCorrections.y( ) ), -std::cos( currentAngleCorrections.y( ) ),
                      0.0, std::cos( currentAngleCorrections.y( ) ), -std::sin( currentAngleCorrections.y( ) ) ).finished( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.z( ), Eigen::Vector3d::UnitZ( ) ) ).toRotationMatrix( ) *
                    polarMotionRotation.toRotationMatrix( );

            partialsOfRotationMatrix +=
                    ( Psi_m * std::sin( correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) +
                                        currentPhaseAngleCorrection ) *
                      std::sin( currentAngleCorrections.y( ) ) *
                      std::cos( correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) +
                                currentPhaseAngleCorrection ) *
                      a_m / ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                      ( a_m * correctionIterator->second.first + freeCoreNutationRate *
                        correctionIterator->second.second * std::sin( angleIAtEpoch ) ) -
                      std::sin( correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) +
                                currentPhaseAngleCorrection )*
                      a_m / ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                      ( a_m * correctionIterator->second.second + freeCoreNutationRate *
                        correctionIterator->second.first / std::sin( angleIAtEpoch ) ) *
                      std::cos( currentAngleCorrections.y( ) ) ) *
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) *
                      Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) ).toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << -std::sin( currentAngleCorrections.z( ) ), -std::cos( currentAngleCorrections.z( ) ), 0.0,
                      std::cos( currentAngleCorrections.z( ) ), -std::sin( currentAngleCorrections.z( ) ), 0.0,
                      0.0, 0.0, 0.0 ).finished( ) *
                    polarMotionRotation.toRotationMatrix( );
        }
    }

    return partialsOfRotationMatrix;

}

//! Function to calculate a partial of rotation matrix derivative from a body-fixed to inertial frame w.r.t.
//! core factor.
Eigen::Matrix3d calculatePartialOfRotationMatrixFromLocalFrameDerivativeWrtCoreFactor(
        const std::shared_ptr< ephemerides::PlanetaryOrientationAngleCalculator > planetaryOrientationCalculator,
        const Eigen::Quaterniond &rotationFromMeanOrbitToIcrf,
        const Eigen::Quaterniond& polarMotionRotation,
        const double ephemerisTime )
{
    double meanMotion = planetaryOrientationCalculator->getBodyMeanMotion( );

    double bodyMeanAnomalyAtEpoch = planetaryOrientationCalculator->getBodyMeanAnomalyAtEpoch( );

    std::map< double, std::pair< double, double > > meanMotionDirectNutationCorrections =
            planetaryOrientationCalculator->getMeanMotionDirectNutationCorrections( );

    std::vector< std::map< double, std::pair< double, double > > > meanMotionTimeDependentPhaseNutationCorrections =
            planetaryOrientationCalculator->getMeanMotionTimeDependentPhaseNutationCorrections( );

    std::vector< std::function< double( const double ) > > phaseAngleCorrectionFunctions =
            planetaryOrientationCalculator->getphaseAngleCorrectionFunctions( );

    double coreFactor = planetaryOrientationCalculator->getCorefactor( );

    double freeCoreNutationRate = planetaryOrientationCalculator->getFreeCoreNutationRate( );

    double angleIAtEpoch = planetaryOrientationCalculator->getAngleIAtEpoch( );

    Eigen::Vector3d currentAngleCorrections = planetaryOrientationCalculator->updateAndGetRotationAngles( ephemerisTime );

    double currentMeanPhiAngleDerivative = planetaryOrientationCalculator->getcurrentMeanPhiAngleDerivative( ephemerisTime );

    Eigen::Matrix3d partialsOfRotationMatrix = Eigen::Matrix3d::Zero( );

    for( std::map< double, std::pair< double, double > >::iterator correctionIterator = meanMotionDirectNutationCorrections.begin( );
        correctionIterator != meanMotionDirectNutationCorrections.end( ); correctionIterator++ )
    {
        double a_m = correctionIterator->first * meanMotion;

        double Psi_m = correctionIterator->second.second + coreFactor * a_m /
                ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                ( a_m * correctionIterator->second.second + freeCoreNutationRate *
                  correctionIterator->second.first / std::sin( angleIAtEpoch ) );

        partialsOfRotationMatrix +=
                currentMeanPhiAngleDerivative *
                std::sin(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) ) *
                a_m / ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                ( a_m * correctionIterator->second.second + freeCoreNutationRate *
                  correctionIterator->second.first / std::sin( angleIAtEpoch ) ) *
                rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                ( Eigen::Matrix3d( ) << -std::sin( currentAngleCorrections.x( ) ), -std::cos( currentAngleCorrections.x( ) ), 0.0,
                  std::cos( currentAngleCorrections.x( ) ), -std::sin( currentAngleCorrections.x( ) ), 0.0,
                  0.0, 0.0, 0.0 ).finished( ) *
                ( Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) ).toRotationMatrix( ) *
                ( Eigen::Matrix3d( ) << -std::sin( currentAngleCorrections.z( ) ), -std::cos( currentAngleCorrections.z( ) ), 0.0,
                  std::cos( currentAngleCorrections.z( ) ), -std::sin( currentAngleCorrections.z( ) ), 0.0,
                  0.0, 0.0, 0.0 ).finished( ) *
                polarMotionRotation.toRotationMatrix( );

        partialsOfRotationMatrix +=
                currentMeanPhiAngleDerivative *
                std::cos(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) ) *
                a_m / ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                ( a_m * correctionIterator->second.first + freeCoreNutationRate *
                  correctionIterator->second.second * std::sin( angleIAtEpoch ) ) *
                rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) ).toRotationMatrix( ) *
                ( Eigen::Matrix3d( ) << 0.0, 0.0, 0.0,
                  0.0, -std::sin( currentAngleCorrections.y( ) ), -std::cos( currentAngleCorrections.y( ) ),
                  0.0, std::cos( currentAngleCorrections.y( ) ), -std::sin( currentAngleCorrections.y( ) ) ).finished( ) *
                ( Eigen::Matrix3d( ) << -std::sin( currentAngleCorrections.z( ) ), -std::cos( currentAngleCorrections.z( ) ), 0.0,
                  std::cos( currentAngleCorrections.z( ) ), -std::sin( currentAngleCorrections.z( ) ), 0.0,
                  0.0, 0.0, 0.0 ).finished( ) *
                polarMotionRotation.toRotationMatrix( );

        partialsOfRotationMatrix +=
                currentMeanPhiAngleDerivative *
                ( Psi_m * std::sin( correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) ) *
                  std::sin( currentAngleCorrections.y( ) ) *
                  std::cos(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) ) *
                  a_m / ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                  ( a_m * correctionIterator->second.first + freeCoreNutationRate *
                    correctionIterator->second.second * std::sin( angleIAtEpoch ) ) -
                  std::sin(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) )*
                  a_m / ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                  ( a_m * correctionIterator->second.second + freeCoreNutationRate *
                    correctionIterator->second.first / std::sin( angleIAtEpoch ) ) *
                  std::cos( currentAngleCorrections.y( ) ) ) *
                rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) *
                  Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) ).toRotationMatrix( ) *
                ( Eigen::Matrix3d( ) << -std::cos( currentAngleCorrections.z( ) ), std::sin( currentAngleCorrections.z( ) ), 0.0,
                  -std::sin( currentAngleCorrections.z( ) ), -std::cos( currentAngleCorrections.z( ) ), 0.0,
                  0.0, 0.0, 0.0 ).finished( ) *
                polarMotionRotation.toRotationMatrix( );
    }

    double currentPhaseAngleCorrection;

    for( unsigned int i = 0; i < meanMotionTimeDependentPhaseNutationCorrections.size( ); i++ )
    {
        currentPhaseAngleCorrection = phaseAngleCorrectionFunctions[ i ]( ephemerisTime );
        for( std::map< double, std::pair< double, double > >::iterator correctionIterator =
            meanMotionTimeDependentPhaseNutationCorrections[ i ].begin( );
            correctionIterator != meanMotionTimeDependentPhaseNutationCorrections[ i ].end( ); correctionIterator++ )
        {
            double a_m = correctionIterator->first * meanMotion;

            double Psi_m = correctionIterator->second.second + coreFactor * a_m /
                    ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                    ( a_m * correctionIterator->second.second + freeCoreNutationRate * correctionIterator->second.first /
                       std::sin( angleIAtEpoch ) );

            partialsOfRotationMatrix +=
                    currentMeanPhiAngleDerivative *
                    std::sin( correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) +
                              currentPhaseAngleCorrection ) *
                    a_m / ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                    ( a_m * correctionIterator->second.second + freeCoreNutationRate *
                      correctionIterator->second.first / std::sin( angleIAtEpoch ) ) *
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << -std::sin( currentAngleCorrections.x( ) ), -std::cos( currentAngleCorrections.x( ) ), 0.0,
                      std::cos( currentAngleCorrections.x( ) ), -std::sin( currentAngleCorrections.x( ) ), 0.0,
                      0.0, 0.0, 0.0 ).finished( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) ).toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << -std::sin( currentAngleCorrections.z( ) ), -std::cos( currentAngleCorrections.z( ) ), 0.0,
                      std::cos( currentAngleCorrections.z( ) ), -std::sin( currentAngleCorrections.z( ) ), 0.0,
                      0.0, 0.0, 0.0 ).finished( ) *
                    polarMotionRotation.toRotationMatrix( );

            partialsOfRotationMatrix +=
                    currentMeanPhiAngleDerivative *
                    std::cos( correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) +
                              currentPhaseAngleCorrection ) *
                    a_m / ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                    ( a_m * correctionIterator->second.first + freeCoreNutationRate *
                      correctionIterator->second.second * std::sin( angleIAtEpoch ) ) *
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) ).toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << 0.0, 0.0, 0.0,
                      0.0, -std::sin( currentAngleCorrections.y( ) ), -std::cos( currentAngleCorrections.y( ) ),
                      0.0, std::cos( currentAngleCorrections.y( ) ), -std::sin( currentAngleCorrections.y( ) ) ).finished( ) *
                    ( Eigen::Matrix3d( ) << -std::sin( currentAngleCorrections.z( ) ), -std::cos( currentAngleCorrections.z( ) ), 0.0,
                      std::cos( currentAngleCorrections.z( ) ), -std::sin( currentAngleCorrections.z( ) ), 0.0,
                      0.0, 0.0, 0.0 ).finished( ) *
                    polarMotionRotation.toRotationMatrix( );

            partialsOfRotationMatrix +=
                    currentMeanPhiAngleDerivative *
                    ( Psi_m * std::sin( correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) +
                                        currentPhaseAngleCorrection ) *
                      std::sin( currentAngleCorrections.y( ) ) *
                      std::cos( correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) +
                                currentPhaseAngleCorrection ) *
                      a_m / ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                      ( a_m * correctionIterator->second.first + freeCoreNutationRate *
                        correctionIterator->second.second * std::sin( angleIAtEpoch ) ) -
                      std::sin( correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) +
                                currentPhaseAngleCorrection )*
                      a_m / ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                      ( a_m * correctionIterator->second.second + freeCoreNutationRate *
                        correctionIterator->second.first / std::sin( angleIAtEpoch ) ) *
                      std::cos( currentAngleCorrections.y( ) ) ) *
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) *
                      Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) ).toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << -std::cos( currentAngleCorrections.z( ) ), std::sin( currentAngleCorrections.z( ) ), 0.0,
                      -std::sin( currentAngleCorrections.z( ) ), -std::cos( currentAngleCorrections.z( ) ), 0.0,
                      0.0, 0.0, 0.0 ).finished( ) *
                    polarMotionRotation.toRotationMatrix( );
        }
    }

        return partialsOfRotationMatrix;

}

//! Function to calculate a partial of rotation matrix from a body-fixed to inertial frame w.r.t.
//! free core nutation rate.
Eigen::Matrix3d calculatePartialOfRotationMatrixFromLocalFrameWrtFreeCoreNutationRate(
        const std::shared_ptr< ephemerides::PlanetaryOrientationAngleCalculator > planetaryOrientationCalculator,
        const Eigen::Quaterniond &rotationFromMeanOrbitToIcrf,
        const Eigen::Quaterniond& polarMotionRotation,
        const double ephemerisTime )
{
    double meanMotion = planetaryOrientationCalculator->getBodyMeanMotion( );

    double bodyMeanAnomalyAtEpoch = planetaryOrientationCalculator->getBodyMeanAnomalyAtEpoch( );

    std::map< double, std::pair< double, double > > meanMotionDirectNutationCorrections =
            planetaryOrientationCalculator->getMeanMotionDirectNutationCorrections( );

    std::vector< std::map< double, std::pair< double, double > > > meanMotionTimeDependentPhaseNutationCorrections =
            planetaryOrientationCalculator->getMeanMotionTimeDependentPhaseNutationCorrections( );

    std::vector< std::function< double( const double ) > > phaseAngleCorrectionFunctions =
            planetaryOrientationCalculator->getphaseAngleCorrectionFunctions( );

    double coreFactor = planetaryOrientationCalculator->getCorefactor( );

    double freeCoreNutationRate = planetaryOrientationCalculator->getFreeCoreNutationRate( );

    double angleIAtEpoch = planetaryOrientationCalculator->getAngleIAtEpoch( );

    Eigen::Vector3d currentAngleCorrections = planetaryOrientationCalculator->updateAndGetRotationAngles( ephemerisTime );

    Eigen::Matrix3d partialsOfRotationMatrix = Eigen::Matrix3d::Zero( );

    for( std::map< double, std::pair< double, double > >::iterator correctionIterator = meanMotionDirectNutationCorrections.begin( );
        correctionIterator != meanMotionDirectNutationCorrections.end( ); correctionIterator++ )
    {
        double a_m = correctionIterator->first * meanMotion;

        double Psi_m = correctionIterator->second.second + coreFactor * a_m /
                ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                ( a_m * correctionIterator->second.second + freeCoreNutationRate *
                  correctionIterator->second.first / std::sin( angleIAtEpoch ) );

        partialsOfRotationMatrix +=
                std::sin(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) ) *
                (coreFactor * a_m * ( correctionIterator->second.first * ( a_m * a_m + freeCoreNutationRate * freeCoreNutationRate ) /
                  std::sin( angleIAtEpoch ) + 2 * freeCoreNutationRate * a_m * correctionIterator->second.second  ) /
                 ( ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) *
                   ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) ) ) *
                rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                ( Eigen::Matrix3d( ) << -std::sin( currentAngleCorrections.x( ) ), -std::cos( currentAngleCorrections.x( ) ), 0.0,
                  std::cos( currentAngleCorrections.x( ) ), -std::sin( currentAngleCorrections.x( ) ), 0.0,
                  0.0, 0.0, 0.0 ).finished( ) *
                ( Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) *
                  Eigen::AngleAxisd( currentAngleCorrections.z( ), Eigen::Vector3d::UnitZ( ) ) ).toRotationMatrix( ) *
                polarMotionRotation.toRotationMatrix( );

        partialsOfRotationMatrix +=
                std::cos(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) ) *
                (coreFactor * a_m * ( 2 * correctionIterator->second.first * freeCoreNutationRate * a_m +
                  std::sin( angleIAtEpoch ) * ( a_m * a_m + freeCoreNutationRate * freeCoreNutationRate ) *
                  correctionIterator->second.second ) / ( ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) *
                  ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) ) ) *
                rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) ).toRotationMatrix( ) *
                ( Eigen::Matrix3d( ) << 0.0, 0.0, 0.0,
                  0.0, -std::sin( currentAngleCorrections.y( ) ), -std::cos( currentAngleCorrections.y( ) ),
                  0.0, std::cos( currentAngleCorrections.y( ) ), -std::sin( currentAngleCorrections.y( ) ) ).finished( ) *
                ( Eigen::AngleAxisd( currentAngleCorrections.z( ), Eigen::Vector3d::UnitZ( ) ) ).toRotationMatrix( ) *
                polarMotionRotation.toRotationMatrix( );

        partialsOfRotationMatrix +=
                ( Psi_m * std::sin( correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) ) *
                  std::sin( currentAngleCorrections.y( ) ) *
                  std::cos(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) ) *
                  (coreFactor * a_m * ( 2 * correctionIterator->second.first * freeCoreNutationRate * a_m +
                    std::sin( angleIAtEpoch ) * ( a_m * a_m + freeCoreNutationRate * freeCoreNutationRate ) *
                    correctionIterator->second.second ) / ( ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) *
                    ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) ) ) -
                  std::sin(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) )*
                  (coreFactor * a_m * ( correctionIterator->second.first * ( a_m * a_m + freeCoreNutationRate * freeCoreNutationRate ) /
                    std::sin( angleIAtEpoch ) + 2 * freeCoreNutationRate * a_m * correctionIterator->second.second  ) /
                   ( ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) *
                     ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) ) ) *
                  std::cos( currentAngleCorrections.y( ) ) ) *
                rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) *
                  Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) ).toRotationMatrix( ) *
                ( Eigen::Matrix3d( ) << -std::sin( currentAngleCorrections.z( ) ), -std::cos( currentAngleCorrections.z( ) ), 0.0,
                  std::cos( currentAngleCorrections.z( ) ), -std::sin( currentAngleCorrections.z( ) ), 0.0,
                  0.0, 0.0, 0.0 ).finished( ) *
                polarMotionRotation.toRotationMatrix( );
    }

    double currentPhaseAngleCorrection;

    for( unsigned int i = 0; i < meanMotionTimeDependentPhaseNutationCorrections.size( ); i++ )
    {
        currentPhaseAngleCorrection = phaseAngleCorrectionFunctions[ i ]( ephemerisTime );
        for( std::map< double, std::pair< double, double > >::iterator correctionIterator =
            meanMotionTimeDependentPhaseNutationCorrections[ i ].begin( );
            correctionIterator != meanMotionTimeDependentPhaseNutationCorrections[ i ].end( ); correctionIterator++ )
        {
            double a_m = correctionIterator->first * meanMotion;

            double Psi_m = correctionIterator->second.second + coreFactor * a_m /
                    ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                    ( a_m * correctionIterator->second.second + freeCoreNutationRate *
                      correctionIterator->second.first / std::sin( angleIAtEpoch ) );

            partialsOfRotationMatrix +=
                    std::sin(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) +
                             currentPhaseAngleCorrection ) *
                    (coreFactor * a_m * ( correctionIterator->second.first *
                                          ( a_m * a_m + freeCoreNutationRate * freeCoreNutationRate ) /
                      std::sin( angleIAtEpoch ) + 2 * freeCoreNutationRate * a_m * correctionIterator->second.second  ) /
                     ( ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) *
                       ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) ) ) *
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << -std::sin( currentAngleCorrections.x( ) ), -std::cos( currentAngleCorrections.x( ) ), 0.0,
                      std::cos( currentAngleCorrections.x( ) ), -std::sin( currentAngleCorrections.x( ) ), 0.0,
                      0.0, 0.0, 0.0 ).finished( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) *
                      Eigen::AngleAxisd( currentAngleCorrections.z( ), Eigen::Vector3d::UnitZ( ) ) ).toRotationMatrix( ) *
                    polarMotionRotation.toRotationMatrix( );

            partialsOfRotationMatrix +=
                    std::cos(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) +
                             currentPhaseAngleCorrection ) *
                    (coreFactor * a_m * ( 2 * correctionIterator->second.first * freeCoreNutationRate * a_m +
                      std::sin( angleIAtEpoch ) * ( a_m * a_m + freeCoreNutationRate * freeCoreNutationRate ) *
                      correctionIterator->second.second ) / ( ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) *
                      ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) ) ) *
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) ).toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << 0.0, 0.0, 0.0,
                      0.0, -std::sin( currentAngleCorrections.y( ) ), -std::cos( currentAngleCorrections.y( ) ),
                      0.0, std::cos( currentAngleCorrections.y( ) ), -std::sin( currentAngleCorrections.y( ) ) ).finished( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.z( ), Eigen::Vector3d::UnitZ( ) ) ).toRotationMatrix( ) *
                    polarMotionRotation.toRotationMatrix( );

            partialsOfRotationMatrix +=
                    ( Psi_m * std::sin( correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) +
                                        currentPhaseAngleCorrection ) *
                      std::sin( currentAngleCorrections.y( ) ) *
                      std::cos(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) +
                               currentPhaseAngleCorrection ) *
                      (coreFactor * a_m * ( 2 * correctionIterator->second.first * freeCoreNutationRate * a_m +
                        std::sin( angleIAtEpoch ) * ( a_m * a_m + freeCoreNutationRate * freeCoreNutationRate ) *
                        correctionIterator->second.second ) / ( ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) *
                        ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) ) ) -
                      std::sin(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) +
                               currentPhaseAngleCorrection ) *
                      (coreFactor * a_m * ( correctionIterator->second.first * ( a_m * a_m + freeCoreNutationRate * freeCoreNutationRate ) /
                        std::sin( angleIAtEpoch ) + 2 * freeCoreNutationRate * a_m * correctionIterator->second.second  ) /
                       ( ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) *
                         ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) ) ) *
                      std::cos( currentAngleCorrections.y( ) ) ) *
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) *
                      Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) ).toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << -std::sin( currentAngleCorrections.z( ) ), -std::cos( currentAngleCorrections.z( ) ), 0.0,
                      std::cos( currentAngleCorrections.z( ) ), -std::sin( currentAngleCorrections.z( ) ), 0.0,
                      0.0, 0.0, 0.0 ).finished( ) *
                    polarMotionRotation.toRotationMatrix( );
        }
    }

    return partialsOfRotationMatrix;

}

//! Function to calculate a partial of rotation matrix derivative from a body-fixed to inertial frame w.r.t.
//! free core nutation rate.
Eigen::Matrix3d calculatePartialOfRotationMatrixFromLocalFrameDerivativeWrtFreeCoreNutationRate(
        const std::shared_ptr< ephemerides::PlanetaryOrientationAngleCalculator > planetaryOrientationCalculator,
        const Eigen::Quaterniond &rotationFromMeanOrbitToIcrf,
        const Eigen::Quaterniond& polarMotionRotation,
        const double ephemerisTime )
{
    double meanMotion = planetaryOrientationCalculator->getBodyMeanMotion( );

    double bodyMeanAnomalyAtEpoch = planetaryOrientationCalculator->getBodyMeanAnomalyAtEpoch( );

    std::map< double, std::pair< double, double > > meanMotionDirectNutationCorrections =
            planetaryOrientationCalculator->getMeanMotionDirectNutationCorrections( );

    std::vector< std::map< double, std::pair< double, double > > > meanMotionTimeDependentPhaseNutationCorrections =
            planetaryOrientationCalculator->getMeanMotionTimeDependentPhaseNutationCorrections( );

    std::vector< std::function< double( const double ) > > phaseAngleCorrectionFunctions =
            planetaryOrientationCalculator->getphaseAngleCorrectionFunctions( );

    double coreFactor = planetaryOrientationCalculator->getCorefactor( );

    double freeCoreNutationRate = planetaryOrientationCalculator->getFreeCoreNutationRate( );

    double angleIAtEpoch = planetaryOrientationCalculator->getAngleIAtEpoch( );

    Eigen::Vector3d currentAngleCorrections = planetaryOrientationCalculator->updateAndGetRotationAngles( ephemerisTime );

    double currentMeanPhiAngleDerivative = planetaryOrientationCalculator->getcurrentMeanPhiAngleDerivative( ephemerisTime );

    Eigen::Matrix3d partialsOfRotationMatrix = Eigen::Matrix3d::Zero( );

    for( std::map< double, std::pair< double, double > >::iterator correctionIterator = meanMotionDirectNutationCorrections.begin( );
        correctionIterator != meanMotionDirectNutationCorrections.end( ); correctionIterator++ )
    {
        double a_m = correctionIterator->first * meanMotion;

        double Psi_m = correctionIterator->second.second + coreFactor * a_m /
                ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                ( a_m * correctionIterator->second.second + freeCoreNutationRate *
                  correctionIterator->second.first / std::sin( angleIAtEpoch ) );

        partialsOfRotationMatrix +=
                currentMeanPhiAngleDerivative *
                std::sin(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) ) *
                (coreFactor * a_m * ( correctionIterator->second.first * ( a_m * a_m + freeCoreNutationRate * freeCoreNutationRate ) /
                  std::sin( angleIAtEpoch ) + 2 * freeCoreNutationRate * a_m * correctionIterator->second.second  ) /
                 ( ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) *
                   ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) ) ) *
                rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                ( Eigen::Matrix3d( ) << -std::sin( currentAngleCorrections.x( ) ), -std::cos( currentAngleCorrections.x( ) ), 0.0,
                  std::cos( currentAngleCorrections.x( ) ), -std::sin( currentAngleCorrections.x( ) ), 0.0,
                  0.0, 0.0, 0.0 ).finished( ) *
                ( Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) ).toRotationMatrix( ) *
                ( Eigen::Matrix3d( ) << -std::sin( currentAngleCorrections.z( ) ), -std::cos( currentAngleCorrections.z( ) ), 0.0,
                  std::cos( currentAngleCorrections.z( ) ), -std::sin( currentAngleCorrections.z( ) ), 0.0,
                  0.0, 0.0, 0.0 ).finished( ) *
                polarMotionRotation.toRotationMatrix( );

        partialsOfRotationMatrix +=
                currentMeanPhiAngleDerivative *
                std::cos(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) ) *
                (coreFactor * a_m * ( 2 * correctionIterator->second.first * freeCoreNutationRate * a_m +
                  std::sin( angleIAtEpoch ) * ( a_m * a_m + freeCoreNutationRate * freeCoreNutationRate ) *
                  correctionIterator->second.second ) / ( ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) *
                  ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) ) ) *
                rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) ).toRotationMatrix( ) *
                ( Eigen::Matrix3d( ) << 0.0, 0.0, 0.0,
                  0.0, -std::sin( currentAngleCorrections.y( ) ), -std::cos( currentAngleCorrections.y( ) ),
                  0.0, std::cos( currentAngleCorrections.y( ) ), -std::sin( currentAngleCorrections.y( ) ) ).finished( ) *
                ( Eigen::Matrix3d( ) << -std::sin( currentAngleCorrections.z( ) ), -std::cos( currentAngleCorrections.z( ) ), 0.0,
                  std::cos( currentAngleCorrections.z( ) ), -std::sin( currentAngleCorrections.z( ) ), 0.0,
                  0.0, 0.0, 0.0 ).finished( ) *
                polarMotionRotation.toRotationMatrix( );

        partialsOfRotationMatrix +=
                currentMeanPhiAngleDerivative *
                ( Psi_m * std::sin( correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) ) *
                  std::sin( currentAngleCorrections.y( ) ) *
                  std::cos(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) ) *
                  (coreFactor * a_m * ( 2 * correctionIterator->second.first * freeCoreNutationRate * a_m +
                    std::sin( angleIAtEpoch ) * ( a_m * a_m + freeCoreNutationRate * freeCoreNutationRate ) *
                    correctionIterator->second.second ) / ( ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) *
                    ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) ) ) -
                  std::sin(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) )*
                  (coreFactor * a_m * ( correctionIterator->second.first * ( a_m * a_m + freeCoreNutationRate * freeCoreNutationRate ) /
                    std::sin( angleIAtEpoch ) + 2 * freeCoreNutationRate * a_m * correctionIterator->second.second  ) /
                   ( ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) *
                     ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) ) ) *
                  std::cos( currentAngleCorrections.y( ) ) ) *
                rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) *
                  Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) ).toRotationMatrix( ) *
                ( Eigen::Matrix3d( ) << -std::cos( currentAngleCorrections.z( ) ), std::sin( currentAngleCorrections.z( ) ), 0.0,
                  -std::sin( currentAngleCorrections.z( ) ), -std::cos( currentAngleCorrections.z( ) ), 0.0,
                  0.0, 0.0, 0.0 ).finished( ) *
                polarMotionRotation.toRotationMatrix( );


    }

    double currentPhaseAngleCorrection;

    for( unsigned int i = 0; i < meanMotionTimeDependentPhaseNutationCorrections.size( ); i++ )
    {
        currentPhaseAngleCorrection = phaseAngleCorrectionFunctions[ i ]( ephemerisTime );
        for( std::map< double, std::pair< double, double > >::iterator correctionIterator =
            meanMotionTimeDependentPhaseNutationCorrections[ i ].begin( );
            correctionIterator != meanMotionTimeDependentPhaseNutationCorrections[ i ].end( ); correctionIterator++ )
        {
            double a_m = correctionIterator->first * meanMotion;

            double Psi_m = correctionIterator->second.second + coreFactor * a_m /
                    ( a_m * a_m - freeCoreNutationRate * freeCoreNutationRate ) *
                    ( a_m * correctionIterator->second.second + freeCoreNutationRate *
                      correctionIterator->second.first / std::sin( angleIAtEpoch ) );

            partialsOfRotationMatrix +=
                    currentMeanPhiAngleDerivative *
                    std::sin(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) +
                             currentPhaseAngleCorrection ) *
                    (coreFactor * a_m * ( correctionIterator->second.first * ( a_m * a_m + freeCoreNutationRate * freeCoreNutationRate ) /
                      std::sin( angleIAtEpoch ) + 2 * freeCoreNutationRate * a_m * correctionIterator->second.second  ) /
                     ( ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) *
                       ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) ) ) *
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << -std::sin( currentAngleCorrections.x( ) ), -std::cos( currentAngleCorrections.x( ) ), 0.0,
                      std::cos( currentAngleCorrections.x( ) ), -std::sin( currentAngleCorrections.x( ) ), 0.0,
                      0.0, 0.0, 0.0 ).finished( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) ).toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << -std::sin( currentAngleCorrections.z( ) ), -std::cos( currentAngleCorrections.z( ) ), 0.0,
                      std::cos( currentAngleCorrections.z( ) ), -std::sin( currentAngleCorrections.z( ) ), 0.0,
                      0.0, 0.0, 0.0 ).finished( ) *
                    polarMotionRotation.toRotationMatrix( );

            partialsOfRotationMatrix +=
                    currentMeanPhiAngleDerivative *
                    std::cos(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) +
                             currentPhaseAngleCorrection ) *
                    (coreFactor * a_m * ( 2 * correctionIterator->second.first * freeCoreNutationRate * a_m +
                      std::sin( angleIAtEpoch ) * ( a_m * a_m + freeCoreNutationRate * freeCoreNutationRate ) *
                      correctionIterator->second.second ) / ( ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) *
                      ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) ) ) *
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) ).toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << 0.0, 0.0, 0.0,
                      0.0, -std::sin( currentAngleCorrections.y( ) ), -std::cos( currentAngleCorrections.y( ) ),
                      0.0, std::cos( currentAngleCorrections.y( ) ), -std::sin( currentAngleCorrections.y( ) ) ).finished( ) *
                    ( Eigen::Matrix3d( ) << -std::sin( currentAngleCorrections.z( ) ), -std::cos( currentAngleCorrections.z( ) ), 0.0,
                      std::cos( currentAngleCorrections.z( ) ), -std::sin( currentAngleCorrections.z( ) ), 0.0,
                      0.0, 0.0, 0.0 ).finished( ) *
                    polarMotionRotation.toRotationMatrix( );

            partialsOfRotationMatrix +=
                    currentMeanPhiAngleDerivative *
                    ( Psi_m * std::sin( correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) +
                                        currentPhaseAngleCorrection ) *
                      std::sin( currentAngleCorrections.y( ) ) *
                      std::cos(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) +
                               currentPhaseAngleCorrection ) *
                      (coreFactor * a_m * ( 2 * correctionIterator->second.first * freeCoreNutationRate * a_m +
                        std::sin( angleIAtEpoch ) * ( a_m * a_m + freeCoreNutationRate * freeCoreNutationRate ) *
                        correctionIterator->second.second ) / ( ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) *
                        ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) ) ) -
                      std::sin(correctionIterator->first * ( meanMotion * ephemerisTime + bodyMeanAnomalyAtEpoch ) +
                               currentPhaseAngleCorrection )*
                      (coreFactor * a_m * ( correctionIterator->second.first * ( a_m * a_m + freeCoreNutationRate * freeCoreNutationRate ) /
                        std::sin( angleIAtEpoch ) + 2 * freeCoreNutationRate * a_m * correctionIterator->second.second  ) /
                       ( ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) *
                         ( freeCoreNutationRate * freeCoreNutationRate - a_m * a_m ) ) ) *
                      std::cos( currentAngleCorrections.y( ) ) ) *
                    rotationFromMeanOrbitToIcrf.toRotationMatrix( ) *
                    ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) *
                      Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) ).toRotationMatrix( ) *
                    ( Eigen::Matrix3d( ) << -std::cos( currentAngleCorrections.z( ) ), std::sin( currentAngleCorrections.z( ) ), 0.0,
                      -std::sin( currentAngleCorrections.z( ) ), -std::cos( currentAngleCorrections.z( ) ), 0.0,
                      0.0, 0.0, 0.0 ).finished( ) *
                    polarMotionRotation.toRotationMatrix( );

        }
    }

        return partialsOfRotationMatrix;

}

//! Function to calculate the partial of the position of a vector, which is given in a body-fixed frame, in the inertial
//! frame wrt a parameter.
Eigen::Matrix< double, 3, Eigen::Dynamic > RotationMatrixPartial::calculatePartialOfInertialPositionWrtParameter(
        const double time,
        const Eigen::Vector3d vectorInLocalFrame )
{
    std::vector< Eigen::Matrix3d > rotationMatrixPartials = calculatePartialOfRotationMatrixToBaseFrameWrParameter( time );
    Eigen::Matrix< double, 3, Eigen::Dynamic > rotatedVectorPartial = Eigen::Matrix< double, 3, Eigen::Dynamic >::Zero(
                3, rotationMatrixPartials.size( ) );

    for( unsigned int i = 0; i < rotationMatrixPartials.size( ); i++ )
    {
        rotatedVectorPartial.block( 0, i, 3, 1 ) = rotationMatrixPartials[ i ] * vectorInLocalFrame;
    }
    return rotatedVectorPartial;
}

//! Function to calculate the partial of the velocity of a vector, which is given in a body-fixed frame, in the inertial
//! frame wrt a parameter.
Eigen::Matrix< double, 3, Eigen::Dynamic > RotationMatrixPartial::calculatePartialOfInertialVelocityWrtParameter(
        const double time,
        const Eigen::Vector3d vectorInLocalFrame )
{
    if( rotationModel_ == nullptr )
    {
        throw std::runtime_error( "Error when caling RotationMatrixPartial::calculatePartialOfInertialVelocityWrtParameter, rotation model is nullptr" );
    }

    // Compute rotation matrix (derivative) partials
    std::vector< Eigen::Matrix3d > rotationMatrixPartials =
            calculatePartialOfRotationMatrixToBaseFrameWrParameter( time );
    std::vector< Eigen::Matrix3d > rotationMatrixDerivativePartials =
            calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter( time );

    // Compute current rotation matrix and derivative
    Eigen::Matrix3d currentRotationToBaseFrame = rotationModel_->getRotationToBaseFrame( time ).toRotationMatrix( );
    Eigen::Matrix3d currentRotationToTargetFrameDerivative = rotationModel_->getDerivativeOfRotationToTargetFrame(
                time );

    Eigen::Matrix< double, 3, Eigen::Dynamic > rotatedVectorPartial = Eigen::Matrix< double, 3, Eigen::Dynamic >::Zero(
                3, rotationMatrixPartials.size( ) );

    // Compute inertial velocity partial
    for( unsigned int i = 0; i < rotationMatrixPartials.size( ); i++ )
    {
        rotatedVectorPartial.block( 0, i, 3, 1 ) =
                -( rotationMatrixPartials.at( i ) * currentRotationToTargetFrameDerivative * currentRotationToBaseFrame +
                 currentRotationToBaseFrame * rotationMatrixDerivativePartials.at( i ).transpose( ) * currentRotationToBaseFrame +
                  currentRotationToBaseFrame * currentRotationToTargetFrameDerivative * rotationMatrixPartials.at( i ) ) * vectorInLocalFrame;
    }
    return rotatedVectorPartial;
}

Eigen::Matrix< double, 1, 6 > calculatePartialOfDirectLibrationAngleWrtCartesianStates(
        const Eigen::Vector6d& currentState,
        const double scaledLibrationAmplitude )
{
    Eigen::Matrix< double, 3, 6 > testPartial;
    testPartial.block( 0, 0, 3, 3 ) = tudat::linear_algebra::getCrossProductMatrix(
                -currentState.segment( 3, 3 ) );
    testPartial.block( 0, 3, 3, 3 ) = tudat::linear_algebra::getCrossProductMatrix(
                currentState.segment( 0, 3 ) );
    Eigen::Vector3d crossProduct =
            ( currentState.segment< 3 >( 0 ).cross( currentState.segment< 3 >( 3 ) ) );
    Eigen::Vector3d preMultiplier = crossProduct / crossProduct.norm( );
    Eigen::Matrix< double, Eigen::Dynamic, 6 > totalPartial = - ( currentState.segment< 3 >( 0 ).dot( currentState.segment< 3 >( 3 ) ) ) *
            preMultiplier.transpose( ) * testPartial / ( crossProduct.norm( ) * crossProduct.norm( ) );
    totalPartial.block( 0, 0, 1, 3 ) += currentState.segment< 3 >( 3 ).transpose( ) / crossProduct.norm( );
    totalPartial.block( 0, 3, 1, 3 ) += currentState.segment< 3 >( 0 ).transpose( ) / crossProduct.norm( );
    return scaledLibrationAmplitude * totalPartial;
//    Eigen::Vector3d positionVector = currentState.segment( 0, 3 );
//    double positionNorm = positionVector.norm( );

//    Eigen::Vector3d velocityVector = currentState.segment( 3, 3 );
//    double velocityNorm = velocityVector.norm( );

//    Eigen::Vector3d crossProduct = positionVector.cross( velocityVector );


//    double crossProductNorm = crossProduct.norm( );
//    double crossProductPartialScaling =
//            positionVector.dot( velocityVector ) / ( crossProductNorm * crossProductNorm );

//    Eigen::Matrix< double, 1, 6 > angleDerivatives;
//    angleDerivatives.setZero( );

//    angleDerivatives.block( 0, 0, 1, 3 ) =
//            velocityVector.transpose( ) + crossProductPartialScaling * positionVector.transpose( ) * (
//                velocityNorm * velocityNorm * Eigen::Matrix3d::Identity( ) - velocityVector * velocityVector.transpose( ) );
//    angleDerivatives.block( 0, 3, 1, 3 ) =
//            positionVector.transpose( ) + crossProductPartialScaling * velocityVector.transpose( ) * (
//                positionNorm * positionNorm * Eigen::Matrix3d::Identity( ) - positionVector * positionVector.transpose( ) );

//    return scaledLibrationAmplitude * angleDerivatives / crossProductNorm;
}

//! Function to compute the required partial derivative of rotation matrix.
std::vector< Eigen::Matrix3d > SynchronousRotationMatrixPartialWrtTranslationalState::
calculatePartialOfRotationMatrixToBaseFrameWrParameter( const double time )
{
    // Get current state/rotation
    Eigen::Matrix3d currentRotationMatrix =
            synchronousRotationaModel_->getFullyLockedRotationToBaseFrame( time );;
    Eigen::Vector6d currentState =
            synchronousRotationaModel_->getCurrentRelativeState( time );
    Eigen::Vector3d positionVector = currentState.segment( 0, 3 );
    Eigen::Vector3d velocityVector = currentState.segment( 3, 3 );
    double positionNorm = positionVector.norm( );

    // Compute r and w vectors
    Eigen::Vector3d rVector = -currentRotationMatrix.block( 0, 0, 3, 1 );
    Eigen::Vector3d wVector = currentRotationMatrix.block( 0, 2, 3, 1 );
    Eigen::Vector3d unnormalizedWVector = positionVector.cross( velocityVector );
    double unnormalizedWVectorNorm = unnormalizedWVector.norm( );

    // Compute r vector partials
    Eigen::Matrix3d rVectorDerivativeWrtPosition =
            Eigen::Matrix3d::Identity( ) / positionNorm - positionVector * positionVector.transpose( ) / (
                positionNorm * positionNorm * positionNorm );

    // Compute partials for unnormalized w vector
    Eigen::Matrix3d unnormalizedWVectorDerivativeWrtPosition =
            -linear_algebra::getCrossProductMatrix( velocityVector );
    Eigen::Matrix3d unnormalizedWVectorDerivativeWrtVelocity =
            linear_algebra::getCrossProductMatrix( positionVector );
    Eigen::Matrix3d wPartialScalingTerm =
            ( Eigen::Matrix3d::Identity( ) / unnormalizedWVectorNorm -
              unnormalizedWVector * unnormalizedWVector.transpose( ) /
              ( unnormalizedWVectorNorm * unnormalizedWVectorNorm * unnormalizedWVectorNorm ) );

    // Compute w vector partials
    Eigen::Matrix3d wVectorDerivativeWrtPosition =
            wPartialScalingTerm * unnormalizedWVectorDerivativeWrtPosition;
    Eigen::Matrix3d wVectorDerivativeWrtVelocity =
            wPartialScalingTerm * unnormalizedWVectorDerivativeWrtVelocity;

    // Compute s vector partials
    Eigen::Matrix3d sVectorDerivativeWrtPosition =
            linear_algebra::getCrossProductMatrix( wVector ) * rVectorDerivativeWrtPosition -
            linear_algebra::getCrossProductMatrix( rVector ) * wVectorDerivativeWrtPosition;
    Eigen::Matrix3d sVectorDerivativeWrtVelocity =
            -linear_algebra::getCrossProductMatrix( rVector ) * wVectorDerivativeWrtVelocity;

    // Build vector of partial derivatives
    std::vector< Eigen::Matrix3d > rotationMatrixPartials;
    rotationMatrixPartials.resize( 6 );
    Eigen::Matrix3d correctionRotation = synchronousRotationaModel_->getLibrationRotation( time );
    for( int i = 0; i < 3; i++ )
    {
        rotationMatrixPartials[ i ].block( 0, 0, 3, 1 ) =
                -rVectorDerivativeWrtPosition.block( 0, i, 3, 1 );
        rotationMatrixPartials[ i ].block( 0, 1, 3, 1 ) =
                -sVectorDerivativeWrtPosition.block( 0, i, 3, 1 );
        rotationMatrixPartials[ i ].block( 0, 2, 3, 1 ) =
                wVectorDerivativeWrtPosition.block( 0, i, 3, 1 );
        rotationMatrixPartials[ i ] = rotationMatrixPartials[ i ] * correctionRotation;

        rotationMatrixPartials[ i + 3 ].block( 0, 0, 3, 1 ).setZero( );
        rotationMatrixPartials[ i + 3 ].block( 0, 1, 3, 1 ) =
                -sVectorDerivativeWrtVelocity.block( 0, i, 3, 1 );
        rotationMatrixPartials[ i + 3 ].block( 0, 2, 3, 1 ) =
                wVectorDerivativeWrtVelocity.block( 0, i, 3, 1 );
        rotationMatrixPartials[ i + 3 ] = rotationMatrixPartials[ i + 3 ] * correctionRotation;

        if( directLongitudeLibrationCalculator_ != nullptr )
        {
            Eigen::Matrix3d lockedRotation =
                    synchronousRotationaModel_->getFullyLockedRotationToBaseFrame( time );
            Eigen::Matrix3d librationRotationDerivativeWrtAngle =
                    reference_frames::getDerivativeOfZAxisRotationWrtAngle( correctionRotation );

            Eigen::Matrix< double, 1, 6 > librationDerivatives =
                    calculatePartialOfDirectLibrationAngleWrtCartesianStates(
                        currentState, directLongitudeLibrationCalculator_->getScaledLibrationAmplitude( ) );
            //        for( int i = 0; i < 6; i++ )
            //        {
            //            rotationMatrixPartials[ i ] += lockedRotation * librationRotationDerivativeWrtAngle * librationDerivatives( i );
            //        }
        }

    }

    return rotationMatrixPartials;

}

std::vector< Eigen::Matrix3d > NumericalRotationMatrixPartialWrtTranslationalState::
calculatePartialOfRotationMatrixToBaseFrameWrParameter( const double time )
{
    std::vector< Eigen::Matrix3d > statePartials;
    Eigen::Matrix3d upperturbedRotationMatrix, downperturbedRotationMatrix;

    Eigen::Vector6d nominalState = getStateFunction_( );
    Eigen::Vector6d upperturbedState, downperturbedState;

    for( int i = 0; i < 6; i++ )
    {
        upperturbedState = nominalState;
        upperturbedState( i ) += statePerturbations_( i );
        setStateFunction_( upperturbedState );
        upperturbedRotationMatrix = rotationMatrixToBaseFrameFunction_( time );

        downperturbedState = nominalState;
        downperturbedState( i ) -= statePerturbations_( i );
        setStateFunction_( downperturbedState );
        downperturbedRotationMatrix = rotationMatrixToBaseFrameFunction_( time );

        statePartials.push_back( -( upperturbedRotationMatrix - downperturbedRotationMatrix ) / (
                                     2.0 * statePerturbations_( i ) ) );

    }
    setStateFunction_( nominalState );

    return statePartials;
}

}

}
