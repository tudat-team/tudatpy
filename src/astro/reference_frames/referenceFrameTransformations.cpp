/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Mooij, E. The Motion of a Vehicle in a Planetary Atmosphere, TU Delft, 1997.
 *      Seidelmann, P. K. (Ed.). (2005). Explanatory supplement to the astronomical almanac.
 *              Univ Science Books.
 *
 *    Notes
 *      Because potential speed improvement it was chosen to use AngleAxisd and quaternions
 *      but to get things working, the rotation angle inputted in angleAxisd need to be inverted.
 *      In the future it might be better to change it to write out the complete transformation for
 *      clarity, or work with directional cosine matrices.
 */

#include <iostream>
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/basic/basicMathematicsFunctions.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
namespace tudat
{

namespace reference_frames
{


//! Function to compute pole right ascension and declination, as well as prime meridian of date, from rotation matrix
Eigen::Vector3d calculateInertialToPlanetFixedRotationAnglesFromMatrix(
        const Eigen::Matrix3d& rotationMatrixFromInertialToPlanetFixedFrame )
{
    Eigen::Vector3d rotationAngles;
    rotationAngles.x( ) = basic_mathematics::computeModulo(
                std::atan2( rotationMatrixFromInertialToPlanetFixedFrame( 2, 0 ),
                            -rotationMatrixFromInertialToPlanetFixedFrame( 2, 1 ) ) - mathematical_constants::PI / 2.0,
                2.0 * mathematical_constants::PI );//right ascension
    rotationAngles.y( ) = -std::acos( rotationMatrixFromInertialToPlanetFixedFrame( 2, 2 ) ) + mathematical_constants::PI / 2.0 ; // declination
    rotationAngles.z( ) = std::atan2( rotationMatrixFromInertialToPlanetFixedFrame( 0, 2 ),
                                      rotationMatrixFromInertialToPlanetFixedFrame( 1, 2 ) ); // longitude of prime meridian
    return rotationAngles;
}

//! Wrapper function to transform a vector to a different frame from a single rotation function.
Eigen::Vector3d transformVectorFromQuaternionFunction(
        const Eigen::Vector3d& originalVector,
        const std::function< Eigen::Quaterniond( ) > rotation )
{
    return rotation( ) * originalVector;
}

//! Wrapper function to transform a vector to a different frame from a single rotation function.
void transformVectorReferenceFromQuaternionFunction(
        Eigen::Vector3d& transformedVector,
        const Eigen::Vector3d& originalVector,
        const std::function< Eigen::Quaterniond( ) > rotation )
{
    transformedVector = rotation( ) * originalVector;
}


//! Wrapper function to transform a vector to a different frame from a single transformation function.
Eigen::Vector3d transformVectorFunctionFromVectorFunctions(
        const std::function< Eigen::Vector3d( ) > originalVector,
        const std::function< Eigen::Vector3d( const Eigen::Vector3d& ) > transformationFunction )
{
    Eigen::Vector3d original = originalVector( );
    return transformationFunction( original );
}

void transformVectorFunctionFromVectorReferenceFunctions(
        Eigen::Vector3d& transformedVector,
        const std::function< Eigen::Vector3d&( ) > originalVector,
        std::function< void( Eigen::Vector3d&, const Eigen::Vector3d& ) > transformationFunction )
{
    transformationFunction( transformedVector, originalVector( ) );
}

//! Wrapper function to transform a vector to a different frame from a list of transformation function.
Eigen::Vector3d transformVectorFromVectorFunctions(
        const Eigen::Vector3d& originalVector,
        const std::vector< std::function< Eigen::Vector3d( const Eigen::Vector3d& ) > >& rotationsList )
{
    Eigen::Vector3d currentVector = originalVector;
    Eigen::Vector3d newVector;

    // Apply each of the required tranformations.
    for( unsigned int i = 0; i < rotationsList.size( ); i++ )
    {
        newVector = rotationsList.at( i )( currentVector );
        currentVector = newVector;
    }
    return currentVector;
}

void transformVectorReferenceFromVectorFunctions(
        Eigen::Vector3d& transformedVector,
        const Eigen::Vector3d& originalVector,
        const std::vector< std::function< Eigen::Vector3d( const Eigen::Vector3d& ) > >& rotationsList )
{
    transformedVector = originalVector;
    Eigen::Vector3d newVector;

    // Apply each of the required tranformations.
    for( unsigned int i = 0; i < rotationsList.size( ); i++ )
    {
        newVector = rotationsList.at( i )( transformedVector );
        transformedVector = newVector;
    }
}


//! Get rotating planetocentric (R) to inertial (I) reference frame transformation matrix.
Eigen::Matrix3d getRotatingPlanetocentricToInertialFrameTransformationMatrix( const double angleFromXItoXR )
{
    // Declare local variables.
    // Declare local matrix.
    Eigen::Matrix3d localMatrix_;

    // Set local matrix.
    localMatrix_ = reference_frames::
            getInertialToPlanetocentricFrameTransformationMatrix( angleFromXItoXR );

    // Return transformation matrix.
    return localMatrix_.transpose( );
}

//! Get rotating planetocentric (R) to inertial (I) reference frame transformation quaternion.
Eigen::Quaterniond getRotatingPlanetocentricToInertialFrameTransformationQuaternion(
        const double angleFromXItoXR )
{
    // Compute transformation quaternion
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd eigenRotationObject = Eigen::AngleAxisd( -1.0 * -angleFromXItoXR,
                                                               Eigen::Vector3d::UnitZ( ) );
    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond( eigenRotationObject );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}

//! Get rotation from planet-fixed to inertial frame.
Eigen::Quaterniond getRotatingPlanetocentricToInertialFrameTransformationQuaternion(
        const double declinationOfPole,
        const double rightAscensionOfPole,
        const double longitudeOfPrimeMeridian )
{
    // Compute transformation quaternion.
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd firstRotationAroundZaxis =
            Eigen::AngleAxisd( longitudeOfPrimeMeridian, Eigen::Vector3d::UnitZ( ) );
    Eigen::AngleAxisd rotationAroundXaxis =
            Eigen::AngleAxisd(
                ( mathematical_constants::PI / 2.0 - declinationOfPole ),
                Eigen::Vector3d::UnitX( ) );
    Eigen::AngleAxisd secondRotationAroundZaxis = Eigen::AngleAxisd(
                rightAscensionOfPole
                + mathematical_constants::PI / 2.0, Eigen::Vector3d::UnitZ( ) );
    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond(
                ( secondRotationAroundZaxis * rotationAroundXaxis *  firstRotationAroundZaxis ) );
    return frameTransformationQuaternion;
}

Eigen::Matrix3d getRotatingPlanetocentricToInertialFrameTransformationMatrix(
        const double declinationOfPole,
        const double rightAscensionOfPole,
        const double longitudeOfPrimeMeridian )
{
    return getRotatingPlanetocentricToInertialFrameTransformationQuaternion(
                declinationOfPole, rightAscensionOfPole, longitudeOfPrimeMeridian ).toRotationMatrix( );
}

//! Get inertial (I) to rotating planetocentric (R) reference frame transformtion matrix.
Eigen::Matrix3d getInertialToPlanetocentricFrameTransformationMatrix(
        const double angleFromXItoXR )
{
    // Compute rotation about Z-Axis.
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd eigenRotationObject = Eigen::AngleAxisd( -1.0 * angleFromXItoXR,
                                                               Eigen::Vector3d::UnitZ( ) );

    // Return transformation matrix.
    return eigenRotationObject.toRotationMatrix( );
}

//! Get rotation from velocity based TNW frame to inertial frame (I) frame.
Eigen::Matrix3d getTnwToInertialRotation(
        const Eigen::Vector6d& vehicleState,
        const Eigen::Vector6d& centralBodyState,
        const bool doesNaxisPointAwayFromCentralBody )
{
    Eigen::Vector3d vehicleVelocity, vehicleRadius;
    vehicleRadius = vehicleState.segment( 0, 3 ) - centralBodyState.segment( 0, 3 );
    vehicleVelocity = vehicleState.segment( 3, 3 ) - centralBodyState.segment( 3, 3 );

    Eigen::Vector3d unitT = vehicleVelocity.normalized( );// / vehicleVelocity.norm( );
    if ( vehicleRadius.cross( vehicleVelocity ).norm( ) == 0.0 )
    {
        std::string errorMessage = "Division by zero: radius and velocity are in the same direction.";
        throw std::runtime_error( errorMessage );
    }

    Eigen::Vector3d unitW =  ( ( ( doesNaxisPointAwayFromCentralBody ) ? -1.0 : 1.0 ) *
                               ( vehicleRadius.cross( vehicleVelocity ) ).normalized( ) );

    Eigen::Vector3d unitN = ( unitW.cross( unitT ) ).normalized( );

    Eigen::Matrix3d transformationMatrix;
    transformationMatrix << unitT( 0 ), unitN( 0 ), unitW( 0 ),
            unitT( 1 ), unitN( 1 ), unitW( 1 ),
            unitT( 2 ), unitN( 2 ), unitW( 2 );

    return transformationMatrix;
}

Eigen::Matrix3d getTnwToInertialRotation(const Eigen::Vector6d& vehicleInertialState,
                                         const bool doesNaxisPointAwayFromCentralBody )
{
    return getTnwToInertialRotation(
                vehicleInertialState, Eigen::Vector6d::Zero( ), doesNaxisPointAwayFromCentralBody );
}

Eigen::Matrix3d getInertialToTnwRotation(const Eigen::Vector6d& vehicleInertialState,
                                         const bool doesNaxisPointAwayFromCentralBody )
{
    return getTnwToInertialRotation(
                vehicleInertialState, Eigen::Vector6d::Zero( ), doesNaxisPointAwayFromCentralBody ).transpose( );
}


//! Get rotation from velocity based TNW frame to inertial frame (I) frame.
Eigen::Matrix3d getTnwToInertialRotationFromFunctions(
        const std::function< Eigen::Vector6d( ) >& vehicleStateFunction,
        const std::function< Eigen::Vector6d( ) >& centralBodyStateFunction,
        const bool doesNaxisPointAwayFromCentralBody )
{
    return getTnwToInertialRotation(
                vehicleStateFunction( ), centralBodyStateFunction( ), doesNaxisPointAwayFromCentralBody );
}

//! Get rotation from velocity based TNW frame to planetocentric frame.
Eigen::Quaterniond getTnwToPlanetocentricRotationKeplerian(
        const Eigen::Matrix< double, 6, 1 > spacecraftKeplerianState )
{
    double eccentricity = spacecraftKeplerianState( 1 );
    double inclination = spacecraftKeplerianState( 2 );
    double argumentOfPeriapsis = spacecraftKeplerianState( 3 );
    double rightAscensionOfAscendingNode = spacecraftKeplerianState( 4 );
    double trueAnomaly = spacecraftKeplerianState( 5 );

    double flightPathAngle = std::atan( ( eccentricity * std::sin( trueAnomaly ) ) /
                                        ( 1.0 + eccentricity * std::cos( trueAnomaly ) ) );

    // Compute first rotation around Z axis.
    Eigen::AngleAxisd firstRotationAroundZaxis(
                -( -mathematical_constants::PI * 0.5 + flightPathAngle - ( trueAnomaly + argumentOfPeriapsis ) ),
                Eigen::Vector3d::UnitZ( ) );

    // Compute rotation around X axis.
    Eigen::AngleAxisd rotationAroundXaxis( inclination, Eigen::Vector3d::UnitX( ) );

    // Compute second rotation around Z axis.
    Eigen::AngleAxisd secondRotationAroundZaxis( rightAscensionOfAscendingNode, Eigen::Vector3d::UnitZ( ) );

    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond(
                ( secondRotationAroundZaxis * rotationAroundXaxis * firstRotationAroundZaxis) );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}

//! Function to compute the rotation matrix to RSW frame, from the frame in which the input state is given.
Eigen::Matrix3d getInertialToRswSatelliteCenteredFrameRotationMatrix(
        const Eigen::Vector6d& bodyState )
{
    Eigen::Vector3d vehicleVelocity, vehicleRadius;
    vehicleRadius = bodyState.segment( 0, 3 );
    vehicleVelocity = bodyState.segment( 3, 3 );

    Eigen::Vector3d unitR = vehicleRadius.normalized( );// / vehicleVelocity.norm( );
    if ( vehicleRadius.cross( vehicleVelocity ).norm( ) == 0.0 )
    {
        std::string errorMessage = "Division by zero: radius and velocity are in the same direction in RSW frame.";
        throw std::runtime_error( errorMessage );
    }

    Eigen::Vector3d unitW = ( vehicleRadius.cross( vehicleVelocity ) ).normalized( );

    Eigen::Vector3d unitS = ( unitW.cross( unitR ) ).normalized( );

    Eigen::Matrix3d transformationMatrix;
    transformationMatrix << unitR( 0 ), unitR( 1 ), unitR( 2 ),
            unitS( 0 ), unitS( 1 ), unitS( 2 ),
            unitW( 0 ), unitW( 1 ), unitW( 2 );
    return transformationMatrix;
}

Eigen::Matrix3d getRswSatelliteCenteredToInertialFrameRotationMatrix(
        const Eigen::Vector6d& bodyState )
{
    return getInertialToRswSatelliteCenteredFrameRotationMatrix( bodyState ).transpose( );
}

Eigen::Matrix3d getPqwPerifocalToInertialRotationMatrix(
        const Eigen::Vector6d& bodyState,
        const double gravitationalParameter )
{
    return getInertialToPqwPerifocalRotationMatrix( bodyState, gravitationalParameter ).transpose( );
}

Eigen::Matrix3d getInertialToPqwPerifocalRotationMatrix(
        const Eigen::Vector6d& bodyState,
        const double gravitationalParameter )
{
    Eigen::Matrix3d rotationInerialToRsw =
            getInertialToRswSatelliteCenteredFrameRotationMatrix( bodyState );
    double trueAnomaly = orbital_element_conversions::convertCartesianToKeplerianElements( bodyState, gravitationalParameter )( 5 );
    return Eigen::AngleAxisd( trueAnomaly,
                              Eigen::Vector3d::UnitZ( ) ).toRotationMatrix( ) * rotationInerialToRsw;
}


Eigen::Matrix3d getEqwEquinoctialToInertialRotationMatrix(
        const Eigen::Vector6d& bodyState,
        const double gravitationalParameter )
{
    return getInertialToEqwEquinoctialRotationMatrix( bodyState, gravitationalParameter ).transpose( );
}

Eigen::Matrix3d getInertialToEqwEquinoctialRotationMatrix(
        const Eigen::Vector6d& bodyState,
        const double gravitationalParameter )
{
    Eigen::Matrix3d rotationInerialToRsw =
            getInertialToRswSatelliteCenteredFrameRotationMatrix( bodyState );
    Eigen::Vector6d keplerElements =
            orbital_element_conversions::convertCartesianToKeplerianElements( bodyState, gravitationalParameter );
    Eigen::Matrix3d intermediateRotation =
            Eigen::AngleAxisd( -keplerElements( orbital_element_conversions::longitudeOfAscendingNodeIndex ),
                               Eigen::Vector3d::UnitZ( ) ).toRotationMatrix( );

    if( keplerElements( orbital_element_conversions::inclinationIndex ) > mathematical_constants::PI / 2.0 )
    {
        return intermediateRotation.transpose( ) *
                Eigen::AngleAxisd( -orbital_element_conversions::inclinationIndex, Eigen::Vector3d::UnitX( ) ).toRotationMatrix( ) *
                intermediateRotation *
                rotationInerialToRsw;
    }
    else
    {
        return intermediateRotation *
                Eigen::AngleAxisd( orbital_element_conversions::inclinationIndex, Eigen::Vector3d::UnitX( ) ).toRotationMatrix( ) *
                intermediateRotation.transpose( ) *
                rotationInerialToRsw;
    }
}
//! Get inertial (I) to rotating planetocentric (R) reference frame transformtion quaternion.
Eigen::Quaterniond getInertialToPlanetocentricFrameTransformationQuaternion(
        const double angleFromXItoXR )
{
    // Compute transformation quaternion.
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd eigenRotationObject = Eigen::AngleAxisd( -1.0 * angleFromXItoXR,
                                                               Eigen::Vector3d::UnitZ( ) );

    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond( eigenRotationObject );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}

//! Get rotation from inertial to planet-fixed frame.
Eigen::Quaterniond getInertialToPlanetocentricFrameTransformationQuaternion(
        const double declinationOfPole,
        const double rightAscensionOfPole,
        const double longitudeOfPrimeMeridian )
{
    // Compute transformation quaternion.
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd secondRotationAroundZaxis =
            Eigen::AngleAxisd( -longitudeOfPrimeMeridian, Eigen::Vector3d::UnitZ( ) );
    Eigen::AngleAxisd rotationAroundXaxis =
            Eigen::AngleAxisd(
                -( mathematical_constants::PI / 2.0 - declinationOfPole ),
                Eigen::Vector3d::UnitX( ) );
    Eigen::AngleAxisd firstRotationAroundZaxis = Eigen::AngleAxisd(
                - ( rightAscensionOfPole + mathematical_constants::PI / 2.0 ),
                Eigen::Vector3d::UnitZ( ) );
    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond(
                ( secondRotationAroundZaxis * rotationAroundXaxis *  firstRotationAroundZaxis ) );
    return frameTransformationQuaternion;
}

Eigen::Matrix3d getInertialToPlanetocentricFrameTransformationMatrix(
        const double declinationOfPole,
        const double rightAscensionOfPole,
        const double longitudeOfPrimeMeridian )
{
    return Eigen::Matrix3d(
                getInertialToPlanetocentricFrameTransformationQuaternion(
                    declinationOfPole, rightAscensionOfPole, longitudeOfPrimeMeridian ) );
}

//! Create a Quaterniond rotation state object from four quaternion values in a Vector4d
Eigen::Quaterniond getQuaternionObjectFromQuaternionValues(
        const Eigen::Vector4d& vectorWithQuaternion )
{
    // Set transformation quaternion.
    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond(
                vectorWithQuaternion( 0 ), vectorWithQuaternion( 1 ),
                vectorWithQuaternion( 2 ), vectorWithQuaternion( 3 ) );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}

//! Get transformation matrix from Planetocentric (R) to the Local vertical (V) frame.
Eigen::Matrix3d getRotatingPlanetocentricToLocalVerticalFrameTransformationMatrix(
        const double longitude, const double latitude )
{
    return getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion(
                longitude, latitude ).toRotationMatrix( );
}

//! Get transformation matrix from local vertical (V) to the Planetocentric frame (R).
Eigen::Matrix3d getLocalVerticalToRotatingPlanetocentricFrameTransformationMatrix(
        const double longitude, const double latitude )
{
    return getRotatingPlanetocentricToLocalVerticalFrameTransformationMatrix(
                longitude, latitude ).transpose( );
}

//! Get transformation matrix from the TA/TG to the V-frame.
Eigen::Matrix3d getTrajectoryToLocalVerticalFrameTransformationMatrix(
        const double flightPathAngle, const double headingAngle )
{
    return getTrajectoryToLocalVerticalFrameTransformationQuaternion(
                flightPathAngle, headingAngle ).toRotationMatrix( );
}

//! Get transformation quaternion from the TA/TG to the V-frame.
Eigen::Quaterniond getTrajectoryToLocalVerticalFrameTransformationQuaternion(
        const double flightPathAngle, const double headingAngle )
{
    // Compute transformation quaternion.
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd rotationAroundZaxis = Eigen::AngleAxisd(
                -1.0 * -headingAngle, Eigen::Vector3d::UnitZ( ) );
    Eigen::AngleAxisd rotationAroundYaxis = Eigen::AngleAxisd(
                -1.0 * -flightPathAngle, Eigen::Vector3d::UnitY( ) );

    return Eigen::Quaterniond( rotationAroundZaxis * rotationAroundYaxis );
}

//! Get transformation matrix from the local V- to TA/TG-frame.
Eigen::Matrix3d getLocalVerticalFrameToTrajectoryTransformationMatrix(
        const double flightPathAngle, const double headingAngle )
{
    return getTrajectoryToLocalVerticalFrameTransformationMatrix(
                flightPathAngle, headingAngle ).transpose( );
}

//! Get transformation quaternion from V- to the TA/TG-frame.
Eigen::Quaterniond getLocalVerticalFrameToTrajectoryTransformationQuaternion(
        const double flightPathAngle, const double headingAngle )
{
    return getTrajectoryToLocalVerticalFrameTransformationQuaternion(
                flightPathAngle, headingAngle ).inverse( );
}

//! Get transformation matrix from the TA- to the AA-frame.
Eigen::Matrix3d getTrajectoryToAerodynamicFrameTransformationMatrix(
        const double bankAngle )
{
    return getTrajectoryToAerodynamicFrameTransformationQuaternion(
                bankAngle ).toRotationMatrix( );
}

//! Get transformation quaternion from the TA- to the AA-frame.
Eigen::Quaterniond getTrajectoryToAerodynamicFrameTransformationQuaternion(
        const double bankAngle )
{
    // Compute transformation quaternion.
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd rotationAroundXaxis
            = Eigen::AngleAxisd( bankAngle, Eigen::Vector3d::UnitX( ) );
    return Eigen::Quaterniond( rotationAroundXaxis );
}

//! Get transformation matrix from the AA- to the TA-frame.
Eigen::Matrix3d getAerodynamicToTrajectoryFrameTransformationMatrix(
        const double bankAngle )
{
    return getTrajectoryToAerodynamicFrameTransformationMatrix( bankAngle ).transpose( );
}

//! Get transformation quaternion from the AA- to the TA-frame.
Eigen::Quaterniond getAerodynamicToTrajectoryFrameTransformationQuaternion(
        const double bankAngle )
{
    return getTrajectoryToAerodynamicFrameTransformationQuaternion( bankAngle ).inverse( );
}

//! Get transformation matrix fom the B- to the AA-frame.
Eigen::Matrix3d getBodyToAirspeedBasedAerodynamicFrameTransformationMatrix(
        const double angleOfAttack, const double angleOfSideslip )
{
    return getBodyToAirspeedBasedAerodynamicFrameTransformationQuaternion(
                angleOfAttack, angleOfSideslip ).toRotationMatrix( );
}

//! Get transformation quaternion fom the B- to the AA-frame.
Eigen::Quaterniond getBodyToAirspeedBasedAerodynamicFrameTransformationQuaternion(
        const double angleOfAttack, const double angleOfSideslip )
{
    // Compute transformation quaternion.
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd rotationAroundZaxis
            = Eigen::AngleAxisd( -1.0 * angleOfSideslip, Eigen::Vector3d::UnitZ( ) );
    Eigen::AngleAxisd rotationAroundYaxis
            = Eigen::AngleAxisd( -1.0 * -angleOfAttack, Eigen::Vector3d::UnitY( ) );

    return Eigen::Quaterniond( rotationAroundZaxis * rotationAroundYaxis );
}

//! Get transformation matrix fom the AA- to the B-frame.
Eigen::Matrix3d getAirspeedBasedAerodynamicToBodyFrameTransformationMatrix(
        const double angleOfAttack, const double angleOfSideslip )
{
    return getBodyToAirspeedBasedAerodynamicFrameTransformationMatrix(
                angleOfAttack, angleOfSideslip ).transpose( );
}

//! Get transformation quaternion fom the AA- to the B-frame.
Eigen::Quaterniond getAirspeedBasedAerodynamicToBodyFrameTransformationQuaternion(
        const double angleOfAttack, const double angleOfSideslip )
{
    return getBodyToAirspeedBasedAerodynamicFrameTransformationQuaternion(
                angleOfAttack, angleOfSideslip ).inverse( );
}

//! Calculate current heading angle.
double calculateHeadingAngle( const Eigen::Vector3d& velocityInVerticalFrame )
{
    return std::atan2( velocityInVerticalFrame( 1 ), velocityInVerticalFrame( 0 ) );
}

//! Calculate current flight path angle.
double calculateFlightPathAngle( const Eigen::Vector3d& velocityInVerticalFrame )
{
    return -std::asin( velocityInVerticalFrame( 2 ) / velocityInVerticalFrame.norm( ) );
}

//! Get transformation quaternion ECEF to ENU V-frame
Eigen::Quaterniond getRotatingPlanetocentricToEnuLocalVerticalFrameTransformationQuaternion(
        double longitude, double latitude )
{
    return getEnuLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion(
                longitude, latitude ).inverse( );
}

//! Get transformation quaternion between V-frame and ECEF
Eigen::Quaterniond getEnuLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion(
        double longitude, double latitude )
{
    // Compute transformation quaternion.
    // source: http://www.navipedia.net/index.php/Transformations_between_ECEF_and_ENU_coordinates
    // Note the sign change (-1.0), because how angleAxisd is defined.
    Eigen::AngleAxisd RotationAroundZaxis = Eigen::AngleAxisd(
                longitude + mathematical_constants::PI / 2.0, Eigen::Vector3d::UnitZ( ) );
    Eigen::AngleAxisd RotationAroundXaxis =
            Eigen::AngleAxisd( ( mathematical_constants::PI / 2.0 - latitude ),
                               Eigen::Vector3d::UnitX( ) );
    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond(
                ( RotationAroundZaxis * RotationAroundXaxis ) );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}

//! Get transformation matrix from J2000 to ECLIPJ2000
Eigen::Matrix3d getJ2000toECLIPJ2000TransformationMatrix ()
{
    return (Eigen::Matrix3d() <<
            1, 0, 0,
            0, 0.9174820620691818, 0.3977771559319137,
            0, -0.3977771559319137, 0.9174820620691818).finished() ;
}

//! Get transformation matrix from ECLIPJ2000 to J2000
Eigen::Matrix3d getECLIPJ2000toJ2000TransformationMatrix ()
    {
        return (Eigen::Matrix3d() <<
                1, 0, 0,
                0, 0.9174820620691818, -0.3977771559319137,
                0, 0.3977771559319137, 0.9174820620691818).finished() ;
    }

//! Function to compute the derivative of a rotation about the x-axis w.r.t. the rotation angle
Eigen::Matrix3d getDerivativeOfXAxisRotationWrtAngle( const double angle )
{
    return ( Eigen::Matrix3d( ) <<
             0.0, 0.0, 0.0,
             0.0, -std::sin( angle ), std::cos( angle ),
             0.0, -std::cos( angle ), -std::sin( angle ) ).finished( );
}

//! Function to compute the derivative of a rotation about the x-axis w.r.t. the rotation angle
Eigen::Matrix3d getDerivativeOfXAxisRotationWrtAngle( const Eigen::Matrix3d& rotationMatrix )
{
    return X_AXIS_ROTATION_MATRIX_DERIVATIVE_PREMULTIPLIER * rotationMatrix;
}

//! Function to compute the derivative of a rotation about the y-axis w.r.t. the rotation angle
Eigen::Matrix3d getDerivativeOfYAxisRotationWrtAngle( const double angle )
{
    return ( Eigen::Matrix3d( ) <<
             -std::sin( angle ), 0.0, -std::cos( angle ),
             0.0, 0.0, 0.0,
             -std::cos( angle ), 0.0, std::sin( angle ) ).finished( );
}

//! Function to compute the derivative of a rotation about the y-axis w.r.t. the rotation angle
Eigen::Matrix3d getDerivativeOfYAxisRotationWrtAngle( const Eigen::Matrix3d& rotationMatrix )
{
    return Y_AXIS_ROTATION_MATRIX_DERIVATIVE_PREMULTIPLIER * rotationMatrix;
}

//! Function to compute the derivative of a rotation about the z-axis w.r.t. the rotation angle
Eigen::Matrix3d getDerivativeOfZAxisRotationWrtAngle( const double angle )
{
    return ( Eigen::Matrix3d( ) <<
             -std::sin( angle ), std::cos( angle ), 0.0 ,
             -std::cos( angle ), -std::sin( angle ), 0.0,
             0.0, 0.0, 0.0 ).finished( );
}

//! Function to compute the derivative of a rotation about the z-axis w.r.t. the rotation angle
Eigen::Matrix3d getDerivativeOfZAxisRotationWrtAngle( const Eigen::Matrix3d& rotationMatrix )
{
    return Z_AXIS_ROTATION_MATRIX_DERIVATIVE_PREMULTIPLIER * rotationMatrix;
}

//! Function to compute a body-fixed relative cartesian position
Eigen::Vector3d getBodyFixedCartesianPosition(
        const std::function< Eigen::Vector3d( ) > positionFunctionOfCentralBody,
        const std::function< Eigen::Vector3d( ) > positionFunctionOfRelativeBody,
        const std::function< Eigen::Quaterniond( ) > orientationFunctionOfCentralBody )
{
    return orientationFunctionOfCentralBody( ) * (
                positionFunctionOfRelativeBody( ) - positionFunctionOfCentralBody( ) );
}

//! Function to compute a body-fixed relative spherical position
Eigen::Vector3d getBodyFixedSphericalPosition(
        const std::function< Eigen::Vector3d( ) > positionFunctionOfCentralBody,
        const std::function< Eigen::Vector3d( ) > positionFunctionOfRelativeBody,
        const std::function< Eigen::Quaterniond( ) > orientationFunctionOfCentralBody )
{
    Eigen::Vector3d sphericalPosition = coordinate_conversions::convertCartesianToSpherical(
                getBodyFixedCartesianPosition( positionFunctionOfCentralBody, positionFunctionOfRelativeBody,
                                               orientationFunctionOfCentralBody ) );
    sphericalPosition( 1 ) = mathematical_constants::PI / 2.0 - sphericalPosition( 1 );
    return sphericalPosition;
}

} // namespace reference_frames

} // namespace tudat
