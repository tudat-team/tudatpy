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
 *      Mooij, E. The Motion of a vehicle in a Planetary Atmosphere, TU Delft, 1997.
 *      Seidelmann, P. K. (Ed.). (2005). Explanatory supplement to the astronomical almanac.
 *              Univ Science Books.
 */

#ifndef TUDAT_REFERENCE_FRAME_TRANSFORMATIONS_H
#define TUDAT_REFERENCE_FRAME_TRANSFORMATIONS_H

#include <cmath>
#include <vector>

#include <functional>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/basic/basicMathematicsFunctions.h"
#include "tudat/basics/basicTypedefs.h"

namespace tudat
{

namespace reference_frames
{


//! Enum to define ids for various reference frames for calculating between inertial and body-fixed
//! frame, using transformation chain via aerodynamic frame.
enum AerodynamicsReferenceFrames
{
    inertial_frame = -1,
    corotating_frame = 0,
    vertical_frame = 1,
    trajectory_frame = 2,
    aerodynamic_frame = 3,
    body_frame = 4
};

// Enum defining identifiers of frames in which a user-specified thrust is defined.
//! @get_docstring(ThrustFrames.__docstring__)
enum SatelliteReferenceFrames
{
    unspecified_reference_frame = -1,
    global_reference_frame = 0,
    tnw_reference_frame = 1,
    rsw_reference_frame = 2
};

//! Function to get a string representing a 'named identification' of a reference frame.
/*!
 * Function to get a string representing a 'named identification' of a reference frame.
 * \param frame Type of reference frame
 * \return String with reference frame id.
 */
std::string getAerodynamicFrameName( const AerodynamicsReferenceFrames frame );


//! Function to compute pole right ascension and declination, as well as prime meridian of date, from rotation matrix
/*!
 *  Function to compute pole right ascension and declination, as well as prime meridian of date, from rotation matrix
 * \param rotationMatrixFromInertialToPlanetFixedFrame Rotation matrix from which Euler angles are to be determined.
 * \return Pole right ascension and declination, and prime meridian of date, from rotation matrix
 */
Eigen::Vector3d calculateInertialToPlanetFixedRotationAnglesFromMatrix(
        const Eigen::Matrix3d& rotationMatrixFromInertialToPlanetFixedFrame );

//! Wrapper function to transform a vector to a different frame from a single rotation function.
/*!
 * Wrapper function to transform a vector to a different frame from a single rotation function.
 * \param originalVector Vector that is to be rotated to a new frame
 * \param rotation Function returning the current rotation to the new frame
 * \return Vector originalVector, transformed to new frame.
 */
Eigen::Vector3d transformVectorFromQuaternionFunction(
        const Eigen::Vector3d& originalVector,
        const std::function< Eigen::Quaterniond( ) > rotation );

//! Wrapper function to transform a vector to a different frame from a single transformation function.
/*!
 * Wrapper function to transform a vector to a different frame from a single transformation function.
 * \param originalVector Vector that is to be transformed to a new frame
 * \param transformationFunction Function transforming a vector to a new frame
 * \return Vector originalVector, transformed to new frame.
 */
Eigen::Vector3d transformVectorFunctionFromVectorFunctions(
        const std::function< Eigen::Vector3d( ) > originalVector,
        const std::function< Eigen::Vector3d( const Eigen::Vector3d& ) > transformationFunction );

void transformVectorFunctionFromVectorReferenceFunctions(
        Eigen::Vector3d& transformedVector,
        const std::function< Eigen::Vector3d&( ) > originalVector,
        std::function< void( Eigen::Vector3d&, const Eigen::Vector3d& ) > transformationFunction );


//! Wrapper function to transform a vector to a different frame from a list of transformation function.
/*!
 * Wrapper function to transform a vector to a different frame from a list of transformation function.
 * \param originalVector Vector that is to be transformed to a new frame
 * \param rotationsList List of transformation function, each of which transforms a vector to a new frame. The functions
 * in this list are called in descending order.
 * \return Vector originalVector, transformed to new frame.
 */
Eigen::Vector3d transformVectorFromVectorFunctions(
        const Eigen::Vector3d& originalVector,
        const std::vector< std::function< Eigen::Vector3d( const Eigen::Vector3d& ) > >& rotationsList );

void transformVectorReferenceFromVectorFunctions(
        Eigen::Vector3d& transformedVector,
        const Eigen::Vector3d& originalVector,
        const std::vector< std::function< Eigen::Vector3d( const Eigen::Vector3d& ) > >& rotationsList );

//! Get rotating planetocentric (R) to inertial (I) reference frame transformation matrix.
/*!
 * Returns tranformation matrix from rotating planetocentric reference frame (R) to inertial
 * reference frame (I).
 * \param angleFromXItoXR Angle between X-axis of the planetocentric reference frame and
 *          X-axis of the inertial reference frame in [rad]. This angle is same as
 *          the rotational rate of the central body [rad/s] times the time from epoch [s].
 * \return Reference frame (R) to inertial reference frame (I) transformation matrix.
 */
Eigen::Matrix3d getRotatingPlanetocentricToInertialFrameTransformationMatrix(
        const double angleFromXItoXR );

//! Get rotating planetocentric (R) to inertial (I) reference frame transformation quaternion.
/*!
 * Returns tranformation quaternion from rotating planetocentric reference frame (R) to inertial
 * reference frame (I). It's an eigen library transformation and can be applied directly to a
 * vector Vector_new = Quaternion * Vector_old.
 * \param angleFromXItoXR Angle between X-axis of the planetocentric reference frame and
 *          X-axis of the inertial reference frame in [rad]. This angle is same as
 *          the rotational rate of the central body [rad/s] times the time from epoch [s].
 * \return Reference frame (R) to inertial reference frame (I) transformation quaternion.
 */
Eigen::Quaterniond getRotatingPlanetocentricToInertialFrameTransformationQuaternion(
        const double angleFromXItoXR );

//! Get rotation from planet-fixed to inertial frame.
/*!
 * Retuns rotation from planet-fixed to inertial frame, assuming that the equatorial plane is not
 * equal to x-y plane of inertial frame. Orientation of body-fixed frame is obtained from right
 * ascension and declination of body's pole and the location of the prime meridian
 * (Seidelmann et al. 2005).
 * \param declinationOfPole Declination of body's pole in inertial frame.
 * \param rightAscensionOfPole Right ascension of body's pole in inertial frame.
 * \param longitudeOfPrimeMeridian Longitude of body prime meridian.
 * \return Rotation quaternion computed.
 */
Eigen::Quaterniond getRotatingPlanetocentricToInertialFrameTransformationQuaternion(
        const double declinationOfPole,
        const double rightAscensionOfPole,
        const double longitudeOfPrimeMeridian );

Eigen::Matrix3d getRotatingPlanetocentricToInertialFrameTransformationMatrix(
        const double declinationOfPole,
        const double rightAscensionOfPole,
        const double longitudeOfPrimeMeridian );

//! Get inertial (I) to rotating planetocentric (R) reference frame transformtion matrix.
/*!
 * Returns transformation matrix from inertial referenceframe (I) to the rotating planetocentric
 * reference frame (R).
 * \param angleFromXItoXR Angle between X-axis of the planetocentric reference frame and
 *          X-axis of the inertial reference frame in [rad]. This angle is same as
 *          the rotational rate of the central body [rad/s] times the time from epoch [s].
 * \return Inertial (I) to planetocentric reference frame (R) transformation matrix.
 */
Eigen::Matrix3d getInertialToPlanetocentricFrameTransformationMatrix(
        const double angleFromXItoXR );

//! Get rotation from velocity based TNW frame to inertial (I) frame.
/*!
 * Returns rotation from inertial (i) to the velocity based TNW frame. The velocity based TNW frame is
 * a right-handed orthogonal frame defined as follows:
 * x-axis tangent to the velocity direction,
 * y-axis in the orbital plane and pointing inwards (if doesNaxisPointAwayFromCentralBody is false),
 * i.e. to the left when looking in velocity-direction,
 * z-axis normal to the orbital plane.
 * \param vehicleState State of the vehicle for which the TNW frame is to be computed.
 * \param centralBodyState State of the central body w.r.t. which the TNW frame is to be computed.
 * \param doesNaxisPointAwayFromCentralBody Boolean denoting whether the local y-axis points away from (if true) or
 * towards (if false) central body.
 * \return Velocity based TNW to inertial (I) frame transformation matrix.
 */
Eigen::Matrix3d getTnwToInertialRotation(const Eigen::Vector6d& vehicleState,
                                                       const Eigen::Vector6d& centralBodyState,
                                                       const bool doesNaxisPointAwayFromCentralBody = true );

Eigen::Matrix3d getTnwToInertialRotation(const Eigen::Vector6d& vehicleInertialState,
                                         const bool doesNaxisPointAwayFromCentralBody = true );

Eigen::Matrix3d getInertialToTnwRotation(const Eigen::Vector6d& vehicleInertialState,
                                         const bool doesNaxisPointAwayFromCentralBody = true );

//! Get rotation from velocity based TNW frame to inertial (I) frame.
/*!
 * Returns rotation from inertial (i) to the velocity based TNW frame. The velocity based TNW frame is
 * a right-handed orthogonal frame defined as follows:
 * x-axis tangent to the velocity direction,
 * y-axis in the orbital plane and pointing inwards (if doesNaxisPointAwayFromCentralBody is false),
 * i.e. to the left when looking in velocity-direction,
 * z-axis normal to the orbital plane.
 * \param vehicleStateFunction Function returning the state of the vehicle for which the TNW frame is to be computed
 * \param centralBodyStateFunction Function returning the state of the central body w.r.t. which the TNW frame is to be
 * computed
 * \param doesNaxisPointAwayFromCentralBody Boolean denoting whether the local y-axis points away from (if true) or
 * towards (if false) central body.
 * \return Velocity based TNW to inertial (I) frame transformation matrix.
 */
Eigen::Matrix3d getTnwToInertialRotationFromFunctions(
        const std::function< Eigen::Vector6d( ) >& vehicleStateFunction,
        const std::function< Eigen::Vector6d( ) >& centralBodyStateFunction,
        bool doesNaxisPointAwayFromCentralBody = true );

//! Get rotation from velocity based TNW frame to planet-fixed frame.
/*!
 * Returns rotation from the velocity based TNW frame to the planet-fixed frame. The velocity based TNW frame is
 * a right-handed orthogonal frame defined as follows:
 * x-axis tangent to the velocity direction,
 * y-axis in the orbital plane and pointing inwards, i.e. to the left when looking in velocity-direction,
 * z-axis normal to the orbital plane.
 * \param spacecraftKeplerianState containging the following elements:
 *          semi-major axis -> not used
 *          eccentricity
 *          inclination
 *          argumentOfPeriapsis
 *          longitudeOfAscendingNode
 *          trueAnomaly
 * \return Computed rotation quaternion.
 */
//! Get rotation from velocity based TNW frame to planetocentric frame.
Eigen::Quaterniond getTnwToPlanetocentricRotationKeplerian(
        const Eigen::Matrix< double, 6, 1 > spacecraftKeplerianState );

//! Function to compute the rotation matrix to RSW frame, from the frame in which the input state is given.
/*!
 * Function to compute the rotation matrix to RSW frame, from the frame in which the input state is given.
 * \param bodyState State for which the RSW frame rotation is to be computed.
 * \return Rotation matrix to RSW frame
 */
Eigen::Matrix3d getInertialToRswSatelliteCenteredFrameRotationMatrix(
        const Eigen::Vector6d& bodyState );

Eigen::Matrix3d getRswSatelliteCenteredToInertialFrameRotationMatrix(
        const Eigen::Vector6d& bodyState );

Eigen::Matrix3d getPqwPerifocalToInertialRotationMatrix(
        const Eigen::Vector6d& bodyState,
        const double gravitationalParameter );

Eigen::Matrix3d getInertialToPqwPerifocalRotationMatrix(
        const Eigen::Vector6d& bodyState,
        const double gravitationalParameter  );

Eigen::Matrix3d getEqwEquinoctialToInertialRotationMatrix(
        const Eigen::Vector6d& bodyState,
        const double gravitationalParameter );

Eigen::Matrix3d getInertialToEqwEquinoctialRotationMatrix(
        const Eigen::Vector6d& bodyState,
        const double gravitationalParameter );

//! Get inertial (I) to rotating planetocentric (R) reference frame transformation quaternion.
/*!
 * Returns transformation quaternion from inertial referenceframe (I) to the rotating
 * planetocentric reference frame (R). It's an eigen library transformation and can be applied
 * directly to a vector Vector_new = Quaternion * Vector_old.
 * \param angleFromXItoXR Angle between X-axis of the planetocentric reference frame and
 *          X-axis of the inertial reference frame in [rad]. This angle is same as
 *          the rotational rate of the central body [rad/s] times the time from epoch [s].
 * \return Inertial (I) to planetocentric reference frame (R) transformation quaternion.
 */
Eigen::Quaterniond getInertialToPlanetocentricFrameTransformationQuaternion(
        const double angleFromXItoXR );

//! Get rotation from inertial to planet-fixed frame.
/*!
 * Returns rotation from inertial to planet-fixed frame, assuming that the equatorial plane is not
 * equal to x-y plane of inertial frame. Orientation of body-fixed frame is obtained from right
 * ascension and declination of body's pole and the location of the prime meridian
 * (Seidelmann et al. 2005).
 * \param declinationOfPole Declination of body's pole in inertial frame.
 * \param rightAscensionOfPole Right ascension of body's pole in inertial frame.
 * \param longitudeOfPrimeMeridian Longitude of body prime meridian.
 * \return Rotation quaternion computed.
 */
Eigen::Quaterniond getInertialToPlanetocentricFrameTransformationQuaternion(
        const double declinationOfPole,
        const double rightAscensionOfPole,
        const double longitudeOfPrimeMeridian );

Eigen::Matrix3d getInertialToPlanetocentricFrameTransformationMatrix(
        const double declinationOfPole,
        const double rightAscensionOfPole,
        const double longitudeOfPrimeMeridian );

//! Create a Quaterniond rotation state object from four quaternion values in vector 4d.
/*!
 * Creates a Quaterniond rotation state object from four quaternion values.
 * This function is not related to any specific rotation matrix, but can be used for general
 * purposes. It's an eigen library transformation and can be applied
 * directly to a vector Vector_new = Quaternion * Vector_old.
 * Note that is also possible to create a quaternion object directly from a Vector4d,
 * but Eigen will rearrange the order of the coefficients ([q2 q3 q4 q1]). This function
 * retreives the individual entries of the Vector4d an uses them as four doubles as input
 * arguments for the constructor of the Quateriond. From the Eigen code documentation:
 * \warning Note the order of the arguments: the real \a w coefficient first,
 *              while internally the coefficients are stored in the following order:
 *              [\c x, \c y, \c z, \c w]
 * \param vectorWithQuaternion A vector containing the quaternions of the rotation state
 * \return Transformation quaternion.
 */
Eigen::Quaterniond getQuaternionObjectFromQuaternionValues(
        const Eigen::Vector4d& vectorWithQuaternion );

//! Get transformation matrix from Planetocentric (R) to the Local vertical (V) frame.
/*!
 * Returns the frame transformation matrix from the Planetocentric (R) to the Local vertical
 * (V) reference frame. The Z-axis is aligned with the local gravity vector. Whether or not,
 * this is in the direction of the center of the central body, depends which kind of latitude
 * is provided (geocentric, geodetic, gravitation latitude).
 * The X-axis is directed to the north.
 * \param longitude The longitude in the planetocentric reference frame in [rad].
 * \param latitude The latitude in the planetocentric reference frame in [rad].
 * \return Transformation matrix from Planetocentric (R) to the local vertical (V) frame.
 */
Eigen::Matrix3d getRotatingPlanetocentricToLocalVerticalFrameTransformationMatrix(
        const double longitude, const double latitude );

//! Get transformation quaternion from Planetocentric (R) to the local vertical (V) frame.
/*!
 * Returns the frame transformation quaternion from the Planetocentric (R) to the Local vertical
 * (V) reference frame. The Z-axis is aligned with the local gravity vector. Whether or not,
 * this is in the direction of the center of the central body, depends which kind of latitude
 * is provided (geocentric, geodetic, gravitation latitude).
 * The X-axis is directed to the north.
 * \param longitude The longitude in the planetocentric reference frame in [rad].
 * \param latitude The latitude in the planetocentric reference frame in [rad].
 * \return Transformation quaternion from Planetocentric (R) to the local vertical (V) frame.
 */
template< typename AngleScalarType = double >
Eigen::Quaternion< AngleScalarType > getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion(
        const AngleScalarType longitude, const AngleScalarType latitude )
{
    // Compute transformation quaternion.
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxis< AngleScalarType > RotationAroundZaxis = Eigen::AngleAxis< AngleScalarType >(
                -longitude, Eigen::Matrix< AngleScalarType, 3, 1 >::UnitZ( ) );
    Eigen::AngleAxis< AngleScalarType > RotationAroundYaxis = Eigen::AngleAxis< AngleScalarType >(
                -( -latitude - mathematical_constants::getPi< AngleScalarType >( ) /
                   mathematical_constants::getFloatingInteger< AngleScalarType >( 2 ) ),
                Eigen::Matrix< AngleScalarType, 3, 1 >::UnitY( ) );
    Eigen::Quaternion< AngleScalarType > frameTransformationQuaternion = Eigen::Quaternion< AngleScalarType >(
                ( RotationAroundYaxis * RotationAroundZaxis ) );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}

//! Get transformation matrix from local vertical (V) to the Planetocentric frame (R).
/*!
 * Returns the frame transformation matrix from the local vertical (V) to the
 * Planetocentric frame (R) reference frame. The Z-axis is aligned with the local gravity vector.
 * Whether or not, this is in the direction of the center of the central body, depends which kind
 * of latitude is provided (geocentric, geodetic, gravitation latitude).
 * The X-axis is directed to the north.
 * \param longitude The longitude in the planetocentric reference frame in [rad].
 * \param latitude The latitude in the planetocentric reference frame in [rad].
 * \return Transformation matrix from local vertical (V) to the Planetocentric (R) frame.
 */
Eigen::Matrix3d getLocalVerticalToRotatingPlanetocentricFrameTransformationMatrix(
        const double longitude, const double latitude );

//! Get transformation quaternion from local vertical (V) to the Planetocentric frame (R).
/*!
 * Returns the frame transformation quaternion from the local vertical (V) to the
 * Planetocentric (R) reference frame. The Z-axis is aligned with the local gravity vector.
 * Whether or not, this is in the direction of the center of the central body, depends which kind
 * of latitude is provided (geocentric, geodetic, gravitation latitude).
 * The X-axis is directed to the north.
 * \param longitude The longitude in the planetocentric reference frame in [rad].
 * \param latitude The latitude in the planetocentric reference frame in [rad].
 * \return Transformation quaternion from local vertical (V) to the Planetocentric (R) frame.
 */
template< typename AngleScalarType = double >
Eigen::Quaternion< AngleScalarType > getLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion(
        const AngleScalarType longitude, const AngleScalarType latitude )
{
    return getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion(
            longitude, latitude ).inverse( );
}


//! Get transformation matrix from the TA/TG to the V-frame.
/*!
 * Returns the frame transformation matrix from the trajectory frame (T) to the local vertical (V)
 * reference frame. Depending on whether the flight-path angle and heading angle express the
 * velocity relative to the rotating planetocentric frame or relative to the local atmosphere,
 * the resulting rotation holds for the groundspeed-based trajectory (TG) or airspeed-based
 * trajectory (TA) frame respectively.
 * \param flightPathAngle The trajectory's flight-path angle, positive upwards, in [rad].
 * \param headingAngle The trajectory's heading angle with respect to the North in [rad].
 * \return Transformation matrix from trajectory frame (T) to the local vertical frame (V).
 */
Eigen::Matrix3d getTrajectoryToLocalVerticalFrameTransformationMatrix(
        const double flightPathAngle, const double headingAngle );

//! Get transformation quaternion from the TA/TG to the V-frame.
/*!
 * Returns the frame transformation quaternion from the trajectory (T) to the local vertical (V)
 * reference frame. Depending on whether the flight-path angle and heading angle express the
 * velocity relative to the rotating planetocentric frame or relative to the local atmosphere,
 * the resulting rotation holds for the groundspeed-based trajectory (TG) or airspeed-based
 * trajectory (TA) frame respectively.
 * \param flightPathAngle The trajectory's flight-path angle, positive upwards, in [rad].
 * \param headingAngle The trajectory's heading angle with respect to the North in [rad].
 * \return Transformation quaternion from trajectory frame (T) to the local vertical frame (V).
 */
Eigen::Quaterniond getTrajectoryToLocalVerticalFrameTransformationQuaternion(
        const double flightPathAngle, const double headingAngle );

//! Get transformation matrix from the local V- to TA/TG-frame.
/*!
 * Returns the frame transformation matrix from the local vertical (V) to the trajectory (T)
 * reference frame. Depending on whether the flight-path angle and heading angle express the
 * velocity relative to the rotating planetocentric frame or relative to the local atmosphere,
 * the resulting rotation holds for the groundspeed-based trajectory (TG) or airspeed-based
 * trajectory (TA) frame respectively.
 * \param flightPathAngle The trajectory's flight-path angle, positive upwards, in [rad].
 * \param headingAngle The trajectory's heading angle with respect to the North in [rad].
 * \return Transformation matrix from the local vertical (V) to the trajectory (T) frame.
 */
Eigen::Matrix3d getLocalVerticalFrameToTrajectoryTransformationMatrix(
        const double flightPathAngle, const double headingAngle );

//! Get transformation quaternion from V- to the TA/TG-frame.
/*!
 * Returns the transformation quaternion from the local vertical (V) to the trajectory (T)
 * reference frame. Depending on whether the flight-path angle and heading angle express the
 * velocity relative to the rotating planetocentric frame or relative to the local atmosphere,
 * the resulting rotation holds for the groundspeed-based trajectory (TG) or airspeed-based
 * trajectory (TA) frame respectively.
 * \param flightPathAngle The trajectory's flight-path angle, positive upwards, in [rad].
 * \param headingAngle The trajectory's heading angle with respect to the North in [rad].
 * \return Transformation quaternion from the local vertical (V) to the trajectory (T) frame.
 */
Eigen::Quaterniond getLocalVerticalFrameToTrajectoryTransformationQuaternion(
        const double flightPathAngle, const double headingAngle );

//! Get transformation matrix from the TA- to the AA-frame.
/*!
 * Returns the transformation matrix from the airspeed-based trajectory (TA) to the airspeed-based
 * aerodynamic frame (AA). These frames differ from each other only by the bank-angle, representing
 * one Euler-rotation.
 * \param bankAngle The object's bank angle in [rad].
 * \return Transformation matrix from the TA- to the AA-frame.
 */
Eigen::Matrix3d getTrajectoryToAerodynamicFrameTransformationMatrix( const double bankAngle );

//! Get transformation quaternion from the TA- to the AA-frame.
/*!
 * Returns the transformation quaternion from the airspeed-based trajectory (TA) to the airspeed-
 * based aerodynamic frame (AA). These frames differ from each other only by the bank-angle,
 * representing one Euler-rotation.
 * \param bankAngle The object's bank angle in [rad].
 * \return Transformation quaternion from the TA- to the AA-frame.
 */
Eigen::Quaterniond getTrajectoryToAerodynamicFrameTransformationQuaternion( 
        const double bankAngle );

//! Get transformation matrix from the AA- to the TA-frame.
/*!
 * Returns the transformation matrix from the airspeed-based aerodynamic (AA) to the airspeed-based
 * trajectory frame (TA). These frames differ from each other only by the bank-angle, representing
 * one Euler-rotation.
 * \param bankAngle The object's bank angle in [rad].
 * \return Transformation matrix from the AA- to the TA-frame.
 */
Eigen::Matrix3d getAerodynamicToTrajectoryFrameTransformationMatrix( const double bankAngle );

//! Get transformation quaternion from the AA- to the TA-frame.
/*!
 * Returns the transformation quaternion from the airspeed-based aerodynamic (AA) to the airspeed-
 * based trajectory frame (TA). These frames differ from each other only by the bank-angle,
 * representing one Euler-rotation.
 * \param bankAngle The object's bank angle in [rad].
 * \return Transformation quaternion from the AA- to the TA-frame.
 */
Eigen::Quaterniond getAerodynamicToTrajectoryFrameTransformationQuaternion(
        const double bankAngle );

//! Get transformation matrix fom the B- to the AA-frame.
/*!
 * Returns the transformation matrix from the body-fixed (B) to the airspeed-based aerodynamic
 * frame (AA).
 * \param angleOfAttack The angle of attack in [rad].
 * \param angleOfSideslip The angle of sideslip in [rad].
 * \return Transformation matrix from the B- to the AA-frame.
 */
Eigen::Matrix3d getBodyToAirspeedBasedAerodynamicFrameTransformationMatrix(
        const double angleOfAttack, const double angleOfSideslip );

//! Get transformation quaternion fom the B- to the AA-frame.
/*!
 * Returns the transformation quaternion from the body-fixed (B) to the airspeed-based aerodynamic
 * frame (AA).
 * \param angleOfAttack The angle of attack in [rad].
 * \param angleOfSideslip The angle of sideslip in [rad].
 * \return Transformation quaternion from the B- to the AA-frame.
 */
Eigen::Quaterniond getBodyToAirspeedBasedAerodynamicFrameTransformationQuaternion(
        const double angleOfAttack, const double angleOfSideslip );

//! Get transformation matrix fom the AA- to the B-frame.
/*!
 * Returns the transformation matrix from the airspeed-based aerodynamic (AA) to the body-fixed
 * frame (B).
 * \param angleOfAttack The angle of attack in [rad].
 * \param angleOfSideslip The angle of sideslip in [rad].
 * \return Transformation matrix from the AA- to the B-frame.
 */
Eigen::Matrix3d getAirspeedBasedAerodynamicToBodyFrameTransformationMatrix(
        const double angleOfAttack, const double angleOfSideslip );

//! Get transformation quaternion fom the AA- to the B-frame.
/*!
 * Returns the transformation quaternion from the airspeed-based aerodynamic (AA) to the body-fixed
 * frame (B).
 * \param angleOfAttack The angle of attack in [rad].
 * \param angleOfSideslip The angle of sideslip in [rad].
 * \return Transformation quaternion from the AA- to the B-frame.
 */
Eigen::Quaterniond getAirspeedBasedAerodynamicToBodyFrameTransformationQuaternion(
        const double angleOfAttack, const double angleOfSideslip );

//! Calculate current heading angle.
/*!
 * Calculate heading angle from velocity in vertical (TNW) frame.
 * \param velocityInVerticalFrame Current Cartesian velocity in vertical frame.
 * \return Current heading angle.
 */
double calculateHeadingAngle( const Eigen::Vector3d& velocityInVerticalFrame );

//! Calculate current flight path angle. Angle is defined positive upwards.
/*!
 *  Calculate flight path angle from velocity in vertical (TNW) frame.
 *  Angle is defined positive upwards.
 *  \param velocityInVerticalFrame Current Cartesian velocity in vertical frame.
 *  \return Current flight path angle.
 */
double calculateFlightPathAngle( const Eigen::Vector3d& velocityInVerticalFrame );

//! Get ECEF to V-frame quaternion
/*!
 *  Get the transformation from the co-rotating planetocentric frame to local vertical.
 *  \param longitude Longitude of position.
 *  \param latitude Latitude of position.
 *  \return Transformation quaternion.
 */
Eigen::Quaterniond getRotatingPlanetocentricToEnuLocalVerticalFrameTransformationQuaternion(
        const double longitude, const double latitude );

//! Get V-frame to ECEF quaternion
/*!
 *  Get the transformation from the local vertical to the co-rotating planetocentric frame.
 *  \sa http://www.navipedia.net/index.php/Transformations_between_ECEF_and_ENU_coordinates
 *  \param longitude Longitude of position.
 *  \param latitude Latitude of position.
 *  \return Transformation quaternion.
 */
Eigen::Quaterniond getEnuLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion(
        const double longitude, const double latitude );

//! Get transformation matrix from J2000 to ECLIPJ2000. Transformation matrix with values from SPICE
/*!
 * @return Transformation matrix
 */
Eigen::Matrix3d getJ2000toECLIPJ2000TransformationMatrix ();

//! Get transformation matrix from ECLIPJ2000 to J2000. Transformation matrix with values from SPICE
/*!
 * @return Transformation matrix
 */
Eigen::Matrix3d getECLIPJ2000toJ2000TransformationMatrix ();

//! Pre-multiplier used to take derivative of rotation matrix about x-axis w.r.t. the rotation angle
static const Eigen::Matrix3d X_AXIS_ROTATION_MATRIX_DERIVATIVE_PREMULTIPLIER =
        ( Eigen::Matrix3d( ) <<
          0.0, 0.0, 0.0,
          0.0, 0.0, 1.0,
          0.0, -1.0, 0.0 ).finished( );

//! Pre-multiplier used to take derivative of rotation matrix about y-axis w.r.t. the rotation angle
static const Eigen::Matrix3d Y_AXIS_ROTATION_MATRIX_DERIVATIVE_PREMULTIPLIER =
        ( Eigen::Matrix3d( ) <<
          0.0, 0.0, -1.0,
          0.0, 0.0, 0.0,
          1.0, 0.0, 0.0 ).finished( );

//! Pre-multiplier used to take derivative of rotation matrix about z-axis w.r.t. the rotation angle
static const Eigen::Matrix3d Z_AXIS_ROTATION_MATRIX_DERIVATIVE_PREMULTIPLIER =
        ( Eigen::Matrix3d( ) <<
          0.0, 1.0, 0.0,
          -1.0, 0.0, 0.0,
          0.0, 0.0, 0.0 ).finished( );

//! Function to compute the derivative of a rotation about the x-axis w.r.t. the rotation angle
/*!
 * Function to compute the derivative of a rotation about the x-axis w.r.t. the rotation angle
 * \param angle Angle about which rotation is taken.
 * \return Derivative of a rotation about the x-axis w.r.t. the rotation angle
 */
Eigen::Matrix3d getDerivativeOfXAxisRotationWrtAngle( const double angle );

//! Function to compute the derivative of a rotation about the x-axis w.r.t. the rotation angle
/*!
 * Function to compute the derivative of a rotation about the x-axis w.r.t. the rotation angle
 * \param rotationMatrix Rotation matrix for which partial is to be computed
 * \return Derivative of a rotation about the x-axis w.r.t. the rotation angle
 */
Eigen::Matrix3d getDerivativeOfXAxisRotationWrtAngle( const Eigen::Matrix3d& rotationMatrix );

//! Function to compute the derivative of a rotation about the y-axis w.r.t. the rotation angle
/*!
 * Function to compute the derivative of a rotation about the y-axis w.r.t. the rotation angle
 * \param angle Angle about which rotation is taken.
 * \return Derivative of a rotation about the x-axis w.r.t. the rotation angle
 */
Eigen::Matrix3d getDerivativeOfYAxisRotationWrtAngle( const double angle );

//! Function to compute the derivative of a rotation about the y-axis w.r.t. the rotation angle
/*!
 * Function to compute the derivative of a rotation about the y-axis w.r.t. the rotation angle
 * \param rotationMatrix Rotation matrix for which partial is to be computed
 * \return Derivative of a rotation about the x-axis w.r.t. the rotation angle
 */
Eigen::Matrix3d getDerivativeOfYAxisRotationWrtAngle( const Eigen::Matrix3d& rotationMatrix );

//! Function to compute the derivative of a rotation about the z-axis w.r.t. the rotation angle
/*!
 * Function to compute the derivative of a rotation about the z-axis w.r.t. the rotation angle
 * \param angle Angle about which rotation is taken.
 * \return Derivative of a rotation about the x-axis w.r.t. the rotation angle
 */
Eigen::Matrix3d getDerivativeOfZAxisRotationWrtAngle( const double angle );

//! Function to compute the derivative of a rotation about the z-axis w.r.t. the rotation angle
/*!
 * Function to compute the derivative of a rotation about the z-axis w.r.t. the rotation angle
 * \param rotationMatrix Rotation matrix for which partial is to be computed
 * \return Derivative of a rotation about the x-axis w.r.t. the rotation angle
 */
Eigen::Matrix3d getDerivativeOfZAxisRotationWrtAngle( const Eigen::Matrix3d& rotationMatrix );

//! Function to compute a body-fixed relative cartesian position
/*!
 * Function to compute a body-fixed relative cartesian position
 * \param positionFunctionOfCentralBody Position function of central body
 * \param positionFunctionOfRelativeBody Position function of point of which body-fixed state is to be computed
 * \param orientationFunctionOfCentralBody Function returning rotation from inertial to body-fixed frame.
 * \return Body-fixed relative cartesian position
 */
Eigen::Vector3d getBodyFixedCartesianPosition(
        const std::function< Eigen::Vector3d( ) > positionFunctionOfCentralBody,
        const std::function< Eigen::Vector3d( ) > positionFunctionOfRelativeBody,
        const std::function< Eigen::Quaterniond( ) > orientationFunctionOfCentralBody );

//! Function to compute a body-fixed relative spherical position
/*!
 * Function to compute a body-fixed relative sphericall position
 * \param positionFunctionOfCentralBody Position function of central body
 * \param positionFunctionOfRelativeBody Position function of point of which body-fixed state is to be computed
 * \param orientationFunctionOfCentralBody Function returning rotation from inertial to body-fixed frame.
 * \return Body-fixed relative sphericall position
 */
Eigen::Vector3d getBodyFixedSphericalPosition(
        const std::function< Eigen::Vector3d( ) > positionFunctionOfCentralBody,
        const std::function< Eigen::Vector3d( ) > positionFunctionOfRelativeBody,
        const std::function< Eigen::Quaterniond( ) > orientationFunctionOfCentralBody );

Eigen::Matrix3d getRotationBetweenSatelliteFrames(
        const Eigen::Vector6d relativeInertialCartesianState,
        const SatelliteReferenceFrames originalFrame,
        const SatelliteReferenceFrames targetFrame );

/*! Compute the rotation matrix between ITRF2014 and another ITRF frame.
 *
 * Compute the rotation matrix between ITRF2014 and another ITRF frame. According to the IERS convetions of 2010, Eq.
 * 4.2, Eq. 4.4, Tab. 4.1.
 *
 * @param targetFrame Desired frame
 * @return Rotation matrix
 */
Eigen::Matrix3d getItrf2014ToArbitraryItrfRotationMatrix( const std::string& targetFrame );

/*! Compute the derivative of the rotation matrix between ITRF2014 and another ITRF frame.
 *
 * Compute the rotation matrix derivative between ITRF2014 and another ITRF frame. According to the IERS convetions of ´
 * 2010, Eq. 4.2, Eq. 4.4, Tab. 4.1.
 *
 * @param targetFrame Desired frame
 * @return Rotation matrix derivative
 */
Eigen::Matrix3d getItrf2014ToArbitraryItrfRotationMatrixDerivative( const std::string& targetFrame );

/*! Compute the translation vector between ITRF2014 and another ITRF frame.
 *
 * Compute the translation vector between ITRF2014 and another ITRF frame. According to the IERS convetions of ´
 * 2010, Eq. 4.2, Eq. 4.4, Tab. 4.1.
 *
 * @param targetFrame Desired frame
 * @return Translation vector
 */
Eigen::Vector6d getItrf2014ToArbitraryItrfTranslation( const std::string& targetFrame );

/*! Function to convert the state of a ground station from ITRF2014 to another ITRS realization.
 *
 * Function to convert the state of a ground station from ITRF2014 to another ITRS realization. Conversion is executed
 * according to the IERS convetions of 2010, Eq. 4.2, Eq. 4.4, Tab. 4.1. The ground station is assumed to have a constant
 * velocity.
 *
 * @param groundStationStateAtEpoch State of the ground station.
 * @param epoch Epoch for which the state of the ground station is valid.
 * @param targetFrame Desired frame
 * @return Ground station state at specified epoch with respect to targetFrame.
 */
Eigen::Vector6d convertGroundStationStateItrf2014ToArbitraryItrf(
        const Eigen::Vector6d& groundStationStateAtEpoch,
        double epoch,
        const std::string& targetFrame );

/*! Function to convert the state of a ground station from some ITRS realization to ITRF2014.
 *
 * Function to convert the state of a ground station from some ITRS realization to ITRF2014. Conversion is executed
 * according to the IERS convetions of 2010, Eq. 4.2, Eq. 4.4, Tab. 4.1. The ground station is assumed to have a constant
 * velocity.
 *
 * @param groundStationStateAtEpoch State of the ground station.
 * @param epoch Epoch for which the state of the ground station is valid.
 * @param baseFrame Frame with resoect to which the specified state is defined.
 * @return Ground station state at specified epoch with respect to ITRF2014.
 */
Eigen::Vector6d convertGroundStationStateArbitraryItrfToItrf2014(
        const Eigen::Vector6d& groundStationStateAtEpoch,
        double epoch,
        const std::string& baseFrame );

/*! Function to convert between two ITRF frames.
 *
 * Function to convert between two ITRF frames. Conversion is executed
 * according to the IERS convetions of 2010, Eq. 4.2, Eq. 4.4, Tab. 4.1. The ground station is assumed to have a constant
 * velocity.
 *
 * @param groundStationStateAtEpoch State of the ground station.
 * @param epoch Epoch for which the state of the ground station is valid.
 * @param baseFrame Frame with resoect to which the input state is defined.
 * @param targetFrame Frame with resoect to which the output state is defined.
 * @return Ground station state at specified epoch with respect to targetFrame.
 */
Eigen::Vector6d convertGroundStationStateBetweenItrfFrames(
        const Eigen::Vector6d& groundStationStateAtEpoch,
        double epoch,
        const std::string& baseFrame,
        const std::string& targetFrame );


} // namespace reference_frames

} // namespace tudat

#endif // TUDAT_REFERENCE_FRAME_TRANSFORMATIONS_H
