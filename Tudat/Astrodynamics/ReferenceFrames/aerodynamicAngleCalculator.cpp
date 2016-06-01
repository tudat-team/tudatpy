/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include "Tudat/Astrodynamics/ReferenceFrames/aerodynamicAngleCalculator.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{

namespace reference_frames
{

//! Function to update the orientation angles to the current state.
void AerodynamicAngleCalculator::update( )
{
    // Clear all current rotation matrices.
    currentRotationMatrices_.clear( );

    // Get current body-fixed state.
    currentBodyFixedState_ = bodyFixedStateFunction_( );
    Eigen::Vector3d sphericalCoordinates = coordinate_conversions::convertCartesianToSpherical(
                currentBodyFixedState_.segment( 0, 3 ) );

    // Calculate latitude and longitude.
    currentAerodynamicAngles_[ latitude_angle ] =
            mathematical_constants::PI / 2.0 - sphericalCoordinates( 1 );
    currentAerodynamicAngles_[ longitude_angle ] = sphericalCoordinates( 2 );

    // Calculate vertical <-> aerodynamic <-> body-fixed angles if neede.
    if( calculateVerticalToAerodynamicFrame_ )
    {
        Eigen::Vector3d verticalFrameVelocity =
                getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion(
                    currentAerodynamicAngles_.at( longitude_angle ),
                    currentAerodynamicAngles_.at( latitude_angle ) ) *
                currentBodyFixedState_.segment( 3, 3 );

        currentAerodynamicAngles_[ heading_angle ] = calculateHeadingAngle( verticalFrameVelocity );
        currentAerodynamicAngles_[ flight_path_angle ] =
                calculateFlightPathAngle( verticalFrameVelocity );

        currentAerodynamicAngles_[ angle_of_attack ] = angleOfAttackFunction_( );
        currentAerodynamicAngles_[ angle_of_sideslip ] = angleOfSideslipFunction_( );
        currentAerodynamicAngles_[ bank_angle ] = bankAngleFunction_( );
    }
}

//! Function to get the rotation quaternion between two frames
Eigen::Quaterniond AerodynamicAngleCalculator::getRotationQuaternionBetweenFrames(
        const AerodynamicsReferenceFrames originalFrame,
        const AerodynamicsReferenceFrames targetFrame )
{
    // Initialize rotation to identity matrix.
    Eigen::Quaterniond rotationToFrame = Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) );

    // Inertial frame is not specificied in current object.
    if( originalFrame == inertial_frame || targetFrame == inertial_frame)
    {
        throw( "Error in AerodynamicAngleCalculator, cannot calculate to/from inertial frame" );
    }

    // Check if update settings are consistent with requested frames.
    if( !calculateVerticalToAerodynamicFrame_ &&
            ( originalFrame > vertical_frame || targetFrame > vertical_frame ) )
    {
        throw( "Error in AerodynamicAngleCalculator, instance ends at vertical frame" );
    }

    // Set current frame pair.
    std::pair< AerodynamicsReferenceFrames, AerodynamicsReferenceFrames > currentRotationPair =
            std::make_pair( originalFrame, targetFrame );

    // Calculate rotation matrix if current rotation is not yet calculated.
    if( currentRotationMatrices_.count( currentRotationPair ) == 0 )
    {
        // Get indices of required frames.
        int currentFrameIndex = static_cast< int >( originalFrame );
        int targetFrameIndex = static_cast< int >( targetFrame );

        // Check if any rotation is needed.
        if( currentFrameIndex != targetFrameIndex )
        {
            // Check 'direction' of transformation through AerodynamicsReferenceFrames list.
            bool isTargetFrameUp;
            if( targetFrameIndex > currentFrameIndex )
            {
                isTargetFrameUp = 1;
            }
            else if( targetFrameIndex < currentFrameIndex )
            {
                isTargetFrameUp = 0;
            }
            else
            {
                throw std::runtime_error(
                            "Error when identifying target frame direction in AerodynamicAngleCalculator." );

            }

            // Add rotation sequence until final frame is reached.
            while( currentFrameIndex != targetFrameIndex )
            {
                switch( currentFrameIndex )
                {
                case static_cast< int >( corotating_frame ):

                    if( isTargetFrameUp )
                    {
                        rotationToFrame =
                                getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion(
                                    currentAerodynamicAngles_.at( longitude_angle ),
                                    currentAerodynamicAngles_.at( latitude_angle ) ) *
                                rotationToFrame;
                    }
                    else
                    {
                        throw std::runtime_error(
                                    "Error, corotating_frame is end frame in AerodynamicAngleCalculator" );
                    }
                    break;
                case static_cast< int >( vertical_frame ):
                    if( isTargetFrameUp )
                    {

                        rotationToFrame =
                                getLocalVerticalFrameToTrajectoryTransformationQuaternion(
                                    currentAerodynamicAngles_.at( flight_path_angle ),
                                    currentAerodynamicAngles_.at( heading_angle ) ) * rotationToFrame;
                    }
                    else
                    {
                        rotationToFrame =
                                getLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion(
                                    currentAerodynamicAngles_.at( longitude_angle ),
                                    currentAerodynamicAngles_.at( latitude_angle ) ) *
                                rotationToFrame;
                    }
                    break;
                case static_cast< int >( trajectory_frame ):
                    if( isTargetFrameUp )
                    {
                        rotationToFrame =
                                getTrajectoryToAerodynamicFrameTransformationQuaternion(
                                    currentAerodynamicAngles_.at( bank_angle ) ) *
                                rotationToFrame;

                    }
                    else
                    {
                        rotationToFrame =
                                getTrajectoryToLocalVerticalFrameTransformationQuaternion(
                                    currentAerodynamicAngles_.at( flight_path_angle ),
                                    currentAerodynamicAngles_.at( heading_angle ) ) *
                                rotationToFrame;
                    }
                    break;
                case static_cast< int >( aerodynamic_frame ):
                    if( isTargetFrameUp )
                    {
                        rotationToFrame =
                                getAirspeedBasedAerodynamicToBodyFrameTransformationQuaternion(
                                    currentAerodynamicAngles_.at( angle_of_attack ),
                                    currentAerodynamicAngles_.at( angle_of_sideslip ) ) *
                                rotationToFrame;
                    }
                    else
                    {
                        rotationToFrame =
                                getAerodynamicToTrajectoryFrameTransformationQuaternion(
                                    currentAerodynamicAngles_.at( bank_angle ) ) *
                                rotationToFrame;
                    }
                    break;
                case static_cast< int >( body_frame ):
                    if( isTargetFrameUp )
                    {

                        throw std::runtime_error(
                                    "Error, body frame is end frame in AerodynamicAngleCalculator." );
                    }
                    else
                    {
                        rotationToFrame =
                                getBodyToAirspeedBasedAerodynamicFrameTransformationQuaternion(
                                    currentAerodynamicAngles_.at( angle_of_attack ),
                                    currentAerodynamicAngles_.at( angle_of_sideslip ) ) *
                                rotationToFrame;
                    }
                    break;
                default:
                    throw std::runtime_error(
                                "Error, index " +
                                boost::lexical_cast< std::string>( currentFrameIndex ) +
                                "not found in AerodynamicAngleCalculator" );

                }

                // Increment/decrement current frame.
                if( isTargetFrameUp )
                {
                    currentFrameIndex++;
                }
                else
                {
                    currentFrameIndex--;
                }
            }
        }

        // Set current rotation (as well as inverse).
        currentRotationMatrices_[ currentRotationPair ] = rotationToFrame;
        currentRotationMatrices_[ std::make_pair( targetFrame, originalFrame ) ] =
                rotationToFrame.inverse( );
    }
    else
    {
        rotationToFrame = currentRotationMatrices_.at( currentRotationPair );
    }
    return rotationToFrame;
}

//! Function to get a single orientation angle.
double AerodynamicAngleCalculator::getAerodynamicAngle(
        const AerodynamicsReferenceFrameAngles angleId )
{
    double angleValue = TUDAT_NAN;
    if( currentAerodynamicAngles_.count( angleId ) == 0 )
    {
        throw std::runtime_error( "Error in AerodynamicAngleCalculator, angleId " +
                                  boost::lexical_cast< std::string >( angleId ) + "not found" );
    }
    else
    {
        angleValue = currentAerodynamicAngles_.at( angleId );
    }
    return angleValue;
}

//! Get a function to transform aerodynamic force from local to propagation frame.
boost::function< Eigen::Vector3d( const Eigen::Vector3d& ) >
getAerodynamicForceTransformationFunction(
        const boost::shared_ptr< AerodynamicAngleCalculator > aerodynamicAngleCalculator,
        const AerodynamicsReferenceFrames accelerationFrame,
        const boost::function< Eigen::Quaterniond( ) > bodyFixedToInertialFrameFunction,
        const AerodynamicsReferenceFrames propagationFrame )
{
    boost::function< Eigen::Vector3d( const Eigen::Vector3d& ) > transformationFunction;

    // If propagation frame is the inertial frame, use bodyFixedToInertialFrameFunction.
    if( propagationFrame == inertial_frame )
    {
        std::vector< boost::function< Eigen::Vector3d( const Eigen::Vector3d& ) > > rotationsList;

        // Get accelerationFrame to corotating frame transformation.
        boost::function< Eigen::Quaterniond( ) > firstRotation =
                boost::bind( &AerodynamicAngleCalculator::getRotationQuaternionBetweenFrames,
                             aerodynamicAngleCalculator, accelerationFrame, corotating_frame );
        rotationsList.push_back(
                    boost::bind(
                        static_cast< Eigen::Vector3d(&)(
                            const Eigen::Vector3d&,
                            const boost::function< Eigen::Quaterniond( ) > ) >( &transformVector ),
                        _1, firstRotation ) );

        // Add corotating to inertial frame.
        rotationsList.push_back(
                    boost::bind(
                        static_cast< Eigen::Vector3d(&)(
                            const Eigen::Vector3d&,
                            const boost::function< Eigen::Quaterniond( ) > ) >( &transformVector ),
                        _1, bodyFixedToInertialFrameFunction ) );

        // Create transformation function.
        transformationFunction = boost::bind(
                    static_cast< Eigen::Vector3d(&)(
                        const Eigen::Vector3d&,
                        const std::vector< boost::function<
                        Eigen::Vector3d( const Eigen::Vector3d& ) > >& ) >( &transformVector ),
                    _1, rotationsList );
    }
    else
    {
        // Get accelerationFrame to propagationFrame frame transformation directly.
        boost::function< Eigen::Quaterniond( ) > rotationFunction =
                boost::bind( &AerodynamicAngleCalculator::getRotationQuaternionBetweenFrames,
                             aerodynamicAngleCalculator, accelerationFrame, propagationFrame );

        // Create transformation function.
        transformationFunction = boost::bind(
                    static_cast< Eigen::Vector3d(&)(
                        const Eigen::Vector3d&,
                        const boost::function< Eigen::Quaterniond( ) > ) >( &transformVector ), _1,
                    rotationFunction );
    }

    return transformationFunction;
}

} // namespace reference_frames

} // namespace tudat



