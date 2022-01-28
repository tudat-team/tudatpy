/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <iostream>

#include <boost/bind/bind.hpp>
using namespace boost::placeholders;

#include <memory>
#include <boost/make_shared.hpp>

#include "tudat/astro/basic_astro/sphericalStateConversions.h"
#include "tudat/astro/reference_frames/aerodynamicAngleCalculator.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/basic/rotationRepresentations.h"

namespace tudat
{

namespace reference_frames
{

//! Function to get a string representing a 'named identification' of a reference frame.
std::string getAerodynamicFrameName( const AerodynamicsReferenceFrames frame )
{
    std::string frameName;
    switch( frame )
    {
    case inertial_frame:
        frameName = "inertial frame ";
        break;
    case corotating_frame:
        frameName = "corotating frame ";
        break;
    case vertical_frame:
        frameName = "vertical frame ";
        break;
    case trajectory_frame:
        frameName = "trajectory frame ";
        break;
    case aerodynamic_frame:
        frameName = "aerodynamic frame ";
        break;
    case body_frame:
        frameName = "body frame ";
        break;
    default:
        std::string errorMessage = "Error, aerodynamic frame type " +
                std::to_string( frame ) +
                "not found when retrieving frame name ";
        throw std::runtime_error( errorMessage );
    }
    return frameName;
}

//! Function to get a string representing a 'named identification' of an aerodynamic angle
std::string getAerodynamicAngleName( const AerodynamicsReferenceFrameAngles angle )
{
    std::string angleName;
    switch( angle )
    {
    case latitude_angle:
        angleName = "latitude angle ";
        break;
    case longitude_angle:
        angleName = "longitude angle ";
        break;
    case heading_angle:
        angleName = "heading angle ";
        break;
    case flight_path_angle:
        angleName = "flight path angle ";
        break;
    case angle_of_attack:
        angleName = "angle of attack ";
        break;
    case angle_of_sideslip:
        angleName = "sideslip angle ";
        break;
    case bank_angle:
        angleName = "bank angle ";
        break;
    default:
        std::string errorMessage = "Error, aerodynamic angle type " +
                std::to_string( angle ) +
                "not found when retrieving angle name ";
        throw std::runtime_error( errorMessage );
    }
    return angleName;
}


//! Function to update the orientation angles to the current state.
void AerodynamicAngleCalculator::update( const double currentTime, const bool updateBodyOrientation )
{
    // Clear all current rotation matrices.
    currentRotationMatrices_.clear( );

    // Get current body-fixed state.
    if( !( currentTime == currentTime_ ) )
    {
        currentBodyFixedGroundSpeedBasedState_ = bodyFixedStateFunction_( );
        currentRotationFromCorotatingToInertialFrame_ = rotationFromCorotatingToInertialFrame_( );

        Eigen::Vector3d sphericalCoordinates = coordinate_conversions::convertCartesianToSpherical< double >(
                    currentBodyFixedGroundSpeedBasedState_.segment( 0, 3 ) );

        // Calculate latitude and longitude.
        currentAerodynamicAngles_[ latitude_angle ] =
                mathematical_constants::PI / 2.0 - sphericalCoordinates( 1 );
        currentAerodynamicAngles_[ longitude_angle ] = sphericalCoordinates( 2 );

        // Compute wind velocity vector
        Eigen::Vector3d localWindVelocity = Eigen::Vector3d::Zero( );
        if( windModel_ != nullptr )
        {
            localWindVelocity = getRotationQuaternionBetweenFrames(
                        windModel_->getAssociatedFrame( ), corotating_frame ) *
                    windModel_->getCurrentBodyFixedCartesianWindVelocity(
                        shapeModel_->getAltitude( currentBodyFixedGroundSpeedBasedState_.segment( 0, 3 ) ),
                        currentAerodynamicAngles_[ longitude_angle ],
                        currentAerodynamicAngles_[ latitude_angle ],
                        currentTime );
        }

        // Compute airspeed-based velocity vector
        currentBodyFixedAirspeedBasedState_ = currentBodyFixedGroundSpeedBasedState_;
        currentBodyFixedAirspeedBasedState_.segment( 3, 3 ) -= localWindVelocity;

        // Calculate vertical <-> aerodynamic <-> body-fixed angles if neede.
        if( calculateVerticalToAerodynamicFrame_ )
        {
            Eigen::Vector3d verticalFrameVelocity =
                    getRotationQuaternionBetweenFrames( corotating_frame, vertical_frame ) *
                    currentBodyFixedAirspeedBasedState_.segment( 3, 3 );

            currentAerodynamicAngles_[ heading_angle ] = calculateHeadingAngle( verticalFrameVelocity );
            currentAerodynamicAngles_[ flight_path_angle ] =
                    calculateFlightPathAngle( verticalFrameVelocity );
        }

        currentTime_ = currentTime;
    }

    if( updateBodyOrientation  && !( currentBodyAngleTime_ == currentTime ) )
    {
        if( !( angleUpdateFunction_ == nullptr ) )
        {
            angleUpdateFunction_( currentTime );
        }

        if( !( angleOfAttackFunction_ == nullptr ) )
        {
            currentAerodynamicAngles_[ angle_of_attack ] = angleOfAttackFunction_( );
        }

        if( !( angleOfSideslipFunction_ == nullptr ) )
        {
            currentAerodynamicAngles_[ angle_of_sideslip ] = angleOfSideslipFunction_( );
        }

        if( !( bankAngleFunction_ == nullptr ) )
        {
            currentAerodynamicAngles_[ bank_angle ] = bankAngleFunction_( );
        }

        currentBodyAngleTime_ = currentTime;
    }
    else if( !( currentBodyAngleTime_ == currentTime ) )
    {
        currentAerodynamicAngles_[ angle_of_attack ] = 0.0;
        currentAerodynamicAngles_[ angle_of_sideslip ] = 0.0;
        currentAerodynamicAngles_[ bank_angle ] = 0.0;
    }
}


//! Function to get the rotation quaternion between two frames
void AerodynamicAngleCalculator::getRotationQuaternionReferenceBetweenFrames(
        Eigen::Quaterniond& rotationToFrame,
        const AerodynamicsReferenceFrames originalFrame,
        const AerodynamicsReferenceFrames targetFrame )
{
    // Initialize rotation to identity matrix.
    rotationToFrame = Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) );

    // Check if update settings are consistent with requested frames.
    if( !calculateVerticalToAerodynamicFrame_ &&
            ( originalFrame > vertical_frame || targetFrame > vertical_frame ) )
    {
        throw std::runtime_error( "Error in AerodynamicAngleCalculator, instance ends at vertical frame" );
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
                case static_cast< int >( inertial_frame ):
                    if( isTargetFrameUp )
                    {
                        rotationToFrame = currentRotationFromCorotatingToInertialFrame_.inverse( ) *
                                rotationToFrame;
                    }
                    else
                    {
                        throw std::runtime_error(
                                    "Error, inertial_frame is end frame in AerodynamicAngleCalculator" );
                    }
                    break;
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
                        rotationToFrame = currentRotationFromCorotatingToInertialFrame_ *
                                rotationToFrame;
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
                        if( ! ( currentAerodynamicAngles_.at( bank_angle ) == 0 ) )
                        {
                            rotationToFrame =
                                    getTrajectoryToAerodynamicFrameTransformationQuaternion(
                                        currentAerodynamicAngles_.at( bank_angle ) ) *
                                    rotationToFrame;
                        }
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
                        if( !( currentAerodynamicAngles_.at( angle_of_attack ) == 0 &&
                               currentAerodynamicAngles_.at( angle_of_sideslip ) == 0 ) )
                        {
                            rotationToFrame =
                                    getAirspeedBasedAerodynamicToBodyFrameTransformationQuaternion(
                                        currentAerodynamicAngles_.at( angle_of_attack ),
                                        currentAerodynamicAngles_.at( angle_of_sideslip ) ) *
                                    rotationToFrame;
                        }
                    }
                    else
                    {
                        if( ! ( currentAerodynamicAngles_.at( bank_angle ) == 0 ) )
                        {
                            rotationToFrame =
                                    getAerodynamicToTrajectoryFrameTransformationQuaternion(
                                        currentAerodynamicAngles_.at( bank_angle ) ) *
                                    rotationToFrame;
                        }
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
                                "Error, index " + std::to_string( currentFrameIndex ) +
                                "not found in AerodynamicAngleCalculator." );
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
}

//! Function to get the rotation quaternion between two frames
Eigen::Quaterniond AerodynamicAngleCalculator::getRotationQuaternionBetweenFrames(
        const AerodynamicsReferenceFrames originalFrame,
        const AerodynamicsReferenceFrames targetFrame )
{
    Eigen::Quaterniond rotationToFrame;
    getRotationQuaternionReferenceBetweenFrames( rotationToFrame, originalFrame, targetFrame );
    return rotationToFrame;
}

//! Function to get a single orientation angle.
double AerodynamicAngleCalculator::getAerodynamicAngle(
        const AerodynamicsReferenceFrameAngles angleId )
{
    return currentAerodynamicAngles_.at( angleId );
}

//! Function to set the trajectory<->body-fixed orientation angles.
void AerodynamicAngleCalculator::setOrientationAngleFunctions(
        const std::function< double( ) > angleOfAttackFunction,
        const std::function< double( ) > angleOfSideslipFunction,
        const std::function< double( ) > bankAngleFunction,
        const std::function< void( const double ) > angleUpdateFunction,
        const bool silenceWarnings )
{
    if( !( angleOfAttackFunction == nullptr ) )
    {
        if( !( angleOfAttackFunction_ == nullptr ) && !silenceWarnings )
        {
            std::cerr << "Warning, overriding existing angle of attack function in AerodynamicAngleCalculator" << std::endl;
        }
        angleOfAttackFunction_ = angleOfAttackFunction;
    }

    if( !( angleOfSideslipFunction == nullptr ) )
    {
        if( !( angleOfSideslipFunction_ == nullptr ) && !silenceWarnings  )
        {
            std::cerr << "Warning, overriding existing angle of sideslip function in AerodynamicAngleCalculator" << std::endl;
        }
        angleOfSideslipFunction_ = angleOfSideslipFunction;
    }

    if( !( bankAngleFunction == nullptr ) )
    {
        if( !( bankAngleFunction_ == nullptr ) && !silenceWarnings  )
        {
            std::cerr << "Warning, overriding existing bank angle function in AerodynamicAngleCalculator" << std::endl;
        }
        bankAngleFunction_ = bankAngleFunction;
    }

    if( !( angleUpdateFunction == nullptr ) )
    {
        if( !( angleUpdateFunction_ == nullptr ) && !silenceWarnings  )
        {
            std::cerr << "Warning, overriding existing aerodynamic angle update function in AerodynamicAngleCalculator" << std::endl;
        }
        angleUpdateFunction_ = angleUpdateFunction;
    }
}

//! Function to set constant trajectory<->body-fixed orientation angles.
void AerodynamicAngleCalculator::setOrientationAngleFunctions(
        const double angleOfAttack,
        const double angleOfSideslip,
        const double bankAngle,
        const bool silenceWarnings )
{
    std::function< double( ) > angleOfAttackFunction =
            ( ( angleOfAttack == angleOfAttack ) ? [ = ]( ){ return angleOfAttack; } : std::function< double( ) >( ) );
    std::function< double( ) > angleOfSideslipFunction =
            ( ( angleOfSideslip == angleOfSideslip ) ? [ = ]( ){ return angleOfSideslip; } : std::function< double( ) >( ) );
    std::function< double( ) > bankAngleFunction =
            ( ( bankAngle == bankAngle ) ? [ = ]( ){ return bankAngle; }: std::function< double( ) >( ) );
    setOrientationAngleFunctions( angleOfAttackFunction, angleOfSideslipFunction, bankAngleFunction, nullptr, silenceWarnings );
}

//! Get a function to transform aerodynamic force from local to propagation frame.
std::function< Eigen::Vector3d( const Eigen::Vector3d& ) >
getAerodynamicForceTransformationFunction(
        const std::shared_ptr< AerodynamicAngleCalculator > aerodynamicAngleCalculator,
        const AerodynamicsReferenceFrames accelerationFrame,
        const std::function< Eigen::Quaterniond( ) > bodyFixedToInertialFrameFunction,
        const AerodynamicsReferenceFrames propagationFrame )
{
    std::function< Eigen::Vector3d( const Eigen::Vector3d& ) > transformationFunction;

    // If propagation frame is the inertial frame, use bodyFixedToInertialFrameFunction.
    if( propagationFrame == inertial_frame )
    {
        std::vector< std::function< Eigen::Vector3d( const Eigen::Vector3d& ) > > rotationsList;

        // Get accelerationFrame to corotating frame transformation.
        std::function< Eigen::Quaterniond( ) > firstRotation =
                std::bind( &AerodynamicAngleCalculator::getRotationQuaternionBetweenFrames,
                           aerodynamicAngleCalculator, accelerationFrame, corotating_frame );
        rotationsList.push_back(
                    std::bind( &transformVectorFromQuaternionFunction,
                               std::placeholders::_1, firstRotation ) );

        // Add corotating to inertial frame.
        rotationsList.push_back(
                    std::bind( &transformVectorFromQuaternionFunction,
                               std::placeholders::_1, bodyFixedToInertialFrameFunction ) );

        // Create transformation function.
        transformationFunction = std::bind( &transformVectorFromVectorFunctions,
                                            std::placeholders::_1, rotationsList );
    }
    else
    {
        // Get accelerationFrame to propagationFrame frame transformation directly.
        std::function< Eigen::Quaterniond( ) > rotationFunction =
                std::bind( &AerodynamicAngleCalculator::getRotationQuaternionBetweenFrames,
                           aerodynamicAngleCalculator, accelerationFrame, propagationFrame );

        // Create transformation function.
        transformationFunction = std::bind( &transformVectorFromQuaternionFunction, std::placeholders::_1,
                                            rotationFunction );
    }

    return transformationFunction;
}

std::function< void( Eigen::Vector3d&, const Eigen::Vector3d& ) >
getAerodynamicForceTransformationReferenceFunction(
        const std::shared_ptr< AerodynamicAngleCalculator > aerodynamicAngleCalculator,
        const AerodynamicsReferenceFrames accelerationFrame,
        const std::function< Eigen::Quaterniond&( ) > bodyFixedToInertialFrameFunction,
        const AerodynamicsReferenceFrames propagationFrame )
{
    std::function< void( Eigen::Vector3d&, const Eigen::Vector3d& ) > transformationFunction;

    // If propagation frame is the inertial frame, use bodyFixedToInertialFrameFunction.
    if( propagationFrame == inertial_frame )
    {
        std::vector< std::function< Eigen::Vector3d( const Eigen::Vector3d& ) > > rotationsList;

        // Get accelerationFrame to corotating frame transformation.
        std::function< Eigen::Quaterniond( ) > firstRotation =
                std::bind( &AerodynamicAngleCalculator::getRotationQuaternionBetweenFrames,
                           aerodynamicAngleCalculator, accelerationFrame, corotating_frame );
        rotationsList.push_back(
                    std::bind( &transformVectorFromQuaternionFunction,
                               std::placeholders::_1, firstRotation ) );

        // Add corotating to inertial frame.
        rotationsList.push_back(
                    std::bind( &transformVectorFromQuaternionFunction,
                               std::placeholders::_1, bodyFixedToInertialFrameFunction ) );

        // Create transformation function.
        transformationFunction = std::bind( &transformVectorReferenceFromVectorFunctions,
                                            std::placeholders::_1, std::placeholders::_2, rotationsList );
    }
    else
    {
        throw std::runtime_error( "TODO: Write error message" );
    }

    return transformationFunction;
}

//! Function to update the aerodynamic angles to current time.
void AerodynamicAnglesClosure::updateAngles( const double currentTime )
{
    // Retrieve rotation matrix that is to be converted to orientation angles.
    currentRotationFromBodyToTrajectoryFrame_ =
            ( ( imposedRotationFromInertialToBodyFixedFrame_( currentTime ).toRotationMatrix( ) *
                aerodynamicAngleCalculator_->getRotationQuaternionBetweenFrames(
                    trajectory_frame, inertial_frame ).toRotationMatrix( ) ) ).transpose( );

    // Compute associated Euler angles and set as orientation angles.
    Eigen::Vector3d eulerAngles = basic_mathematics::get132EulerAnglesFromRotationMatrix(
                currentRotationFromBodyToTrajectoryFrame_ );
    currentBankAngle_ = eulerAngles( 0 );
    currentAngleOfSideslip_ = eulerAngles( 1 );
    currentAngleOfAttack_ = -eulerAngles( 2 );
}

//! Function to make aerodynamic angle computation consistent with imposed body-fixed to inertial rotation.
void setAerodynamicDependentOrientationCalculatorClosure(
        const std::function< Eigen::Quaterniond( const double ) > imposedRotationFromInertialToBodyFixedFrame,
        std::shared_ptr< AerodynamicAngleCalculator > aerodynamicAngleCalculator )
{
    std::shared_ptr< AerodynamicAnglesClosure > aerodynamicAnglesClosure =
            std::make_shared< AerodynamicAnglesClosure >(
                imposedRotationFromInertialToBodyFixedFrame, aerodynamicAngleCalculator );
    aerodynamicAngleCalculator->setOrientationAngleFunctions(
                std::bind( &AerodynamicAnglesClosure::getCurrentAngleOfAttack, aerodynamicAnglesClosure ),
                std::bind( &AerodynamicAnglesClosure::getCurrentAngleOfSideslip, aerodynamicAnglesClosure ),
                std::bind( &AerodynamicAnglesClosure::getCurrentBankAngle, aerodynamicAnglesClosure ),
                std::bind( &AerodynamicAnglesClosure::updateAngles, aerodynamicAnglesClosure, std::placeholders::_1 ) );
}

//! Function to make aerodynamic angle computation consistent with existing DependentOrientationCalculator
void setAerodynamicDependentOrientationCalculatorClosure(
        std::shared_ptr< DependentOrientationCalculator > dependentOrientationCalculator,
        std::shared_ptr< AerodynamicAngleCalculator > aerodynamicAngleCalculator )
{
    std::function< Eigen::Quaterniond( const double ) > imposedRotationFromInertialToBodyFixedFrame =
            std::bind( &DependentOrientationCalculator::computeAndGetRotationToLocalFrame, dependentOrientationCalculator, std::placeholders::_1 );
    setAerodynamicDependentOrientationCalculatorClosure(
                imposedRotationFromInertialToBodyFixedFrame,
                aerodynamicAngleCalculator );
}

//! Function to make aerodynamic angle computation consistent with existing rotational ephemeris
void setAerodynamicDependentOrientationCalculatorClosure(
        std::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris,
        std::shared_ptr< AerodynamicAngleCalculator > aerodynamicAngleCalculator )
{
    setAerodynamicDependentOrientationCalculatorClosure(
                std::bind( &ephemerides::RotationalEphemeris::getRotationToTargetFrame,
                           rotationalEphemeris, std::placeholders::_1 ),
                aerodynamicAngleCalculator );
}

} // namespace reference_frames

} // namespace tudat



