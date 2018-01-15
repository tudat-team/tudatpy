/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <iostream>

#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalStateConversions.h"
#include "Tudat/Astrodynamics/ReferenceFrames/aerodynamicAngleCalculator.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

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
        if( windModel_ != NULL )
        {
            localWindVelocity = windModel_->getCurrentWindVelocity(
                        shapeModel_->getAltitude( currentBodyFixedGroundSpeedBasedState_.segment( 0, 3 ) ),
                        currentAerodynamicAngles_[ longitude_angle ],
                        currentAerodynamicAngles_[ latitude_angle ],
                        currentTime );
        }

        // Compute airspeed-based velocity vector
        currentBodyFixedAirspeedBasedState_ = currentBodyFixedGroundSpeedBasedState_;
        currentBodyFixedAirspeedBasedState_.segment( 3, 3 ) += localWindVelocity;

        // Calculate vertical <-> aerodynamic <-> body-fixed angles if neede.
        if( calculateVerticalToAerodynamicFrame_ )
        {
            Eigen::Vector3d verticalFrameVelocity =
                    getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion(
                        currentAerodynamicAngles_.at( longitude_angle ),
                        currentAerodynamicAngles_.at( latitude_angle ) ) *
                    currentBodyFixedAirspeedBasedState_.segment( 3, 3 );

            currentAerodynamicAngles_[ heading_angle ] = calculateHeadingAngle( verticalFrameVelocity );
            currentAerodynamicAngles_[ flight_path_angle ] =
                    calculateFlightPathAngle( verticalFrameVelocity );
        }

        currentTime_ = currentTime;
    }

    if( updateBodyOrientation  && !( currentBodyAngleTime_ == currentTime ) )
    {
        if( !angleUpdateFunction_.empty( ) )
        {
            angleUpdateFunction_( currentTime );
        }

        if( !angleOfAttackFunction_.empty( ) )
        {
            currentAerodynamicAngles_[ angle_of_attack ] = angleOfAttackFunction_( );
        }

        if( !angleOfSideslipFunction_.empty( ) )
        {
            currentAerodynamicAngles_[ angle_of_sideslip ] = angleOfSideslipFunction_( );
        }

        if( !bankAngleFunction_.empty( ) )
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
Eigen::Quaterniond AerodynamicAngleCalculator::getRotationQuaternionBetweenFrames(
        const AerodynamicsReferenceFrames originalFrame,
        const AerodynamicsReferenceFrames targetFrame )
{
    // Initialize rotation to identity matrix.
    Eigen::Quaterniond rotationToFrame = Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) );


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
                                std::to_string( currentFrameIndex ) +
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
                                  std::to_string( angleId ) + "not found" );
    }
    else
    {
        angleValue = currentAerodynamicAngles_.at( angleId );
    }
    return angleValue;
}

//! Function to set the trajectory<->body-fixed orientation angles.
void AerodynamicAngleCalculator::setOrientationAngleFunctions(
        const boost::function< double( ) > angleOfAttackFunction,
        const boost::function< double( ) > angleOfSideslipFunction,
        const boost::function< double( ) > bankAngleFunction,
        const boost::function< void( const double ) > angleUpdateFunction )
{
    if( !angleOfAttackFunction.empty( ) )
    {
        if( !angleOfAttackFunction_.empty( ) )
        {
            std::cerr << "Warning, overriding existing angle of attack function in AerodynamicAngleCalculator" << std::endl;
        }
        angleOfAttackFunction_ = angleOfAttackFunction;
    }

    if( !angleOfSideslipFunction.empty( ) )
    {
        if( !angleOfSideslipFunction_.empty( ) )
        {
            std::cerr << "Warning, overriding existing angle of sideslip function in AerodynamicAngleCalculator" << std::endl;
        }
        angleOfSideslipFunction_ = angleOfSideslipFunction;
    }

    if( !bankAngleFunction.empty( ) )
    {
        if( !bankAngleFunction_.empty( ) )
        {
            std::cerr << "Warning, overriding existing bank angle function in AerodynamicAngleCalculator" << std::endl;
        }
        bankAngleFunction_ = bankAngleFunction;
    }

    if( !angleUpdateFunction.empty( ) )
    {
        if( !angleUpdateFunction_.empty( ) )
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
        const double bankAngle )
{
    boost::function< double( ) > angleOfAttackFunction =
            ( ( angleOfAttack == angleOfAttack ) ? boost::lambda::constant( angleOfAttack ) : boost::function< double( ) >( ) );
    boost::function< double( ) > angleOfSideslipFunction =
            ( ( angleOfSideslip == angleOfSideslip ) ? boost::lambda::constant( angleOfSideslip ) : boost::function< double( ) >( ) );
    boost::function< double( ) > bankAngleFunction =
            ( ( bankAngle == bankAngle ) ? boost::lambda::constant( bankAngle ): boost::function< double( ) >( ) );
    setOrientationAngleFunctions( angleOfAttackFunction, angleOfSideslipFunction, bankAngleFunction );
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
                    boost::bind( &transformVectorFromQuaternionFunction,
                                 _1, firstRotation ) );

        // Add corotating to inertial frame.
        rotationsList.push_back(
                    boost::bind( &transformVectorFromQuaternionFunction,
                                 _1, bodyFixedToInertialFrameFunction ) );

        // Create transformation function.
        transformationFunction = boost::bind( &transformVectorFromVectorFunctions,
                                              _1, rotationsList );
    }
    else
    {
        // Get accelerationFrame to propagationFrame frame transformation directly.
        boost::function< Eigen::Quaterniond( ) > rotationFunction =
                boost::bind( &AerodynamicAngleCalculator::getRotationQuaternionBetweenFrames,
                             aerodynamicAngleCalculator, accelerationFrame, propagationFrame );

        // Create transformation function.
        transformationFunction = boost::bind( &transformVectorFromQuaternionFunction, _1,
                                              rotationFunction );
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
    Eigen::Vector3d eulerAngles = reference_frames::get132EulerAnglesFromRotationMatrix(
                currentRotationFromBodyToTrajectoryFrame_ );
    currentBankAngle_ = eulerAngles.x( );
    currentAngleOfSideslip_ = eulerAngles.y( );
    currentAngleOfAttack_ = -eulerAngles.z( );
}

//! Function to make aerodynamic angle computation consistent with imposed body-fixed to inertial rotation.
void setAerodynamicDependentOrientationCalculatorClosure(
        const boost::function< Eigen::Quaterniond( const double ) > imposedRotationFromInertialToBodyFixedFrame,
        boost::shared_ptr< AerodynamicAngleCalculator > aerodynamicAngleCalculator )
{
    boost::shared_ptr< AerodynamicAnglesClosure > aerodynamicAnglesClosure =
            boost::make_shared< AerodynamicAnglesClosure >(
                imposedRotationFromInertialToBodyFixedFrame, aerodynamicAngleCalculator );
    aerodynamicAngleCalculator->setOrientationAngleFunctions(
                boost::bind( &AerodynamicAnglesClosure::getCurrentAngleOfAttack, aerodynamicAnglesClosure ),
                boost::bind( &AerodynamicAnglesClosure::getCurrentAngleOfSideslip, aerodynamicAnglesClosure ),
                boost::bind( &AerodynamicAnglesClosure::getCurrentBankAngle, aerodynamicAnglesClosure ),
                boost::bind( &AerodynamicAnglesClosure::updateAngles, aerodynamicAnglesClosure, _1 ) );
}

//! Function to make aerodynamic angle computation consistent with existing DependentOrientationCalculator
void setAerodynamicDependentOrientationCalculatorClosure(
        boost::shared_ptr< DependentOrientationCalculator > dependentOrientationCalculator,
        boost::shared_ptr< AerodynamicAngleCalculator > aerodynamicAngleCalculator )
{
    setAerodynamicDependentOrientationCalculatorClosure(
                boost::bind( &DependentOrientationCalculator::getRotationToLocalFrame, dependentOrientationCalculator, _1 ),
                aerodynamicAngleCalculator );
}

//! Function to make aerodynamic angle computation consistent with existing rotational ephemeris
void setAerodynamicDependentOrientationCalculatorClosure(
        boost::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris,
        boost::shared_ptr< AerodynamicAngleCalculator > aerodynamicAngleCalculator )
{
    setAerodynamicDependentOrientationCalculatorClosure(
                boost::bind( &ephemerides::RotationalEphemeris::getRotationToTargetFrame,
                             rotationalEphemeris, _1 ),
                aerodynamicAngleCalculator );
}

} // namespace reference_frames

} // namespace tudat



