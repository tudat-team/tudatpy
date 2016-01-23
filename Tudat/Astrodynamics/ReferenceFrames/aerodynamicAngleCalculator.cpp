#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/ReferenceFrames/aerodynamicAngleCalculator.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{

namespace reference_frames
{

void AerodynamicAngleCalculator::update( )
{
    currentRotationMatrices_.clear( );

    currentBodyFixedState_ = bodyFixedStateFunction_( );
    Eigen::Vector3d sphericalCoordinates = coordinate_conversions::convertCartesianToSpherical(
                currentBodyFixedState_.segment( 0, 3 ) );

    currentAerodynamicAngles_[ latitude_angle ] = mathematical_constants::PI / 2.0 - sphericalCoordinates( 1 );
    currentAerodynamicAngles_[ longitude_angle ] = sphericalCoordinates( 2 );

    if( calculateVerticalToAerodynamicFrame_ )
    {
        Eigen::Vector3d verticalFrameVelocity =
                getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion(
                    currentAerodynamicAngles_.at( longitude_angle ),
                    currentAerodynamicAngles_.at( latitude_angle ) ) * currentBodyFixedState_.segment(
                    3, 3 );

        currentAerodynamicAngles_[ heading_angle ] = std::atan2(
                    verticalFrameVelocity( 1 ), verticalFrameVelocity( 0 ) );
        currentAerodynamicAngles_[ flight_path_angle ] = std::asin(
                    verticalFrameVelocity( 2 ) / verticalFrameVelocity.norm( ) );

        currentAerodynamicAngles_[ angle_of_attack ] = angleOfAttackFunction_( );
        currentAerodynamicAngles_[ angle_of_sideslip ] = angleOfSideslipFunction_( );
        currentAerodynamicAngles_[ bank_angle ] = bankAngleFunction_( );
    }
}

Eigen::Quaterniond AerodynamicAngleCalculator::getRotationQuaternionBetweenFrames(
        const AerodynamicsReferenceFrames originalFrame,
        const AerodynamicsReferenceFrames targetFrame )
{
    Eigen::Quaterniond rotationToFrame = Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) );

    if( originalFrame == inertial_frame || targetFrame == inertial_frame)
    {
        throw( "" );
    }

    // Check if update settings are consistent with requested frames.
    if( !calculateVerticalToAerodynamicFrame_ &&
            ( originalFrame > vertical_frame || targetFrame > vertical_frame ) )
    {
        throw( "" );
    }

    // Set current frame pair.
    std::pair< AerodynamicsReferenceFrames, AerodynamicsReferenceFrames > currentRotationPair =
            std::make_pair( originalFrame, targetFrame );

    // Calculate rotation matrix if current rotation is not yet calculated.
    if( currentRotationMatrices_.count( currentRotationPair ) == 0 )
    {
        int currentFrameIndex = static_cast< int >( originalFrame );
        int targetFrameIndex = static_cast< int >( targetFrame );

        if( currentFrameIndex != targetFrameIndex )
        {
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
                throw( "" );

            }
            while( currentFrameIndex != targetFrameIndex )
            {
                switch( currentFrameIndex )
                {
                case static_cast< int >( corotating_frame ):
                {
                    if( isTargetFrameUp )
                    {
                        rotationToFrame =
                                getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion(
                                    currentAerodynamicAngles_.at( longitude_angle ),
                                    currentAerodynamicAngles_.at( latitude_angle ) ) * rotationToFrame;
                    }
                    else
                    {
                        throw( "" );
                    }
                }
                case static_cast< int >( vertical_frame ):
                    if( isTargetFrameUp )
                    {

                        rotationToFrame =
                                getLocalVerticalFrameToTrajectoryTransformationQuaternion(
                                    flight_path_angle, heading_angle ) * rotationToFrame;
                    }
                    else
                    {
                        rotationToFrame =
                                getLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion(
                                    currentAerodynamicAngles_.at( longitude_angle ),
                                    currentAerodynamicAngles_.at( latitude_angle ) ) * rotationToFrame;
                    }
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
                                    currentAerodynamicAngles_.at( heading_angle ) ) * rotationToFrame;
                    }
                case static_cast< int >( aerodynamic_frame ):
                    if( isTargetFrameUp )
                    {
                        rotationToFrame =
                                getAirspeedBasedAerodynamicToBodyFrameTransformationQuaternion(
                                    currentAerodynamicAngles_.at( angle_of_attack ),
                                    currentAerodynamicAngles_.at( angle_of_sideslip ) ) * rotationToFrame;
                    }
                    else
                    {
                        rotationToFrame =
                                getAerodynamicToTrajectoryFrameTransformationQuaternion(
                                    currentAerodynamicAngles_.at( bank_angle ) ) *
                                rotationToFrame;
                    }
                case static_cast< int >( body_frame ):
                    if( isTargetFrameUp )
                    {
                        throw( "" );
                    }
                    else
                    {
                        rotationToFrame =
                                getBodyToAirspeedBasedAerodynamicFrameTransformationQuaternion(
                                    currentAerodynamicAngles_.at( angle_of_attack ),
                                    currentAerodynamicAngles_.at( angle_of_sideslip ) ) * rotationToFrame;
                    }
                default:
                    throw( "" );
                }
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

double AerodynamicAngleCalculator::getAerodynamicAngle( const AerodynamicsReferenceFrameAngles angleId )
{
    double angleValue = TUDAT_NAN;
    if( currentAerodynamicAngles_.count( angleId ) == 0 )
    {
        throw( "" );
    }
    else
    {
        angleValue = currentAerodynamicAngles_.at( angleId );
    }
    return angleValue;
}

boost::function< Eigen::Vector3d( const Eigen::Vector3d& ) >
getAerodynamicForceTransformationFunction(
        const boost::shared_ptr< AerodynamicAngleCalculator > aerodynamicAngleCalculator,
        const AerodynamicsReferenceFrames accelerationFrame,
        const boost::function< Eigen::Quaterniond( ) > bodyFixedToInertialFrameFunction,
        const AerodynamicsReferenceFrames propagationFrame )
{
    boost::function< Eigen::Vector3d( const Eigen::Vector3d& ) > transformationFunction;
    AerodynamicsReferenceFrames aerodynamicBaseFrame = propagationFrame;

    if( aerodynamicBaseFrame == inertial_frame )
    {
        std::vector< boost::function< Eigen::Vector3d( const Eigen::Vector3d& ) > > rotationsList;

        boost::function< Eigen::Quaterniond( ) > firstRotation =
                boost::bind( &AerodynamicAngleCalculator::getRotationQuaternionBetweenFrames,
                             aerodynamicAngleCalculator, accelerationFrame, corotating_frame );
        rotationsList.push_back(
                    boost::bind(
                        static_cast< Eigen::Vector3d(&)(
                            const Eigen::Vector3d&,
                            const boost::function< Eigen::Quaterniond( ) > ) >( &transformVector ),
                        _1, firstRotation ) );
        rotationsList.push_back(
                    boost::bind(
                        static_cast< Eigen::Vector3d(&)(
                            const Eigen::Vector3d&,
                            const boost::function< Eigen::Quaterniond( ) > ) >( &transformVector ),
                        _1, bodyFixedToInertialFrameFunction ) );
        transformationFunction = boost::bind(
                    static_cast< Eigen::Vector3d(&)(
                        const Eigen::Vector3d&,
                        const std::vector< boost::function< Eigen::Vector3d( const Eigen::Vector3d& ) > >& ) >( &transformVector ),
                    _1, rotationsList );
    }
    else
    {

        boost::function< Eigen::Quaterniond( ) > rotationFunction =
                boost::bind( &AerodynamicAngleCalculator::getRotationQuaternionBetweenFrames,
                             aerodynamicAngleCalculator, accelerationFrame, aerodynamicBaseFrame );
        transformationFunction = boost::bind(
                    static_cast< Eigen::Vector3d(&)(
                        const Eigen::Vector3d&,
                        const boost::function< Eigen::Quaterniond( ) > ) >( &transformVector ), _1,
                     rotationFunction );
    }

    return transformationFunction;
}

}

}



