#include <iostream>
#include <iomanip>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/ephemerides/synchronousRotationalEphemeris.h"


namespace tudat
{

namespace ephemerides
{

double DirectLongitudeLibrationCalculator::getLibrationAngleWrtFullySynchronousRotationFromScaledLibration(
        const Eigen::Vector6d& relativeState,
        const double time,
        const double scaledLibrationAmplitude )
{
    double eccentricitySineEccentricAnomaly =
            ( ( relativeState.segment< 3 >( 0 ) ).dot( relativeState.segment< 3 >( 3 ) ) ) /
            (( relativeState.segment< 3 >( 0 ) ).cross( relativeState.segment< 3 >( 3 ) ) ).norm( );
    return scaledLibrationAmplitude * eccentricitySineEccentricAnomaly;
}

//! Calculate rotation quaternion from target frame to base frame.
Eigen::Matrix3d SynchronousRotationalEphemeris::getFullyLockedRotationToBaseFrame(
        const Eigen::Vector6d& relativeState,
        const double currentTime )
{
    // Get rotation to RSW frame
    Eigen::Matrix3d rotationToBaseFrame = reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrix(
                relativeState ).transpose( );
    rotationToBaseFrame.block( 0, 0, 3, 2 ) *= -1.0;

    return rotationToBaseFrame;
}

//! Calculate rotation quaternion from target frame to base frame.
Eigen::Matrix3d SynchronousRotationalEphemeris::getFullyLockedRotationToBaseFrame( const double currentTime )
{
    Eigen::Vector6d relativeState = relativeStateFunction_( currentTime, isBodyInPropagation_ );
    return getFullyLockedRotationToBaseFrame( relativeState, currentTime );
}

Eigen::Matrix3d SynchronousRotationalEphemeris::getLibrationRotation(
        const Eigen::Vector6d& relativeState,
        const double currentTime )
{
    if( isLibrationOn_ )
    {
        double librationAngle = longitudeLibrationCalculator_->getLibrationAngleWrtFullySynchronousRotation(
                    relativeState, currentTime );
//        std::cout<<"Lib: "<<librationAngle<<std::endl;
        return Eigen::AngleAxisd( -librationAngle, Eigen::Vector3d::UnitZ( ) ).toRotationMatrix( );

    }
    else
    {
//        std::cout<<"No lib: "<<std::endl;
        return Eigen::Matrix3d::Identity( );
    }
}

Eigen::Matrix3d SynchronousRotationalEphemeris::getLibrationRotation(
        const double currentTime )
{
    Eigen::Vector6d relativeState = relativeStateFunction_( currentTime, isBodyInPropagation_ );
    return getLibrationRotation( relativeState, currentTime );
}

//! Calculate rotation quaternion from target frame to base frame.
Eigen::Quaterniond SynchronousRotationalEphemeris::getRotationToBaseFrame( const double currentTime )
{
    Eigen::Vector6d relativeState = relativeStateFunction_( currentTime, isBodyInPropagation_ );
    return Eigen::Quaterniond(
                getFullyLockedRotationToBaseFrame( relativeState, currentTime ) *
                getLibrationRotation( relativeState, currentTime ) );


}

//! Function to calculate the derivative of the rotation matrix from target frame to base frame.
Eigen::Matrix3d SynchronousRotationalEphemeris::getDerivativeOfRotationToBaseFrame( const double currentTime )
{
    if( !warningPrinted_ )
    {
        std::cerr<<"Warning, time-derivative of synchronous rotation matrix not yet implemented (using zero matrix)"<<std::endl;
        warningPrinted_ = true;
    }
    return Eigen::Matrix3d::Zero( );
}

}

}

