#include <iostream>
#include <iomanip>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/tidallyLockedRotationalEphemeris.h"


namespace tudat
{

namespace ephemerides
{

Eigen::Vector6d getStateFromSelectedStateFunction(
        const double currentTime,
        const bool useFirstFunction,
        const std::function< Eigen::Vector6d( const double ) > stateFunction1,
        const std::function< Eigen::Vector6d( const double ) > stateFunction2 )
{
    return ( useFirstFunction ) ? ( stateFunction1( currentTime ) ) : ( stateFunction2( currentTime ) );
}


Eigen::Quaterniond TidallyLockedRotationalEphemeris::getRotationToBaseFrame( const double currentTime )
{
    Eigen::Vector6d relativeState = relativeStateFunction_( currentTime, isBodyInPropagation_ );

    Eigen::Matrix3d rotationToBaseFrame = reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrix(
                relativeState ).transpose( );
    rotationToBaseFrame.block( 0, 0, 3, 2 ) *= -1.0;

    return Eigen::Quaterniond( rotationToBaseFrame );
}

Eigen::Matrix3d TidallyLockedRotationalEphemeris::getDerivativeOfRotationToBaseFrame( const double currentTime )
{
    throw std::runtime_error( "Error, time-derivative of tidally-locked rotation matrix not yet implemented" );
}

}

}

