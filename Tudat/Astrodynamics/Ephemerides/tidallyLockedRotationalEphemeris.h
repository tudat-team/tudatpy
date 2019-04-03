#ifndef TIDALLYLOCKEDROTATIONALEPHEMERIS_H
#define TIDALLYLOCKEDROTATIONALEPHEMERIS_H

#include <functional>

#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

namespace tudat
{

namespace ephemerides
{

Eigen::Vector6d getStateFromSelectedStateFunction(
        const double currentTime,
        const bool useFirstFunction,
        const std::function< Eigen::Vector6d( const double ) > stateFunction1,
        const std::function< Eigen::Vector6d( const double ) > stateFunction2 );

class TidallyLockedRotationalEphemeris: public RotationalEphemeris
{
public:
    TidallyLockedRotationalEphemeris(
            const std::function< Eigen::Vector6d( const double, bool ) > relativeStateFunction,
            const std::string& centralBodyName,
            const std::string& baseFrameOrientation,
            const std::string& targetFrameOrientation ):
        RotationalEphemeris( baseFrameOrientation, targetFrameOrientation ),
        relativeStateFunction_( relativeStateFunction ),
        isBodyInPropagation_( 0 ),
        centralBodyName_( centralBodyName )
    { }

    ~TidallyLockedRotationalEphemeris( ){ }

    Eigen::Quaterniond getRotationToBaseFrame( const double currentTime );

    Eigen::Quaterniond getRotationToTargetFrame(
            const double currentTime )
    {
        return getRotationToBaseFrame( currentTime ).inverse( );
    }

    Eigen::Matrix3d getDerivativeOfRotationToBaseFrame( const double currentTime );

    Eigen::Matrix3d getDerivativeOfRotationToTargetFrame(
            const double currentTime )
    {
        return getDerivativeOfRotationToBaseFrame( currentTime ).transpose( );
    }

    void setIsBodyInPropagation( const bool isBodyInPropagation )
    {
        isBodyInPropagation_ = isBodyInPropagation;
    }

private:


    const std::function< Eigen::Vector6d( const double, bool ) > relativeStateFunction_;

    bool isBodyInPropagation_;

    std::string centralBodyName_;
};

}

}


#endif // TIDALLYLOCKEDROTATIONALEPHEMERIS_H
