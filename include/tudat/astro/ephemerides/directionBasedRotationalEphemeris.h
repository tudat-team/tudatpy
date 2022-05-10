#ifndef TUDAT_DIRECTION_BASED_ROTATIONAL_EPHEMERIS_H
#define TUDAT_DIRECTION_BASED_ROTATIONAL_EPHEMERIS_H

#include "tudat/astro/ephemerides/rotationalEphemeris.h"

namespace tudat
{

namespace ephemerides
{

class DirectionBasedRotationalEphemeris: public ephemerides::RotationalEphemeris
{
public:

    //! Constructor.
    /*!
     * Constructor, sets frames between which rotation is determined.
     * \param baseFrameOrientation Base frame identifier.
     * \param targetFrameOrientation Target frame identifier.
     */
    DirectionBasedRotationalEphemeris(
            const std::function< Eigen::Vector3d( const double ) > inertialBodyAxisDirectionFunction,
            const Eigen::Vector3d associatedBodyFixedDirection,
            const std::string& baseFrameOrientation,
            const std::string& targetFrameOrientation )
        : RotationalEphemeris( baseFrameOrientation, targetFrameOrientation ),
          inertialBodyAxisDirectionFunction_( inertialBodyAxisDirectionFunction ),
          associatedBodyFixedDirection_( associatedBodyFixedDirection ),
          currentTime_( TUDAT_NAN )
    { }

    //! Virtual destructor.
    /*!
     * Virtual destructor.
     */
    virtual ~DirectionBasedRotationalEphemeris( ) { }

    virtual Eigen::Quaterniond getRotationToBaseFrame(
            const double currentTime )
    {
        resetCurrentTime( currentTime );
    }

    virtual Eigen::Quaterniond getRotationToTargetFrame(
            const double currentTime )
    {
        return getRotationToBaseFrame( currentTime ).inverse( );
    }

    virtual Eigen::Matrix3d getDerivativeOfRotationToBaseFrame(
            const double currentTime )
    {
        return Eigen::Matrix3d::Constant( TUDAT_NAN );
    }

    virtual Eigen::Matrix3d getDerivativeOfRotationToTargetFrame(
            const double currentTime )

    {
        return Eigen::Matrix3d::Constant( TUDAT_NAN );
    }


    virtual void resetCurrentTime( const double currentTime = TUDAT_NAN )
    {
        if( !( currentTime == currentTime_ ) )
        {
            currentTime_ = currentTime;
            if( currentTime_ == currentTime_ )
            {
                currentInertialDirection_ = inertialBodyAxisDirectionFunction_( currentTime_ );
            }
            else
            {
                currentInertialDirection_.setConstant( TUDAT_NAN );
            }
        }
    }

    Eigen::Vector3d getAssociatedBodyFixedDirection( )
    {
        return associatedBodyFixedDirection_;
    }

    Eigen::Vector3d getCurrentInertialDirection( const double currentTime )
    {
        if( currentTime != currentTime_ )
        {
            resetCurrentTime( currentTime );
        }
        return currentInertialDirection_;
    }


protected:

    std::function< Eigen::Vector3d( const double ) > inertialBodyAxisDirectionFunction_;

    Eigen::Vector3d associatedBodyFixedDirection_;

    Eigen::Vector3d currentInertialDirection_;

    double currentTime_;

};


} // namespace ephemerides

} // namespace tudat

#endif // TUDAT_DIRECTION_BASED_ROTATIONAL_EPHEMERIS_H
