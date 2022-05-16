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
            const std::string& targetFrameOrientation,
            const std::function< double( const double ) > freeRotationAngleFunction = nullptr )
        : RotationalEphemeris( baseFrameOrientation, targetFrameOrientation ),
          inertialBodyAxisDirectionFunction_( inertialBodyAxisDirectionFunction ),
          associatedBodyFixedDirection_( associatedBodyFixedDirection ),
          currentTime_( TUDAT_NAN ),
          freeRotationAngleFunction_( freeRotationAngleFunction )

    {
        if( associatedBodyFixedDirection != Eigen::Vector3d::UnitX( ) )
        {
            throw std::runtime_error( "Error in DirectionBasedRotationalEphemeris, only x-axis body-fixed direction is currenly supported" );
        }
    }

    //! Virtual destructor.
    /*!
     * Virtual destructor.
     */
    virtual ~DirectionBasedRotationalEphemeris( ) { }

    virtual Eigen::Quaterniond getRotationToBaseFrame(
            const double currentTime )
    {
        update( currentTime );
        return Eigen::Quaterniond( );
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

    virtual void update( const double currentTime )
    {
        if( !( currentTime == currentTime_ ) )
        {
            currentTime_ = currentTime;
            if( currentTime_ == currentTime_ )
            {
                currentInertialDirection_ = inertialBodyAxisDirectionFunction_( currentTime_ ).normalized( );
//                calculateEulerAngles( );
            }
            else
            {
                currentInertialDirection_.setConstant( TUDAT_NAN );
                eulerAngles_.setConstant( TUDAT_NAN );
            }
        }
    }


    virtual void resetCurrentTime(  )
    {
        update( TUDAT_NAN );
    }

    Eigen::Vector3d getAssociatedBodyFixedDirection( )
    {
        return associatedBodyFixedDirection_;
    }

    Eigen::Vector3d getCurrentInertialDirection( const double currentTime )
    {
        if( currentTime != currentTime_ )
        {
            update( currentTime );
        }
        return currentInertialDirection_;
    }


protected:

    void calculateEulerAngles( )
    {
        eulerAngles_( 0 ) = std::atan2( currentInertialDirection_.y( ), currentInertialDirection_.x( ) );
        eulerAngles_( 1 ) = std::atan2( currentInertialDirection_.z( ), currentInertialDirection_.segment( 0, 2 ).norm( ) );
        if( freeRotationAngleFunction_ != nullptr )
        {
            eulerAngles_( 2 ) = 0.0;
        }
        else
        {
            eulerAngles_( 2 ) = freeRotationAngleFunction_( currentTime_ );
        }
    }

    std::function< Eigen::Vector3d( const double ) > inertialBodyAxisDirectionFunction_;

    Eigen::Vector3d associatedBodyFixedDirection_;

    Eigen::Vector3d currentInertialDirection_;

    double currentTime_;

    std::function< double( const double ) > freeRotationAngleFunction_;

    Eigen::Vector3d eulerAngles_;

};


} // namespace ephemerides

} // namespace tudat

#endif // TUDAT_DIRECTION_BASED_ROTATIONAL_EPHEMERIS_H
