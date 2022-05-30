#ifndef TUDAT_DIRECTION_BASED_ROTATIONAL_EPHEMERIS_H
#define TUDAT_DIRECTION_BASED_ROTATIONAL_EPHEMERIS_H

#include "tudat/astro/ephemerides/rotationalEphemeris.h"

namespace tudat
{

namespace ephemerides
{

enum SatelliteBasedFrames
{
    unspecified_satellite_based_frame = -1,
    inertial_satellite_based_frame = 0,
    tnw_satellite_based_frame = 1
};

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
          freeRotationAngleFunction_( freeRotationAngleFunction ),
          currentTime_( TUDAT_NAN ),
          currentEulerAnglesTime_( TUDAT_NAN )
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
            const double currentTime );

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

    virtual void update( const double currentTime );

    virtual void resetCurrentTime( );

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

    void setInertialBodyAxisDirectionFunction(
            const std::function< Eigen::Vector3d( const double ) >& inertialBodyAxisDirectionFunction )
    {
        inertialBodyAxisDirectionFunction_ = inertialBodyAxisDirectionFunction;
    }

    void setFreeRotationAngleFunction( const std::function< double( const double ) >& freeRotationAngleFunction )
    {
        freeRotationAngleFunction_= freeRotationAngleFunction;
    }



protected:

    Eigen::Vector3d getEulerAngles( const double currentTime )
    {
        if( currentTime_ != currentTime )
        {
            update( currentTime );
            calculateEulerAngles( );
        }
        else if( currentEulerAnglesTime_ != currentTime )
        {
            calculateEulerAngles( );
        }
        return eulerAngles_;
    }

    void calculateEulerAngles( );

    std::function< Eigen::Vector3d( const double ) > inertialBodyAxisDirectionFunction_;

    Eigen::Vector3d associatedBodyFixedDirection_;

    Eigen::Vector3d currentInertialDirection_;

    std::function< double( const double ) > freeRotationAngleFunction_;

    double currentTime_;

    double currentEulerAnglesTime_;

    Eigen::Vector3d eulerAngles_;

};


} // namespace ephemerides

} // namespace tudat

#endif // TUDAT_DIRECTION_BASED_ROTATIONAL_EPHEMERIS_H
