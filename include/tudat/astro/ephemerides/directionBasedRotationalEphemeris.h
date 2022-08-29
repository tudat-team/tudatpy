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
    tnw_satellite_based_frame = 1,
    rsw_satellite_based_frame = 2
};

/*!
 *  Base class for computation of inertial direction algorithm to be used in DirectionBasedRotationalEphemeris
 *  rotation model. The rotationMatrixToPropagationFrame input allows the frame in which the direction is defined
 *  to be non-inertial (but RSW, TNW, etc.)
 */
class InertialBodyFixedDirectionCalculator
{
public:
    InertialBodyFixedDirectionCalculator(
            const std::function< Eigen::Matrix3d( const double ) > rotationMatrixToPropagationFrame = nullptr ):
        rotationMatrixToPropagationFrame_( rotationMatrixToPropagationFrame ){ }

    virtual ~InertialBodyFixedDirectionCalculator( ){ }

    virtual Eigen::Vector3d getDirection( const double time ) = 0;

    virtual void resetCurrentTime( ) = 0;

    virtual void update( const double time ) = 0;

    Eigen::Vector3d getInertialDirection( const double time )
    {
        if( rotationMatrixToPropagationFrame_ != nullptr )
        {
            return rotationMatrixToPropagationFrame_( time ) * getDirection( time );
        }
        else
        {
            return getDirection( time );
        }
    }

protected:

    std::function< Eigen::Matrix3d( const double ) > rotationMatrixToPropagationFrame_;




};


/*!
 *  Custom class for computation of inertial direction algorithm to be used in DirectionBasedRotationalEphemeris
 *  rotation model. This class takes an arbitrary direction function as input.
 *  The rotationMatrixToPropagationFrame input allows the frame in which the direction is defined
 *  to be non-inertial (but RSW, TNW, etc.)
 */
class CustomBodyFixedDirectionCalculator: public InertialBodyFixedDirectionCalculator
{
public:
    CustomBodyFixedDirectionCalculator(
            const std::function< Eigen::Vector3d( const double ) > inertialBodyAxisDirectionFunction,
            const std::function< Eigen::Matrix3d( const double ) > rotationMatrixToPropagationFrame = nullptr ):
        InertialBodyFixedDirectionCalculator( rotationMatrixToPropagationFrame ),
        inertialBodyAxisDirectionFunction_( inertialBodyAxisDirectionFunction ),
        currentTime_( TUDAT_NAN )
        { }

    ~CustomBodyFixedDirectionCalculator( ){ }

    Eigen::Vector3d getDirection( const double time )
    {
        update( time );
        return currentBodyAxisDirection_;
    }

    virtual void resetCurrentTime( )
    {
        currentBodyAxisDirection_.setConstant( TUDAT_NAN );
        currentTime_ = TUDAT_NAN;
        try { inertialBodyAxisDirectionFunction_( TUDAT_NAN ); }
        catch( ... ) { }
    }

    virtual void update( const double time );

protected:

    std::function< Eigen::Vector3d( const double ) > inertialBodyAxisDirectionFunction_;

    double currentTime_;

    Eigen::Vector3d currentBodyAxisDirection_;

};

/*!
 *  Class for computation of inertial direction algorithm to be used in DirectionBasedRotationalEphemeris
 *  rotation model, with the direction based on current position/velocity.
 *  This class takes an arbitrary direction function as input.
 */
class StateBasedBodyFixedDirectionCalculator: public InertialBodyFixedDirectionCalculator
{
public:
    StateBasedBodyFixedDirectionCalculator(
            const std::string& centralBody,
            const bool isColinearWithVelocity,
            const bool directionIsOppositeToVector,
            const std::function< void( Eigen::Vector6d& ) > relativeStateFunction ):
        InertialBodyFixedDirectionCalculator( nullptr ),
        centralBody_( centralBody ),
        isColinearWithVelocity_( isColinearWithVelocity ),
        directionIsOppositeToVector_( directionIsOppositeToVector ),
        relativeStateFunction_( relativeStateFunction ),
        currentTime_( TUDAT_NAN ){ }

    ~StateBasedBodyFixedDirectionCalculator( ){ }

    Eigen::Vector3d getDirection( const double time )
    {
        update( time );
        return currentDirection;
    }

    virtual void resetCurrentTime( )
    {
        currentRelativeState_.setConstant( TUDAT_NAN );
        currentTime_ = TUDAT_NAN;
    }

    virtual void update( const double time )
    {
        if( currentTime_ != time )
        {
            relativeStateFunction_( currentRelativeState_ );
            if( isColinearWithVelocity_ )
            {
                currentDirection = ( directionIsOppositeToVector_ ? -1.0 : 1.0 ) * currentRelativeState_.segment( 3, 3 );
            }
            else
            {
                currentDirection = ( directionIsOppositeToVector_ ? -1.0 : 1.0 ) * currentRelativeState_.segment( 0, 3 );
            }
            currentTime_ = time;
        }
    }

protected:

    std::string centralBody_;

    bool isColinearWithVelocity_;

    bool directionIsOppositeToVector_;

    std::function< void( Eigen::Vector6d& ) > relativeStateFunction_;

    Eigen::Vector6d currentRelativeState_;

    double currentTime_;

    Eigen::Vector3d currentDirection;


};

/*!
 *  Class to define a vehicle rotation model based on aligning a body-fixed axis with a given
 *  (time-variable) inertial axis. This aligning of vectors, plus the definition of a
 *  free rotation angle, defines the body-fixed <-> inertial frame.
 */
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
            const std::shared_ptr< InertialBodyFixedDirectionCalculator > directionCalculator,
            const Eigen::Vector3d& associatedBodyFixedDirection,
            const std::string& baseFrameOrientation,
            const std::string& targetFrameOrientation,
            const std::function< double( const double ) > freeRotationAngleFunction = nullptr )
        : RotationalEphemeris( baseFrameOrientation, targetFrameOrientation ),
          directionCalculator_( directionCalculator ),
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
        update( currentTime );
        return currentInertialDirection_;
    }

    void setInertialBodyAxisDirectionCalculator(
            const std::shared_ptr< InertialBodyFixedDirectionCalculator > directionCalculator )
    {
        directionCalculator_ = directionCalculator;
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

    std::shared_ptr< InertialBodyFixedDirectionCalculator > directionCalculator_;

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
