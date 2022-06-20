#ifndef TUDAT_AERODYNAMIC_ANGLE_ROTATIONAL_EPHEMERIS_H
#define TUDAT_AERODYNAMIC_ANGLE_ROTATIONAL_EPHEMERIS_H

#include "tudat/astro/aerodynamics/trimOrientation.h"
#include "tudat/astro/ephemerides/rotationalEphemeris.h"
#include "tudat/astro/reference_frames/aerodynamicAngleCalculator.h"

namespace tudat
{

namespace ephemerides
{

/*!
 * Rotation model that calculates the body's orientation based on the current angles of
 * attack, sideslip angle and bank angle. These angles cann be defined as a custom function
 * and are 0 by default. These angles transform from body-fixed to trajectory frame
 * The transformation to inertial frame is handled through the AerodynamicAngleCalculator
 * class, which requires the current state of the vehicle.
 */
class AerodynamicAngleRotationalEphemeris: public ephemerides::RotationalEphemeris
{
public:

    AerodynamicAngleRotationalEphemeris(
            const std::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator,
            const std::string& baseFrameOrientation,
            const std::string& targetFrameOrientation,
            const std::function< Eigen::Vector3d( const double ) > aerodynamicAngleFunction = nullptr )
        : RotationalEphemeris( baseFrameOrientation, targetFrameOrientation ),
          aerodynamicAngleCalculator_( aerodynamicAngleCalculator ),
          aerodynamicAngleFunction_( aerodynamicAngleFunction ),
          currentTime_( TUDAT_NAN ),
          isBodyInPropagation_( false )
    {
        aerodynamicAngleCalculator->setAerodynamicAngleClosureIsIncomplete( );

    }

    //! Virtual destructor.
    /*!
     * Virtual destructor.
     */
    virtual ~AerodynamicAngleRotationalEphemeris( ) { }

    virtual Eigen::Quaterniond getRotationToBaseFrame(
            const double currentTime )
    {
        update( currentTime );
        return Eigen::Quaterniond( aerodynamicAngleCalculator_->getRotationMatrixBetweenFrames(
                                       reference_frames::body_frame, reference_frames::inertial_frame ) );
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

    void update( const double currentTime );

    void resetCurrentTime(  );

    Eigen::Vector3d getBodyAngles( const double currentTime )
    {
        update( currentTime );
        return currentBodyAngles_;
    }

    std::shared_ptr< reference_frames::AerodynamicAngleCalculator > getAerodynamicAngleCalculator( )
    {
        return aerodynamicAngleCalculator_;
    }

    void setAerodynamicAngleFunction( const std::function< Eigen::Vector3d( const double ) > aerodynamicAngleFunction )
    {
        aerodynamicAngleFunction_ = aerodynamicAngleFunction;
    }

    void addSideslipBankAngleFunctions( std::function< Eigen::Vector2d( const double ) > sideslipAndBankAngleFunction )
    {
        if( aerodynamicAngleFunction_ == nullptr )
        {
            aerodynamicAngleFunction_ = [=](const double time)
            {
                return ( Eigen::Vector3d( ) <<0.0, sideslipAndBankAngleFunction( time ) ).finished( );
            };
        }
        else
        {
            std::function< Eigen::Vector3d( const double ) > oldAerodynamicAngleFunction_ = aerodynamicAngleFunction_;

            aerodynamicAngleFunction_ = [=](const double time)
            {
                return ( Eigen::Vector3d( ) <<oldAerodynamicAngleFunction_( time )( 0 ), sideslipAndBankAngleFunction( time ) ).finished( );
            };
        }
    }

    void setIsBodyInPropagation( const bool isBodyInPropagation )
    {
        isBodyInPropagation_ = isBodyInPropagation;
    }

protected:

    void updateBodyAngles( );

    std::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator_;

    std::function< Eigen::Vector3d( const double ) > aerodynamicAngleFunction_;

    Eigen::Vector3d currentBodyAngles_;

    double currentTime_;

    bool isBodyInPropagation_;

};


} // namespace ephemerides


namespace reference_frames
{


Eigen::Vector3d computeBodyFixedAeroAngles(
        const Eigen::Matrix3d& inertialToBodyFixedFrame,
        const Eigen::Matrix3d& trajectoryToInertialFrame );

/*!
 * Class to communicate aerodynamic angles (attack, sidelip, bank) to the AerodynamicAngleCalculator
 * class from an arbitrary ephemeris object
 */
class FromGenericEphemerisAerodynamicAngleInterface: public BodyFixedAerodynamicAngleInterface
{
public:
    FromGenericEphemerisAerodynamicAngleInterface(
            const std::shared_ptr< ephemerides::RotationalEphemeris > ephemeris ):
        BodyFixedAerodynamicAngleInterface( body_fixed_angles_from_generic_ephemeris ),
        ephemeris_( ephemeris ){ }

    virtual ~FromGenericEphemerisAerodynamicAngleInterface( ){ }

    Eigen::Vector3d getAngles( const double time,
                               const Eigen::Matrix3d& trajectoryToInertialFrame );

    void resetCurrentTime( );
private:

    std::shared_ptr< ephemerides::RotationalEphemeris > ephemeris_;

};

/*!
 * Class to communicate aerodynamic angles (attack, sidelip, bank) to the AerodynamicAngleCalculator
 * class from an AerodynamicAngleRotationalEphemeris object
 */
class FromAeroEphemerisAerodynamicAngleInterface: public BodyFixedAerodynamicAngleInterface
{
public:
    FromAeroEphemerisAerodynamicAngleInterface(
            const std::shared_ptr< ephemerides::AerodynamicAngleRotationalEphemeris > ephemeris ):
        BodyFixedAerodynamicAngleInterface( body_fixed_angles_from_aero_based_ephemeris ),
        ephemeris_( ephemeris ){ }

    virtual ~FromAeroEphemerisAerodynamicAngleInterface( ){ }

    Eigen::Vector3d getAngles( const double time,
                               const Eigen::Matrix3d& trajectoryToInertialFrame );

    void resetCurrentTime( );

private:

    std::shared_ptr< ephemerides::AerodynamicAngleRotationalEphemeris > ephemeris_;

};

}

} // namespace tudat

#endif // TUDAT_AERODYNAMIC_ANGLE_ROTATIONAL_EPHEMERIS_H
