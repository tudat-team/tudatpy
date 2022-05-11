#ifndef TUDAT_AERODYNAMIC_ANGLE_ROTATIONAL_EPHEMERIS_H
#define TUDAT_AERODYNAMIC_ANGLE_ROTATIONAL_EPHEMERIS_H

#include "tudat/astro/aerodynamics/trimOrientation.h"
#include "tudat/astro/ephemerides/rotationalEphemeris.h"
#include "tudat/astro/reference_frames/aerodynamicAngleCalculator.h"

namespace tudat
{

namespace ephemerides
{

class AerodynamicAngleRotationalEphemeris: public ephemerides::RotationalEphemeris
{
public:

    //! Constructor.
    /*!
     * Constructor, sets frames between which rotation is determined.
     * \param baseFrameOrientation Base frame identifier.
     * \param targetFrameOrientation Target frame identifier.
     */
    AerodynamicAngleRotationalEphemeris(
            const std::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator,
            const std::string& baseFrameOrientation,
            const std::string& targetFrameOrientation,
            const std::function< Eigen::Vector3d( const double ) > aerodynamicAngleFunction = nullptr )
        : RotationalEphemeris( baseFrameOrientation, targetFrameOrientation ),
          aerodynamicAngleCalculator_( aerodynamicAngleCalculator ),
          aerodynamicAngleFunction_( aerodynamicAngleFunction ),
          currentTime_( TUDAT_NAN )
    {
        aerodynamicAngleCalculator->setAerodynamicAngleClosureIsIncomplete( );
//        aerodynamicAngleCalculator->setBodyFixedAngleInterface(
//                    std::make_shared< reference_frames::FromAeroEphemerisAerodynamicAngleInterface >( this ) );
    }

    //! Virtual destructor.
    /*!
     * Virtual destructor.
     */
    virtual ~AerodynamicAngleRotationalEphemeris( ) { }

    virtual Eigen::Quaterniond getRotationToBaseFrame(
            const double currentTime )
    {
        resetCurrentTime( currentTime );
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

    void updateBodyAngles( const double currentTime )
    {
        if( aerodynamicAngleFunction_ != nullptr )
        {
            currentBodyAngles_ = aerodynamicAngleFunction_( currentTime_ );
        }
        else
        {
            currentBodyAngles_.setZero( );
        }
    }
    virtual void resetCurrentTime( const double currentTime = TUDAT_NAN )
    {
        if( !( currentTime == currentTime_ ) )
        {
            currentTime_ = currentTime;
            if( currentTime_ == currentTime_ )
            {
                aerodynamicAngleCalculator_->update( currentTime, false );
                updateBodyAngles( currentTime );
                aerodynamicAngleCalculator_->update( currentTime, true );
            }
            else
            {
                aerodynamicAngleCalculator_->resetCurrentTime( TUDAT_NAN );
            }
        }
    }

    Eigen::Vector3d getBodyAngles( const double currentTime )
    {
        if( ( currentTime != currentTime_ ) )
        {
            updateBodyAngles( currentTime );
//            std::cout<<currentTime<<" "<<currentTime_<<std::endl;
//            throw std::runtime_error( "Error when getting body angles from AerodynamicAngleRotationalEphemeris, times are inconsistent" );
        }
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


protected:

    std::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator_;

    std::function< Eigen::Vector3d( const double ) > aerodynamicAngleFunction_;

    Eigen::Vector3d currentBodyAngles_;

    double currentTime_;

};


////! Function to make aerodynamic angle computation consistent with imposed body-fixed to inertial rotation.
///*!
// * Function to make aerodynamic angle computation consistent with imposed body-fixed to inertial rotation.
// * \param imposedRotationFromInertialToBodyFixedFrame Inertial to body-fixed frame rotation to which the
// * aerodynamicAngleCalculator object is to be made consistent
// * \param aerodynamicAngleCalculator Object from which the aerodynamic angles are computed.
// */
//void setAerodynamicDependentOrientationCalculatorClosure(
//        const std::function< Eigen::Quaterniond( const double ) > imposedRotationFromInertialToBodyFixedFrame,
//        std::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator );

////! Function to make aerodynamic angle computation consistent with existing rotational ephemeris
///*!
// * Function to make aerodynamic angle computation consistent with existing  rotational ephemeris
// * \param rotationalEphemeris Object computing the current orientation of the body. Aerodynamic angles are to be computed
// * from output given by this class.
// * \param aerodynamicAngleCalculator Object from which the aerodynamic angles are computed.
// */
//void setAerodynamicDependentOrientationCalculatorClosure(
//        std::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris,
//        std::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator );

void verifyAerodynamicDependentOrientationCalculatorClosure(
        std::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris,
        std::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator );


} // namespace ephemerides


namespace reference_frames
{


Eigen::Vector3d computeBodyFixedAeroAngles(
        const Eigen::Matrix3d& inertialToBodyFixedFrame,
        const Eigen::Matrix3d& trajectoryToInertialFrame );

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
private:

    std::shared_ptr< ephemerides::RotationalEphemeris > ephemeris_;

};

class FromAeroEphemerisAerodynamicAngleInterface: public BodyFixedAerodynamicAngleInterface
{
public:
    FromAeroEphemerisAerodynamicAngleInterface(
            const std::shared_ptr< ephemerides::AerodynamicAngleRotationalEphemeris > ephemeris ):
        BodyFixedAerodynamicAngleInterface( body_fixed_angles_from_aero_based_ephemeris ),
        ephemeris_( ephemeris ){ }

    virtual ~FromAeroEphemerisAerodynamicAngleInterface( ){ }

    Eigen::Vector3d getAngles( const double time,
                               const Eigen::Matrix3d& trajectoryToInertialFrame )
    {
        ephemeris_->resetCurrentTime( time );
        return ephemeris_->getBodyAngles( time );
    }

private:

    std::shared_ptr< ephemerides::AerodynamicAngleRotationalEphemeris > ephemeris_;

};

}

} // namespace tudat

#endif // TUDAT_AERODYNAMIC_ANGLE_ROTATIONAL_EPHEMERIS_H
