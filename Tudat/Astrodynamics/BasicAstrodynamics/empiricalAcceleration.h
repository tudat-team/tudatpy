#ifndef TUDAT_EMPIRICALACCELERATION_H
#define TUDAT_EMPIRICALACCELERATION_H

#include <iomanip>

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

namespace tudat
{

namespace basic_astrodynamics
{

enum EmpiricalAccelerationFunctionalShapes
{
    constant_empirical = 0,
    sine_empirical = 1,
    cosine_empirical = 2
};


enum EmpiricalAccelerationComponents
{
    radial_empicial_acceleration_component = 0,
    along_track_empicial_acceleration_component = 1,
    across_track_empicial_acceleration_component = 2
};


//! Class for calculating an empirical acceleration, based on a once per orbit model, (Montenbruck and Gill, 2000)
/*!
 *  Class for calculating an empirical acceleration, based on a once per orbit model, (Montenbruck and Gill, 2000). Empirical
 *  accelerations are typically used in orbit determination to abosrb unmodelled small accelerations. The magnitude of the
 *  acceleration is determined by three vectors, a constant acceleration and two accelerations which scale with the sine and
 *  cosine fo the true anomaly of the  accelerated body, respectively. The acceleration is evaluated in the satellite's velocity
 *  frame and then transformed to the inertial frame.
 */
class EmpiricalAcceleration: public AccelerationModel< Eigen::Vector3d >
{
public:
    EmpiricalAcceleration(
            const Eigen::Vector3d constantAcceleration,
            const Eigen::Vector3d sineAcceleration,
            const Eigen::Vector3d cosineAcceleration,
            const boost::function< Eigen::Vector6d( ) > bodyStateFunction,
            const boost::function< double( ) > centralBodyGravitationalParameterFunction,
            const boost::function< Eigen::Vector6d( ) > centralBodyStateFunction =
            boost::lambda::constant( Eigen::Vector6d::Zero( ) ) ):
        bodyStateFunction_( bodyStateFunction ), centralBodyStateFunction_( centralBodyStateFunction ),
        centralBodyGravitationalParameterFunction_( centralBodyGravitationalParameterFunction )
    {
        Eigen::Matrix3d accelerationComponents;
        accelerationComponents.block( 0, 0, 3, 1 ) = constantAcceleration;
        accelerationComponents.block( 0, 1, 3, 1 ) = sineAcceleration;
        accelerationComponents.block( 0, 2, 3, 1 ) = cosineAcceleration;

        accelerationComponentsFunction_ = boost::lambda::constant( accelerationComponents );

        updateAccelerationComponents( 0.0 );

        areAccelerationComponentsTimeDependent_ = 0;
    }

    ~EmpiricalAcceleration( ){ }

    Eigen::Vector3d getAcceleration( )
    {
        return currentToInertialFrame_ * currentLocalAcclereration_ ;
    }

    void updateAccelerationComponents( const double currentTime = 0.0 )
    {
        Eigen::Matrix3d accelerationComponents = accelerationComponentsFunction_( currentTime );
        currentConstantAcceleration_ = accelerationComponents.block( 0, 0, 3, 1 );
        currentSineAcceleration_ = accelerationComponents.block( 0, 1, 3, 1 );
        currentCosineAcceleration_ = accelerationComponents.block( 0, 2, 3, 1 );
    }

    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            // Calulate current relative state of accelerated body
            currentState_ = bodyStateFunction_( ) - centralBodyStateFunction_( );

            // Calculate current body-fixed state of accelerated body.
            currentToInertialFrame_ = Eigen::Quaterniond(
                        reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrx( currentState_ ) );

            // Calculate current true anomaly of accelerated body.
            currentTrueAnomaly_ = orbital_element_conversions::convertCartesianToKeplerianElements(
                        currentState_, centralBodyGravitationalParameterFunction_( ) )( 5 );

            updateAccelerationComponents( currentTime );

            // Calculate acceleration
            currentLocalAcclereration_ = ( currentConstantAcceleration_ +
                                           currentSineAcceleration_ * std::sin( currentTrueAnomaly_ ) +
                                           currentCosineAcceleration_ * std::cos( currentTrueAnomaly_ ) );

            // Perform sanity check.
            if( currentLocalAcclereration_ != currentLocalAcclereration_ )
            {
                throw std::runtime_error( "Error when computing emprical acceleration, result is NaN" );
            }
            this->currentTime_ = currentTime;
        }
    }

    Eigen::Matrix3d getAccelerationComponents( )
    {
        if( areAccelerationComponentsTimeDependent_ )
        {
            std::cerr<<"Warning when getting time-invariant empirical acceleration components; original componets are time-varying"<<std::endl;
        }
        return accelerationComponentsFunction_( 0.0 );
    }

    Eigen::Vector3d getAccelerationComponent( const EmpiricalAccelerationFunctionalShapes functionalShape )
    {
        Eigen::Vector3d accelerationComponent;

        switch( functionalShape )
        {
        case constant_empirical:
            accelerationComponent = currentConstantAcceleration_;
            break;
        case sine_empirical:
            accelerationComponent = currentSineAcceleration_;
            break;
        case cosine_empirical:
            accelerationComponent = currentCosineAcceleration_;
            break;
        }

        return accelerationComponent;
    }

    void resetAccelerationComponents( const Eigen::Matrix3d& newAccelerationComponents )
    {
        if( areAccelerationComponentsTimeDependent_ )
        {
            std::cerr<<"Warning when resetting time-invariant empirical acceleration components; original componets are time-varying"<<std::endl;
        }

        areAccelerationComponentsTimeDependent_ = 0;
        accelerationComponentsFunction_ = boost::lambda::constant( newAccelerationComponents );

    }

    void resetAccelerationComponentsFunction(
            const boost::function< Eigen::Matrix3d( const double ) > & accelerationComponentsFunction )
    {
        areAccelerationComponentsTimeDependent_ = 1;
        accelerationComponentsFunction_ = accelerationComponentsFunction;

    }


    Eigen::Vector6d getCurrentState( ){
        return currentState_;
    }

    Eigen::Quaterniond getCurrentToInertialFrame( )
    {
        return currentToInertialFrame_;
    }

    Eigen::Vector3d getCurrentLocalAcceleration( )
    {
        return currentLocalAcclereration_;
    }

    double getCurrentTrueAnomaly( )
    {
        return currentTrueAnomaly_;
    }

    double getCurrentGravitationalParameter( )
    {
        return centralBodyGravitationalParameterFunction_( );
    }

private:


    boost::function< Eigen::Matrix3d( const double ) > accelerationComponentsFunction_;

    bool areAccelerationComponentsTimeDependent_;

    Eigen::Vector3d  currentConstantAcceleration_;

    Eigen::Vector3d  currentSineAcceleration_;

    Eigen::Vector3d  currentCosineAcceleration_;


    //! State function of the body that is undergoing the empirical acceleration.
    boost::function< Eigen::Vector6d( ) > bodyStateFunction_;

    //! State function of the body that is being orbited
    boost::function< Eigen::Vector6d( ) > centralBodyStateFunction_;

    //! Current state of the body that is undergoing the empirical acceleration, relative to central body, in global frame.
    Eigen::Vector6d currentState_;

    Eigen::Quaterniond currentToInertialFrame_;

    //! Current true anomaly of accelerated body in its orbit about the central body.
    double currentTrueAnomaly_;

    //! Current empirical acceleration in body-fixed frame.
    Eigen::Vector3d currentLocalAcclereration_;

    //! Function returning the gravitational parameter of the central body (ofr calculation of Kepler elements)
    boost::function< double( ) > centralBodyGravitationalParameterFunction_;
};

}

}

#endif // TUDAT_EMPIRICALACCELERATION_H
