/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Montebruck O, Gill E. Satellite Orbits, Springer, 2000.
 *
 */

#ifndef TUDAT_EMPIRICALACCELERATION_H
#define TUDAT_EMPIRICALACCELERATION_H


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

//! Enum defining shape of empirical accelerations
enum EmpiricalAccelerationFunctionalShapes
{
    constant_empirical = 0,
    sine_empirical = 1,
    cosine_empirical = 2
};

inline std::string getEmpiricalAccelerationFunctionalShapeString(
        const EmpiricalAccelerationFunctionalShapes functionalShape )
{
    std::string parameterDescription = "";
    switch( functionalShape )
    {
    case constant_empirical:
        parameterDescription = "constant";
        break;
    case sine_empirical:
        parameterDescription = "sine of true anomaly";
        break;
    case cosine_empirical:
        parameterDescription = "cosine of true anomaly";
        break;
    default:
        throw std::runtime_error( "Error when getting functional shape string, type not recognized" );
    }
    return parameterDescription;
}


//! Enum defining component of empirical accelerations
enum EmpiricalAccelerationComponents
{
    radial_empirical_acceleration_component = 0,
    along_track_empirical_acceleration_component = 1,
    across_track_empirical_acceleration_component = 2
};


//! Class for calculating an empirical acceleration, based on a once per orbit model, (Montenbruck and Gill, 2000)
/*!
 *  Class for calculating an empirical acceleration, based on a once per orbit model, (Montenbruck and Gill, 2000). Empirical
 *  accelerations are typically used in orbit determination to abosrb unmodelled small accelerations. The magnitude of the
 *  acceleration is determined by three vectors, a constant acceleration and two accelerations which scale with the sine and
 *  cosine fo the true anomaly of the  accelerated body, respectively. The acceleration is evaluated in the satellite's RSW
 *  frame and then transformed to the inertial frame.
 */
class EmpiricalAcceleration: public AccelerationModel< Eigen::Vector3d >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param constantAcceleration Constant empirical acceleration in RSW frame.
     * \param sineAcceleration Empirical acceleration that is scaled by sine of true anomaly in RSW frame.
     * \param cosineAcceleration Empirical acceleration that is scaled by cosine of true anomaly in RSW frame.
     * \param bodyStateFunction Functon that returns the state of the body undergoing acceleration
     * \param centralBodyGravitationalParameterFunction Functon that returns the gravitational parameter of central body.
     * \param centralBodyStateFunction Functon that returns the state of central body.
     */
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
        // Set empirical acceleration components
        Eigen::Matrix3d accelerationComponents;
        accelerationComponents.block( 0, 0, 3, 1 ) = constantAcceleration;
        accelerationComponents.block( 0, 1, 3, 1 ) = sineAcceleration;
        accelerationComponents.block( 0, 2, 3, 1 ) = cosineAcceleration;
        accelerationComponentsFunction_ = boost::lambda::constant( accelerationComponents );

        updateAccelerationComponents( 0.0 );
        areAccelerationComponentsTimeDependent_ = 0;
    }

    //! Destructor
    ~EmpiricalAcceleration( ){ }

    //! Function to retrieve current empirical acceleration
    /*!
     * Function to retrieve current empirical acceleration
     * \return Current empirical acceleration, expressed in propagation frame
     */
    Eigen::Vector3d getAcceleration( )
    {
        return currentRotationFromRswToInertialFrame_ * currentLocalAcclereration_ ;
    }

    //! Function to update constituent elements of empirical acceleration to current time
    /*!
     * Function to update constituent elements of empirical acceleration to current time
     * \param currentTime Time to which empirical acceleration elements are to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            // Calulate current relative state of accelerated body
            currentState_ = bodyStateFunction_( ) - centralBodyStateFunction_( );

            // Calculate current body-fixed state of accelerated body.
            currentRotationFromRswToInertialFrame_ = Eigen::Quaterniond(
                        reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrx( currentState_ ) ).inverse( );

            // Calculate current true anomaly of accelerated body.
            currentTrueAnomaly_ = orbital_element_conversions::convertCartesianToKeplerianElements(
                        currentState_, centralBodyGravitationalParameterFunction_( ) )( 5 );

            // Calculate acceleration
            updateAccelerationComponents( currentTime );
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

    //! Function to retrieve empirical acceleration components in RSW frame at a given time.
    /*!
     *  Function to retrieve empirical acceleration components in RSW frame at a given time. If components are not time-dependent,
     *  input timemay be NaN
     *  \param currentTime Time at which acceleration components are to be retrieved
     *  \return Empirical acceleration components in RSW frame. Constant, sine and cosine terms are given in first, second and
     *  third column of return matrix, respectively.
     */
    Eigen::Matrix3d getAccelerationComponents( const double currentTime = TUDAT_NAN )
    {
        if( areAccelerationComponentsTimeDependent_ && !( currentTime == currentTime ) )
        {
            throw std::runtime_error(
                        "Error when getting time-invariant empirical acceleration components; original componets are time-varying and input time is NaN" );
        }
        return accelerationComponentsFunction_( 0.0 );
    }

    //! Function to retrieve empirical acceleration component in RSW frame, as set by last call to updateMembers function
    /*!
     *  Function to retrieve empirical acceleration component in RSW frame, as set by last call to updateMembers function
     *  \param functionalShape Component of acceleration that is to be retrieved
     *  \return Current empirical acceleration component in RSW frame.
     */
    Eigen::Vector3d getCurrentAccelerationComponent( const EmpiricalAccelerationFunctionalShapes functionalShape )
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

    //! Function to reset time-independent empirical acceleration components
    /*!
     *  Function to reset time-independent empirical acceleration components
     *  \param newAccelerationComponents Time-independent empirical acceleration components. Constant, sine and cosine terms are
     *  given in first, second and third column of return matrix, respectively.
     */
    void resetAccelerationComponents( const Eigen::Matrix3d& newAccelerationComponents )
    {
        if( areAccelerationComponentsTimeDependent_ )
        {
            std::cerr << "Warning when resetting time-invariant empirical acceleration components; original componets are time-varying" << std::endl;
        }

        areAccelerationComponentsTimeDependent_ = 0;
        accelerationComponentsFunction_ = boost::lambda::constant( newAccelerationComponents );

    }

    //! Function to reset time-dependent empirical acceleration components
    /*!
     *  Function to reset time-dependent empirical acceleration components
     *  \param accelerationComponentsFunction Function returning empirical acceleration components as a function of time.
     *  Constant, sine and cosine terms are given in first, second and third column of return matrix, respectively.
     */
    void resetAccelerationComponentsFunction(
            const boost::function< Eigen::Matrix3d( const double ) > & accelerationComponentsFunction )
    {
        areAccelerationComponentsTimeDependent_ = 1;
        accelerationComponentsFunction_ = accelerationComponentsFunction;

    }

    //! Function to retrieve current state of the body that is undergoing the empirical acceleration, relative to central body
    /*!
     *  Function to retrieve Current state of the body that is undergoing the empirical acceleration, relative to central body,
     *  in propagation frame.
     *  \return Current state of the body that is undergoing the empirical acceleration, relative to central body, in global frame.
     */
    Eigen::Vector6d getCurrentState( )
    {
        return currentState_;
    }

    //! Function to retrieve quaternion defining the rotation from RSW to inertial frame.
    /*!
     * Function to retrieve quaternion defining the rotation from RSW to inertial frame.
     * \return Quaternion defining the rotation from RSW to inertial frame.
     */
    Eigen::Quaterniond getCurrentToInertialFrame( )
    {
        return currentRotationFromRswToInertialFrame_;
    }

    //! Function to retrieve current empirical acceleration in RSW frame.
    /*!
     * Function to retrieve current empirical acceleration in RSW frame.
     * \return Current empirical acceleration in RSW frame.
     */
    Eigen::Vector3d getCurrentLocalAcceleration( )
    {
        return currentLocalAcclereration_;
    }

    //! Function to retrieve current true anomaly of accelerated body in its orbit about the central body.
    /*!
     * Function to retrieve current true anomaly of accelerated body in its orbit about the central body.
     * \return Current true anomaly of accelerated body in its orbit about the central body.
     */
    double getCurrentTrueAnomaly( )
    {
        return currentTrueAnomaly_;
    }

    //! Function to retrieve gravitational parameter of the central body
    /*!
     * Function to retrieve gravitational parameter of the central body
     * \return Gravitational parameter of the central body
     */
    double getCurrentGravitationalParameter( )
    {
        return centralBodyGravitationalParameterFunction_( );
    }

private:

    //! Function to update components of empirical accelerations to current time
    /*!
     * Function to update components of empirical accelerations to current time
     * \param currentTime Time to which empirical accelerations are to be updated.
     */
    void updateAccelerationComponents( const double currentTime = 0.0 )
    {
        Eigen::Matrix3d accelerationComponents = accelerationComponentsFunction_( currentTime );

        currentConstantAcceleration_ = accelerationComponents.block( 0, 0, 3, 1 );
        currentSineAcceleration_ = accelerationComponents.block( 0, 1, 3, 1 );
        currentCosineAcceleration_ = accelerationComponents.block( 0, 2, 3, 1 );
    }

    //! Function returning empirical acceleration components as a function of time.
    /*!
     *  Function returning empirical acceleration components as a function of time.
     *  Constant, sine and cosine terms are given in first, second and third column of return matrix, respectively.
     */
    boost::function< Eigen::Matrix3d( const double ) > accelerationComponentsFunction_;

    //! Boolean denoting whether empirical accelerations are time-dependent.
    bool areAccelerationComponentsTimeDependent_;


    //! Value of constant empirical acceleration, in RSW frame, as computed by last call to updateMembers function.
    Eigen::Vector3d currentConstantAcceleration_;

    //! Value of sine empirical acceleration, in RSW frame, as computed by last call to updateMembers function.
    Eigen::Vector3d currentSineAcceleration_;

    //! Value of cosine empirical acceleration, in RSW frame, as computed by last call to updateMembers function.
    Eigen::Vector3d currentCosineAcceleration_;


    //! State function of the body that is undergoing the empirical acceleration.
    boost::function< Eigen::Vector6d( ) > bodyStateFunction_;

    //! State function of the body that is being orbited
    boost::function< Eigen::Vector6d( ) > centralBodyStateFunction_;

    //! Current state of the body that is undergoing the empirical acceleration, relative to central body, in global frame.
    Eigen::Vector6d currentState_;

    //! Quaternion defining the rotation from RSW to inertial frame.
    Eigen::Quaterniond currentRotationFromRswToInertialFrame_;

    //! Current true anomaly of accelerated body in its orbit about the central body.
    double currentTrueAnomaly_;

    //! Current empirical acceleration in RSW frame.
    Eigen::Vector3d currentLocalAcclereration_;

    //! Function returning the gravitational parameter of the central body (for calculation of Kepler elements)
    boost::function< double( ) > centralBodyGravitationalParameterFunction_;
};

}

}

#endif // TUDAT_EMPIRICALACCELERATION_H
