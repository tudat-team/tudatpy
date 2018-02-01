/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_ROTATIONALMOTIONSTATEDERIVATIVE_H
#define TUDAT_ROTATIONALMOTIONSTATEDERIVATIVE_H


#include <vector>
#include <map>
#include <string>

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModel.h"

#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"

namespace tudat
{

namespace propagators
{

//! Function to evaluated the classical rotational equations of motion (Euler equations)
/*!
 * Function to evaluated the classical rotational equations of motion (Euler equations). The function returns the time-derivative
 * of a body's angular velocity vector, expressed in its body-fixed frame
 * \param inertiaTensor Inertia tensor of body, expressed in its body-fixed frame
 * \param totalTorque Total torque acting on body, expressed in its body-fixed frame
 * \param angularVelocityVector Current angular velocity vector of body, expressed in its body-fixed frame
 * \param inertiaTensorTimeDerivative Time derivative of inertiaTensor (default zero)
 * \return Time-derivative of a body's angular velocity vector, expressed in its body-fixed frame
 */
Eigen::Vector3d evaluateRotationalEquationsOfMotion(
        const Eigen::Matrix3d& inertiaTensor, const Eigen::Vector3d& totalTorque,
        const Eigen::Vector3d& angularVelocityVector,
        const Eigen::Matrix3d& inertiaTensorTimeDerivative = Eigen::Matrix3d::Zero( ) );

//! Function to obtain the matrix by which a quaternion vector is to be pre-multiplied to obtain this quaternion's time-derivative
/*!
 * Function to obtain the matrix by which a quaternion vector (representing body-fixed to inertial frame rotation) is to be
 * pre-multiplied to obtain this quaternion's time-derivative
 * \param angularVelocityVectorInBodyFixedFrame  Current angular velocity vector of body, expressed in its body-fixed frame
 * \return Matrix by which a quaternion vector (representing body-fixed to inertial frame rotation) is to be
 * pre-multiplied to obtain this quaternion's time-derivative
 */
Eigen::Matrix4d getQuaterionToQuaternionRateMatrix( const Eigen::Vector3d& angularVelocityVectorInBodyFixedFrame );

//! Function to obtain the time derivative of a quaternion (in vector representation) of body-fixed to inertial frame
/*!
 * Function to obtain the time derivative of a quaternion (in vector representation) of body-fixed to inertial frame
 * \param currentQuaternionToBaseFrame Quaternion (in vector representation) that defined the rotation from body-fixed to inertial
 * frame.
 * \param angularVelocityVectorInBodyFixedFrame Current angular velocity vector of body, expressed in its body-fixed frame
 * \return Time derivative of a quaternion (in vector representation) of body-fixed to inertial frame
 */
Eigen::Vector4d calculateQuaternionDerivative(
        const Eigen::Vector4d& currentQuaternionToBaseFrame, const Eigen::Vector3d& angularVelocityVectorInBodyFixedFrame );

//! Class for computing the state derivative for rotational dynamics of N bodies.
/*!
 *  Class for computing the state derivative for rotational dynamics of N bodies., using quaternion from body-fixed to inertial
 *  frame (in quaternion format) and angular velocity-vector of body expressed in body-fixed frame as the rotational state of a
 *  single body
 */
template< typename StateScalarType = double, typename TimeType = double >
class RotationalMotionStateDerivative: public propagators::SingleStateTypeDerivative< StateScalarType, TimeType >
{
public:

    using propagators::SingleStateTypeDerivative< StateScalarType, TimeType >::calculateSystemStateDerivative;

    //! Constructor
    /*!
     * Constructor
     * \param torqueModelsPerBody List of torque models (first map key body undergoing acceleration, second map key body exerting
     * acceleration)
     * \param bodiesToPropagate List of names of bodies for which rotational state is to be propagated
     * \param bodyInertiaTensorFunctions List of functions returning inertia tensors of bodiesToPropagate (in same order)
     * \param bodyInertiaTensorTimeDerivativeFunctions List of functions returning time derivatives of inertia tensors of
     *  bodiesToPropagate (in same order). Default empty, denoting time-invariant inertia tensors.
     */
    RotationalMotionStateDerivative(
            const basic_astrodynamics::TorqueModelMap& torqueModelsPerBody,
            const std::vector< std::string >& bodiesToPropagate,
            std::vector< boost::function< Eigen::Matrix3d( ) > > bodyInertiaTensorFunctions,
            std::vector< boost::function< Eigen::Matrix3d( ) > > bodyInertiaTensorTimeDerivativeFunctions =
            std::vector< boost::function< Eigen::Matrix3d( ) > >( ) ):
        propagators::SingleStateTypeDerivative< StateScalarType, TimeType >(
            propagators::rotational_state ),
        torqueModelsPerBody_( torqueModelsPerBody ),
        bodiesToPropagate_( bodiesToPropagate ),
        bodyInertiaTensorFunctions_( bodyInertiaTensorFunctions ),
        bodyInertiaTensorTimeDerivativeFunctions_( bodyInertiaTensorTimeDerivativeFunctions )
    {
        // Check input consistency
        if( bodiesToPropagate_.size( ) != bodyInertiaTensorFunctions_.size( ) )
        {
            throw std::runtime_error(
                        "Error when making rotational state derivative model, inertia tensor list is of incompatible size" );
        }

        if( bodyInertiaTensorTimeDerivativeFunctions_.size( ) == 0 )
        {
            for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
            {
                bodyInertiaTensorTimeDerivativeFunctions_.push_back( boost::lambda::constant( Eigen::Matrix3d::Zero( ) ) );
            }
        }
        else if( bodiesToPropagate_.size( ) != bodyInertiaTensorTimeDerivativeFunctions_.size( ) )
        {
            throw std::runtime_error(
                        "Error when making rotational state derivative model, inertia tensor time derivative list is of incompatible size" );
        }
    }

    //! Destructor
    ~RotationalMotionStateDerivative( ){ }

    //! Function to update the state derivative model to the current time.
    /*!
     * Function to update the state derivative model (i.e. torque models) to the
     * current time. Note that this function only updates the state derivative model itself, the
     * environment models must be updated before calling this function.
     * \param currentTime Time at which state derivative is to be calculated
     */
    void updateStateDerivativeModel( const TimeType currentTime )
    {
        for( torqueModelMapIterator = torqueModelsPerBody_.begin( );
             torqueModelMapIterator != torqueModelsPerBody_.end( ); torqueModelMapIterator++ )
        {
            for( innerTorqueIterator  = torqueModelMapIterator->second.begin( ); innerTorqueIterator !=
                 torqueModelMapIterator->second.end( ); innerTorqueIterator++ )
            {
                for( unsigned int j = 0; j < innerTorqueIterator->second.size( ); j++ )
                {
                    innerTorqueIterator->second[ j ]->updateMembers( currentTime );
                }
            }
        }
    }


    //! Function to convert the propagator-specific form of the state to the conventional form in the global frame.
    /*!
     * Function to convert the propagator-specific form of the state to the conventional form in the
     * global frame.  The conventional form for translational dynamics this is the quaternion from body-fixed to propagation frame
     * (in vector form) and the body's angular velocity vector in body-fixed frame. For this class the 'conventional' and
     * 'propagator-specific' forms are equal.
     * \param internalSolution State in propagator-specific form (i.e. form that is used in
     * numerical integration).
     * \param time Current time at which the state is valid.
     * \param currentRotationalState State (internalSolution), converted to the rotational state in 'conventional' format (equal
     * in this case).
     */
    void convertCurrentStateToGlobalRepresentation(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& internalSolution, const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentRotationalState )
    {
        currentRotationalState = internalSolution;
        for( unsigned int i = 0; i < bodiesToPropagate_.size( ); i++ )
        {
            currentRotationalState.block( 7 * i, 0, 4, 1 ) =
                    internalSolution.block( 7 * i, 0, 4, 1 ).normalized( );
        }
    }

    //! Function to clear reference/cached values of rotational state derivative model
    /*!
     *  Function to clear reference/cached values of rotational state derivative model. For each derived class, this
     *  entails resetting the current time in the acceleration models to NaN.
     */
    void clearStateDerivativeModel( )
    {
        for( torqueModelMapIterator = torqueModelsPerBody_.begin( );
             torqueModelMapIterator != torqueModelsPerBody_.end( ); torqueModelMapIterator++ )
        {
            for( innerTorqueIterator  = torqueModelMapIterator->second.begin( ); innerTorqueIterator !=
                 torqueModelMapIterator->second.end( ); innerTorqueIterator++ )
            {
                for( unsigned int j = 0; j < innerTorqueIterator->second.size( ); j++ )
                {
                    innerTorqueIterator->second[ j ]->resetTime( TUDAT_NAN );
                }
            }
        }
    }

    //! Calculates the state derivative of the rotational motion of the system.
    /*!
     *  Calculates the state derivative of the rotational motion of the system at the given time and rotational state.
     *  \param time Time (seconds since reference epoch) at which the system is to be updated.
     *  \param stateOfSystemToBeIntegrated List of 7 * bodiesToPropagate_.size( ), containing rotation quaternion/angular
     *  velocity of the bodies being propagated. The order of the values is defined by the order of bodies in
     *  bodiesToPropagate_
     *  \param stateDerivative Current state derivative (quaternion rate+angular acceleration) of system of bodies
     *  integrated numerically (returned by reference).
     */
    void calculateSystemStateDerivative(
            const TimeType time,
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        stateDerivative = Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( stateOfSystemToBeIntegrated.rows( ), 1 );
        std::vector< Eigen::Vector3d > torquesActingOnBodies = sumTorquesPerBody( );

        for( unsigned int i = 0; i < torquesActingOnBodies.size( ); i++ )
        {
            Eigen::Matrix< StateScalarType, 4, 1 > currentQuaternion = ( stateOfSystemToBeIntegrated.block( 7 * i, 0, 4, 1 ) ).normalized( );
            Eigen::Matrix< StateScalarType, 3, 1 > currentBodyFixedRotationRate = stateOfSystemToBeIntegrated.block( 7 * i + 4, 0, 3, 1 );

            stateDerivative.block( 7 * i, 0, 4, 1 ) = calculateQuaternionDerivative(
                        currentQuaternion.template cast< double >( ), currentBodyFixedRotationRate.template cast< double >( ) ).
                    template cast< StateScalarType >( );
            stateDerivative.block( 7 * i + 4, 0, 3, 1 ) = evaluateRotationalEquationsOfMotion(
                        bodyInertiaTensorFunctions_.at( i )( ), torquesActingOnBodies.at( i ),
                        currentBodyFixedRotationRate.template cast< double >( ),
                        bodyInertiaTensorTimeDerivativeFunctions_.at( i )( ) ).template cast< StateScalarType >( );

        }
    }

    //! Function to convert the state in the conventional form to the propagator-specific form.
    /*!
     * Function to convert the state in the conventional form to the propagator-specific form. For this propagator,
     * the two are equivalent, and this function returns the input state.
     * \param outputSolution State in 'conventional form'
     * \param time Current time at which the state is valid (not used in this class).
     * \return State (outputSolution), converted to the 'propagator-specific form' (which is equal to outputSolution).
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& outputSolution, const TimeType& time )
    {
        return outputSolution;
    }

    //! Function to convert the propagator-specific form of the state to the conventional form.
    /*!
     * Function to convert the propagator-specific form of the state to the conventional form. For the this propagator,
     * the two are equivalent, and this function returns the input state.
     * \param internalSolution State in propagator-specific form (which is equal to outputSolution to conventional form for
     * this propagator)
     * \param time Current time at which the state is valid (not used in this class).
     * \param currentLocalSoluton State (internalSolution), converted to the 'conventional form',
     *  which is equal to outputSolution for this class (returned by reference).
     */
    void convertToOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution, const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentLocalSoluton )
    {
        currentLocalSoluton = internalSolution;
    }
    //! Function to get list of names of bodies that are to be integrated numerically.
    /*!
     * Function to get list of names of bodies that are to be integrated numerically.
     * \return List of names of bodies that are to be integrated numerically.
     */
    std::vector< std::string > getBodiesToBeIntegratedNumerically( )
    {
        return bodiesToPropagate_;
    }

    //! Function to return the size of the state handled by the object
    /*!
     * Function to return the size of the state handled by the object
     * \return Size of the state under consideration (7 times the number if integrated bodies).
     */
    int getStateSize( )
    {
        return 7 * bodiesToPropagate_.size( );
    }

    Eigen::Vector3d getTotalTorqueForBody(
            const std::string& bodyName )
    {
        // Check if body is propagated.
        Eigen::Vector3d totalTorque = Eigen::Vector3d::Zero( );
        if( std::find( bodiesToPropagate_.begin( ),
                       bodiesToPropagate_.end( ),
                       bodyName ) == bodiesToPropagate_.end( ) )
        {
            std::string errorMessage = "Error when getting total torque for body " + bodyName +
                    ", no such torque is found";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            if( torqueModelsPerBody_.count( bodyName ) != 0 )
            {
                basic_astrodynamics::SingleBodyTorqueModelMap torquesOnBody =
                        torqueModelsPerBody_.at( bodyName );

                // Iterate over all torques acting on body
                for( innerTorqueIterator  = torquesOnBody.begin( );
                     innerTorqueIterator != torquesOnBody.end( );
                     innerTorqueIterator++ )
                {
                    for( unsigned int j = 0; j < innerTorqueIterator->second.size( ); j++ )
                    {
                        // Calculate torque and add to state derivative.
                        totalTorque += innerTorqueIterator->second[ j ]->getTorque( );
                    }
                }
            }
        }
        return totalTorque;
    }

    basic_astrodynamics::TorqueModelMap getTorquesMap( )
    {
        return torqueModelsPerBody_;
    }

protected:

    //! Function to get the total torques acting on each body, expressed in the body-fixed frames
    /*!
     * Function to get the total torques acting on each body, expressed in the body-fixed frames. The environment
     * and torque models must have been updated to the current state before calling this
     * function.
     * \return Total torques acting on each body, expressed in the body-fixed frames (in order of bodiesToPropagate_).
     */
    std::vector< Eigen::Vector3d > sumTorquesPerBody( )
    {
        using namespace basic_astrodynamics;

        std::vector< Eigen::Vector3d > torques;
        torques.resize( bodiesToPropagate_.size( ) );

        // Iterate over all bodies with torques
        for( unsigned int i = 0; i < bodiesToPropagate_.size( ); i++ )
        {
            torques[ i ].setZero( );

            if( torqueModelsPerBody_.count( bodiesToPropagate_[ i ] ) != 0 )
            {
                for( innerTorqueIterator  = torqueModelsPerBody_[ bodiesToPropagate_[ i ] ].begin( );
                     innerTorqueIterator != torqueModelsPerBody_[ bodiesToPropagate_[ i ] ].end( );
                     innerTorqueIterator++ )
                {
                    for( unsigned int j = 0; j < innerTorqueIterator->second.size( ); j++ )
                    {
                        torques[ i ] += ( innerTorqueIterator->second[ j ]->getTorque( ) );
                    }
                }
            }
        }

        return torques;
    }

    //! A map containing the list of torques acting on each body,
    /*!
     * A map containing the list of torques acting on each body, identifying the body being
     * acted on and the body acted on by an acceleration. The map has as key a string denoting the
     * name of the body the list of torques, provided as the value corresponding to a key, is
     * acting on.  This map-value is again a map with string as key, denoting the body exerting the
     * acceleration, and as value a pointer to an acceleration model.
     */
    basic_astrodynamics::TorqueModelMap torqueModelsPerBody_;

    //! List of names of bodies for which rotational state is to be propagated
    std::vector< std::string > bodiesToPropagate_;

    //! List of functions returning inertia tensors of bodiesToPropagate (in same order)
    std::vector< boost::function< Eigen::Matrix3d( ) > > bodyInertiaTensorFunctions_;

    //!  List of functions returning time derivatives of inertia tensors of bodiesToPropagate (in same order)
    std::vector< boost::function< Eigen::Matrix3d( ) > > bodyInertiaTensorTimeDerivativeFunctions_;


    //! Predefined iterator to save (de-)allocation time.
    basic_astrodynamics::TorqueModelMap::iterator torqueModelMapIterator;

    //! Predefined iterator to save (de-)allocation time.
    basic_astrodynamics::SingleBodyTorqueModelMap::iterator innerTorqueIterator;

};

}

}

#endif // TUDAT_ROTATIONALMOTIONSTATEDERIVATIVE_H
