/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_ROTATIONAL_MOTION_STATE_DERIVATIVE_H
#define TUDAT_ROTATIONAL_MOTION_STATE_DERIVATIVE_H

#include <vector>
#include <map>
#include <string>

#include <memory>
#include <functional>

#include "tudat/astro/basic_astro/torqueModel.h"

#include "tudat/astro/propagators/singleStateTypeDerivative.h"
#include "tudat/simulation/environment_setup/body.h"

namespace tudat
{

namespace propagators
{

// Enum listing propagator types for rotational dynamics that can be used.
//! @get_docstring(AvailableAcceleration.__docstring__)
enum RotationalPropagatorType
{
    undefined_rotational_propagator = -1,
    quaternions = 0,
    modified_rodrigues_parameters = 1,
    exponential_map = 2
};

// Function to evaluated the classical rotational equations of motion (Euler equations)
/*
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

// Class for computing the state derivative for rotational dynamics of N bodies.
/*
 *  Class for computing the state derivative for rotational dynamics of N bodies, using quaternion from body-fixed to inertial
 *  frame (in quaternion format) and angular velocity-vector of body expressed in body-fixed frame as the rotational state of a
 *  single body
 */
template< typename StateScalarType = double, typename TimeType = double >
class RotationalMotionStateDerivative: public propagators::SingleStateTypeDerivative< StateScalarType, TimeType >
{
public:

    using propagators::SingleStateTypeDerivative< StateScalarType, TimeType >::calculateSystemStateDerivative;

    // Constructor.
    /*
     * Constructor.
     * \param torqueModelsPerBody List of torque models (first map key body undergoing acceleration, second map key body exerting
     * acceleration).
     * \param propagatorType Type of propagator that is to be used (i.e. quaternions, etc.)
     * \param bodiesToPropagate List of names of bodies for which rotational state is to be propagated.
     * \param bodyInertiaTensorFunctions List of functions returning inertia tensors of bodiesToPropagate (in same order).
     * \param bodyInertiaTensorTimeDerivativeFunctions List of functions returning time derivatives of inertia tensors of
     * bodiesToPropagate (in same order). Default empty, denoting time-invariant inertia tensors.
     */
    RotationalMotionStateDerivative(
            const basic_astrodynamics::TorqueModelMap& torqueModelsPerBody,
            const RotationalPropagatorType propagatorType,
            const std::vector< std::string >& bodiesToPropagate,
            std::vector< std::function< Eigen::Matrix3d( ) > > bodyInertiaTensorFunctions,
            std::vector< std::function< Eigen::Matrix3d( ) > > bodyInertiaTensorTimeDerivativeFunctions =
            std::vector< std::function< Eigen::Matrix3d( ) > >( ) ):
        propagators::SingleStateTypeDerivative< StateScalarType, TimeType >(
            propagators::rotational_state ),
        torqueModelsPerBody_( torqueModelsPerBody ),
        propagatorType_( propagatorType ),
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
                bodyInertiaTensorTimeDerivativeFunctions_.push_back( [ ]( ){ return Eigen::Matrix3d::Zero( ); } );
            }
        }
        else if( bodiesToPropagate_.size( ) != bodyInertiaTensorTimeDerivativeFunctions_.size( ) )
        {
            throw std::runtime_error(
                        "Error when making rotational state derivative model, inertia tensor time derivative list is of incompatible size" );
        }

        for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
        {
            if( torqueModelsPerBody_.count( bodiesToPropagate.at( i ) ) == 0 )
            {
                torqueModelsPerBody_[ bodiesToPropagate.at( i ) ][ bodiesToPropagate.at( i ) ] =
                        std::vector< std::shared_ptr< basic_astrodynamics::TorqueModel > >( );
            }
        }

        verifyInput( );
    }

    // Destructor
    virtual ~RotationalMotionStateDerivative( ){ }

    // Function to clear any reference/cached values of state derivative model
    /*
     * Function to clear any reference/cached values of state derivative model, in addition to those performed in the
     * clearTranslationalStateDerivativeModel function. Default implementation is empty.
     */
    virtual void clearDerivedRotationalStateDerivativeModel( ){ }

    // Function to clear reference/cached values of acceleration models
    /*
     * Function to clear reference/cached values of acceleration models, to ensure that they are all recalculated.
     */
    void clearRotationalStateDerivativeModel( )
    {
        for( torqueModelMapIterator = torqueModelsPerBody_.begin( );
             torqueModelMapIterator != torqueModelsPerBody_.end( ); torqueModelMapIterator++ )
        {
            for( innerTorqueIterator = torqueModelMapIterator->second.begin( ); innerTorqueIterator !=
                 torqueModelMapIterator->second.end( ); innerTorqueIterator++ )
            {
                for( unsigned int j = 0; j < innerTorqueIterator->second.size( ); j++ )
                {
                    innerTorqueIterator->second[ j ]->resetTime( TUDAT_NAN );
                }
            }
        }
    }

    // Function to clear reference/cached values of translational state derivative model
    /*
     * Function to clear reference/cached values of translational state derivative model. For each derived class, this
     * entails resetting the current time in the acceleration models to NaN (see clearRotationalStateDerivativeModel).
     * Every derived class requiring additional values to be cleared should implement the
     * clearDerivedRotationalStateDerivativeModel function.
     */
    void clearStateDerivativeModel(  )
    {
        clearRotationalStateDerivativeModel( );
        clearDerivedRotationalStateDerivativeModel( );
    }

    // Function to update the state derivative model to the current time.
    /*
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
            for( innerTorqueIterator = torqueModelMapIterator->second.begin( ); innerTorqueIterator !=
                 torqueModelMapIterator->second.end( ); innerTorqueIterator++ )
            {
                for( unsigned int j = 0; j < innerTorqueIterator->second.size( ); j++ )
                {
                    innerTorqueIterator->second[ j ]->updateMembers( currentTime );
                }
            }
        }
    }

    // Function to convert the propagator-specific form of the state to the conventional form in the global frame.
    /*
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
        this->convertToOutputSolution( internalSolution, time, currentRotationalState );
    }

    // Function to get list of names of bodies that are to be integrated numerically.
    /*
     * Function to get list of names of bodies that are to be integrated numerically.
     * \return List of names of bodies that are to be integrated numerically.
     */
    std::vector< std::string > getBodiesToBeIntegratedNumerically( )
    {
        return bodiesToPropagate_;
    }

    // Function to return the size of the state handled by the object.
    /*
     * Function to return the size of the state handled by the object.
     * \return Size of the state under consideration (7 times the number of integrated bodies).
     */
    int getConventionalStateSize( )
    {
        return 7 * bodiesToPropagate_.size( );
    }

    // Function to get total torque acting on specific body.
    /*
     * Function to get total torque acting on specific body.
     * \return bodyName Name of body for which total torque needs to be retrieved.
     */
    Eigen::Vector3d getTotalTorqueForBody( const std::string& bodyName )
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

    // Function to get map of torques acting on each body.
    /*
     * Function to get map of torques acting on each body.
     * \return TorqueModelMap Map of torques acting on each body.
     */
    basic_astrodynamics::TorqueModelMap getTorquesMap( )
    {
        return torqueModelsPerBody_;
    }

    // Function to get type of propagator that is to be used (i.e., quaternions, etc.).
    /*
     * Function to type of propagator that is to be used (i.e., quaternions, etc.).
     * \return Type of propagator that is to be used (i.e., quaternions, etc.).
     */
    RotationalPropagatorType getRotationalPropagatorType( )
    {
        return propagatorType_;
    }

protected:

    void verifyInput( )
    {
        for( unsigned int i = 0; i < bodiesToPropagate_.size( ); i++ )
        {
            if( torqueModelsPerBody_.count( bodiesToPropagate_.at( i ) ) == 0 )
            {
                throw std::runtime_error( "Error, requested propagation of rotational dynamics of body " +
                                          bodiesToPropagate_.at( i ) +
                                          ", but no torque models provided" );
            }
        }

        for( auto it : torqueModelsPerBody_ )
        {
            if( std::find( bodiesToPropagate_.begin( ),
                           bodiesToPropagate_.end( ),
                           it.first ) == bodiesToPropagate_.end( ) )
            {
                throw std::runtime_error( "Error, provided torque models for body " +
                                          it.first +
                                          ", but this body is not included in list of bodies for which rotation is to be propagated." );
            }
        }
    }

    // Function to get the total torques acting on each body, expressed in the body-fixed frames
    /*
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

    // A map containing the list of torques acting on each body,
    /*
     * A map containing the list of torques acting on each body, identifying the body being
     * acted on and the body acted on by an acceleration. The map has as key a string denoting the
     * name of the body the list of torques, provided as the value corresponding to a key, is
     * acting on. This map-value is again a map with string as key, denoting the body exerting the
     * acceleration, and as value a pointer to an acceleration model.
     */
    basic_astrodynamics::TorqueModelMap torqueModelsPerBody_;

    // Type of propagator that is to be used (i.e., quaternions, etc.)
    RotationalPropagatorType propagatorType_;

    // List of names of bodies for which rotational state is to be propagated
    std::vector< std::string > bodiesToPropagate_;

    // List of functions returning inertia tensors of bodiesToPropagate (in same order)
    std::vector< std::function< Eigen::Matrix3d( ) > > bodyInertiaTensorFunctions_;

    //  List of functions returning time derivatives of inertia tensors of bodiesToPropagate (in same order)
    std::vector< std::function< Eigen::Matrix3d( ) > > bodyInertiaTensorTimeDerivativeFunctions_;

    // Predefined iterator to save (de-)allocation time.
    basic_astrodynamics::TorqueModelMap::iterator torqueModelMapIterator;

    // Predefined iterator to save (de-)allocation time.
    basic_astrodynamics::SingleBodyTorqueModelMap::iterator innerTorqueIterator;

};


extern template class RotationalMotionStateDerivative< double, double >;

#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
extern template class RotationalMotionStateDerivative< long double, double >;
extern template class RotationalMotionStateDerivative< double, Time >;
extern template class RotationalMotionStateDerivative< long double, Time >;
#endif

} // namespace propagators

} // namespace tudat

#endif // TUDAT_ROTATIONAL_MOTION_STATE_DERIVATIVE_H
