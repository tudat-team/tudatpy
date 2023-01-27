/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ENVIRONMENTUPDATER_H
#define TUDAT_ENVIRONMENTUPDATER_H

#include <vector>
#include <string>
#include <map>




#include <functional>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/astro/propagators/environmentUpdateTypes.h"

namespace tudat
{

namespace propagators
{

//! Class used to update the environment during numerical integration.
/*!
 *  Class used to update the environment during numerical integration. The class ensures that the
 *  current state of the numerical integration is properly set, and that all the environment models
 *  that are used during the numerical integration are updated to the current time and state in the
 *  correct order.
 */
template< typename StateScalarType, typename TimeType >
class EnvironmentUpdater
{
public:

    //! Constructor
    /*!
     * Constructor, provides the required settings for updating the environment.
     * \param bodyList List of body objects, this list encompasses all environment object in the
     * simulation.
     * \param updateSettings List of updates of the environment models that are required.
     * The list defines per model type (key) the bodies for which this environment model should be
     * updated (values)
     * \param integratedStates This map provides the list of identifiers for the numerically
     * integrated states. The type of integrated state (key) is defined for each reference
     * (value). Note that for the numerical integration of translational motion, the entry
     * in the pair will have a second entry that is empty (""), with the first entry defining
     * the body that is integrated.
     */
    EnvironmentUpdater(
            const simulation_setup::SystemOfBodies& bodyList,
            const std::map< EnvironmentModelsToUpdate, std::vector< std::string > >& updateSettings,
            const std::map< IntegratedStateType,
            std::vector< std::tuple< std::string, std::string, PropagatorType > > >& integratedStates =
            ( std::map< IntegratedStateType,
            std::vector< std::tuple< std::string, std::string, PropagatorType > > >( ) ) ):
        bodyList_( bodyList ), integratedStates_( integratedStates )
    {
        // Set update function to be evaluated as dependent variables of state and time during each
        // integration time step.
        setUpdateFunctions( updateSettings );
    }

    //! Function to update the environment to the current state and time.
    /*!
     * Function to update the environment to the current state and time. This function calls the
     * dependent variable functions set by the setUpdateFunctions function. By default, the
     * numerically integrated states are set in the environment first. This may be overridden by
     * using the setIntegratedStatesFromEnvironment variable, which forces the function to ignore
     * specific integrated states and update them from the existing environment modelsad.
     * \param currentTime Current time.
     * \param integratedStatesToSet Current list of integrated states, with specific integrated
     * states defined by integratedStates_ member variable. Note that these states must have been
     * converted to the global state before input to this function, as is done by the
     * convertCurrentStateToGlobalRepresentationPerType function of the DynamicsStateDerivativeModel.
     * \param setIntegratedStatesFromEnvironment Integrated state types which are not to be used for
     * updating the environment, but which are to be set from existing environment models instead.
     */
    void updateEnvironment(
            const TimeType currentTime,
            const std::unordered_map< IntegratedStateType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >&
            integratedStatesToSet,
            const std::vector< IntegratedStateType >& setIntegratedStatesFromEnvironment =
            std::vector< IntegratedStateType >( ) )
    {
        // Check consistency of input.
        if( ( integratedStatesToSet.size( ) + setIntegratedStatesFromEnvironment.size( ) ) != integratedStates_.size( ) )
        {
            throw std::runtime_error( "Error when updating environment, input size is inconsistent " +
                                      std::to_string( integratedStatesToSet.size( ) ) + " " +
                                      std::to_string( setIntegratedStatesFromEnvironment.size( ) ) + " " +
                                      std::to_string( integratedStates_.size( ) ) );
        }

        for( unsigned int i = 0; i < resetFunctionVector_.size( ); i++ )
        {
            resetFunctionVector_.at( i ).template get< 2 >( )( );
        }

        // Set integrated state variables in environment.
        setIntegratedStatesInEnvironment( integratedStatesToSet );

        // Set current state from environment for override settings setIntegratedStatesFromEnvironment
        setStatesFromEnvironment( setIntegratedStatesFromEnvironment, currentTime );

        // Evaluate time-dependent update functions (dependent variables of state and time)
        // determined by setUpdateFunctions
        for( unsigned int i = 0; i < updateFunctionVector_.size( ); i++ )
        {
            updateFunctionVector_.at( i ).template get< 2 >( )( currentTime );
        }
    }

private:

    //! Function to set numerically integrated states in environment.
    /*!
     * Function to set numerically integrated states in environment.  Note that these states must
     * have been converted to the global state before input to this function, as is done by the
     * convertCurrentStateToGlobalRepresentationPerType function of the
     * DynamicsStateDerivativeModel.
     * \param integratedStatesToSet Integrated states which are to be set in environment.
     */
    void setIntegratedStatesInEnvironment(
            const std::unordered_map< IntegratedStateType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >&
            integratedStatesToSet )
    {
        // Iterate over state types and set states in environment
        for( integratedStateIterator_ = integratedStatesToSet.begin( );
             integratedStateIterator_ != integratedStatesToSet.end( );
             integratedStateIterator_++ )
        {
            switch( integratedStateIterator_->first )
            {
            case translational_state:
            {
                // Set translational states for bodies provided as input.
                for( unsigned int i = 0; i < integratedStates_[ translational_state ].size( ); i++ )
                {
                    bodyList_.at( std::get< 0 >( integratedStates_[ translational_state ][ i ] ) )->template
                            setTemplatedState< StateScalarType >(
                                integratedStateIterator_->second.segment( i * 6, 6 ) );
                }
                break;
            }
            case rotational_state:
            {
                std::vector< std::tuple< std::string, std::string, PropagatorType > > bodiesWithIntegratedStates =
                        integratedStates_.at( rotational_state );
                for( unsigned int i = 0; i < bodiesWithIntegratedStates.size( ); i++ )
                {
                    bodyList_.at( std::get< 0 >( bodiesWithIntegratedStates[ i ] ) )->setCurrentRotationalStateToLocalFrame(
                                integratedStateIterator_->second.segment( i * 7, 7 ).template cast< double >( ) );
                }
                break;
            }
            case body_mass_state:
            {
                // Set mass for bodies provided as input.
                std::vector< std::tuple< std::string, std::string, PropagatorType > > bodiesWithIntegratedMass =
                        integratedStates_.at( body_mass_state );

                for( unsigned int i = 0; i < bodiesWithIntegratedMass.size( ); i++ )
                {
                    bodyList_.at( std::get< 0 >( bodiesWithIntegratedMass[ i ] ) )
                            ->setConstantBodyMass( integratedStateIterator_->second( i ) );
                }
                break;
            }
            case custom_state:
            {
                break;
            }
            default:
                throw std::runtime_error( "Error, could not find integrated state settings for " +
                                          std::to_string( integratedStateIterator_->first ) );
            }
        }
    }
    //! Function to explicitly use existing environment models to update current states of integrated bodies
    /*!
     * Function to explicitly use existing environment models to update current states of integrated
     * bodies, overriding the numerically integrated states.
     * \param statesToSet Integrated state types which are not to be used for updating the environment,
     * but which are to be set from existing environment models instead.
     * \param currentTime Time to which environment is to be updated.
     */
    void setStatesFromEnvironment(
            const std::vector< IntegratedStateType >& statesToSet,
            const TimeType currentTime )
    {
        // Iterate over selected state types.
        for( unsigned int i = 0; i < statesToSet.size( ); i++ )
        {
            switch( statesToSet.at( i ) )
            {
            case translational_state:
            {
                // Iterate over all integrated translational states.
                std::vector< std::tuple< std::string, std::string, PropagatorType > > bodiesWithIntegratedStates =
                        integratedStates_[ translational_state ];
                for( unsigned int i = 0; i < bodiesWithIntegratedStates.size( ); i++ )
                {
                    bodyList_.at( std::get< 0 >( bodiesWithIntegratedStates[ i ] ) )->
                            template setStateFromEphemeris< StateScalarType, TimeType >( currentTime );

                }
                break;
            }
            case rotational_state:
            {
                std::vector< std::tuple< std::string, std::string, PropagatorType > > bodiesWithIntegratedStates =
                        integratedStates_.at( rotational_state );
                for( unsigned int i = 0; i < bodiesWithIntegratedStates.size( ); i++ )
                {
                    bodyList_.at( std::get< 0 >( bodiesWithIntegratedStates[ i ] ) )->template setCurrentRotationalStateToLocalFrameFromEphemeris< TimeType >(
                                currentTime );
                }
                break;
            }
            case body_mass_state:
            {
                // Iterate over all integrated masses.
                std::vector< std::tuple< std::string, std::string, PropagatorType > > bodiesWithIntegratedStates =
                        integratedStates_.at( body_mass_state );
                for( unsigned int i = 0; i < bodiesWithIntegratedStates.size( ); i++ )
                {
                    bodyList_.at( std::get< 0 >( bodiesWithIntegratedStates[ i ] ) )->
                            updateMass( currentTime );

                }
                break;
            }
            default:
                throw std::runtime_error( "Error, could not find  state settings for " +
                                          std::to_string( statesToSet.at( i ) ) );
            }
        }
    }

    //! Function to set the order in which the updateFunctionVector_ is to be updated.
    /*!
     *  Function to set the order in which the updateFunctionVector_ is to be updated. Order is determined recursively in
     *  this function,  the variable iterationNumber keeps track of the number of nested calls to the function.
     *  \param iterationNumber Number of subsequent calls to this funtion
     */
    void setUpdateFunctionOrder( const int iterationNumber = 0 )
    {
        bool rerunUpdateOrder = 0;

        // Iterate over all update functions
        for( unsigned int i = 0; i < updateFunctionVector_.size( ); i++ )
        {
            // Check if environment model is rotational state.
            if( updateFunctionVector_.at( i ).template get< 0 >( ) == body_rotational_state_update )
            {
                // Check id body has no rotational ephemeris (i.e. if rotation comes from iterationNumber ).
                std::shared_ptr< ephemerides::AerodynamicAngleRotationalEphemeris > angleBasedRotationModel =
                        std::dynamic_pointer_cast< ephemerides::AerodynamicAngleRotationalEphemeris >(
                            bodyList_.at( updateFunctionVector_.at( i ).template get< 1 >( ) )->getRotationalEphemeris( ) );
                if( angleBasedRotationModel != nullptr )
                {
                    std::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator =
                            angleBasedRotationModel->getAerodynamicAngleCalculator( );

                    unsigned int centralTranslationalUpdateIndex = 0;
                    unsigned int centralRotationalUpdateIndex = 0;
                    unsigned int vehicleTranslationalUpdateIndex = 0;
                    unsigned int vehicleRotationalUpdateIndex = i;
                    unsigned int flightCoditionsUpdateIndex = 0;

                    bool centralTranslationalUpdateIndexSet = false;
                    bool centralRotationalUpdateIndexSet = false;
                    bool vehicleTranslationalUpdateIndexSet = false;
                    bool flightConditionsUpdateIndexSet = false;

                    // Check if the state or orientation of the central body of AerodynamicAngleCalculator is updated.
                    for( unsigned int j = 0; j < updateFunctionVector_.size( ); j++ )
                    {
                        if( ( updateFunctionVector_.at( j ).template get< 0 >( ) == body_translational_state_update ) &&
                                ( updateFunctionVector_.at( j ).template get< 1 >( ) ==
                                  aerodynamicAngleCalculator->getCentralBodyName( ) ) )
                        {
                            if( updateFunctionVector_.at( i ).template get< 0 >( ) == body_rotational_state_update )
                            {
                                // Check id body has no rotational ephemeris (i.e. if rotation comes from iterationNumber ).
                                if( bodyList_.at( updateFunctionVector_.at( i ).template get< 1 >( ) )->getRotationalEphemeris( ) == nullptr )
                                {
//                                    // Check if DependentOrientationCalculator is an AerodynamicAngleCalculator.
//                                    std::shared_ptr< reference_frames::DependentOrientationCalculator > dependentOrientationCalculator =
//                                            bodyList_.at( updateFunctionVector_.at( i ).template get< 1 >( ) )->
//                                            getDependentOrientationCalculator( );
//                                    std::shared_ptr< reference_frames::AerodynamicAngleCalculator > aerodynamicAngleCalculator =
//                                            std::dynamic_pointer_cast< reference_frames::AerodynamicAngleCalculator >(
//                                                dependentOrientationCalculator );

//                                    // Check if properties of AerodynamicAngleCalculator are such that a different update order is warranted.
//                                    if( std::dynamic_pointer_cast< reference_frames::AerodynamicAngleCalculator >(
//                                                dependentOrientationCalculator ) != nullptr )
                                    {
                                        unsigned int centralTranslationalUpdateIndex = 0;
                                        unsigned int centralRotationalUpdateIndex = 0;
                                        unsigned int vehicleTranslationalUpdateIndex = 0;
                                        unsigned int vehicleRotationalUpdateIndex = i;
                                        unsigned int flightCoditionsUpdateIndex = 0;

                                        bool centralTranslationalUpdateIndexSet = false;
                                        bool centralRotationalUpdateIndexSet = false;
                                        bool vehicleTranslationalUpdateIndexSet = false;
                                        bool flightConditionsUpdateIndexSet = false;

                                        // Check if the state or orientation of the central body of AerodynamicAngleCalculator is updated.
                                        for( unsigned int j = 0; j < updateFunctionVector_.size( ); j++ )
                                        {
                                            if( ( updateFunctionVector_.at( j ).template get< 0 >( ) == body_translational_state_update ) &&
                                                    ( updateFunctionVector_.at( j ).template get< 1 >( ) ==
                                                      aerodynamicAngleCalculator->getCentralBodyName( ) ) )
                                            {
                                                centralTranslationalUpdateIndex = j;
                                                centralTranslationalUpdateIndexSet = true;
                                            }

                                            if( ( updateFunctionVector_.at( j ).template get< 0 >( ) == body_rotational_state_update ) &&
                                                    ( updateFunctionVector_.at( j ).template get< 1 >( ) ==
                                                      aerodynamicAngleCalculator->getCentralBodyName( ) ) )
                                            {
                                                centralRotationalUpdateIndex = j;
                                                centralRotationalUpdateIndexSet = true;
                                            }

                                            if( ( updateFunctionVector_.at( j ).template get< 0 >( ) == body_translational_state_update ) &&
                                                    ( updateFunctionVector_.at( j ).template get< 1 >( ) ==
                                                      updateFunctionVector_.at( i ).template get< 1 >( ) ) )
                                            {
                                                vehicleTranslationalUpdateIndex = j;
                                                vehicleTranslationalUpdateIndexSet = true;
                                            }

                                            if( ( updateFunctionVector_.at( j ).template get< 0 >( ) == vehicle_flight_conditions_update ) &&
                                                    ( updateFunctionVector_.at( j ).template get< 1 >( ) ==
                                                      updateFunctionVector_.at( i ).template get< 1 >( ) ) )
                                            {
                                                flightCoditionsUpdateIndex = j;
                                                flightConditionsUpdateIndexSet = true;
                                            }
                                        }


                                        std::vector< int > indices;
                                        std::vector< boost::tuple< EnvironmentModelsToUpdate, std::string, std::function< void( const double ) > > > updatesToMove;

                                        if( centralTranslationalUpdateIndexSet )
                                        {
                                            indices.push_back( centralTranslationalUpdateIndex );
                                            updatesToMove.push_back( updateFunctionVector_.at( centralTranslationalUpdateIndex ) );
                                        }

                                        if( centralRotationalUpdateIndexSet )
                                        {
                                            indices.push_back( centralRotationalUpdateIndex );
                                            updatesToMove.push_back( updateFunctionVector_.at( centralRotationalUpdateIndex ) );
                                        }

                                        if( vehicleTranslationalUpdateIndexSet )
                                        {
                                            indices.push_back( vehicleTranslationalUpdateIndex );
                                            updatesToMove.push_back( updateFunctionVector_.at( vehicleTranslationalUpdateIndex ) );
                                        }

                                        if( flightConditionsUpdateIndexSet )
                                        {
                                            indices.push_back( flightCoditionsUpdateIndex );
                                            updatesToMove.push_back( updateFunctionVector_.at( flightCoditionsUpdateIndex ) );
                                        }

                                        indices.push_back( vehicleRotationalUpdateIndex );
                                        updatesToMove.push_back( updateFunctionVector_.at( vehicleRotationalUpdateIndex ) );
                                        std::vector< int > unorderedIndices = indices;
                                        std::sort( indices.begin( ), indices.end( ) );

                                        if( indices != unorderedIndices )
                                        {

                                            for( unsigned int k = 0; k < updatesToMove.size( ); k++ )
                                            {
                                                updateFunctionVector_[ indices.at( k ) ] = updatesToMove.at( k );
                                            }
                                            rerunUpdateOrder = true;
                                            break;
                                        }
                                    }
                                }
                            }
                            centralTranslationalUpdateIndex = j;
                            centralTranslationalUpdateIndexSet = true;
                        }

                        if( ( updateFunctionVector_.at( j ).template get< 0 >( ) == body_rotational_state_update ) &&
                                ( updateFunctionVector_.at( j ).template get< 1 >( ) ==
                                  aerodynamicAngleCalculator->getCentralBodyName( ) ) )
                        {
                            centralRotationalUpdateIndex = j;
                            centralRotationalUpdateIndexSet = true;
                        }

                        if( ( updateFunctionVector_.at( j ).template get< 0 >( ) == body_translational_state_update ) &&
                                ( updateFunctionVector_.at( j ).template get< 1 >( ) ==
                                  updateFunctionVector_.at( i ).template get< 1 >( ) ) )
                        {
                            vehicleTranslationalUpdateIndex = j;
                            vehicleTranslationalUpdateIndexSet = true;
                        }

                        if( ( updateFunctionVector_.at( j ).template get< 0 >( ) == vehicle_flight_conditions_update ) &&
                                ( updateFunctionVector_.at( j ).template get< 1 >( ) ==
                                  updateFunctionVector_.at( i ).template get< 1 >( ) ) )
                        {
                            flightCoditionsUpdateIndex = j;
                            flightConditionsUpdateIndexSet = true;
                        }
                    }


                    std::vector< int > indices;
                    std::vector< boost::tuple< EnvironmentModelsToUpdate, std::string, std::function< void( const double ) > > > updatesToMove;

                    if( centralTranslationalUpdateIndexSet )
                    {
                        indices.push_back( centralTranslationalUpdateIndex );
                        updatesToMove.push_back( updateFunctionVector_.at( centralTranslationalUpdateIndex ) );
                    }

                    if( centralRotationalUpdateIndexSet )
                    {
                        indices.push_back( centralRotationalUpdateIndex );
                        updatesToMove.push_back( updateFunctionVector_.at( centralRotationalUpdateIndex ) );
                    }

                    if( vehicleTranslationalUpdateIndexSet )
                    {
                        indices.push_back( vehicleTranslationalUpdateIndex );
                        updatesToMove.push_back( updateFunctionVector_.at( vehicleTranslationalUpdateIndex ) );
                    }

                    if( flightConditionsUpdateIndexSet )
                    {
                        indices.push_back( flightCoditionsUpdateIndex );
                        updatesToMove.push_back( updateFunctionVector_.at( flightCoditionsUpdateIndex ) );
                    }

                    indices.push_back( vehicleRotationalUpdateIndex );
                    updatesToMove.push_back( updateFunctionVector_.at( vehicleRotationalUpdateIndex ) );
                    std::vector< int > unorderedIndices = indices;
                    std::sort( indices.begin( ), indices.end( ) );

                    if( indices != unorderedIndices )
                    {

                        for( unsigned int k = 0; k < updatesToMove.size( ); k++ )
                        {
                            updateFunctionVector_[ indices.at( k ) ] = updatesToMove.at( k );
                        }
                        rerunUpdateOrder = true;
                        break;
                    }
                }
            }
        }

        // Define escape condition in case of infinite loop.
        if( iterationNumber > 10000 )
        {
            throw std::runtime_error( "Error when finding update order; stuck in infinite loop" );
        }

        // Rerun function if needed.
        if( rerunUpdateOrder )
        {
            setUpdateFunctionOrder( iterationNumber + 1 );
        }
    }

    //! Function to set the update functions for the environment from the required update settings.
    /*!
     * Function to set the update functions for the environment from the required update settings.
     * \param updateSettings Settings for the environment updates.
     */
    void setUpdateFunctions( const std::map< EnvironmentModelsToUpdate,
                             std::vector< std::string > >& updateSettings )
    {
        std::map< EnvironmentModelsToUpdate,
                std::vector< std::pair< std::string, std::function< void( const double ) > > > > updateTimeFunctionList;

        // Iterate over all required updates and set associated update function in lists
        for( std::map< EnvironmentModelsToUpdate,
             std::vector< std::string > >::const_iterator updateIterator =
             updateSettings.begin( ); updateIterator != updateSettings.end( ); updateIterator++ )
        {
            // Get list of bodies for which current environment type is to be updated.
            std::vector< std::string > currentBodies = updateIterator->second;
            for( unsigned int i = 0; i < currentBodies.size( ); i++ )
            {
                if( currentBodies.at( i ) != "" )
                {
                    // Check whether body exists
                    if( bodyList_.count( currentBodies.at( i ) ) == 0 )
                    {
                        throw std::runtime_error(
                                    "Error when setting environment update functions, could not find body " +
                                    currentBodies.at( i ) );
                    }

                    // Find type of environment.
                    switch( updateIterator->first )
                    {

                    // If requested body is not propagated, add to list.
                    case body_translational_state_update:
                    {
                        bool addUpdate = 1;

                        // Check if mass is propagated
                        if( integratedStates_.count( translational_state ) > 0 )
                        {
                            // Check if current body is propagated
                            std::string bodyToCheck = currentBodies.at( i );
                            std::vector< std::string > integratedTranslationalStates =
                                    utilities::getFirstTupleEntryVector( integratedStates_.at( translational_state ) );
                            if( std::find( integratedTranslationalStates.begin( ),
                                           integratedTranslationalStates.end( ),
                                           bodyToCheck ) != integratedTranslationalStates.end( ) )
                            {
                                addUpdate = 0;
                            }
                        }

                        // Add state update function to list.
                        if( addUpdate == 1 )
                        {
                            std::function< void( const TimeType ) > stateSetFunction =
                                    std::bind(
                                        &simulation_setup::Body
                                        ::setStateFromEphemeris< StateScalarType, TimeType >,
                                        bodyList_.at( currentBodies.at( i ) ), std::placeholders::_1 );

                            updateTimeFunctionList[ body_translational_state_update ].push_back(
                                        std::make_pair( currentBodies.at( i ), stateSetFunction ) );

                            resetFunctionVector_.push_back(
                                        boost::make_tuple(
                                            body_translational_state_update, currentBodies.at( i ),
                                            std::bind( &simulation_setup::Body::recomputeStateOnNextCall,
                                                       bodyList_.at( currentBodies.at( i ) ) ) ) );
                        }
                        break;
                    }
                    case body_rotational_state_update:
                    {

                        bool addUpdate = 1;
                        if( integratedStates_.count( rotational_state ) > 0 )
                        {
                            std::string bodyToCheck = currentBodies.at( i );
                            std::vector< std::string > integratedRotationalStates =
                                    utilities::getFirstTupleEntryVector( integratedStates_.at( rotational_state ) );
                            if( std::find( integratedRotationalStates.begin( ), integratedRotationalStates.end( ), bodyToCheck ) !=
                                    integratedRotationalStates.end( ) )
                            {
                                addUpdate = 0;
                            }
                        }

                        if( addUpdate == 1 )
                        {

                            // Check if rotational ephemeris exists
                            if(  ( bodyList_.at( currentBodies.at( i ) )->getRotationalEphemeris( ) != nullptr )
                                 //                                 ||
                                 //                                 ( bodyList_.at( currentBodies.at( i ) )->getDependentOrientationCalculator( ) != nullptr )
                                 )
                            {
                                std::function< void( const TimeType ) > rotationalStateSetFunction =
                                        std::bind( &simulation_setup::Body
                                                   ::setCurrentRotationalStateToLocalFrameFromEphemeris< TimeType >,
                                                   bodyList_.at( currentBodies.at( i ) ), std::placeholders::_1 );
                                updateTimeFunctionList[ body_rotational_state_update ].push_back(
                                            std::make_pair( currentBodies.at( i ), rotationalStateSetFunction ) );
                                if( bodyList_.at( currentBodies.at( i ) )->getRotationalEphemeris( ) != nullptr )
                                {
                                    resetFunctionVector_.push_back(
                                                boost::make_tuple(
                                                    body_rotational_state_update, currentBodies.at( i ),
                                                    std::bind( &ephemerides::RotationalEphemeris::
                                                               resetCurrentTime, bodyList_.at( currentBodies.at( i ) )->
                                                               getRotationalEphemeris( ) ) ) );
                                }

                                //                                if( bodyList_.at( currentBodies.at( i ) )->getRotationalEphemeris( ) == nullptr )
                                //                                {
                                //                                    resetFunctionVector_.push_back(
                                //                                                boost::make_tuple(
                                //                                                    body_rotational_state_update, currentBodies.at( i ),
                                //                                                    std::bind( &reference_frames::DependentOrientationCalculator::
                                //                                                                 resetCurrentTime, bodyList_.at( currentBodies.at( i ) )->
                                //                                                                 getDependentOrientationCalculator( ), TUDAT_NAN ) ) );
                                //                                }
                            }
                            else
                            {
                                throw std::runtime_error(
                                            "Request rotation update of " + currentBodies.at( i ) +
                                            ", but body has no rotational ephemeris" );
                            }
                        }

                        break;

                    }
                    case body_mass_update:
                    {
                        bool addUpdate = 1;

                        // Check if translational state is propagated
                        if( integratedStates_.count( body_mass_state ) > 0 )
                        {
                            // Check if current body is propagated
                            std::string bodyToCheck = currentBodies.at( i );
                            std::vector< std::string > integratedBodyMasses =
                                    utilities::getFirstTupleEntryVector( integratedStates_.at( body_mass_state ) );
                            if( std::find( integratedBodyMasses.begin( ),
                                           integratedBodyMasses.end( ),
                                           bodyToCheck ) != integratedBodyMasses.end( ) )
                            {
                                addUpdate = 0;
                            }
                        }

                        if( addUpdate )
                        {
                            updateTimeFunctionList[ body_mass_update ].push_back(
                                        std::make_pair( currentBodies.at( i ),
                                                        std::bind( &simulation_setup::Body::updateMass,
                                                                   bodyList_.at( currentBodies.at( i ) ), std::placeholders::_1  ) ) );
                        }
                        break;
                    }
                    case spherical_harmonic_gravity_field_update:
                    {

                        // Check if body has time-dependent sh field
                        std::shared_ptr< gravitation::TimeDependentSphericalHarmonicsGravityField >
                                gravityField = std::dynamic_pointer_cast
                                < gravitation::TimeDependentSphericalHarmonicsGravityField >
                                (  bodyList_.at( currentBodies.at( i ) )->getGravityFieldModel( ) );
                        if( gravityField != nullptr )
                        {
                            updateTimeFunctionList[ spherical_harmonic_gravity_field_update ].push_back(
                                        std::make_pair(
                                            currentBodies.at( i ),
                                            std::bind( &gravitation
                                                       ::TimeDependentSphericalHarmonicsGravityField
                                                       ::update,
                                                       gravityField, std::placeholders::_1 ) ) );
                        }
                        // If no sh field at all, throw eeror.
                        else if( std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >
                                 (  bodyList_.at( currentBodies.at( i ) )->getGravityFieldModel( ) ) == nullptr )
                        {
                            throw std::runtime_error( "Request sh update of " + currentBodies.at( i ) +
                                                      ", but body has no sh model" );
                        }

                        break;
                    }
                    case vehicle_flight_conditions_update:
                    {
                        // Check if current body has flight conditions set.
                        if( bodyList_.at( currentBodies.at( i ) )->getFlightConditions( ) != nullptr )
                        {
                            // If vehicle has flight conditions, add flight conditions update
                            // function to update list.
                            updateTimeFunctionList[ vehicle_flight_conditions_update ].push_back(
                                        std::make_pair(
                                            currentBodies.at( i ), std::bind(
                                                &aerodynamics::FlightConditions::updateConditions,
                                                bodyList_.at( currentBodies.at( i ) )
                                                ->getFlightConditions( ), std::placeholders::_1 ) ) );

                            resetFunctionVector_.push_back(
                                        boost::make_tuple(
                                            vehicle_flight_conditions_update, currentBodies.at( i ),
                                            std::bind( &aerodynamics::FlightConditions::
                                                       resetCurrentTime, bodyList_.at( currentBodies.at( i ) )->
                                                       getFlightConditions( ) ) ) );
                        }
                        else
                        {
                            throw std::runtime_error(
                                        "Request flight condition update of " + currentBodies.at( i ) +
                                        ", but body has no flight conditions" );
                        }
                        break;
                    }
                    case radiation_source_model_update:
                    {
                        // Check if current body has radiation source model set.
                        if( bodyList_.at( currentBodies.at( i ) )->getRadiationSourceModel() == nullptr )
                        {
                            throw std::runtime_error(
                                        "Request radiation source model update of " + currentBodies.at( i ) +
                                        ", but body has no radiation source model" );
                        }
                        // If vehicle has radiation source model, add its update function to update list.
                        updateTimeFunctionList[ radiation_source_model_update ].push_back(
                                    std::make_pair(
                                        currentBodies.at( i ), std::bind(
                                            &electromagnetism::RadiationSourceModel::updateMembers,
                                            bodyList_.at( currentBodies.at( i ) )
                                            ->getRadiationSourceModel(), std::placeholders::_1 ) ) );
                        break;
                    }
                    case radiation_pressure_target_model_update:
                    {
                        // Check if current body has radiation pressure target model set.
                        if( bodyList_.at( currentBodies.at( i ) )->getRadiationPressureTargetModel() == nullptr )
                        {
                            throw std::runtime_error(
                                        "Request radiation pressure target model update of " + currentBodies.at( i ) +
                                        ", but body has no radiation pressure target model" );
                        }
                        // If vehicle has radiation pressure target model, add its update function to update list.
                        updateTimeFunctionList[ radiation_pressure_target_model_update ].push_back(
                                    std::make_pair(
                                        currentBodies.at( i ), std::bind(
                                            &electromagnetism::RadiationPressureTargetModel::updateMembers,
                                            bodyList_.at( currentBodies.at( i ) )
                                            ->getRadiationPressureTargetModel(), std::placeholders::_1 ) ) );
                        break;
                    }
                    case radiation_pressure_interface_update:
                    {
                        // RP-OLD
                        // Get body radiation pressure interface(s) (one per source)
                        std::map< std::string, std::shared_ptr< electromagnetism
                                ::RadiationPressureInterface > >
                                radiationPressureInterfaces =
                                bodyList_.at( currentBodies.at( i ) )->getRadiationPressureInterfaces( );

                        if( radiationPressureInterfaces.size( ) == 0 )
                        {
                            throw std::runtime_error(
                                        "Request radiation pressure update of " + currentBodies.at( i ) +
                                        ", but body has no radiation pressure interfaces" );
                        }
                        else if( radiationPressureInterfaces.size( ) > 1 )
                        {
                            std::cerr << "Warning, requested radiation pressure update of " << currentBodies.at( i ) <<
                                         ", but body has multiple radiation pressure interfaces: updating all." << std::endl;
                        }

                        // Add each interface update function to update list.
                        for( std::map< std::string,
                             std::shared_ptr< electromagnetism::RadiationPressureInterface > >
                             ::iterator iterator = radiationPressureInterfaces.begin( );
                             iterator != radiationPressureInterfaces.end( ); iterator++ )
                        {
                            updateTimeFunctionList[ radiation_pressure_interface_update ].push_back(
                                        std::make_pair( currentBodies.at( i ),
                                                        std::bind(
                                                            &electromagnetism
                                                            ::RadiationPressureInterface
                                                            ::updateInterface,
                                                            iterator->second, std::placeholders::_1 ) ) );
                        }
                        break;
                    }
                    }
                }
            }
        }

        // Create list of update functions.
        for( std::map< EnvironmentModelsToUpdate, std::vector< std::pair< std::string,
             std::function< void( const double ) > > > >::iterator updateTimeIterator  = updateTimeFunctionList.begin( );
             updateTimeIterator != updateTimeFunctionList.end( ); updateTimeIterator++ )
        {
            for( unsigned int i = 0; i < updateTimeIterator->second.size( ); i++ )
            {
                updateFunctionVector_.push_back(
                            boost::make_tuple( updateTimeIterator->first,  updateTimeIterator->second.at( i ).first,
                                               updateTimeIterator->second.at( i ).second ) );
            }
        }

        // Set update order of functions.
        setUpdateFunctionOrder( );
    }

    //! List of body objects, this list encompasses all environment object in the simulation.
    simulation_setup::SystemOfBodies bodyList_;


    //! list of identifiers for the numerically integrated states
    /*!
     * This map provides the list of identifiers for the numerically
     * integrated states. The type of integrated state (key) is defined for each reference
     * (value). Note that for the numerical integration of translational motion, the entry
     * in the pair will have a second entry that is empty (""), with the first entry defining
     * the body that is integrated.
     */
    std::map< IntegratedStateType, std::vector< std::tuple< std::string, std::string, PropagatorType > > >
    integratedStates_;

    //! List of time-dependent functions to call to update the environment.
    std::vector< boost::tuple< EnvironmentModelsToUpdate, std::string, std::function< void( const double ) > > >
    updateFunctionVector_;

    //! List of time-dependent functions to call to reset the time of the environment (to NaN signal recomputation for next
    //! time step).
    std::vector< boost::tuple< EnvironmentModelsToUpdate, std::string, std::function< void( ) > > > resetFunctionVector_;




    //! Predefined state history iterator for computational efficiency.
    typename std::unordered_map< IntegratedStateType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >::const_iterator
    integratedStateIterator_;


};

extern template class EnvironmentUpdater< double, double >;


} // namespace propagators

} // namespace tudat

#endif // TUDAT_ENVIRONMENTUPDATER_H
