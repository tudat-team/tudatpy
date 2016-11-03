/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROPAGATIONSETTINGS_H
#define TUDAT_PROPAGATIONSETTINGS_H

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <unordered_map>

#include <boost/make_shared.hpp>
#include <boost/lexical_cast.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/massRateModel.h"
#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"
#include "Tudat/Astrodynamics/Propagators/nBodyStateDerivative.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationOutputSettings.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationTerminationSettings.h"

namespace tudat
{

namespace propagators
{

//! Base class for defining setting of a propagator
/*!
 *  Base class for defining setting of a propagator. This class is non-functional, and each state type requires its
 *  own derived class (which may have multiple derived classes of its own).
 */
template< typename StateScalarType >
class PropagatorSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param stateType Type of state being propagated
     * \param initialBodyStates Initial state used as input for numerical integration
     * \param terminationSettings Settings for creating the object that checks whether the propagation is finished.
     * \param dependentVariablesToSave Settings for the dependent variables that are to be saved during propagation
     * (default none).
     * \param printInterval Variable indicating how often (once per printInterval_ seconds or propagation independenty
     * variable) the current state and time are to be printed to console (default never).
     */
    PropagatorSettings( const IntegratedStateType stateType,
                        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialBodyStates,
                        const boost::shared_ptr< PropagationTerminationSettings > terminationSettings,
                        const boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            boost::shared_ptr< DependentVariableSaveSettings >( ),
                        const double printInterval = TUDAT_NAN ):
        stateType_( stateType ), initialStates_( initialBodyStates ), stateSize_( initialBodyStates.rows( ) ),
        terminationSettings_( terminationSettings ), dependentVariablesToSave_( dependentVariablesToSave ),
        printInterval_( printInterval){ }

    //! Virtual destructor.
    virtual ~PropagatorSettings( ){ }

    //!T ype of state being propagated
    IntegratedStateType stateType_;

    //! Function to retrieve the initial state used as input for numerical integration
    /*!
     * Function to retrieve the initial state used as input for numerical integration
     * \return Initial state used as input for numerical integration
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialStates( )
    {
        return initialStates_;
    }

    //! Function to reset the initial state used as input for numerical integration
    /*!
     * Function to reset the initial state used as input for numerical integration
     * \param initialBodyStates New initial state used as input for numerical integration
     */
    virtual void resetInitialStates( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates )
    {
        initialStates_ = initialBodyStates;
        stateSize_ = initialStates_.rows( );
    }

    //! Get total size of the propagated state.
    /*!
     * Get total size of the propagated state.
     * \return Total size of the propagated state.
     */
    int getStateSize( )
    {
        return stateSize_;
    }

    //! Function to retrieve settings for creating the object that checks whether the propagation is finished.
    /*!
     * Function to retrieve settings for creating the object that checks whether the propagation is finished.
     * \return Settings for creating the object that checks whether the propagation is finished.
     */
    boost::shared_ptr< PropagationTerminationSettings > getTerminationSettings( )
    {
        return terminationSettings_;
    }

    //! Function to retrieve settings for the dependent variables that are to be saved during propagation (default none).
    /*!
     * Function to retrieve settings for the dependent variables that are to be saved during propagation (default none).
     * \return Settings for the dependent variables that are to be saved during propagation (default none).
     */
    boost::shared_ptr< DependentVariableSaveSettings > getDependentVariablesToSave( )
    {
        return dependentVariablesToSave_;
    }

    //! Function to retrieve how often the current state and time are to be printed to console
    /*!
     * Function to retrieve how often the current state and time are to be printed to console
     * \return Time intercal with which the current state and time are to be printed to console (default NaN, meaning never).
     */
    double getPrintInterval( )
    {
        return printInterval_;
    }


protected:

    //!  Initial state used as input for numerical integration
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialStates_;

    //! Total size of the propagated state.
    int stateSize_;

    //! Settings for creating the object that checks whether the propagation is finished.
    boost::shared_ptr< PropagationTerminationSettings > terminationSettings_;

    //! Settings for the dependent variables that are to be saved during propagation (default none).
    boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave_;

    //! Variable indicating how often (once per printInterval_ seconds or propagation independenty variable) the
    //! current state and time are to be printed to console (default never).
    double printInterval_;

};

//! Class for defining settings for propagating translational dynamics.
/*!
 *  Class for defining settings for propagating translational dynamics. The propagator defines the form of the equations of
 *  motion (i.e. Cowell, Encke, Gauss etc.). This base class can be used for Cowell propagator.
 *  Other propagators have dedicated derived class.
 */
template< typename StateScalarType = double >
class TranslationalStatePropagatorSettings: public PropagatorSettings< StateScalarType >
{
public:

    //! Constructor for generic stopping conditions.
    /*!
     * Constructor for generic stopping conditions.
     * \param centralBodies List of bodies w.r.t. which the bodies in bodiesToIntegrate_ are propagated.
     * \param accelerationsMap A map containing the list of accelerations acting on each body, identifying
     *  the body being acted on and the body acted on by an acceleration. The map has as key a string denoting
     *  the name of the body the list of accelerations, provided as the value corresponding to a key, is acting on.
     *  This map-value is again a map with string as key, denoting the body exerting the acceleration, and as value
     *  a pointer to an acceleration model.
     * \param bodiesToIntegrate List of bodies for which the translational state is to be propagated.
     * \param initialBodyStates Initial state used as input for numerical integration
     * \param terminationSettings Settings for creating the object that checks whether the propagation is finished.
     * \param propagator Type of translational state propagator to be used
     * \param dependentVariablesToSave Settings for the dependent variables that are to be saved during propagation
     * (default none).
     * \param printInterval Variable indicating how often (once per printInterval_ seconds or propagation independenty
     * variable) the current state and time are to be printed to console (default never).
     */
    TranslationalStatePropagatorSettings( const std::vector< std::string >& centralBodies,
                                          const basic_astrodynamics::AccelerationMap& accelerationsMap,
                                          const std::vector< std::string >& bodiesToIntegrate,
                                          const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
                                          const boost::shared_ptr< PropagationTerminationSettings > terminationSettings,
                                          const TranslationalPropagatorType propagator = cowell,
                                          const boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            boost::shared_ptr< DependentVariableSaveSettings >( ),
                                          const double printInterval = TUDAT_NAN ):
        PropagatorSettings< StateScalarType >( transational_state, initialBodyStates, terminationSettings,
                                               dependentVariablesToSave, printInterval ),
        centralBodies_( centralBodies ),
        accelerationsMap_( accelerationsMap ), bodiesToIntegrate_( bodiesToIntegrate ),
        propagator_( propagator ){ }

    //! Constructor for fixed propagation time stopping conditions.
    /*!
     * Constructor for fixed propagation time stopping conditions.
     * \param centralBodies List of bodies w.r.t. which the bodies in bodiesToIntegrate_ are propagated.
     * \param accelerationsMap A map containing the list of accelerations acting on each body, identifying
     *  the body being acted on and the body acted on by an acceleration. The map has as key a string denoting
     *  the name of the body the list of accelerations, provided as the value corresponding to a key, is acting on.
     *  This map-value is again a map with string as key, denoting the body exerting the acceleration, and as value
     *  a pointer to an acceleration model.
     * \param bodiesToIntegrate List of bodies for which the translational state is to be propagated.
     * \param initialBodyStates Initial state used as input for numerical integration
     * \param endTime Time at which to stop the numerical propagation
     * \param propagator Type of translational state propagator to be used
     * \param dependentVariablesToSave Settings for the dependent variables that are to be saved during propagation
     * (default none).
     * \param printInterval Variable indicating how often (once per printInterval_ seconds or propagation independenty
     * variable) the current state and time are to be printed to console (default never).
     */
    TranslationalStatePropagatorSettings( const std::vector< std::string >& centralBodies,
                                          const basic_astrodynamics::AccelerationMap& accelerationsMap,
                                          const std::vector< std::string >& bodiesToIntegrate,
                                          const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates,
                                          const double endTime,
                                          const TranslationalPropagatorType propagator = cowell,
                                          const boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            boost::shared_ptr< DependentVariableSaveSettings >( ),
                                          const double printInterval = TUDAT_NAN ):
        PropagatorSettings< StateScalarType >(
            transational_state, initialBodyStates,  boost::make_shared< PropagationTimeTerminationSettings >( endTime ),
            dependentVariablesToSave, printInterval ),
        centralBodies_( centralBodies ),
        accelerationsMap_( accelerationsMap ), bodiesToIntegrate_( bodiesToIntegrate ),
        propagator_( propagator ){ }

    //! Destructor
    ~TranslationalStatePropagatorSettings( ){ }

    //! List of bodies w.r.t. which the bodies in bodiesToIntegrate_ are propagated.
    std::vector< std::string > centralBodies_;

    //! A map containing the list of accelerations acting on each body
    /*!
     *  A map containing the list of accelerations acting on each body, identifying
     *  the body being acted on and the body acted on by an acceleration. The map has as key a string denoting
     *  the name of the body the list of accelerations, provided as the value corresponding to a key, is acting on.
     *  This map-value is again a map with string as key, denoting the body exerting the acceleration, and as value
     *  a pointer to an acceleration model.
     */
    basic_astrodynamics::AccelerationMap accelerationsMap_;

    //! List of bodies for which the translational state is to be propagated.
    std::vector< std::string > bodiesToIntegrate_;

    //! Type of translational state propagator to be used
    TranslationalPropagatorType propagator_;

};


//! Class for defining settings for propagating the mass of a body
/*!
 *  Class for defining settings for propagating the mass of a body. The body masses are propagated in their natural
 *  form (i.e. no choice of equations of motion as is the case for translational dynamics)l
 */
template< typename StateScalarType >
class MassPropagatorSettings: public PropagatorSettings< StateScalarType >
{
public:

    //! Constructor of mass state propagator settings, with single mass rate model per body.
    /*!
     * Constructor  of mass state propagator settings, with single mass rate model per body.
     * \param bodiesWithMassToPropagate List of bodies for which the mass is to be propagated.
     * \param massRateModels List of mass rate models per propagated body.
     * \param initialBodyMasses Initial masses used as input for numerical integration.
     * \param terminationSettings Settings for creating the object that checks whether the propagation is finished.
     * \param dependentVariablesToSave Settings for the dependent variables that are to be saved during propagation
     * (default none).
     * \param printInterval Variable indicating how often (once per printInterval_ seconds or propagation independenty
     * variable) the current state and time are to be printed to console (default never).
     */
    MassPropagatorSettings(
            const std::vector< std::string > bodiesWithMassToPropagate,
            const std::map< std::string, boost::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels,
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyMasses,
            const boost::shared_ptr< PropagationTerminationSettings > terminationSettings,
            const boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            boost::shared_ptr< DependentVariableSaveSettings >( ),
            const double printInterval = TUDAT_NAN ):
        PropagatorSettings< StateScalarType >( body_mass_state, initialBodyMasses, terminationSettings,
                                               dependentVariablesToSave, printInterval ),
        bodiesWithMassToPropagate_( bodiesWithMassToPropagate )
    {
        for( std::map< std::string, boost::shared_ptr< basic_astrodynamics::MassRateModel > >::const_iterator
             massRateIterator = massRateModels.begin( ); massRateIterator != massRateModels.end( ); massRateIterator++ )
        {
            massRateModels_[ massRateIterator->first ].push_back( massRateIterator->second );
        }
    }

    //! Constructor of mass state propagator settings
    /*!
     * Constructor  of mass state propagator settings
     * \param bodiesWithMassToPropagate List of bodies for which the mass is to be propagated.
     * \param massRateModels List of mass rate models per propagated body.
     * \param initialBodyMasses Initial masses used as input for numerical integration.
     * \param terminationSettings Settings for creating the object that checks whether the propagation is finished.
     * \param dependentVariablesToSave Settings for the dependent variables that are to be saved during propagation
     * (default none).
     * \param printInterval Variable indicating how often (once per printInterval_ seconds or propagation independenty
     * variable) the current state and time are to be printed to console (default never).
     */
    MassPropagatorSettings(
            const std::vector< std::string > bodiesWithMassToPropagate,
            const std::map< std::string, std::vector< boost::shared_ptr< basic_astrodynamics::MassRateModel > > >
            massRateModels,
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyMasses,
            const boost::shared_ptr< PropagationTerminationSettings > terminationSettings,
            const boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            boost::shared_ptr< DependentVariableSaveSettings >( ),
            const double printInterval = TUDAT_NAN ):
        PropagatorSettings< StateScalarType >( body_mass_state, initialBodyMasses, terminationSettings,
                                               dependentVariablesToSave, printInterval ),
        bodiesWithMassToPropagate_( bodiesWithMassToPropagate ), massRateModels_( massRateModels )
    { }

    //! List of bodies for which the mass is to be propagated.
    std::vector< std::string > bodiesWithMassToPropagate_;

    //! List of mass rate models per propagated body.
    std::map< std::string, std::vector< boost::shared_ptr< basic_astrodynamics::MassRateModel > > > massRateModels_;
};

//! Function to evaluate a floating point state-derivative function as though it was a vector state function
/*!
 *  Function to evaluate a floating point state-derivative function as though it was a vector state function.
 *  \param stateDerivativeFunction Function to compute the state derivative, as a function of current time and state.
 *  \param currentTime Time at which to evaluate the state derivative function
 *  \param currentStateVector Vector of size 1 containing the current state
 *  \return Current state derivative (as vector of size 1).
 */
template< typename StateScalarType = double, typename TimeType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > convertScalarToVectorStateFunction(
        const boost::function< StateScalarType( const TimeType, const StateScalarType ) > stateDerivativeFunction,
        const TimeType currentTime,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& currentStateVector )
{
    if( currentStateVector.rows( ) != 1 )
    {
        throw std::runtime_error( "Error, expected vector of size one when converting scalar to vector state function" );
    }
    return ( Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >( 1 )
             << stateDerivativeFunction( currentTime, currentStateVector( 0 ) ) ).finished( );

}

//! Class used to provide settings for a custom state derivative model
/*!
 *  Class used to provide settings for a custom state derivative model. The custom state derivative function has to be
 *  defined by the user and provided as input here. It may depend on the current state, and any dependent variables may
 *  be computed for other state derivative models.
 */
template< typename StateScalarType = double, typename TimeType = double >
class CustomStatePropagatorSettings: public PropagatorSettings< StateScalarType >
{
public:

    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > StateVectorType;

    //! Constructor for scalar custom state
    /*!
     * Constructor for scalar custom state
     * \param stateDerivativeFunction Function to compute the state derivative, as a function of current time and state.
     * \param initialState Initial state value of custom state
     * \param terminationSettings Settings for creating the object that checks whether the propagation is finished.
     * \param dependentVariablesToSave Settings for the dependent variables that are to be saved during propagation
     * (default none).
     * \param printInterval Variable indicating how often (once per printInterval_ seconds or propagation independenty
     * variable) the current state and time are to be printed to console (default never).
     */
    CustomStatePropagatorSettings(
            const boost::function< StateScalarType( const TimeType, const StateScalarType ) > stateDerivativeFunction,
            const StateScalarType initialState,
            const boost::shared_ptr< PropagationTerminationSettings > terminationSettings,
            const boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            boost::shared_ptr< DependentVariableSaveSettings >( ),
            const double printInterval = TUDAT_NAN ):
        PropagatorSettings< StateScalarType >(
            custom_state, ( StateVectorType( 1 ) << initialState ).finished( ), terminationSettings,
            dependentVariablesToSave, printInterval ),
        stateDerivativeFunction_( boost::bind( &convertScalarToVectorStateFunction< StateScalarType, TimeType >,
                                               stateDerivativeFunction, _1, _2 ) ), stateSize_( 1 )
    { }

    //! Constructor for vector custom state
    /*!
     * Constructor for vector custom state
     * \param stateDerivativeFunction Function to compute the state derivative, as a function of current time and state.
     * \param initialState Initial state value of custom state
     * \param terminationSettings Settings for creating the object that checks whether the propagation is finished.
     * \param dependentVariablesToSave Settings for the dependent variables that are to be saved during propagation
     * (default none).
     * \param printInterval Variable indicating how often (once per printInterval_ seconds or propagation independenty
     * variable) the current state and time are to be printed to console (default never).
     */
    CustomStatePropagatorSettings(
            const boost::function< StateVectorType( const TimeType, const StateVectorType& ) > stateDerivativeFunction,
            const Eigen::VectorXd initialState,
            const boost::shared_ptr< PropagationTerminationSettings > terminationSettings,
            const boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            boost::shared_ptr< DependentVariableSaveSettings >( ),
            const double printInterval = TUDAT_NAN ):
        PropagatorSettings< StateScalarType >( custom_state, initialState, terminationSettings,
                                               dependentVariablesToSave, printInterval ),
        stateDerivativeFunction_( stateDerivativeFunction ), stateSize_( initialState.rows( ) ){ }

    //! Destructor
    ~CustomStatePropagatorSettings( ){ }

    //! Function to compute the state derivative, as a function of current time and state.
    boost::function< StateVectorType( const TimeType, const StateVectorType& ) > stateDerivativeFunction_;

    //! Size of the state that is propagated.
    int stateSize_;

};

//! Function to retrieve the state size for a list of propagator settings.
/*!
 *  Function to retrieve the initial state for a list of propagator settings.
 *  \param propagatorSettingsList List of propagator settings (sorted by type as key). Map value provides list
 *  of propagator settings for given type.
 *  \return Vector of initial states, sorted in order of IntegratedStateType, and then in the order of the
 *  vector of PropagatorSettings of given type.
 */
template< typename StateScalarType >
int getMultiTypePropagatorStateSize(
        const std::map< IntegratedStateType, std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >&
        propagatorSettingsList )
{
    int stateSize = 0;

    // Iterate over all propagation settings and add size to list
    for( typename std::map< IntegratedStateType,
         std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >::const_iterator
         typeIterator = propagatorSettingsList.begin( ); typeIterator != propagatorSettingsList.end( ); typeIterator++ )
    {
        for( unsigned int i = 0; i < typeIterator->second.size( ); i++ )
        {
            stateSize += typeIterator->second.at( i )->getStateSize( );
        }
    }
    return stateSize;
}

//! Function to retrieve the initial state for a list of propagator settings.
/*!
 *  Function to retrieve the initial state for a list of propagator settings.
 *  \param propagatorSettingsList List of propagator settings (sorted by type as key). Map value provides list
 *  of propagator settings for given type.
 *  \return Vector of initial states, sorted in order of IntegratedStateType, and then in the order of the
 *  vector of PropagatorSettings of given type.
 */
template< typename StateScalarType >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > createCombinedInitialState(
        const std::map< IntegratedStateType, std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >&
        propagatorSettingsList )
{
    // Get total size of propagated state
    int totalSize = getMultiTypePropagatorStateSize( propagatorSettingsList );

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > combinedInitialState =
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( totalSize, 1 );

    // Iterate over all propagation settings and add to total list
    int currentIndex = 0;
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentInitialState;
    for( typename std::map< IntegratedStateType,
         std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >::const_iterator
         typeIterator = propagatorSettingsList.begin( ); typeIterator != propagatorSettingsList.end( ); typeIterator++ )
    {
        for( unsigned int i = 0; i < typeIterator->second.size( ); i++ )
        {
            currentInitialState = typeIterator->second.at( i )->getInitialStates( );
            combinedInitialState.segment( currentIndex, currentInitialState.rows( ) ) = currentInitialState;
            currentIndex += currentInitialState.rows( );
        }
    }

    return combinedInitialState;
}

//! Class for defining settings for propagating multiple types of dynamics concurrently.
template< typename StateScalarType >
class MultiTypePropagatorSettings: public PropagatorSettings< StateScalarType >
{
public:

    //! Constructor.
    /*!
     * Constructor
     * \param propagatorSettingsMap List of propagator settings to use (state type as key). List of propagator settigns
     * per type given as vector in map value.
     * \param terminationSettings Settings for creating the object that checks whether the propagation is finished.
     * \param dependentVariablesToSave Settings for the dependent variables that are to be saved during propagation
     * (default none).
     * \param printInterval Variable indicating how often (once per printInterval_ seconds or propagation independenty
     * variable) the current state and time are to be printed to console (default never).
     */
    MultiTypePropagatorSettings(
            const std::map< IntegratedStateType,
            std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > > propagatorSettingsMap,
            const boost::shared_ptr< PropagationTerminationSettings > terminationSettings,
            const boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            boost::shared_ptr< DependentVariableSaveSettings >( ),
            const double printInterval = TUDAT_NAN ):
        PropagatorSettings< StateScalarType >(
            hybrid, createCombinedInitialState< StateScalarType >( propagatorSettingsMap ),
            terminationSettings, dependentVariablesToSave, printInterval ),
        propagatorSettingsMap_( propagatorSettingsMap )
    {

    }

    //! Constructor.
    /*!
     * Constructor
     * \param propagatorSettingsVector Vector of propagator settings to use.
     * \param terminationSettings Settings for creating the object that checks whether the propagation is finished.
     * \param dependentVariablesToSave Settings for the dependent variables that are to be saved during propagation
     * (default none).
     * \param printInterval Variable indicating how often (once per printInterval_ seconds or propagation independenty
     * variable) the current state and time are to be printed to console (default never).
     */
    MultiTypePropagatorSettings(
            const std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > propagatorSettingsVector,
            const boost::shared_ptr< PropagationTerminationSettings > terminationSettings,
            const boost::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            boost::shared_ptr< DependentVariableSaveSettings >( ),
            const double printInterval = TUDAT_NAN ):
        PropagatorSettings< StateScalarType >(
            hybrid, Eigen::VectorXd::Zero( 0 ),
            terminationSettings, dependentVariablesToSave, printInterval )
    {
        for( unsigned int i = 0; i < propagatorSettingsVector.size( ); i++ )
        {
            propagatorSettingsMap_[ propagatorSettingsVector.at( i )->stateType_ ].push_back(
                        propagatorSettingsVector.at( i ) );
        }

        this->initialStates_ = createCombinedInitialState< StateScalarType >( propagatorSettingsMap_ );
    }

    //! Destructor
    ~MultiTypePropagatorSettings( ){ }

    //! Function to reset the initial state used as input for numerical integration
    /*!
     * Function to reset the initial state used as input for numerical integration
     * \param initialBodyStates New initial state used as input for numerical integration, sorted in order of
     * IntegratedStateType, and then in the order of the vector of PropagatorSettings of given type.
     */
    void resetInitialStates( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialBodyStates )
    {
        // Iterate over all propagator settings.
        int currentStartIndex = 0;
        for( typename std::map< IntegratedStateType,
             std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >::iterator
             propagatorIterator = propagatorSettingsMap_.begin( ); propagatorIterator != propagatorSettingsMap_.end( );
             propagatorIterator++ )
        {
            for( unsigned int i = 0; i < propagatorIterator->second.size( ); i++ )
            {
                // Get current state size
                int currentParameterSize = propagatorIterator->second.at( i )->getInitialStates( ).rows( );

                // Check consistency
                if( currentParameterSize + currentStartIndex > initialBodyStates.rows( ) )
                {
                    throw std::runtime_error(
                                "Error when resetting multi-type state, sizes are incompatible " );
                }

                // Reset state for current settings
                propagatorIterator->second.at( i )->resetInitialStates(
                            initialBodyStates.block( currentStartIndex, 0, currentParameterSize, 1 ) );
                currentStartIndex += currentParameterSize;

            }
        }

        // Check consistency
        if( currentStartIndex != initialBodyStates.rows( ) )
        {
            std::string errorMessage = "Error when resetting multi-type state, total size is incompatible "+
                    boost::lexical_cast< std::string >( currentStartIndex ) +
                    boost::lexical_cast< std::string >( initialBodyStates.rows( ) );
            throw std::runtime_error( errorMessage );
        }
    }

    //! List of propagator settings to use
    /*!
     * List of propagator settings to use (state type as key). List of propagator settigns
     * per type given as vector in map value.
     */

    std::map< IntegratedStateType, std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >
    propagatorSettingsMap_;

};

template< typename StateScalarType >
//! Function to retrieve the list of integrated state types and reference ids
/*!
* Function to retrieve the list of integrated state types and reference ids. For translational and rotational dynamics,
* the id refers only to the body being propagated (and the second entry of the pair is empty: ""). For proper time
* propagation, a body and a reference point may be provided, resulting in non-empty first and second pair entries.
* \param propagatorSettings Settings that are to be used for the propagation.
* \return List of integrated state types and reference ids
*/
std::map< IntegratedStateType, std::vector< std::pair< std::string, std::string > > > getIntegratedTypeAndBodyList(
        const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings )
{
    std::map< IntegratedStateType, std::vector< std::pair< std::string, std::string > > > integratedStateList;

    // Identify propagator type
    switch( propagatorSettings->stateType_ )
    {
    case hybrid:
    {
        boost::shared_ptr< MultiTypePropagatorSettings< StateScalarType > > multiTypePropagatorSettings =
                boost::dynamic_pointer_cast< MultiTypePropagatorSettings< StateScalarType > >( propagatorSettings );

        std::map< IntegratedStateType, std::vector< std::pair< std::string, std::string > > > singleTypeIntegratedStateList;


        for( typename std::map< IntegratedStateType,
             std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >::const_iterator
             typeIterator = multiTypePropagatorSettings->propagatorSettingsMap_.begin( );
             typeIterator != multiTypePropagatorSettings->propagatorSettingsMap_.end( ); typeIterator++ )
        {
            if( typeIterator->first != hybrid )
            {
                for( unsigned int i = 0; i < typeIterator->second.size( ); i++ )
                {
                    singleTypeIntegratedStateList = getIntegratedTypeAndBodyList< StateScalarType >(
                                typeIterator->second.at( i ) );

                    if( singleTypeIntegratedStateList.begin( )->first != typeIterator->first
                            || singleTypeIntegratedStateList.size( ) != 1 )
                    {
                        std::string errorMessage = "Error when making integrated state list for hybrid propagator, inconsistency encountered " +
                                boost::lexical_cast< std::string >( singleTypeIntegratedStateList.begin( )->first ) + " " +
                                boost::lexical_cast< std::string >( typeIterator->first ) + " " +
                                boost::lexical_cast< std::string >( singleTypeIntegratedStateList.size( ) ) + " " +
                                boost::lexical_cast< std::string >( singleTypeIntegratedStateList.begin( )->second.size( ) );
                        throw std::runtime_error( errorMessage );
                    }
                    else
                    {
                        for( unsigned int j = 0; j < singleTypeIntegratedStateList[ typeIterator->first ].size( ); j++ )
                        {
                            integratedStateList[ typeIterator->first ].push_back(
                                        singleTypeIntegratedStateList.begin( )->second.at( j ) );
                        }

                    }
                }
            }
            else
            {
                throw std::runtime_error( "Error when making integrated state list, cannot handle hybrid propagator inside hybrid propagator" );
            }
        }
        break;
    }
    case transational_state:
    {
        boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > >
                translationalPropagatorSettings = boost::dynamic_pointer_cast<
                TranslationalStatePropagatorSettings< StateScalarType > >( propagatorSettings );

        if( translationalPropagatorSettings == NULL )
        {
            throw std::runtime_error( "Error getting integrated state type list, translational state input inconsistent" );
        }

        // Retrieve list of integrated bodies in correct formatting.
        std::vector< std::pair< std::string, std::string > > integratedBodies;
        for( unsigned int i = 0; i < translationalPropagatorSettings->bodiesToIntegrate_.size( ); i++ )
        {
            integratedBodies.push_back( std::make_pair( translationalPropagatorSettings->bodiesToIntegrate_.at( i ), "" ) );
        }
        integratedStateList[ transational_state ] = integratedBodies;

        break;
    }
    case body_mass_state:
    {
        boost::shared_ptr< MassPropagatorSettings< StateScalarType > >
                massPropagatorSettings = boost::dynamic_pointer_cast<
                MassPropagatorSettings< StateScalarType > >( propagatorSettings );
        if( massPropagatorSettings == NULL )
        {
            throw std::runtime_error( "Error getting integrated state type list, mass state input inconsistent" );
        }

        // Retrieve list of integrated bodies in correct formatting.
        std::vector< std::pair< std::string, std::string > > integratedBodies;
        for( unsigned int i = 0; i < massPropagatorSettings->bodiesWithMassToPropagate_.size( ); i++ )
        {
            integratedBodies.push_back( std::make_pair(
                                            massPropagatorSettings->bodiesWithMassToPropagate_.at( i ), "" ) );
        }
        integratedStateList[ body_mass_state ] = integratedBodies;

        break;
    }
    case custom_state:
    {
        std::vector< std::pair< std::string, std::string > > customList;
        customList.push_back( std::make_pair( "N/A", "N/A" ) );
        integratedStateList[ custom_state ] = customList;
        break;
    }
    default:
        throw std::runtime_error( "Error, could not process integrated state type in getIntegratedTypeAndBodyList " +
                                  boost::lexical_cast< std::string >( propagatorSettings->stateType_ ) );
    }

    return integratedStateList;
}

} // namespace propagators

} // namespace tudat

namespace std
{

//! Hash for IntegratedStateType enum.
template< >
struct hash< tudat::propagators::IntegratedStateType >
{
    typedef tudat::propagators::IntegratedStateType argument_type;
    typedef size_t result_type;

    result_type operator () (const argument_type& x) const
    {
        using type = typename std::underlying_type<argument_type>::type;
        return std::hash< type >( )( static_cast< type >( x ) );
    }
};

} // namespace std

#endif // TUDAT_PROPAGATIONSETTINGS_H
