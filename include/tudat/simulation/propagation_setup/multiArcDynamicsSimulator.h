/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MULTIARCDYNAMICSSIMULATOR_H
#define TUDAT_MULTIARCDYNAMICSSIMULATOR_H

#include <type_traits>

#include "tudat/simulation/propagation_setup/singleArcDynamicsSimulator.h"

namespace tudat
{

namespace propagators
{
//! Function to get a vector of initial states from a vector of propagator settings
/*!
 *  Function to get a vector of initial states from a vector of propagator settings.
 *  \param propagatorSettings List of propagator settings
 *  \return List of initial states, as retrieved from propagatorSettings list.
 */
template< typename StateScalarType = double >
std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > getInitialStatesPerArc(
        const std::vector< std::shared_ptr< PropagatorSettings< StateScalarType > > > propagatorSettings )
{
    std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > initialStatesList;
    for( unsigned int i = 0; i < propagatorSettings.size( ); i++ )
    {
        initialStatesList.push_back( propagatorSettings.at( i )->getInitialStates( ) );
    }

    return initialStatesList;
}

//! Function to get the initial state of a translational state arc from the previous state's numerical solution
/*!
 *  Function to get the initial state of a translational state arc from the previous state's numerical solution
 *  \param previousArcDynamicsSolution Numerical solution of previous arc
 *  \param currentArcInitialTime Start time of current arc
 *  \return Interpolated initial state of current arc
 */
template< typename StateScalarType = double, typename TimeType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getArcInitialStateFromPreviousArcResult(
        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& previousArcDynamicsSolution,
        const double currentArcInitialTime )
{
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentArcInitialState;
    {
        // Check if overlap exists
        if( previousArcDynamicsSolution.rbegin( )->first < currentArcInitialTime )
        {
            throw std::runtime_error( "Error when getting initial arc state from previous arc: no arc overlap" );
        }
        else
        {
            int currentIndex = 0;
            int initialTimeIndex = -1;

            std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > initialStateInterpolationMap;

            // Set sub-part of previous arc to interpolate for current arc
            for( typename std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >::const_reverse_iterator
                         previousArcIterator = previousArcDynamicsSolution.rbegin( );
                 previousArcIterator != previousArcDynamicsSolution.rend( );
                 previousArcIterator++ )
            {
                initialStateInterpolationMap[ previousArcIterator->first ] = previousArcIterator->second;
                if( initialTimeIndex < 0 )
                {
                    if( previousArcIterator->first < currentArcInitialTime )
                    {
                        initialTimeIndex = currentIndex;
                    }
                }
                else
                {
                    if( currentIndex - initialTimeIndex > 5 )
                    {
                        break;
                    }
                }
                currentIndex++;
            }

            // Interpolate to obtain initial state of current arc
            try
            {
                currentArcInitialState = std::make_shared< interpolators::LagrangeInterpolator<
                        TimeType,
                        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >,
                        std::conditional_t< std::is_same_v< StateScalarType, HighPrecisionStateScalar >,
                                            HighPrecisionStateScalar,
                                            long double > > >( initialStateInterpolationMap, 8 )
                                                 ->interpolate( currentArcInitialTime );
            }
            catch( std::runtime_error& caughtException )
            {
                throw std::runtime_error( "Error in arc initial state interpolation.\nOriginal error: " +
                                          std::string( caughtException.what( ) ) );
            }
        }
    }
    return currentArcInitialState;
}

template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< MultiArcPropagatorSettings< StateScalarType, TimeType > > validateDeprecatedMultiArcSettings(
        const std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > integratorSettings,
        const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
        const bool clearNumericalSolutions = false,
        const bool setIntegratedResult = true,
        const bool updateDependentVariableInterpolator = false )
{
    std::shared_ptr< MultiArcPropagatorSettings< StateScalarType, TimeType > > multiArcPropagatorSettings =
            std::dynamic_pointer_cast< MultiArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings );
    if( multiArcPropagatorSettings == nullptr )
    {
        throw std::runtime_error( "Error in dynamics simulator (deprecated), input must be multi-arc." );
    }

    std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > integratorSettingsList;
    if( integratorSettings.size( ) == 1 && multiArcPropagatorSettings->getSingleArcSettings( ).size( ) > 1 )
    {
        integratorSettingsList = std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > >(
                multiArcPropagatorSettings->getSingleArcSettings( ).size( ), integratorSettings.at( 0 ) );
    }
    else
    {
        integratorSettingsList = integratorSettings;
    }
    std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > independentIntegratorSettings =
            utilities::cloneDuplicatePointers( integratorSettingsList );

    if( multiArcPropagatorSettings->getSingleArcSettings( ).size( ) != independentIntegratorSettings.size( ) )
    {
        throw std::runtime_error( "Error in multi-arc dynamics simulator (deprecated), number of integrator settings is inconsistent." );
    }
    else
    {
        for( unsigned int i = 0; i < multiArcPropagatorSettings->getSingleArcSettings( ).size( ); i++ )
        {
            if( multiArcPropagatorSettings->getSingleArcSettings( ).at( i )->getIntegratorSettings( ) != nullptr &&
                independentIntegratorSettings.at( i ) != nullptr )
            {
                std::cerr << "Warning, multi-arc integrator settings, defined independently, and in propagator settings" << std::endl;
                break;
            }
            multiArcPropagatorSettings->getSingleArcSettings( ).at( i )->setIntegratorSettings( independentIntegratorSettings.at( i ) );
            if( multiArcPropagatorSettings->getSingleArcSettings( ).at( i )->getInitialTime( ) !=
                multiArcPropagatorSettings->getSingleArcSettings( ).at( i )->getInitialTime( ) )
            {
                multiArcPropagatorSettings->getSingleArcSettings( ).at( i )->resetInitialTime(
                        independentIntegratorSettings.at( i )->initialTimeDeprecated_ );
            }
        }
    }
    multiArcPropagatorSettings->getOutputSettings( )->setClearNumericalSolutions( clearNumericalSolutions );
    multiArcPropagatorSettings->getOutputSettings( )->setIntegratedResult( setIntegratedResult );
    multiArcPropagatorSettings->getOutputSettings( )->setUpdateDependentVariableInterpolator( updateDependentVariableInterpolator );

    return multiArcPropagatorSettings;
}

template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< MultiArcPropagatorSettings< StateScalarType, TimeType > > validateDeprecatedMultiArcSettings(
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
        const std::vector< double > propagationStartTimes,
        const bool clearNumericalSolutions = false,
        const bool setIntegratedResult = true,
        const bool updateDependentVariableInterpolator = false )
{
    std::shared_ptr< MultiArcPropagatorSettings< StateScalarType, TimeType > > multiArcPropagatorSettings =
            std::dynamic_pointer_cast< MultiArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings );
    if( multiArcPropagatorSettings == nullptr )
    {
        throw std::runtime_error( "Error in dynamics simulator (deprecated), input must be multi-arc." );
    }

    std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > integratorSettingsList(
            propagationStartTimes.size( ), integratorSettings );

    std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > independentIntegratorSettingsList =
            utilities::cloneDuplicatePointers( integratorSettingsList );

    for( unsigned int i = 0; i < propagationStartTimes.size( ); i++ )
    {
        multiArcPropagatorSettings->getSingleArcSettings( ).at( i )->resetInitialTime( propagationStartTimes.at( i ) );
    }

    return validateDeprecatedMultiArcSettings( independentIntegratorSettingsList,
                                               propagatorSettings,
                                               clearNumericalSolutions,
                                               setIntegratedResult,
                                               updateDependentVariableInterpolator );
}

template< typename StateScalarType = double >
class MultiArcInitialStateProvider
{
public:
    MultiArcInitialStateProvider(
            const std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& initialStatesList,
            const std::vector< std::pair< int, int > > variationalEquationsSize = std::vector< std::pair< int, int > >( ) ):
        initialStatesList_( initialStatesList ), variationalEquationsSize_( variationalEquationsSize ), updateInitialStates_( false )
    {
        useVariationalEquations_ = variationalEquationsSize_.size( ) == 0 ? false : true;
    }

    void restartPropagation( )
    {
        updateInitialStates_ = false;
    }

    bool getUpdateInitialStates( )
    {
        return updateInitialStates_;
    }

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > getArcInitialState( const int arcIndex,
                                                                                         bool& initialStateFromPreviousArc )
    {
        if( arcIndex >= static_cast< int >( initialStatesList_.size( ) ) )
        {
            throw std::runtime_error( "Error whenn getting arc initial state for arc " + std::to_string( arcIndex ) +
                                      ", index exceeds available initial states " );
        }
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > arcInitialStateFromList = initialStatesList_.at( arcIndex );
        if( linear_algebra::doesMatrixHaveNanEntries( arcInitialStateFromList ) )
        {
            initialStateFromPreviousArc = true;
            updateInitialStates_ = true;
        }
        else
        {
            initialStateFromPreviousArc = false;
        }

        Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > initialState;
        if( useVariationalEquations_ )
        {
            int numberOfRows = variationalEquationsSize_.at( arcIndex ).first;
            int numberOfColumns = variationalEquationsSize_.at( arcIndex ).second + 1;

            initialState = Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >::Zero( numberOfRows, numberOfColumns );
            initialState.block( 0, 0, numberOfRows, numberOfRows ).setIdentity( );
            initialState.block( 0, numberOfColumns - 1, numberOfRows, 1 ) = arcInitialStateFromList;
        }
        else
        {
            initialState = arcInitialStateFromList;
        }
        return initialState;
    }

private:
    const std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > initialStatesList_;

    const std::vector< std::pair< int, int > > variationalEquationsSize_;

    bool updateInitialStates_;

    bool useVariationalEquations_;
};

template< typename StateScalarType, typename TimeType, typename SimulationResults >
void checkPropagationResultsObjectConsistency(
        const std::shared_ptr< MultiArcSimulationResults< SingleArcSimulationResults, StateScalarType, TimeType > >
                originalPropagationResults,
        const std::shared_ptr< SimulationResults > comparePropagationResults )
{}

template< typename StateScalarType, typename TimeType >
void checkPropagationResultsObjectConsistency(
        const std::shared_ptr< MultiArcSimulationResults< SingleArcSimulationResults, StateScalarType, TimeType > >
                originalPropagationResults,
        const std::shared_ptr< MultiArcSimulationResults< SingleArcVariationalSimulationResults, StateScalarType, TimeType > >
                comparePropagationResults )
{
    if( originalPropagationResults->getSingleArcResults( ).size( ) != comparePropagationResults->getSingleArcResults( ).size( ) )
    {
        throw std::runtime_error(
                "Error when checking consistency of multi-arc dynamics results with variational input; results objects number of single "
                "arcs is inconsistent" );
    }

    for( unsigned int i = 0; i < originalPropagationResults->getSingleArcResults( ).size( ); i++ )
    {
        if( originalPropagationResults->getSingleArcResults( ) !=
            comparePropagationResults->getSingleArcResults( )->getSingleArcResults( ) )
        {
            throw std::runtime_error(
                    "Error when checking consistency of multi-arc dynamics results with variational input; results objects are "
                    "incosistent" );
        }
    }
}
template< typename StateScalarType, typename TimeType >
void checkPropagationResultsObjectConsistency(
        const std::shared_ptr< MultiArcSimulationResults< SingleArcSimulationResults, StateScalarType, TimeType > >
                originalPropagationResults,
        const std::shared_ptr< MultiArcSimulationResults< SingleArcSimulationResults, StateScalarType, TimeType > >
                comparePropagationResults )
{
    if( originalPropagationResults != comparePropagationResults )
    {
        if( originalPropagationResults->getSingleArcResults( ).size( ) != comparePropagationResults->getSingleArcResults( ).size( ) )
        {
            throw std::runtime_error(
                    "Error when checking consistency of multi-arc dynamics results with dynamics-only input; results objects number of "
                    "single arcs is inconsistent" );
        }
        for( unsigned int i = 0; i < originalPropagationResults->getSingleArcResults( ).size( ); i++ )
        {
            if( originalPropagationResults->getSingleArcResults( ) != comparePropagationResults->getSingleArcResults( ) )
            {
                throw std::runtime_error(
                        "Error when checking consistency of multi-arc dynamics results with dynamics-only input; results objects are "
                        "incosistent" );
            }
        }
    }
}

//! Class for performing full numerical integration of a dynamical system over multiple arcs.
/*!
 *  Class for performing full numerical integration of a dynamical system over multiple arcs, equations of motion are set up
 *  for each arc (and need not be equal for each arc). In this class, the governing equations are set once,
 *  but can be re-integrated for different initial conditions using the same instance of the class.
 */
template< typename StateScalarType = double, typename TimeType = double >
class MultiArcDynamicsSimulator : public DynamicsSimulator< StateScalarType, TimeType >
{
public:
    typedef MultiArcSimulationResults< SingleArcSimulationResults, StateScalarType, TimeType > MultiArcResults;
    using DynamicsSimulator< StateScalarType, TimeType >::bodies_;

    MultiArcDynamicsSimulator( const simulation_setup::SystemOfBodies& bodies,
                               const std::shared_ptr< MultiArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings,
                               const bool areEquationsOfMotionToBeIntegrated = true ):
        DynamicsSimulator< StateScalarType, TimeType >( bodies, propagatorSettings ), multiArcPropagatorSettings_( propagatorSettings )
    {
        if( multiArcPropagatorSettings_ == nullptr )
        {
            throw std::runtime_error( "Error when creating multi-arc dynamics simulator, input is not multi arc" );
        }
        else
        {
            std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType, TimeType > > > singleArcSettings =
                    multiArcPropagatorSettings_->getSingleArcSettings( );

            // Create dynamics simulators
            std::vector< std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType > > > singleArcResults;
            for( unsigned int i = 0; i < singleArcSettings.size( ); i++ )
            {
                singleArcDynamicsSimulators_.push_back( std::make_shared< SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                        bodies,
                        singleArcSettings.at( i ),
                        false,
                        PredefinedSingleArcStateDerivativeModels< StateScalarType, TimeType >( ),
                        true ) );
                singleArcResults.push_back( singleArcDynamicsSimulators_.at( i )->getSingleArcPropagationResults( ) );
                if( singleArcDynamicsSimulators_.at( i )->getPropagatorSettings( )->getOutputSettings( )->getSetIntegratedResult( ) )
                {
                    singleArcDynamicsSimulators_.at( i )->createAndSetIntegratedStateProcessors( );
                }
            }

            std::vector< std::shared_ptr< SingleArcDependentVariablesInterface< TimeType > > > singleArcInterfaces;
            for( unsigned int i = 0; i < singleArcSettings.size( ); i++ )
            {
                singleArcInterfaces.push_back( singleArcDynamicsSimulators_.at( i )
                                                       ->getSingleArcPropagationResults( )
                                                       ->getSingleArcDependentVariablesInterface( ) );
            }

            std::shared_ptr< MultiArcDependentVariablesInterface< TimeType > > dependentVariableInterface =
                    std::make_shared< MultiArcDependentVariablesInterface< TimeType > >(
                            singleArcInterfaces, std::vector< double >( ), std::vector< double >( ) );

            propagationResults_ = std::make_shared< MultiArcResults >( singleArcResults, dependentVariableInterface );

            // Integrate equations of motion if required.
            if( areEquationsOfMotionToBeIntegrated )
            {
                integrateEquationsOfMotion( multiArcPropagatorSettings_->getInitialStates( ) );
            }
        }
    }

    //! Constructor of multi-arc simulator for same integration settings per arc.
    /*!
     *  Constructor of multi-arc simulator for same integration settings per arc.
     *  \param bodies Map of bodies (with names) of all bodies in integration.
     *  \param integratorSettings Integrator settings for numerical integrator, used for all arcs.
     *  \param propagatorSettings Propagator settings for dynamics (must be of multi arc type)
     *  \param arcStartTimes Times at which the separate arcs start
     *  \param areEquationsOfMotionToBeIntegrated Boolean to denote whether equations of motion should be integrated at
     *  the end of the contructor or not.
     *  \param clearNumericalSolutions Boolean to determine whether to clear the raw numerical solution member variables
     *  after propagation and resetting ephemerides (default true).
     *  \param setIntegratedResult Boolean to determine whether to automatically use the integrated results to set
     *  ephemerides (default true).
     */
    MultiArcDynamicsSimulator( const simulation_setup::SystemOfBodies& bodies,
                               const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
                               const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
                               const std::vector< double > arcStartTimes,
                               const bool areEquationsOfMotionToBeIntegrated = true,
                               const bool clearNumericalSolutions = true,
                               const bool setIntegratedResult = true,
                               const bool updateDependentVariableInterpolator = false ):
        MultiArcDynamicsSimulator( bodies,
                                   validateDeprecatedMultiArcSettings( integratorSettings,
                                                                       propagatorSettings,
                                                                       arcStartTimes,
                                                                       clearNumericalSolutions,
                                                                       setIntegratedResult,
                                                                       updateDependentVariableInterpolator ),
                                   areEquationsOfMotionToBeIntegrated )
    {}

    //! Constructor of multi-arc simulator for different integration settings per arc.
    /*!
     *  Constructor of multi-arc simulator for different integration settings per arc.
     *  \param bodies Map of bodies (with names) of all bodies in integration.
     *  \param integratorSettings List of integrator settings for numerical integrator, defined per arc.
     *  \param propagatorSettings Propagator settings for dynamics (must be of multi arc type)
     *  \param areEquationsOfMotionToBeIntegrated Boolean to denote whether equations of motion should be integrated at
     *  the end of the contructor or not.
     *  \param clearNumericalSolutions Boolean to determine whether to clear the raw numerical solution member variables
     *  after propagation and resetting ephemerides (default true).
     *  \param setIntegratedResult Boolean to determine whether to automatically use the integrated results to set
     *  ephemerides (default true).
     */
    MultiArcDynamicsSimulator(
            const simulation_setup::SystemOfBodies& bodies,
            const std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > integratorSettings,
            const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const bool areEquationsOfMotionToBeIntegrated = true,
            const bool clearNumericalSolutions = true,
            const bool setIntegratedResult = true,
            const bool updateDependentVariableInterpolator = false ):
        MultiArcDynamicsSimulator( bodies,
                                   validateDeprecatedMultiArcSettings( integratorSettings,
                                                                       propagatorSettings,
                                                                       clearNumericalSolutions,
                                                                       setIntegratedResult,
                                                                       updateDependentVariableInterpolator ),
                                   areEquationsOfMotionToBeIntegrated )
    {}

    //! Destructor
    ~MultiArcDynamicsSimulator( ) {}

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > getArcInitialState(
            const int arcIndex,
            const std::shared_ptr< MultiArcInitialStateProvider< StateScalarType > > stateProvider )
    {
        bool initialStateFromPreviousArc = false;
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > currentArcInitialState =
                stateProvider->getArcInitialState( arcIndex, initialStateFromPreviousArc );
        if( initialStateFromPreviousArc )
        {
            currentArcInitialState.block( 0, currentArcInitialState.cols( ) - 1, currentArcInitialState.rows( ), 1 ) =
                    getArcInitialStateFromPreviousArcResult(
                            propagationResults_->getSingleArcResults( ).at( arcIndex - 1 )->getEquationsOfMotionNumericalSolution( ),
                            singleArcDynamicsSimulators_.at( arcIndex )->getInitialPropagationTime( ) );
        }
        return currentArcInitialState;
    }

    //! This function numerically (re-)integrates the equations of motion, using concatenated states for all arcs
    /*!
     *  This function numerically (re-)integrates the equations of motion, using the settings set through the constructor
     *  and a new initial state vector provided here. The raw results are set in the equationsOfMotionNumericalSolution_
     *  \param concatenatedInitialStates Initial state vector that is to be used for numerical integration. Note that this state
     *  should be in the correct frame (i.e. corresponding to centralBodies in propagatorSettings_), but not in the propagator-
     *  specific form (i.e Encke, Gauss, etc. for translational dynamics). The states for all arcs must be concatenated in
     *  order into a single Eigen Vector.
     */
    void integrateEquationsOfMotion( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& concatenatedInitialStates )
    {
        std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > splitInitialState;

        int currentIndex = 0;
        for( unsigned int i = 0; i < singleArcDynamicsSimulators_.size( ); i++ )
        {
            int currentSize = singleArcDynamicsSimulators_.at( i )->getPropagatorSettings( )->getConventionalStateSize( );
            splitInitialState.push_back( concatenatedInitialStates.block( currentIndex, 0, currentSize, 1 ) );
            currentIndex += currentSize;
        }

        if( currentIndex != concatenatedInitialStates.rows( ) )
        {
            throw std::runtime_error( "Error when doing multi-arc integration, input state vector size is incompatible with settings" );
        }

        integrateEquationsOfMotion( splitInitialState );
    }

    void integrate( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& concatenatedInitialStates )
    {
        integrateEquationsOfMotion( concatenatedInitialStates );
    }

    //! This function numerically (re-)integrates the equations of motion, using separate states for all arcs
    /*!
     *  This function numerically (re-)integrates the equations of motion, using the settings set through the constructor
     *  and a new initial state vector provided here. The raw results are set in the equationsOfMotionNumericalSolution_
     *  \param initialStatesList Initial state vector that is to be used for numerical integration. Note that this state should
     *  be in the correct frame (i.e. corresponding to centralBodies in propagatorSettings_), but not in the propagator-
     *  specific form (i.e Encke, Gauss, etc. for translational dynamics). The states for all stored, in order, in the input
     *  std vector.
     */
    void integrateEquationsOfMotion( const std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& initialStatesList )
    {
        integrateEquationsOfMotion< MultiArcResults >(
                propagationResults_, std::make_shared< MultiArcInitialStateProvider< StateScalarType > >( initialStatesList ) );
    }

    template< typename MultiArcSimulationResults >
    void integrateEquationsOfMotion( const std::shared_ptr< MultiArcSimulationResults > propagationResults,
                                     const std::shared_ptr< MultiArcInitialStateProvider< StateScalarType > > initialStateProvider )
    {
        checkPropagationResultsObjectConsistency< StateScalarType, TimeType, MultiArcSimulationResults >( propagationResults_,
                                                                                                          propagationResults );

        Eigen::Matrix< StateScalarType, Eigen::Dynamic, MultiArcSimulationResults::single_arc_type::number_of_columns >
                currentArcInitialState;
        std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, MultiArcSimulationResults::single_arc_type::number_of_columns > >
                arcInitialStateList;

        initialStateProvider->restartPropagation( );
        propagationResults->restartPropagation( );
        propagationResults_->restartPropagation( );

        printPrePropagationMessages( );

        // Propagate dynamics for each arc
        for( unsigned int i = 0; i < singleArcDynamicsSimulators_.size( ); i++ )
        {
            currentArcInitialState = getArcInitialState( i, initialStateProvider );
            arcInitialStateList.push_back( currentArcInitialState );

            singleArcDynamicsSimulators_.at( i )
                    ->template integrateEquationsOfMotion< typename MultiArcSimulationResults::single_arc_type >(
                            currentArcInitialState, propagationResults->getSingleArcResults( ).at( i ) );
        }

        printPostPropagationMessages( );
        propagationResults->setPropagationIsPerformed( );
        propagationResults_->setPropagationIsPerformed( );

        if( initialStateProvider->getUpdateInitialStates( ) )
        {
            std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > newInitialStates;
            for( unsigned int i = 0; i < arcInitialStateList.size( ); i++ )
            {
                newInitialStates.push_back( arcInitialStateList.at( i ).block(
                        0, arcInitialStateList.at( i ).cols( ) - 1, arcInitialStateList.at( i ).rows( ), 1 ) );
            }

            multiArcPropagatorSettings_->resetInitialStatesList( newInitialStates );
        }

        processNumericalEquationsOfMotionSolution( );
    }

    //! Function to get the list of DynamicsStateDerivativeModel objects used for each arc
    /*!
     * Function to get the list of DynamicsStateDerivativeModel objects used for each arc
     * \return List of DynamicsStateDerivativeModel objects used for each arc
     */
    std::vector< std::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > > getDynamicsStateDerivative( )
    {
        std::vector< std::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > > dynamicsStateDerivatives;
        for( unsigned int i = 0; i < singleArcDynamicsSimulators_.size( ); i++ )
        {
            dynamicsStateDerivatives.push_back( singleArcDynamicsSimulators_.at( i )->getDynamicsStateDerivative( ) );
        }
        return dynamicsStateDerivatives;
    }

    //! Function to get the list of DynamicsSimulator objects used for each arc
    /*!
     * Function to get the list of DynamicsSimulator objects used for each arc
     * \return List of DynamicsSimulator objects used for each arc
     */
    std::vector< std::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > > getSingleArcDynamicsSimulators( )
    {
        return singleArcDynamicsSimulators_;
    }

    //! This function updates the environment with the numerical solution of the propagation.
    /*!
     *  This function updates the environment with the numerical solution of the propagation. It sets
     *  the propagated dynamics solution as the new input for e.g., the ephemeris object of the boies that were
     *  propagated (for translational states).
     */
    void processNumericalEquationsOfMotionSolution( )
    {
        if( multiArcPropagatorSettings_->getOutputSettings( )->getSetIntegratedResult( ) )
        {
            try
            {
                std::map< IntegratedStateType,
                          std::vector< std::shared_ptr< SingleArcIntegratedStateProcessor< TimeType, StateScalarType > > > >
                        singleArcIntegratedStatesProcessors;

                for( unsigned int i = 0; i < singleArcDynamicsSimulators_.size( ); i++ )
                {
                    std::map< IntegratedStateType, std::shared_ptr< SingleArcIntegratedStateProcessor< TimeType, StateScalarType > > >
                            currentArcStateProcessors = singleArcDynamicsSimulators_.at( i )->getIntegratedStateProcessors( );
                    if( currentArcStateProcessors.size( ) == 0 )
                    {
                        throw std::runtime_error( "Error when resetting multi-arc states, no state processors found" );
                    }
                    for( auto itr : currentArcStateProcessors )
                    {
                        singleArcIntegratedStatesProcessors[ itr.first ].push_back( itr.second );
                    }
                }

                std::map< IntegratedStateType, std::shared_ptr< MultiArcIntegratedStateProcessor< TimeType, StateScalarType > > >
                        multiArcStateProcessors = createMultiArcIntegratedStateProcessors(
                                bodies_, propagationResults_->getArcStartTimes( ), singleArcIntegratedStatesProcessors );
                for( auto itr : multiArcStateProcessors )
                {
                    itr.second->processIntegratedMultiArcStates(
                            propagationResults_->getConcatenatedEquationsOfMotionResults(
                                    multiArcPropagatorSettings_->getOutputSettings( )->getClearNumericalSolutions( ) ),
                            propagationResults_->getArcStartTimes( ) );
                }
            }
            catch( const std::exception& caughtException )
            {
                std::cerr << "Error occured when post-processing mulyi-arc integration results, and seting integrated states in "
                             "environment, caught error is: "
                          << std::endl
                          << std::endl;
                std::cerr << caughtException.what( ) << std::endl << std::endl;
                std::cerr << "The problem may be that there is an insufficient number of data points (epochs) at which propagation results "
                             "are produced for one or more arcs"
                          << std::endl;
                if( multiArcPropagatorSettings_->getOutputSettings( )->getClearNumericalSolutions( ) )
                {
                    propagationResults_->clearSolutionMaps( );
                }
            }
        }
        else if( multiArcPropagatorSettings_->getOutputSettings( )->getClearNumericalSolutions( ) )
        {
            propagationResults_->clearSolutionMaps( );
        }

        // Reset dependent variables interface
        if( multiArcPropagatorSettings_->getOutputSettings( )->getUpdateDependentVariableInterpolator( ) )
        {
            propagationResults_->updateDependentVariableInterface( );
        }
    }

    void printPrePropagationMessages( )
    {
        if( multiArcPropagatorSettings_->getOutputSettings( )->printAnyOutput( ) )
        {
            std::cout << multiArcPropagatorSettings_->getOutputSettings( )->getPropagationStartHeader( ) << std::endl << std::endl;
        }
    }

    void printPostPropagationMessages( )
    {
        if( multiArcPropagatorSettings_->getOutputSettings( )->printAnyOutput( ) )
        {
            std::cout << multiArcPropagatorSettings_->getOutputSettings( )->getPropagationEndHeader( ) << std::endl << std::endl;
        }
    }

    std::shared_ptr< SimulationResults< StateScalarType, TimeType > > getPropagationResults( )
    {
        return propagationResults_;
    }

    std::shared_ptr< MultiArcResults > getMultiArcPropagationResults( )
    {
        return propagationResults_;
    }

    ///////////////////////////////////////////////////
    //////////////// DEPRECATED ///////////////////////
    ///////////////////////////////////////////////////

    //! Function to retrieve the current state and end times of the arcs
    /*!
     * Function to retrieve the current state and end times of the arcs
     * \return The current state and end times of the arcs
     */
    std::vector< double > getArcStartTimes( )
    {
        return propagationResults_->getArcStartTimes( );
    }

    std::vector< double > getArcEndTimes( )
    {
        return propagationResults_->getArcEndTimes( );
    }

    //! Function to return the numerical solution to the equations of motion.
    /*!
     *  Function to return the numerical solution to the equations of motion for last numerical integration. Each vector entry
     *  denotes one arc. Key of map denotes time, values are full propagated state vectors.
     *  \return List of maps of history of numerically integrated states.
     */
    std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > getEquationsOfMotionNumericalSolution( )
    {
        return propagationResults_->getConcatenatedEquationsOfMotionResults( );
    }

    //! Function to return the numerical solution of the dependent variables
    /*!
     *  Function to return the numerical solution of the dependent variables for last numerical integration. Each vector entry
     *  denotes one arc. Key of map denotes time, values are dependent variable vectors
     *  \return List of maps of dependent variable history
     */
    std::vector< std::map< TimeType, Eigen::VectorXd > > getDependentVariableHistory( )
    {
        return propagationResults_->getConcatenatedDependentVariableResults( );
    }

    std::vector< std::map< TimeType, double > > getCumulativeComputationTimeHistory( )
    {
        return propagationResults_->getConcatenatedCumulativeComputationTimeHistory( );
    }

    //! Function to return the numerical solution to the equations of motion (base class interface).
    /*!
     *  Function to return the numerical solution to the equations of motion for last numerical integration. Each vector entry
     *  denotes one arc. Key of map denotes time, values are full propagated state vectors.
     *  \return List of maps of history of numerically integrated states.
     */
    std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > getEquationsOfMotionNumericalSolutionBase( )
    {
        return getEquationsOfMotionNumericalSolution( );
    }

    //! Function to return the numerical solution of the dependent variables (base class interface)
    /*!
     *  Function to return the numerical solution of the dependent variables for last numerical integration. Each vector entry
     *  denotes one arc. Key of map denotes time, values are dependent variable vectors
     *  \return List of maps of dependent variable history
     */
    std::vector< std::map< TimeType, Eigen::VectorXd > > getDependentVariableNumericalSolutionBase( )
    {
        return getDependentVariableHistory( );
    }

    std::vector< std::map< TimeType, double > > getCumulativeComputationTimeHistoryBase( )
    {
        return getCumulativeComputationTimeHistory( );
    }

    //! Get whether the integration was completed successfully.
    /*!
     * @copybrief integrationCompletedSuccessfully
     * \return Whether the integration was completed successfully by reaching the termination condition.
     */
    virtual bool integrationCompletedSuccessfully( ) const
    {
        return propagationResults_->integrationCompletedSuccessfully( );
    }

    std::vector< std::shared_ptr< PropagationTerminationDetails > > getPropagationTerminationReasons( )
    {
        return propagationResults_->getConcatenatedTerminationReasons( );
    }

    ///////////////////////////////////////////////////
    //////////////// END DEPRECATED ///////////////////
    ///////////////////////////////////////////////////

protected:
    //! Objects used to compute the dynamics of the sepatrate arcs
    std::vector< std::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > > singleArcDynamicsSimulators_;

    //! Propagator settings used by this objec
    std::shared_ptr< MultiArcPropagatorSettings< StateScalarType, TimeType > > multiArcPropagatorSettings_;

    std::shared_ptr< MultiArcResults > propagationResults_;
};

template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< HybridArcPropagatorSettings< StateScalarType, TimeType > > validateDeprecatedHybridArcSettings(
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > singleArcIntegratorSettings,
        const std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > multiArcIntegratorSettings,
        const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
        const bool clearNumericalSolutions = true,
        const bool setIntegratedResult = true,
        const bool updateDependentVariableInterpolator = false )
{
    std::shared_ptr< HybridArcPropagatorSettings< StateScalarType, TimeType > > hybridArcPropagatorSettings =
            std::dynamic_pointer_cast< HybridArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings );
    if( hybridArcPropagatorSettings == nullptr )
    {
        throw std::runtime_error( "Error in dynamics simulator (deprecated), input must be hybrid-arc." );
    }

    validateDeprecatedSingleArcSettings< StateScalarType, TimeType >( singleArcIntegratorSettings,
                                                                      hybridArcPropagatorSettings->getSingleArcPropagatorSettings( ) );
    validateDeprecatedMultiArcSettings< StateScalarType, TimeType >( multiArcIntegratorSettings,
                                                                     hybridArcPropagatorSettings->getMultiArcPropagatorSettings( ),
                                                                     clearNumericalSolutions,
                                                                     setIntegratedResult,
                                                                     updateDependentVariableInterpolator );

    hybridArcPropagatorSettings->getOutputSettingsWithCheck( )->setClearNumericalSolutions( clearNumericalSolutions );
    hybridArcPropagatorSettings->getOutputSettingsWithCheck( )->setIntegratedResult( setIntegratedResult );

    return hybridArcPropagatorSettings;
}

template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< HybridArcPropagatorSettings< StateScalarType, TimeType > > validateDeprecatedHybridArcSettings(
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > singleArcIntegratorSettings,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > multiArcIntegratorSettings,
        const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
        const std::vector< double > arcStartTimes,
        const bool clearNumericalSolutions = true,
        const bool setIntegratedResult = true,
        const bool updateDependentVariableInterpolator = false )
{
    std::shared_ptr< HybridArcPropagatorSettings< StateScalarType, TimeType > > hybridArcPropagatorSettings =
            std::dynamic_pointer_cast< HybridArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings );
    if( hybridArcPropagatorSettings == nullptr )
    {
        throw std::runtime_error( "Error in dynamics simulator (deprecated), input must be hybrid-arc." );
    }
    std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > integratorSettingsList(
            arcStartTimes.size( ), multiArcIntegratorSettings );

    std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > independentIntegratorSettingsList =
            utilities::cloneDuplicatePointers( integratorSettingsList );
    for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
    {
        hybridArcPropagatorSettings->getMultiArcPropagatorSettings( )->getSingleArcSettings( ).at( i )->resetInitialTime(
                arcStartTimes.at( i ) );
    }
    hybridArcPropagatorSettings->getSingleArcPropagatorSettings( )->resetInitialTime( singleArcIntegratorSettings->initialTimeDeprecated_ );
    return validateDeprecatedHybridArcSettings( singleArcIntegratorSettings,
                                                independentIntegratorSettingsList,
                                                propagatorSettings,
                                                clearNumericalSolutions,
                                                setIntegratedResult,
                                                updateDependentVariableInterpolator );
}

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS
extern template class MultiArcDynamicsSimulator< double, double >;
#endif

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_MULTIARCDYNAMICSSIMULATOR_H
