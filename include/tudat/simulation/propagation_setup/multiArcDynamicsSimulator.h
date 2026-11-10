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

#include <exception>
#include <thread>

#include "tudat/simulation/environment_setup/copySystemOfBodies.h"
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
                currentArcInitialState =
                        std::make_shared< interpolators::LagrangeInterpolator< TimeType,
                                                                               Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >,
                                                                               long double > >( initialStateInterpolationMap, 8 )
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

    //! Construct a multi-arc simulator from fully configured propagator settings.
    /*!
     * \param bodies Bodies used by the propagation.
     * \param propagatorSettings Multi-arc settings containing each arc's epoch, integrator, and output settings.
     * \param areEquationsOfMotionToBeIntegrated Whether to integrate immediately after construction.
     */
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

            // For multi-threaded propagation, create a copy of the bodies for each thread, and the state derivative models (e.g.
            // accelerations) of each arc on the copy of the thread that propagates it
            if( multiArcPropagatorSettings_->getNumberOfThreads( ) > 1 )
            {
                checkMultiArcSettingsForMultiThreadedPropagation( multiArcPropagatorSettings_ );
                unsigned int numberOfThreads = std::min( multiArcPropagatorSettings_->getNumberOfThreads( ),
                                                         static_cast< unsigned int >( singleArcSettings.size( ) ) );
                for( unsigned int i = 0; i < numberOfThreads; i++ )
                {
                    threadBodies_.push_back( simulation_setup::createSystemOfBodiesCopy< StateScalarType, TimeType >( bodies, i == 0 ) );
                }
                // Arc settings are cloned, so that the models in the settings provided to this object remain linked to the
                // original bodies
                for( unsigned int i = 0; i < singleArcSettings.size( ); i++ )
                {
                    singleArcSettings.at( i ) = singleArcSettings.at( i )->clone( );
                    singleArcSettings.at( i )->resetIntegratedStateModels( threadBodies_.at( getArcThreadIndex( i ) ) );
                }
            }

            // Create dynamics simulators
            std::vector< std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType > > > singleArcResults;
            for( unsigned int i = 0; i < singleArcSettings.size( ); i++ )
            {
                singleArcDynamicsSimulators_.push_back( std::make_shared< SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                        threadBodies_.size( ) == 0 ? bodies : threadBodies_.at( getArcThreadIndex( i ) ),
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

        unsigned int numberOfThreads = threadBodies_.size( );
        if( numberOfThreads == 0 )
        {
            // Propagate dynamics for each arc
            for( unsigned int i = 0; i < singleArcDynamicsSimulators_.size( ); i++ )
            {
                currentArcInitialState = getArcInitialState( i, initialStateProvider );
                arcInitialStateList.push_back( currentArcInitialState );

                singleArcDynamicsSimulators_.at( i )
                        ->template integrateEquationsOfMotion< typename MultiArcSimulationResults::single_arc_type >(
                                currentArcInitialState, propagationResults->getSingleArcResults( ).at( i ) );
            }
        }
        else
        {
            for( unsigned int i = 0; i < singleArcDynamicsSimulators_.size( ); i++ )
            {
                bool initialStateFromPreviousArc = false;
                arcInitialStateList.push_back( initialStateProvider->getArcInitialState( i, initialStateFromPreviousArc ) );
                if( initialStateFromPreviousArc )
                {
                    throw std::runtime_error(
                            "Error in multi-threaded multi-arc propagation, arc initial states cannot be taken from the previous arc" );
                }
            }

            // Propagate arcs (arcIndex modulo numberOfThreads) on thread with index threadIndex
            std::vector< std::exception_ptr > threadExceptions( numberOfThreads );
            std::vector< std::thread > threads;
            for( unsigned int threadIndex = 0; threadIndex < numberOfThreads; threadIndex++ )
            {
                threads.emplace_back( &MultiArcDynamicsSimulator< StateScalarType,
                                                                  TimeType >::template integrateArcsOfThread< MultiArcSimulationResults >,
                                      this,
                                      threadIndex,
                                      numberOfThreads,
                                      std::cref( arcInitialStateList ),
                                      propagationResults,
                                      std::ref( threadExceptions.at( threadIndex ) ) );
            }
            for( unsigned int threadIndex = 0; threadIndex < numberOfThreads; threadIndex++ )
            {
                threads.at( threadIndex ).join( );
            }
            for( unsigned int threadIndex = 0; threadIndex < numberOfThreads; threadIndex++ )
            {
                if( threadExceptions.at( threadIndex ) != nullptr )
                {
                    std::rethrow_exception( threadExceptions.at( threadIndex ) );
                }
            }
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

    //! Function to retrieve the copies of the bodies used by each thread (empty for single-threaded propagation)
    std::vector< simulation_setup::SystemOfBodies > getThreadBodies( )
    {
        return threadBodies_;
    }

    //! Function to retrieve the index of the thread (and of its copy of the bodies in getThreadBodies) that propagates an arc (0 for
    //! single-threaded propagation)
    unsigned int getArcThreadIndex( const unsigned int arcIndex )
    {
        return threadBodies_.size( ) == 0 ? 0 : arcIndex % threadBodies_.size( );
    }

    //! Function to set whether the integrated states are also set in the bodies of each thread after propagation (default false)
    /*!
     * Function to set whether the integrated states are also set in the bodies of each thread after propagation (default false), in
     * addition to the bodies provided to this object (which are always reset, if this is requested in the output settings).
     * \param resetIntegratedStatesInThreadBodies Boolean denoting whether the integrated states are also set in the bodies of each
     * thread
     */
    void setResetIntegratedStatesInThreadBodies( const bool resetIntegratedStatesInThreadBodies )
    {
        resetIntegratedStatesInThreadBodies_ = resetIntegratedStatesInThreadBodies;
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

                // Reset the integrated states in the (master) bodies, and if requested in the bodies of each thread
                std::vector< simulation_setup::SystemOfBodies > bodiesToReset = { bodies_ };
                if( resetIntegratedStatesInThreadBodies_ )
                {
                    bodiesToReset.insert( bodiesToReset.end( ), threadBodies_.begin( ), threadBodies_.end( ) );
                }
                std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > concatenatedResults =
                        propagationResults_->getConcatenatedEquationsOfMotionResults(
                                multiArcPropagatorSettings_->getOutputSettings( )->getClearNumericalSolutions( ) );
                for( const simulation_setup::SystemOfBodies& currentBodiesToReset : bodiesToReset )
                {
                    std::map< IntegratedStateType, std::shared_ptr< MultiArcIntegratedStateProcessor< TimeType, StateScalarType > > >
                            multiArcStateProcessors = createMultiArcIntegratedStateProcessors(
                                    currentBodiesToReset, propagationResults_->getArcStartTimes( ), singleArcIntegratedStatesProcessors );
                    for( auto itr : multiArcStateProcessors )
                    {
                        itr.second->setInterpolatorSettings(
                                multiArcPropagatorSettings_->getOutputSettings( )->getInterpolatorSettings( ) );
                        itr.second->processIntegratedMultiArcStates( concatenatedResults, propagationResults_->getArcStartTimes( ) );
                    }
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
    //! Function to propagate the arcs that are assigned to a single thread in multi-threaded propagation
    /*!
     * Function to propagate the arcs that are assigned to a single thread in multi-threaded propagation, which are the arcs with
     * index threadIndex + k * numberOfThreads (k = 0, 1, ...)
     * \param threadIndex Index of the thread
     * \param numberOfThreads Total number of threads
     * \param arcInitialStateList Initial states of all arcs
     * \param propagationResults Object in which the propagation results are stored
     * \param threadException Exception that is thrown during the propagation (none if no exception is thrown)
     */
    template< typename MultiArcSimulationResults >
    void integrateArcsOfThread(
            const unsigned int threadIndex,
            const unsigned int numberOfThreads,
            const std::vector<
                    Eigen::Matrix< StateScalarType, Eigen::Dynamic, MultiArcSimulationResults::single_arc_type::number_of_columns > >&
                    arcInitialStateList,
            const std::shared_ptr< MultiArcSimulationResults > propagationResults,
            std::exception_ptr& threadException )
    {
        try
        {
            for( unsigned int i = threadIndex; i < singleArcDynamicsSimulators_.size( ); i += numberOfThreads )
            {
                singleArcDynamicsSimulators_.at( i )
                        ->template integrateEquationsOfMotion< typename MultiArcSimulationResults::single_arc_type >(
                                arcInitialStateList.at( i ), propagationResults->getSingleArcResults( ).at( i ) );
            }
        }
        catch( ... )
        {
            threadException = std::current_exception( );
        }
    }

    //! Objects used to compute the dynamics of the sepatrate arcs
    std::vector< std::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > > singleArcDynamicsSimulators_;

    //! Propagator settings used by this objec
    std::shared_ptr< MultiArcPropagatorSettings< StateScalarType, TimeType > > multiArcPropagatorSettings_;

    std::shared_ptr< MultiArcResults > propagationResults_;

    //! Copies of the bodies used by each thread, for multi-threaded propagation (empty for single-threaded propagation)
    std::vector< simulation_setup::SystemOfBodies > threadBodies_;

    //! Boolean denoting whether the integrated states are also set in the bodies of each thread after propagation
    bool resetIntegratedStatesInThreadBodies_ = false;
};

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS
extern template class MultiArcDynamicsSimulator< double, double >;
#endif

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_MULTIARCDYNAMICSSIMULATOR_H
