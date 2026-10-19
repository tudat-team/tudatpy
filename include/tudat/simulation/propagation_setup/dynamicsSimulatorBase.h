/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DYNAMICSSIMULATOR_BASE_H
#define TUDAT_DYNAMICSSIMULATOR_BASE_H

#include <vector>
#include <string>
#include <chrono>

#include "tudat/basics/tudatTypeTraits.h"
#include "tudat/basics/utilities.h"
#include "tudat/astro/propagators/nBodyStateDerivative.h"
#include "tudat/astro/ephemerides/frameManager.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameterSet.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/simulation/propagation_setup/setNumericallyIntegratedStates.h"
#include "tudat/astro/propagators/integrateEquations.h"
#include "tudat/simulation/propagation_setup/createStateDerivativeModel.h"
#include "tudat/simulation/propagation_setup/propagationResults.h"
#include "tudat/simulation/propagation_setup/createEnvironmentUpdater.h"
#include "tudat/simulation/propagation_setup/propagationTermination.h"
#include "tudat/astro/propagators/dynamicsStateDerivativeModel.h"
#include "tudat/math/interpolators/lagrangeInterpolator.h"
#include "tudat/simulation/propagation_setup/dependentVariablesInterface.h"

namespace tudat
{

namespace propagators
{
//! Function to get the states of a set of bodies, w.r.t. some set of central bodies, at the requested time.
/*!
 * Function to get the states of a set of bodies, w.r.t. some set of central bodies, at the requested time.
 * \param bodiesToIntegrate List of bodies for which to retrieve state.
 * \param centralBodies Origins w.r.t. which to retrieve states of bodiesToIntegrate.
 * \param bodies List of bodies to use in simulations.
 * \param initialTime Time at which to retrieve states.
 * \param frameManager OBject with which to calculate frame origin translations.
 * \return Initial state vector (with 6 Cartesian elements per body, in order of bodiesToIntegrate vector).
 */
template< typename TimeType = double, typename StateScalarType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialStatesOfBodiesFromFrameManager(
        const std::vector< std::string >& bodiesToIntegrate,
        const std::vector< std::string >& centralBodies,
        const simulation_setup::SystemOfBodies& bodies,
        const TimeType initialTime,
        const std::shared_ptr< ephemerides::ReferenceFrameManager > frameManager )
{
    // Set initial states of bodies to integrate.
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > systemInitialState =
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( bodiesToIntegrate.size( ) * 6, 1 );
    std::shared_ptr< ephemerides::Ephemeris > ephemerisOfCurrentBody;

    // Iterate over all bodies.
    for( unsigned int i = 0; i < bodiesToIntegrate.size( ); i++ )
    {
        ephemerisOfCurrentBody = bodies.at( bodiesToIntegrate.at( i ) )->getEphemeris( );

        if( !ephemerisOfCurrentBody )
        {
            throw std::runtime_error( "Could not determine initial state for body " + bodiesToIntegrate.at( i ) +
                                      " because it does not have a valid Ephemeris object." );
        }

        // Get body initial state from ephemeris
        systemInitialState.segment( i * 6, 6 ) =
                ephemerisOfCurrentBody->getTemplatedStateFromEphemeris< StateScalarType, TimeType >( initialTime );

        // Correct initial state if integration origin and ephemeris origin are not equal.
        if( centralBodies.at( i ) != ephemerisOfCurrentBody->getReferenceFrameOrigin( ) )
        {
            std::shared_ptr< ephemerides::Ephemeris > correctionEphemeris =
                    frameManager->getEphemeris( ephemerisOfCurrentBody->getReferenceFrameOrigin( ), centralBodies.at( i ) );
            systemInitialState.segment( i * 6, 6 ) -=
                    correctionEphemeris->getTemplatedStateFromEphemeris< StateScalarType, TimeType >( initialTime );
        }
    }
    return systemInitialState;
}

//! Function to get the rotational states states of a set of bodies, at the requested time.
/*!
 * Function to get the rotational states states of a set of bodies, at the requested time.
 * \param bodiesToIntegrate List of bodies for which to retrieve rotational state.
 * \param baseOrientations Reference base frame orientation
 * \param bodies List of bodies to use in simulations.
 * \param initialTime Time at which to retrieve states.
 * \return Initial rotational state vector (with 7 elements: 4 for quaternion; 3 for angular velocity) per body,
 * in order of bodiesToIntegrate vector).
 */
template< typename TimeType = double, typename StateScalarType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialRotationalStatesOfBodies( const std::vector< std::string >& bodiesToIntegrate,
                                                                                        const std::vector< std::string >& baseOrientations,
                                                                                        const simulation_setup::SystemOfBodies& bodies,
                                                                                        const TimeType initialTime )
{
    // Set initial states of bodies to integrate.
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > systemInitialState =
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( bodiesToIntegrate.size( ) * 7, 1 );
    std::shared_ptr< ephemerides::RotationalEphemeris > rotationModelOfCurrentBody;

    // Iterate over all bodies.
    for( unsigned int i = 0; i < bodiesToIntegrate.size( ); i++ )
    {
        rotationModelOfCurrentBody = bodies.at( bodiesToIntegrate.at( i ) )->getRotationalEphemeris( );

        if( !rotationModelOfCurrentBody )
        {
            throw std::runtime_error( "Could not determine initial state for body " + bodiesToIntegrate.at( i ) +
                                      " because it does not have a valid RotationalEphemeris object." );
        }

        // Get body initial state from ephemeris
        systemInitialState.segment( i * 7, 7 ) =
                rotationModelOfCurrentBody->getRotationStateVector( initialTime ).template cast< StateScalarType >( );

        // Correct initial state if integration origin and rotation model origin are not equal.
        if( baseOrientations.at( i ) != rotationModelOfCurrentBody->getBaseFrameOrientation( ) )
        {
            throw std::runtime_error( "Error, cannot get initial rotational state w.r.t. non-base frame" );
        }
    }
    return systemInitialState;
}

//! Function to get the states of a set of bodies, w.r.t. some set of central bodies, at the requested time.
/*!
 * Function to get the states of a set of bodies, w.r.t. some set of central bodies, at the requested time, creates
 * frameManager from input data.
 * \param bodiesToIntegrate List of bodies for which to retrieve state.
 * \param centralBodies Origins w.r.t. which to retrieve states of bodiesToIntegrate.
 * \param bodies List of bodies to use in simulations.
 * \param initialTime Time at which to retrieve states.
 * \return Initial state vector (with 6 Cartesian elements per body, in order of bodiesToIntegrate vector).
 */
template< typename TimeType = double, typename StateScalarType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialStatesOfBodies( const std::vector< std::string >& bodiesToIntegrate,
                                                                              const std::vector< std::string >& centralBodies,
                                                                              const simulation_setup::SystemOfBodies& bodies,
                                                                              const TimeType initialTime )
{
    // Create ReferenceFrameManager and call overloaded function.
    return getInitialStatesOfBodiesFromFrameManager< TimeType, StateScalarType >(
            bodiesToIntegrate, centralBodies, bodies, initialTime, simulation_setup::createFrameManager( bodies.getMap( ) ) );
}

//! Function to get the states of single body, w.r.t. some central body, at the requested time.
/*!
 * Function to get the states of  single body, w.r.t. some central body, at the requested time. This function creates
 * frameManager from input data to perform all required conversions.
 * \param bodyToIntegrate Body for which to retrieve state
 * \param centralBody Origin w.r.t. which to retrieve state of bodyToIntegrate.
 * \param bodies List of bodies to use in simulations.
 * \param initialTime Time at which to retrieve state.
 * \return Initial state vector of bodyToIntegrate
 */
template< typename TimeType = double, typename StateScalarType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialStateOfBody( const std::string& bodyToIntegrate,
                                                                           const std::string& centralBody,
                                                                           const simulation_setup::SystemOfBodies& bodies,
                                                                           const TimeType initialTime )
{
    return getInitialStatesOfBodies< TimeType, StateScalarType >( { bodyToIntegrate }, { centralBody }, bodies, initialTime );
}

//! Function to get the rotational states state of a body, at the requested time.
/*!
 * Function to get the rotational states state of a body, at the requested time..
 * \param bodyToIntegrate Body for which to retrieve rotational state.
 * \param baseOrientation Reference base frame orientation
 * \param bodies List of bodies to use in simulations.
 * \param initialTime Time at which to retrieve states.
 * \return Initial rotational state vector (with 7 elements: 4 for quaternion; 3 for angular velocity)
 */
template< typename TimeType = double, typename StateScalarType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialRotationalStateOfBody( const std::string& bodyToIntegrate,
                                                                                     const std::string& baseOrientation,
                                                                                     const simulation_setup::SystemOfBodies& bodies,
                                                                                     const TimeType initialTime )
{
    return getInitialRotationalStatesOfBodies< TimeType, StateScalarType >(
            std::vector< std::string >{ bodyToIntegrate }, std::vector< std::string >{ baseOrientation }, bodies, initialTime );
}

//! Function to get the state of single body, w.r.t. some central body, at a set of requested times, concatanated into one vector.
/*!
 * Function to get the states of  single body, w.r.t. some central body, at a set of requested times, concatanated into one vector.
 * This function creates frameManager from input data to perform all required conversions.
 * \param bodyToIntegrate Body for which to retrieve state
 * \param centralBodies Origins w.r.t. which to retrieve state of bodyToIntegrate (per arc).
 * \param bodies List of bodies to use in simulations.
 * \param arcStartTimes List of times at which to retrieve states.
 * \return Initial state vectosr of bodyToIntegrate at requested times.
 */
template< typename TimeType = double, typename StateScalarType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialArcWiseStateOfBody( const std::string& bodyToIntegrate,
                                                                                  const std::vector< std::string >& centralBodies,
                                                                                  const simulation_setup::SystemOfBodies& bodies,
                                                                                  const std::vector< TimeType > arcStartTimes )
{
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialStates =
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( 6 * arcStartTimes.size( ), 1 );
    for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
    {
        initialStates.block( 6 * i, 0, 6, 1 ) =
                getInitialStateOfBody< double, StateScalarType >( bodyToIntegrate, centralBodies.at( i ), bodies, arcStartTimes.at( i ) );
    }
    return initialStates;
}

//! Function to print what is inside the propagated state vector
template< typename StateScalarType = double >
void printStateVectorContent( const std::map< std::pair< int, int >, std::string > stateDescriptions, const std::string stateDescription )
{
    std::cout << "[" << stateDescription << " entries], content description" << std::endl;

    // Loop trough propagated state types and body names
    for( auto it : stateDescriptions )
    {
        int startIndex = it.first.first;
        int variableSize = it.first.second;

        // Print index at which given state type of body can be accessed
        if( variableSize == 1 )
        {
            std::cout << "[" << startIndex << "], ";
        }
        else
        {
            std::cout << "[" << startIndex << ":" << startIndex + variableSize - 1 << "], ";
        }
        std::cout << it.second << std::endl;
    }
    std::cout << std::endl;
}

//! Function to print what is inside the propagated state vector
template< typename StateScalarType = double >
void printPropagatedDependentVariableContent( std::map< std::pair< int, int >, std::string > dependentVariableIds )
{
    std::cout << "DEPENDENT VARIABLE VECTOR CONTENTS: " << std::endl;
    if( dependentVariableIds.size( ) > 0 )
    {
        std::cout << "[Vector entries], content description" << std::endl;
        for( auto it : dependentVariableIds )
        {
            int startIndex = it.first.first;
            int variableSize = it.first.second;

            // Print index at which given state type of body can be accessed
            if( variableSize == 1 )
            {
                std::cout << "[" << startIndex << "], ";
            }
            else
            {
                std::cout << "[" << startIndex << ":" << startIndex + variableSize - 1 << "], ";
            }
            std::cout << it.second << std::endl;
        }
    }
    else
    {
        std::cout << "No dependent variables have been selected." << std::endl;
    }
    std::cout << std::endl;
}

template< typename StateScalarType = double, typename TimeType = double >
static void printGenericSingleArcPostPropagationMessages(
        const std::shared_ptr< PropagationPrintSettings > printSettings,
        const std::string& propagationEndHeader,
        const std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType > > propagationResults )
{
    // Retrieve and print number of total function evaluations
    if( printSettings->printPostPropagation( ) )
    {
        std::cout << "PROPAGATION FINISHED." << std::endl;
        if( printSettings->getPrintNumberOfFunctionEvaluations( ) )
        {
            std::cout << "Total Number of Function Evaluations: " << propagationResults->getTotalNumberOfFunctionEvaluations( )
                      << std::endl;
        }
        if( printSettings->getPrintPropagationTime( ) )
        {
            std::cout << "Total propagation clock time: " << propagationResults->getTotalComputationRuntime( ) << " seconds" << std::endl;
        }
        if( printSettings->getPrintTerminationReason( ) )
        {
            std::cout << "Termination reason: " << propagationResults->getPropagationTerminationReason( )->getTerminationReasonString( )
                      << std::endl;
        }
        if( printSettings->getPrintProcessedStateData( ) )
        {
            if( propagationResults->isPropagatedAndProcessedStateEqual( ) )
            {
                std::cout << "Processed state: all state entries are propagated using default propagators, and the processed and "
                             "propagated states are identical."
                          << std::endl;
            }
            else
            {
                std::cout
                        << "Processed state: one or more state blocks are propagated using non-default propagators, the processed state is:"
                        << std::endl;
            }
            printStateVectorContent( propagationResults->getProcessedStateIds( ), "Processed state vector" );
        }
        std::cout << std::endl;
    }
    if( printSettings->printAnyOutput( ) )
    {
        std::cout << propagationEndHeader << std::endl << std::endl;
    }
}

template< typename SimulationResults, typename StateScalarType = double, typename TimeType = double >
class PropagationPrintingInterface
{
public:
    static void printSingleArcPrePropagationMessages( const std::shared_ptr< PropagationPrintSettings > printSettings,
                                                      const std::string& propagationStartHeader,
                                                      const std::shared_ptr< SimulationResults > propagationResults );

    static void printSingleArcPostPropagationMessages(
            const std::shared_ptr< PropagationPrintSettings > printSettings,
            const std::string& propagationEndHeader,
            const std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType > > propagationResults );
};

template< typename StateScalarType, typename TimeType >
class PropagationPrintingInterface< SingleArcSimulationResults< StateScalarType, TimeType >, StateScalarType, TimeType >
{
public:
    static void printSingleArcPrePropagationMessages(
            const std::shared_ptr< PropagationPrintSettings > printSettings,
            const std::string& propagationStartHeader,
            const std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType > > propagationResults )
    {
        if( printSettings->printAnyOutput( ) )
        {
            std::cout << propagationStartHeader << std::endl << std::endl;
        }
        if( printSettings->getPrintPropagatedStateData( ) )
        {
            std::cout << "PROPAGATED STATE DETAILS:" << std::endl;
            std::cout << "Propagating state vector y only, size [" << std::to_string( propagationResults->getPropagatedStateSize( ) )
                      << " x 1]" << std::endl
                      << std::endl;
            printStateVectorContent( propagationResults->getPropagatedStateIds( ), "Propagated state" );
        }
        if( printSettings->getPrintDependentVariableData( ) )
        {
            printPropagatedDependentVariableContent( propagationResults->getDependentVariableId( ) );
        }
    }

    static void printSingleArcPostPropagationMessages(
            const std::shared_ptr< PropagationPrintSettings > printSettings,
            const std::string& propagationEndHeader,
            const std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType > > propagationResults )
    {
        printGenericSingleArcPostPropagationMessages( printSettings, propagationEndHeader, propagationResults );
    }
};

template< typename StateScalarType, typename TimeType >
class PropagationPrintingInterface< SingleArcVariationalSimulationResults< StateScalarType, TimeType >, StateScalarType, TimeType >
{
public:
    static void printSingleArcPrePropagationMessages(
            const std::shared_ptr< PropagationPrintSettings > printSettings,
            const std::string& propagationStartHeader,
            const std::shared_ptr< SingleArcVariationalSimulationResults< StateScalarType, TimeType > > propagationResults )
    {
        if( printSettings->printAnyOutput( ) )
        {
            std::cout << propagationStartHeader << std::endl << std::endl;
        }
        if( printSettings->getPrintPropagatedStateData( ) )
        {
            int totalNumberOfColumns =
                    propagationResults->getStateTransitionMatrixSize( ) + propagationResults->getSensitivityMatrixSize( ) + 1;
            std::cout << "PROPAGATED STATE DETAILS:" << std::endl;
            std::cout << "Propagating state transition matrix Phi(=dx/dx0), Sensitivity matrix S(=dx/dp), and state vector y as single "
                         "matrix [Phi┊S┊y]], total size ["
                      << std::to_string( propagationResults->getDynamicsResults( )->getPropagatedStateSize( ) ) << " x "
                      << std::to_string( totalNumberOfColumns ) << "], " << std::endl;
            if( propagationResults->getDynamicsResults( )->isPropagatedAndProcessedStateEqual( ) )
            {
                std::cout << "all state entries are propagated using default propagators, and the vectors x and y are identical in your "
                             "propagation."
                          << std::endl;
            }
            else
            {
                std::cout << "note that one or more state blocks of y are propagated using non-default propagators, and the vectors x "
                             "(default formulation) and y (propagated formulation) are different in your propagation."
                          << std::endl;
            }
            printStateVectorContent( propagationResults->getDynamicsResults( )->getPropagatedStateIds( ), "State vector y" );
        }
        if( printSettings->getPrintDependentVariableData( ) )
        {
            printPropagatedDependentVariableContent( propagationResults->getDynamicsResults( )->getDependentVariableId( ) );
        }
    }

    static void printSingleArcPostPropagationMessages(
            const std::shared_ptr< PropagationPrintSettings > printSettings,
            const std::string& propagationEndHeader,
            const std::shared_ptr< SingleArcVariationalSimulationResults< StateScalarType, TimeType > > propagationResults )
    {
        printGenericSingleArcPostPropagationMessages(
                printSettings,
                propagationEndHeader,
                SingleArcResultsRetriever< SingleArcVariationalSimulationResults< StateScalarType, TimeType >, StateScalarType, TimeType >::
                        getSingleArcSimulationResults( propagationResults ) );
    }
};

//! Base class for performing full numerical integration of a dynamical system.
/*!
 *  Base class for performing full numerical integration of a dynamical system. Governing equations are set once,
 *  but can be re-integrated for different initial conditions using the same instance of the class.
 *  Derived classes define the specific kind of integration that is performed
 *  (single-arc/multi-arc/etc.)
 */
template< typename StateScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< StateScalarType, TimeType >::value, int >::type = 0 >
class DynamicsSimulator
{
public:
    //! Constructor of simulator.
    /*!
     *  Constructor of simulator, constructs integrator and object for calculating each time step of integration.
     *  \param bodies Map of bodies (with names) of all bodies in integration.
     *  \param clearNumericalSolutions Boolean to determine whether to clear the raw numerical solution member variables
     *  after propagation and resetting ephemerides (default true).
     *  \param setIntegratedResult Boolean to determine whether to automatically use the integrated results to set
     *  ephemerides (default true).
     */
    DynamicsSimulator( const simulation_setup::SystemOfBodies& bodies,
                       const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings ):
        bodies_( bodies ), propagatorSettingsBase_( propagatorSettings )
    {}

    //! Virtual destructor
    virtual ~DynamicsSimulator( ) {}

    //! This function numerically (re-)integrates the equations of motion.
    /*!
     *  This function numerically (re-)integrates the equations of motion, using the settings set through the constructor
     *  and a new initial state vector provided here. The raw results are set in the equationsOfMotionNumericalSolution_
     *  \param initialGlobalStates Initial state vector that is to be used for numerical integration.
     *  Note that this state should be in the correct frame (i.e. corresponding to centralBodies in propagatorSettings_),
     *  but not in the propagator-specific form (i.e Encke, Gauss, etc. for translational dynamics)
     * \sa SingleStateTypeDerivative::convertToOutputSolution
     */
    virtual void integrateEquationsOfMotion(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& initialGlobalStates ) = 0;

    //! Get whether the integration was completed successfully.
    /*!
     * @copybrief integrationCompletedSuccessfully
     * \return Whether the integration was completed successfully by reaching the termination condition.
     */
    virtual bool integrationCompletedSuccessfully( ) const = 0;

    //! Pure virtual function that returns the numerical result of the state propagation
    /*!
     * Pure virtual function that returns the numerical result of the state propagation.
     * \return Numerical result of the state propagation. See derived class documentation for precise contents structure.
     */
    virtual std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
    getEquationsOfMotionNumericalSolutionBase( ) = 0;

    //! Pure virtual function that returns the numerical result of the dependent variable history
    /*!
     * Pure virtual function that returns the numerical result of the dependent variable history
     * \return Numerical result of the  dependent variable history. See derived class documentation for precise contents
     *  structure.
     */
    virtual std::vector< std::map< TimeType, Eigen::VectorXd > > getDependentVariableNumericalSolutionBase( ) = 0;

    virtual std::vector< std::map< TimeType, double > > getCumulativeComputationTimeHistoryBase( ) = 0;

    //! Function to get the map of named bodies involved in simulation.
    /*!
     *  Function to get the map of named bodies involved in simulation.
     *  \return Map of named bodies involved in simulation.
     */
    simulation_setup::SystemOfBodies getSystemOfBodies( )
    {
        return bodies_;
    }

    //! Function to reset the named system of bodies.
    /*!
     *  Function to reset the named system of bodies.
     *  \param bodies The new named system of bodies.
     */
    void resetSystemOfBodies( const simulation_setup::SystemOfBodies& bodies )
    {
        bodies_ = bodies;
    }

    bool getSetIntegratedResult( )
    {
        return propagatorSettingsBase_->getOutputSettingsBase( )->getSetIntegratedResult( );
    }

    //! This function updates the environment with the numerical solution of the propagation.
    /*!
     *  This function updates the environment with the numerical solution of the propagation. For instance, it sets
     *  the propagated translational dynamics solution as the new input for the Ephemeris object of the body that was
     *  propagated. This function is pure virtual and must be implemented in the derived class.
     */
    virtual void processNumericalEquationsOfMotionSolution( ) = 0;

    virtual std::shared_ptr< SimulationResults< StateScalarType, TimeType > > getPropagationResults( ) = 0;

    std::shared_ptr< DependentVariablesInterface< TimeType > > getDependentVariablesInterface( )
    {
        return getPropagationResults( )->getDependentVariablesInterface( );
    }

protected:
    //!  Map of bodies (with names) of all bodies in integration.
    simulation_setup::SystemOfBodies bodies_;

    std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettingsBase_;
};

template< typename StateScalarType = double, typename TimeType = double >
struct PredefinedSingleArcStateDerivativeModels {
public:
    PredefinedSingleArcStateDerivativeModels(
            const std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > >& stateDerivativeModels,
            const std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap >& stateDerivativePartials,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > >& parametersToEstimate = nullptr ):
        stateDerivativeModels_( stateDerivativeModels ), stateDerivativePartials_( stateDerivativePartials ),
        parametersToEstimate_( parametersToEstimate )
    {}

    PredefinedSingleArcStateDerivativeModels( ) {}

    std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > stateDerivativeModels_;

    std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap > stateDerivativePartials_;

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate_;
};

template< typename StateScalarType, typename TimeType, int NumberOfColumns >
struct PostProcessingFunctionProvider {
    static std::function< void( Eigen::Matrix< StateScalarType, Eigen::Dynamic, NumberOfColumns >& ) > getPostProcessingFunction(
            const std::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > stateDerivateModel )
    {
        throw std::runtime_error( "Error, post-processing function can only be retrieved for single-column or dynamic size" );
        return nullptr;
    }
};

template< typename StateScalarType, typename TimeType >
struct PostProcessingFunctionProvider< StateScalarType, TimeType, 1 > {
    static std::function< void( Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& ) > getPostProcessingFunction(
            const std::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > stateDerivateModel )
    {
        return std::bind(
                &DynamicsStateDerivativeModel< TimeType, StateScalarType >::postProcessState, stateDerivateModel, std::placeholders::_1 );
        ;
    }
};

template< typename StateScalarType, typename TimeType >
struct PostProcessingFunctionProvider< StateScalarType, TimeType, Eigen::Dynamic > {
    static std::function< void( Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& ) > getPostProcessingFunction(
            const std::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > stateDerivateModel )
    {
        return std::bind( &DynamicsStateDerivativeModel< TimeType, StateScalarType >::postProcessStateAndVariationalEquations,
                          stateDerivateModel,
                          std::placeholders::_1 );
    }
};

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS
extern template class DynamicsSimulator< double, double >;
#endif

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_DYNAMICSSIMULATOR_BASE_H
