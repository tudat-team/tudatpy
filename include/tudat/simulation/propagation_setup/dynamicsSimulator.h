/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DYNAMICSSIMULATOR_H
#define TUDAT_DYNAMICSSIMULATOR_H

#include <vector>
#include <string>
#include <chrono>



#include "tudat/basics/tudatTypeTraits.h"
#include "tudat/basics/utilities.h"
#include "tudat/astro/propagators/nBodyStateDerivative.h"
#include "tudat/astro/ephemerides/frameManager.h"
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
    for( unsigned int i = 0; i < bodiesToIntegrate.size( ) ; i++ )
    {
        ephemerisOfCurrentBody = bodies.at( bodiesToIntegrate.at( i ) )->getEphemeris( );

        if ( !ephemerisOfCurrentBody )
        {
            throw std::runtime_error( "Could not determine initial state for body " + bodiesToIntegrate.at( i ) +
                                      " because it does not have a valid Ephemeris object." );
        }

        // Get body initial state from ephemeris
        systemInitialState.segment( i * 6 , 6 ) = ephemerisOfCurrentBody->getTemplatedStateFromEphemeris<
                StateScalarType, TimeType >( initialTime );

        // Correct initial state if integration origin and ephemeris origin are not equal.
        if( centralBodies.at( i ) != ephemerisOfCurrentBody->getReferenceFrameOrigin( ) )
        {
            std::shared_ptr< ephemerides::Ephemeris > correctionEphemeris =
                    frameManager->getEphemeris( ephemerisOfCurrentBody->getReferenceFrameOrigin( ), centralBodies.at( i ) );
            systemInitialState.segment( i * 6 , 6 ) -= correctionEphemeris->getTemplatedStateFromEphemeris<
                    StateScalarType, TimeType >( initialTime );
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
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialRotationalStatesOfBodies(
        const std::vector< std::string >& bodiesToIntegrate,
        const std::vector< std::string >& baseOrientations,
        const simulation_setup::SystemOfBodies& bodies,
        const TimeType initialTime )
{
    // Set initial states of bodies to integrate.
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > systemInitialState =
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( bodiesToIntegrate.size( ) * 7, 1 );
    std::shared_ptr< ephemerides::RotationalEphemeris > rotationModelOfCurrentBody;

    // Iterate over all bodies.
    for( unsigned int i = 0; i < bodiesToIntegrate.size( ) ; i++ )
    {
        rotationModelOfCurrentBody = bodies.at( bodiesToIntegrate.at( i ) )->getRotationalEphemeris( );

        if ( ! rotationModelOfCurrentBody )
        {
            throw std::runtime_error( "Could not determine initial state for body " + bodiesToIntegrate.at( i ) +
                                      " because it does not have a valid RotationalEphemeris object." );
        }

        // Get body initial state from ephemeris
        systemInitialState.segment( i * 7 , 7 ) = rotationModelOfCurrentBody->getRotationStateVector(
                    initialTime ).template cast< StateScalarType >( );

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
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialStatesOfBodies(
        const std::vector< std::string >& bodiesToIntegrate,
        const std::vector< std::string >& centralBodies,
        const simulation_setup::SystemOfBodies& bodies,
        const TimeType initialTime )
{
    // Create ReferenceFrameManager and call overloaded function.
    return getInitialStatesOfBodiesFromFrameManager< TimeType, StateScalarType >(
                bodiesToIntegrate, centralBodies, bodies, initialTime,
                simulation_setup::createFrameManager( bodies.getMap( ) ) );
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
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialStateOfBody(
        const std::string& bodyToIntegrate,
        const std::string& centralBody,
        const simulation_setup::SystemOfBodies& bodies,
        const TimeType initialTime )
{
    return getInitialStatesOfBodies< TimeType, StateScalarType >(
    { bodyToIntegrate }, { centralBody }, bodies, initialTime );
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
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialRotationalStateOfBody(
        const std::string& bodyToIntegrate,
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
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialArcWiseStateOfBody(
        const std::string& bodyToIntegrate,
        const std::vector< std::string >& centralBodies,
        const simulation_setup::SystemOfBodies& bodies,
        const std::vector< TimeType > arcStartTimes )
{
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialStates = Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero(
                6 * arcStartTimes.size( ), 1 );
    for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
    {
        initialStates.block( 6 * i, 0, 6, 1 ) = getInitialStateOfBody< double, StateScalarType >(
                    bodyToIntegrate, centralBodies.at( i ), bodies, arcStartTimes.at( i ) );
    }
    return initialStates;
}


//! Function to print what is inside the propagated state vector
template< typename StateScalarType = double >
void printStateVectorContent(
        const std::map< std::pair< int, int >, std::string > stateDescriptions,
        const std::string stateDescription )
{
    std::cout <<"["<<stateDescription<<" entries], content description" << std::endl;
    
    // Loop trough propagated state types and body names
    for ( auto it : stateDescriptions )
    {
        int startIndex = it.first.first;
        int variableSize = it.first.second;

        // Print index at which given state type of body can be accessed
        if (variableSize == 1)
        {
            std::cout << "[" << startIndex << "], ";
        }
        else
        {
            std::cout << "[" << startIndex << ":" << startIndex+variableSize-1<< "], ";
        }
        std::cout<<it.second<<std::endl;
    }
    std::cout<<std::endl;
}


//! Function to print what is inside the propagated state vector
template< typename StateScalarType = double >
void printPropagatedDependentVariableContent (
        std::map< std::pair< int, int >, std::string > dependentVariableIds )
{
    std::cout << "DEPENDENT VARIABLE VECTOR CONTENTS: " << std::endl;
    if( dependentVariableIds.size( ) > 0 )
    {
        std::cout<< "[Vector entries], content description" << std::endl;
        for ( auto it : dependentVariableIds )
        {
            int startIndex = it.first.first;
            int variableSize = it.first.second;

            // Print index at which given state type of body can be accessed
            if (variableSize == 1)
            {
                std::cout << "[" << startIndex << "], ";
            }
            else
            {
                std::cout << "[" << startIndex << ":" << startIndex+variableSize-1<< "], ";
            }
            std::cout<<it.second<<std::endl;
        }
    }
    else
    {
        std::cout<<"No dependent variables have been selected."<<std::endl;
    }
    std::cout<<std::endl;
}

        template< typename StateScalarType = double, typename TimeType = double >
        static void printGenericSingleArcPostPropagationMessages(
                const std::shared_ptr< PropagationPrintSettings > printSettings,
                const std::string& propagationEndHeader,
                const std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType > > propagationResults )
        {
            // Retrieve and print number of total function evaluations
            if ( printSettings->printPostPropagation( ) )
            {
                std::cout << "PROPAGATION FINISHED."<<std::endl;
                if( printSettings->getPrintNumberOfFunctionEvaluations( ) )
                {
                    std::cout << "Total Number of Function Evaluations: "
                              << propagationResults->getTotalNumberOfFunctionEvaluations( ) << std::endl;
                }
                if( printSettings->getPrintPropagationTime( ) )
                {
                    std::cout << "Total propagation clock time: "
                              << propagationResults->getTotalComputationRuntime( )<<" seconds"<<std::endl;
                }
                if( printSettings->getPrintTerminationReason( ) )
                {
                    std::cout << "Termination reason: "<<propagationResults->getPropagationTerminationReason()->getTerminationReasonString( )<<std::endl;
                }
                if( printSettings->getPrintProcessedStateData( ) )
                {
                    if( propagationResults->isPropagatedAndProcessedStateEqual( ) )
                    {
                        std::cout<<"Processed state: all state entries are propagated using default propagators, and the processed and propagated states are identical."<<std::endl;
                    }
                    else
                    {
                        std::cout<<"Processed state: one or more state blocks are propagated using non-default propagators, the processed state is:"<<std::endl;

                    }
                    printStateVectorContent( propagationResults->getProcessedStateIds( ),"Processed state vector" );
                }
                std::cout<<std::endl;
            }
            if( printSettings->printAnyOutput( ) )
            {
                std::cout<<propagationEndHeader<<std::endl<<std::endl;
            }
        }

        template< typename SimulationResults, typename StateScalarType = double, typename TimeType = double >
        class PropagationPrintingInterface
        {
        public:

            static void printSingleArcPrePropagationMessages(
                    const std::shared_ptr< PropagationPrintSettings > printSettings,
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
                    std::cout<<propagationStartHeader<<std::endl<<std::endl;
                }
                if( printSettings->getPrintPropagatedStateData( ) )
                {
                    std::cout<<"PROPAGATED STATE DETAILS:"<<std::endl;
                    std::cout<<"Propagating state vector y only, size ["<<std::to_string( propagationResults->getPropagatedStateSize( ) )<<" x 1]"<<std::endl<<std::endl;
                    printStateVectorContent( propagationResults->getPropagatedStateIds( ),"Propagated state" );
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
                    std::cout<<propagationStartHeader<<std::endl<<std::endl;
                }
                if( printSettings->getPrintPropagatedStateData( ) )
                {
                    int totalNumberOfColumns = propagationResults->getStateTransitionMatrixSize( ) + propagationResults->getSensitivityMatrixSize( ) + 1;
                    std::cout<<"PROPAGATED STATE DETAILS:"<<std::endl;
                    std::cout<<"Propagating state transition matrix Phi(=dx/dx0), Sensitivity matrix S(=dx/dp), and state vector y as single matrix [Phi┊S┊y]], total size ["<<
                        std::to_string( propagationResults->getDynamicsResults( )->getPropagatedStateSize( ) )<<" x "<<std::to_string( totalNumberOfColumns )<< "], "<<std::endl;
                    if( propagationResults->getDynamicsResults( )->isPropagatedAndProcessedStateEqual( ) )
                    {
                        std::cout<<"all state entries are propagated using default propagators, and the vectors x and y are identical in your propagation."<<std::endl;
                    }
                    else
                    {
                        std::cout<<"note that one or more state blocks of y are propagated using non-default propagators, and the vectors x (default formulation) and y (propagated formulation) are different in your propagation."<<std::endl;
                    }
                    printStateVectorContent( propagationResults->getDynamicsResults( )->getPropagatedStateIds( ),"State vector y" );
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
                        printSettings, propagationEndHeader,
                        SingleArcResultsRetriever< SingleArcVariationalSimulationResults< StateScalarType, TimeType >, StateScalarType, TimeType >::getSingleArcSimulationResults( propagationResults )  );
            }
        };



//! Base class for performing full numerical integration of a dynamical system.
/*!
 *  Base class for performing full numerical integration of a dynamical system. Governing equations are set once,
 *  but can be re-integrated for different initial conditions using the same instance of the class.
 *  Derived classes define the specific kind of integration that is performed
 *  (single-arc/multi-arc/etc.)
 */
template< typename StateScalarType = double, typename TimeType = double,
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
    DynamicsSimulator(
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings ):
        bodies_( bodies ),
        propagatorSettingsBase_( propagatorSettings )
    {
    }

    //! Virtual destructor
    virtual ~DynamicsSimulator( ) { }

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
std::shared_ptr< SingleArcPropagatorSettings< StateScalarType, TimeType > > validateDeprecatedSingleArcSettings(
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
        const bool clearNumericalSolutions = false,
        const bool setIntegratedResult = false,
        const bool printNumberOfFunctionEvaluations = false,
        const bool printDependentVariableData = false,
        const bool printStateData = false,
        const bool updateDependentVariableInterpolator = false )
{
    std::shared_ptr< SingleArcPropagatorSettings< StateScalarType, TimeType > > singleArcPropagatorSettings =
            std::dynamic_pointer_cast< SingleArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings );
    if( singleArcPropagatorSettings == nullptr )
    {
        throw std::runtime_error( "Error in dynamics simulator (deprecated), input must be single-arc." );
    }
    else
    {
        if( integratorSettings != nullptr && singleArcPropagatorSettings->getIntegratorSettings( ) != nullptr )
        {
            std::cerr<<"Warning when processing deprecated propagator settitngs, integrator settings defined independently, and in propagator settings"<<std::endl;
        }

        if( integratorSettings != nullptr )
        {
            singleArcPropagatorSettings->resetInitialTime( integratorSettings->initialTimeDeprecated_ );
            if( singleArcPropagatorSettings->getInitialTime( ) != singleArcPropagatorSettings->getInitialTime( ) )
            {
                std::cerr<<"Warning when processing deprecated propagator settitngs, initial propagation time is NaN"<<std::endl;
            }
        }

        singleArcPropagatorSettings->getOutputSettings( )->setClearNumericalSolutions( clearNumericalSolutions );
        singleArcPropagatorSettings->getOutputSettings( )->setIntegratedResult( setIntegratedResult );
        singleArcPropagatorSettings->getOutputSettings( )->setUpdateDependentVariableInterpolator( updateDependentVariableInterpolator );
        singleArcPropagatorSettings->getOutputSettings( )->getPrintSettings( )->reset(
                    printNumberOfFunctionEvaluations, printDependentVariableData,
                    singleArcPropagatorSettings->getOutputSettings( )->getPrintSettings( )->getResultsPrintFrequencyInSeconds( ), 0,
                    false, false, printStateData, false, false, false );

        singleArcPropagatorSettings->setIntegratorSettings( integratorSettings );
    }


    return singleArcPropagatorSettings;
}

template< typename StateScalarType = double, typename TimeType = double >
struct PredefinedSingleArcStateDerivativeModels
{
public:
    PredefinedSingleArcStateDerivativeModels(
            const std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > >& stateDerivativeModels,
            const std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap >& stateDerivativePartials ):
        stateDerivativeModels_( stateDerivativeModels ), stateDerivativePartials_( stateDerivativePartials ){ }

    PredefinedSingleArcStateDerivativeModels( ){ }

    std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > > stateDerivativeModels_;

    std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap > stateDerivativePartials_;
};


template< typename StateScalarType, typename TimeType, int NumberOfColumns >
struct PostProcessingFunctionProvider
{
    static std::function< void( Eigen::Matrix< StateScalarType, Eigen::Dynamic, NumberOfColumns >& ) > getPostProcessingFunction(
            const std::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > stateDerivateModel )
    {
        throw std::runtime_error( "Error, post-processing function can only be retrieved for single-column or dynamic size" );
        return nullptr;
    }
};

template< typename StateScalarType, typename TimeType >
struct PostProcessingFunctionProvider< StateScalarType, TimeType, 1 >
{
    static std::function< void( Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& ) > getPostProcessingFunction(
            const std::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > stateDerivateModel )
    {
        return std::bind(
                &DynamicsStateDerivativeModel< TimeType, StateScalarType >::postProcessState,
                stateDerivateModel, std::placeholders::_1 );;
    }
};

template< typename StateScalarType, typename TimeType >
struct PostProcessingFunctionProvider< StateScalarType, TimeType, Eigen::Dynamic >
{
    static std::function< void( Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& ) > getPostProcessingFunction(
            const std::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > stateDerivateModel )
    {
        return std::bind(
                &DynamicsStateDerivativeModel< TimeType, StateScalarType >::postProcessStateAndVariationalEquations,
                stateDerivateModel, std::placeholders::_1 );
    }
};

//!Class for performing full numerical integration of a dynamical system in a single arc.
/*!
 *  Class for performing full numerical integration of a dynamical system in a single arc, i.e. the equations of motion
 *  have a single initial time, and are propagated once for the full prescribed time interval. This is in contrast to
 *  multi-arc dynamics, where the time interval si cut into pieces. In this class, the governing equations are set once,
 *  but can be re-integrated for different initial conditions using the same instance of the class.
 */
template< typename StateScalarType = double, typename TimeType = double >
class SingleArcDynamicsSimulator: public DynamicsSimulator< StateScalarType, TimeType >
{
public:


    SingleArcDynamicsSimulator(
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< SingleArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings,
            const bool areEquationsOfMotionToBeIntegrated = true,
            const PredefinedSingleArcStateDerivativeModels< StateScalarType, TimeType >& predefinedStateDerivativeModels =
            PredefinedSingleArcStateDerivativeModels< StateScalarType, TimeType >( ) ):
        DynamicsSimulator< StateScalarType, TimeType >(
            bodies, propagatorSettings ),
        propagatorSettings_( propagatorSettings )
    {
        // Check consistency of input settings
        if( propagatorSettings == nullptr )
        {
            throw std::runtime_error( "Error in dynamics simulator, propagator settings not defined." );
        }
        else if( std::dynamic_pointer_cast< SingleArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) == nullptr )
        {
            throw std::runtime_error( "Error in dynamics simulator, input must be single-arc." );
        }
        else
        {
            // Retrieve output and integrator settings TODO: no need to set as member variables; can just retrieve from propagatorSettings_
            outputSettings_ = propagatorSettings_->getOutputSettingsWithCheck( );
            integratorSettings_ = propagatorSettings_->getIntegratorSettings( );
        }
        if( integratorSettings_ == nullptr )
        {
            throw std::runtime_error( "Error in dynamics simulator, integrator settings not defined." );
        }
        checkPropagatedStatesFeasibility( propagatorSettings_, bodies_ );

        // Create objects that reset the environment (e.g. ephemerides) after propagation is required
        if( propagatorSettings_->getOutputSettings( )->getSetIntegratedResult( ) )
        {
            createAndSetIntegratedStateProcessors( );
        }

        // Create object that updates the environment during propagation
        try
        {
            environmentUpdater_ = createEnvironmentUpdaterForDynamicalEquations< StateScalarType, TimeType >(
                    propagatorSettings_, bodies_ );
        }
        catch( const std::runtime_error& error )
        {
            throw std::runtime_error( "Error when creating environment updater: "  + std::string( error.what( ) ) );
        }

        // Create object that calculates the complete state derivatives
        if( predefinedStateDerivativeModels.stateDerivativeModels_.size( ) == 0 )
        {
            dynamicsStateDerivative_ = std::make_shared< DynamicsStateDerivativeModel< TimeType, StateScalarType > >(
                        createStateDerivativeModels< StateScalarType, TimeType >(
                            propagatorSettings_, bodies_, propagatorSettings_->getInitialTime( ) ),
                        std::bind( &EnvironmentUpdater< StateScalarType, TimeType >::updateEnvironment,
                                     environmentUpdater_, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3 ) );
        }
        else
        {
            dynamicsStateDerivative_ = std::make_shared< DynamicsStateDerivativeModel< TimeType, StateScalarType > >(
                        predefinedStateDerivativeModels.stateDerivativeModels_,
                        std::bind( &EnvironmentUpdater< StateScalarType, TimeType >::updateEnvironment,
                                     environmentUpdater_, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3 ) );
        }
        stateDerivativeFunction_ =
                std::bind( &DynamicsStateDerivativeModel< TimeType, StateScalarType >::computeStateDerivative,
                           dynamicsStateDerivative_, std::placeholders::_1, std::placeholders::_2 );

        // Create object that determines if the propagation is to be terminated
        propagationTerminationCondition_ = createPropagationTerminationConditions(
                    propagatorSettings_->getTerminationSettings( ), bodies_,
                    integratorSettings_->initialTimeStep_, dynamicsStateDerivative_->getStateDerivativeModels( ),
                    predefinedStateDerivativeModels.stateDerivativePartials_ );

        sequentialPropagation_ = true;
        if ( propagationTerminationCondition_->getTerminationType( ) == non_sequential_stopping_condition )
        {
            sequentialPropagation_ = false;
            if ( integratorSettings_->initialTimeStep_ < 0.0 )
            {
                throw std::runtime_error( "Error when using non-sequential propagation, the initial integrator time step must be positive (first provided for forward leg, "
                                          "conversion to negative time step for backward leg is automatic)." );
            }
        }

        std::map< IntegratedStateType, std::vector< std::tuple< std::string, std::string, PropagatorType > > > integratedStateAndBodyList =
                getIntegratedTypeAndBodyList( propagatorSettings_ );

        std::map< std::pair< int, int >, std::string > dependentVariableIds_;
        std::map< std::pair< int, int >, std::shared_ptr< SingleDependentVariableSaveSettings > > orderedDependentVariableSettings_;

        // Create functions that compute the dependent variables
        if( propagatorSettings_->getDependentVariablesToSave( ).size( ) > 0 )
        {
            std::pair< std::function< Eigen::VectorXd( ) >, std::map< std::pair< int, int >, std::string > > dependentVariableData =
                    createDependentVariableListFunction< TimeType, StateScalarType >(
                        propagatorSettings_->getDependentVariablesToSave( ), bodies_,
                        orderedDependentVariableSettings_,
                        dynamicsStateDerivative_->getStateDerivativeModels( ),
                        predefinedStateDerivativeModels.stateDerivativePartials_ );
            dependentVariablesFunctions_ = dependentVariableData.first;
            dependentVariableIds_ = dependentVariableData.second;
        }

        // Create object that will contain and process the propagation results
        std::shared_ptr< SingleArcDependentVariablesInterface< TimeType > > dependentVariableInterface =
            std::make_shared< SingleArcDependentVariablesInterface< TimeType > >(
                std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::VectorXd > >( ),
                propagatorSettings_->getDependentVariablesToSave( ),
                dependentVariableIds_,
                orderedDependentVariableSettings_ );

        propagationResults_= std::make_shared< SingleArcSimulationResults< StateScalarType, TimeType > >(
                    integratedStateAndBodyList, propagatorSettings_->getOutputSettingsWithCheck( ),
                    std::bind( &DynamicsStateDerivativeModel< TimeType, StateScalarType >::convertNumericalStateSolutionsToOutputSolutions,
                               dynamicsStateDerivative_,
                               std::placeholders::_1, std::placeholders::_2 ), dependentVariableInterface, sequentialPropagation_ ) ;

        // Integrate equations of motion if required.
        if( areEquationsOfMotionToBeIntegrated )
        {
            integrateEquationsOfMotion( propagatorSettings_->getInitialStates( ) );
        }
    }

    using DynamicsSimulator< StateScalarType, TimeType >::bodies_;

    //! Constructor of simulator.
    /*!
     *  Constructor of simulator, constructs integrator and object for calculating each time step of integration.
     *  \param bodies Map of bodies (with names) of all bodies in integration.
     *  \param integratorSettings Settings for numerical integrator.
     *  \param propagatorSettings Settings for propagator.
     *  \param areEquationsOfMotionToBeIntegrated Boolean to denote whether equations of motion should be integrated
     *  immediately at the end of the contructor or not (default true).
     *  \param clearNumericalSolutions Boolean to determine whether to clear the raw numerical solution member variables
     *  of this class, after propagation and resetting ephemerides (default false).
     *  \param setIntegratedResult Boolean to determine whether to automatically use the integrated results to set
     *  ephemerides (default false).
     *  \param printNumberOfFunctionEvaluations Boolean denoting whether the number of function evaluations should be printed
     *  at the end of propagation.
     *  \param initialClockTime Initial clock time from which to determine cumulative computation time.
     *  By default now( ), i.e. the moment at which this function is called.
     *  \param stateDerivativeModels List of state derivative models used in the simulation.
     */
    SingleArcDynamicsSimulator(
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const PredefinedSingleArcStateDerivativeModels< StateScalarType, TimeType >& predefinedStateDerivativeModels,
            const bool areEquationsOfMotionToBeIntegrated = true,
            const bool clearNumericalSolutions = false,
            const bool setIntegratedResult = false,
            const bool printNumberOfFunctionEvaluations = false,
            const std::chrono::steady_clock::time_point initialClockTime = std::chrono::steady_clock::now( ),
            const bool printDependentVariableData = false,
            const bool printStateData = false,
            const bool updateDependentVariableInterpolator = false ):
        SingleArcDynamicsSimulator( bodies, validateDeprecatedSingleArcSettings(
                                        integratorSettings, propagatorSettings,
                                        clearNumericalSolutions, setIntegratedResult, printNumberOfFunctionEvaluations,
                                        printDependentVariableData, printStateData, updateDependentVariableInterpolator ),
                                    areEquationsOfMotionToBeIntegrated,
                                    predefinedStateDerivativeModels ){ }

    SingleArcDynamicsSimulator(
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const bool areEquationsOfMotionToBeIntegrated = true,
            const bool clearNumericalSolutions = false,
            const bool setIntegratedResult = false,
            const bool printNumberOfFunctionEvaluations = false,
            const bool printDependentVariableData = false,
            const bool printStateData = false,
            const bool updateDependentVariableInterpolator = false ):
        SingleArcDynamicsSimulator(  bodies, integratorSettings,  propagatorSettings,
                                     PredefinedSingleArcStateDerivativeModels< StateScalarType, TimeType >( ),
                                     areEquationsOfMotionToBeIntegrated,
                                     clearNumericalSolutions,
                                     setIntegratedResult,
                                     printNumberOfFunctionEvaluations,
                                     std::chrono::steady_clock::now( ),
                                     printDependentVariableData,
                                     printStateData,
                                     updateDependentVariableInterpolator ){ }

    //! Destructor
    ~SingleArcDynamicsSimulator( ) { }

    //! This function numerically (re-)integrates the equations of motion.
    /*!
     *  This function numerically (re-)integrates the equations of motion, using the settings set through the constructor
     *  and a new initial state vector provided here. The raw results are set in the equationsOfMotionNumericalSolution_
     *  \param initialStates Initial state vector that is to be used for numerical integration. Note that this state should
     *  be in the correct frame (i.e. corresponding to centralBodies in propagatorSettings_), but not in the propagator-
     *  specific form (i.e Encke, Gauss, etc. for translational dynamics)
     * \sa SingleStateTypeDerivative::convertToOutputSolution
     */
    void integrateEquationsOfMotion(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& initialStates )
    {
        integrateEquationsOfMotion< SingleArcSimulationResults< StateScalarType, TimeType > >(
            dynamicsStateDerivative_->convertFromOutputSolution( initialStates, propagatorSettings_->getInitialTime( ) ),
                propagationResults_ );
    }

    void integrate(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& initialStates )
    {
        integrateEquationsOfMotion( initialStates );
    }

    template< typename SimulationResults >
    void integrateEquationsOfMotion(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& processedInitialState,
            const std::shared_ptr< SimulationResults > propagationResults )
    {
        performPropagationPreProcessingSteps( propagationResults );
        propagateDynamics< SimulationResults >( processedInitialState,
                           propagationResults,
                           PostProcessingFunctionProvider< StateScalarType, TimeType, SimulationResults::number_of_columns >::
                                   getPostProcessingFunction( dynamicsStateDerivative_ ) );
        performPropagationPostProcessingSteps( propagationResults );
    }

    //! Function to return the map of state history of numerically integrated bodies (base class interface).
    /*!
     * Function to return the map of state history of numerically integrated bodies (base class interface).
     * \return Vector is size 1, with entry: map of state history of numerically integrated bodies.
     */
    std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > getEquationsOfMotionNumericalSolutionBase( )
    {
        return std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >({ getEquationsOfMotionNumericalSolution( ) } );
    }

    //! Function to return the map of dependent variable history that was saved during numerical propagation (base class interface)
    /*!
     * Function to return the map of dependent variable history that was saved during numerical propagation (base class interface)
     * \return Vector is size 1, with entry: map of dependent variable history that was saved during numerical propagation.
     */
    std::vector< std::map< TimeType, Eigen::VectorXd > > getDependentVariableNumericalSolutionBase( )
    {
        return std::vector< std::map< TimeType, Eigen::VectorXd > >( { getDependentVariableHistory( ) } );
    }

    //! Function to return the map of cumulative computation time history that was saved during numerical propagation.
    /*!
     * Function to return the map of cumulative computation time history that was saved during numerical propagation (base class interface).
     * \return Vector is size 1, with entry: map of cumulative computation time history that was saved during numerical propagation.
     */
    std::vector< std::map< TimeType, double > > getCumulativeComputationTimeHistoryBase( )
    {
        return std::vector< std::map< TimeType, double > >( { getCumulativeComputationTimeHistory( ) } );
    }

    //! Function to get the settings for the numerical integrator.
    /*!
     * Function to get the settings for the numerical integrator.
     * \return The settings for the numerical integrator.
     */
    std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > getIntegratorSettings( )
    {
        return integratorSettings_;
    }

    //! Function to get the function that performs a single state derivative function evaluation.
    /*!
     * Function to get the function that performs a single state derivative function evaluation.
     * \return Function that performs a single state derivative function evaluation.
     */
    std::function< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >
    ( const TimeType, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >&) >
    getStateDerivativeFunction( )
    {
        return stateDerivativeFunction_;
    }

    //! Function to get the settings for the propagator.
    /*!
     * Function to get the settings for the propagator.
     * \return The settings for the propagator.
     */
    std::shared_ptr< SingleArcPropagatorSettings< StateScalarType, TimeType > > getPropagatorSettings( )
    {
        return propagatorSettings_;
    }

    //! Function to get the object that updates the environment.
    /*!
     * Function to get the object responsible for updating the environment based on the current state and time.
     * \return Object responsible for updating the environment based on the current state and time.
     */
    std::shared_ptr< EnvironmentUpdater< StateScalarType, TimeType > > getEnvironmentUpdater( )
    {
        return environmentUpdater_;
    }

    //! Function to get the object that updates and returns state derivative
    /*!
     * Function to get the object that updates current environment and returns state derivative from single function call
     * \return Object that updates current environment and returns state derivative from single function call
     */
    std::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > getDynamicsStateDerivative( )
    {
        return dynamicsStateDerivative_;
    }

    //! Function to retrieve the object defining when the propagation is to be terminated.
    /*!
     * Function to retrieve the object defining when the propagation is to be terminated.
     * \return Object defining when the propagation is to be terminated.
     */
    std::shared_ptr< PropagationTerminationCondition > getPropagationTerminationCondition( )
    {
        return propagationTerminationCondition_;
    }

    //! Function to retrieve the list of object that process the integrated numerical solution by updating the environment
    /*!
     * Function to retrieve the List of object (per dynamics type) that process the integrated numerical solution by
     * updating the environment
     * \return List of object (per dynamics type) that process the integrated numerical solution by updating the environment
     */
    std::map< IntegratedStateType,
    std::shared_ptr< SingleArcIntegratedStateProcessor< TimeType, StateScalarType > > > getIntegratedStateProcessors( )
    {
        return integratedStateProcessors_;
    }

    //! Function to retrieve initial time of propagation
    /*!
     * Function to retrieve initial time of propagation
     * \return Initial time of propagation
     */
    double getInitialPropagationTime( )
    {
        return propagatorSettings_->getInitialTime( );
    }

    //! Function to reset initial propagation time
    /*!
     * Function to reset initial propagation time
     * \param initialPropagationTime New initial propagation time
     */
    void resetInitialPropagationTime( const double initialPropagationTime )
    {
        propagatorSettings_->resetInitialTime( initialPropagationTime );
    }

    //! Function to retrieve the functions that compute the dependent variables at each time step
    /*!
     * Function to retrieve the functions that compute the dependent variables at each time step
     * \return Functions that compute the dependent variables at each time step
     */
    std::function< Eigen::VectorXd( ) > getDependentVariablesFunctions( )
    {
        return dependentVariablesFunctions_;
    }

    //! Function to reset the object that checks whether the simulation has finished from
    //! (newly defined) propagation settings.
    /*!
     *  Function to reset the object that checks whether the simulation has finished from
     *  (newly defined) propagation settings.
     */
    void resetPropagationTerminationConditions( )
    {
        propagationTerminationCondition_ = createPropagationTerminationConditions(
                    propagatorSettings_->getTerminationSettings( ), bodies_,
                    integratorSettings_->initialTimeStep_,
                    dynamicsStateDerivative_->getStateDerivativeModels( ) );
    }

    //! This function updates the environment with the numerical solution of the propagation.
    /*!
     *  This function updates the environment with the numerical solution of the propagation. It sets
     *  the propagated translational dynamics solution as the new input for the Ephemeris object of the body that was
     *  propagated.
     */
    void processNumericalEquationsOfMotionSolution( )
    {
        if( outputSettings_->getSetIntegratedResult( ) )
        {
            try {
                // Create and set interpolators for ephemerides
                resetIntegratedStates( propagationResults_->equationsOfMotionNumericalSolution_,
                                       integratedStateProcessors_ );
            }
            catch ( const std::exception &caughtException ) {
                std::cerr
                        << "Error occured when post-processing single-arc integration results, and seting integrated states in environment, caught error is: "
                        << std::endl << std::endl;
                std::cerr << caughtException.what( ) << std::endl << std::endl;
                std::cerr <<
                          "The problem may be that there is an insufficient number of data points (epochs) at which propagation results are produced. Integrated results are given at" +
                          std::to_string( propagationResults_->equationsOfMotionNumericalSolution_.size( )) + " epochs"
                          << std::endl;
            }

            // Clear numerical solution if so required.
            if ( propagatorSettings_->getOutputSettings( )->getClearNumericalSolutions( ))
            {
                propagationResults_->clearSolutionMaps( );
            }

            for ( auto bodyIterator: bodies_.getMap( )) {
                bodyIterator.second->updateConstantEphemerisDependentMemberQuantities( );
            }
        }
        else if ( propagatorSettings_->getOutputSettings( )->getClearNumericalSolutions( ))
        {
            propagationResults_->clearSolutionMaps( );
        }

        if( propagatorSettings_->getOutputSettings( )->getUpdateDependentVariableInterpolator( ) )
        {
            propagationResults_->updateDependentVariableInterface( );
        }
    }

    void suppressDependentVariableDataPrinting( )
    {
        outputSettings_->getPrintSettings( )->setPrintDependentVariableData( false );
    }

    void enableDependentVariableDataPrinting( )
    {
        outputSettings_->getPrintSettings( )->setPrintDependentVariableData( true );
    }

    void createAndSetIntegratedStateProcessors( )
    {
        frameManager_ = simulation_setup::createFrameManager( bodies_.getMap( ) );
        integratedStateProcessors_ = createIntegratedStateProcessors< TimeType, StateScalarType >(
                    propagatorSettings_, bodies_, frameManager_ );
    }

    std::shared_ptr< SimulationResults< StateScalarType, TimeType > > getPropagationResults( )
    {
        return propagationResults_;
    }

    std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType > > getSingleArcPropagationResults( )
    {
        return propagationResults_;
    }

///////////////////////////////////////////////////
//////////////// DEPRECATED ///////////////////////
///////////////////////////////////////////////////

    //! Function to return the map of state history of numerically integrated bodies.
    /*!
     * Function to return the map of state history of numerically integrated bodies.
     * \return Map of state history of numerically integrated bodies.
     */
    const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& getEquationsOfMotionNumericalSolution( )
    {
        return propagationResults_->equationsOfMotionNumericalSolution_;
    }

    //! Function to return the map of state history of numerically integrated bodies, in propagation coordinates.
    /*!
     * Function to return the map of state history of numerically integrated bodies, in propagation coordinates.
     * \return Map of state history of numerically integrated bodies, in propagation coordinates.
     */
    const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& getEquationsOfMotionNumericalSolutionRaw( )
    {
        return propagationResults_->equationsOfMotionNumericalSolutionRaw_;
    }

    //! Function to return the map of dependent variable history that was saved during numerical propagation.
    /*!
     * Function to return the map of dependent variable history that was saved during numerical propagation.
     * \return Map of dependent variable history that was saved during numerical propagation.
     */
    const std::map< TimeType, Eigen::VectorXd >& getDependentVariableHistory( )
    {
        return propagationResults_->dependentVariableHistory_;
    }

    //! Function to return the map of cumulative computation time history that was saved during numerical propagation.
    /*!
     * Function to return the map of cumulative computation time history that was saved during numerical propagation.
     * \return Map of cumulative computation time history that was saved during numerical propagation.
     */
    std::map< TimeType, double > getCumulativeComputationTimeHistory( )
    {
        return propagationResults_->cumulativeComputationTimeHistory_;
    }

    //! Function to return the map of number of cumulative function evaluations that was saved during numerical propagation.
    /*!
     * Function to return the map of cumulative number of function evaluations that was saved during numerical propagation.
     * \return Map of cumulative number of function evaluations that was saved during numerical propagation.
     */
    std::map< TimeType, unsigned int > getCumulativeNumberOfFunctionEvaluations( )
    {
        return propagationResults_->cumulativeNumberOfFunctionEvaluations_;
    }


    //! Function to retrieve the event that triggered the termination of the last propagation
    /*!
     * Function to retrieve the event that triggered the termination of the last propagation
     * \return Event that triggered the termination of the last propagation
     */
    std::shared_ptr< PropagationTerminationDetails > getPropagationTerminationReason( )
    {
        return propagationResults_->propagationTerminationReason_;
    }

    void setPropagationTerminationReason( const std::shared_ptr< PropagationTerminationDetails > propagationTerminationReason )
    {
        propagationResults_->propagationTerminationReason_ = propagationTerminationReason;
    }

    //! Get whether the integration was completed successfully.
    /*!
     * Get whether the integration was completed successfully.
     * \return Whether the integration was completed successfully by reaching the termination condition.
     */
    virtual bool integrationCompletedSuccessfully( ) const
    {
        return ( propagationResults_->propagationTerminationReason_->getPropagationTerminationReason( ) == termination_condition_reached );
    }

    //! Function to retrieve the dependent variables IDs
    /*!
     * Function to retrieve the dependent variables IDs
     * \return Map listing starting entry of dependent variables in output vector, along with associated ID
     */
    std::map< std::pair< int, int >, std::string > getDependentVariableIds( )
    {
        return propagationResults_->getDependentVariableId( );
    }

    //! Function return whether the propagation is sequential or not (forward and backward leg).
    bool isPropagationSequential( ) const
    {
        return sequentialPropagation_;
    }

///////////////////////////////////////////////////
//////////////// END DEPRECATED ///////////////////
///////////////////////////////////////////////////


protected:

    //! List of object (per dynamics type) that process the integrated numerical solution by updating the environment
    std::map< IntegratedStateType,
            std::shared_ptr< SingleArcIntegratedStateProcessor< TimeType, StateScalarType > > > integratedStateProcessors_;

    //! Object responsible for updating the environment based on the current state and time.
    /*!
     *  Object responsible for updating the environment based on the current state and time. Calling the updateEnvironment
     * function automatically updates all dependent variables that are needed to calulate the state derivative.
     */
    std::shared_ptr< EnvironmentUpdater< StateScalarType, TimeType > > environmentUpdater_;

    //! Interface object that updates current environment and returns state derivative from single function call.
    std::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > dynamicsStateDerivative_;

    //! Function that performs a single state derivative function evaluation.
    /*!
     *  Function that performs a single state derivative function evaluation, will typically be set to
     *  DynamicsStateDerivativeModel< TimeType, StateScalarType >::computeStateDerivative function.
     *  Calling this function will first update the environment (using environmentUpdater_) and then calculate the
     *  full system state derivative.
     */
    std::function< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >
    ( const TimeType, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& ) > stateDerivativeFunction_;

    //! Settings for numerical integrator.
    std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings_;

    //! Settings for propagator.
    std::shared_ptr< SingleArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings_;

    //! Object defining when the propagation is to be terminated.
    std::shared_ptr< PropagationTerminationCondition > propagationTerminationCondition_;

    //! Function returning dependent variables (during numerical propagation)
    std::function< Eigen::VectorXd( ) > dependentVariablesFunctions_;

//    std::map< std::pair< int, int >, std::string > dependentVariableIds_;
//
//    std::map< std::pair< int, int >, std::string > processedStateIds_;
//
//    std::map< std::pair< int, int >, std::string > propagatedStateIds_;

    //! Object for retrieving ephemerides for transformation of reference frame (origins)
    std::shared_ptr< ephemerides::ReferenceFrameManager > frameManager_;

    std::shared_ptr< SingleArcPropagatorProcessingSettings > outputSettings_;

    std::shared_ptr< SingleArcSimulationResults< StateScalarType, TimeType > > propagationResults_;

    //! Boolean denoting whether the propagation is performing sequentially, or both forward and backward (default = true).
    bool sequentialPropagation_;


private:

    //! Function that propagates the dynamics and (if requested) variational equations.
    /*
    *  Function that propagates the dynamics and (if requested) variational equations. Whether the variational
     *  equations are propagated is defined by the choice of SimulationResults template argument (if
     *  SingleArcSimulationResults< StateScalarType, TimeType >: dynamics only;
     *  if SingleArcVariationalSimulationResults< StateScalarType, TimeType >: dynamics and variational equations)
     *  NOTE: This function requires the performPropagationPreProcessingSteps
     *  and performPropagationPostProcessingSteps to be called before/after it. This is done automatically by the
     *  integrateEquationsOfMotion function.
     */
    template< typename SimulationResults >
    void propagateDynamics(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, SimulationResults::number_of_columns >& processedInitialState,
            const std::shared_ptr< SimulationResults > propagationResults,
            const std::function< void( Eigen::Matrix< StateScalarType, Eigen::Dynamic, SimulationResults::number_of_columns >& ) > statePostProcessingFunction )
    {
        // Integrate equations of motion numerically.
        simulation_setup::setAreBodiesInPropagation( bodies_, true );
        dynamicsStateDerivative_->updateStateDerivativeModelSettings( processedInitialState.block(
                0, processedInitialState.cols( ) - 1, processedInitialState.rows(), 1  ) );

        if ( sequentialPropagation_ )
        {
            integrateEquations< SimulationResults, Eigen::Matrix< StateScalarType, Eigen::Dynamic, SimulationResults::number_of_columns >, TimeType >(
                    stateDerivativeFunction_,
                    processedInitialState ,
                    propagatorSettings_->getInitialTime( ),
                    integratorSettings_,
                    propagationTerminationCondition_,
                    propagationResults,
                    dependentVariablesFunctions_,
                    statePostProcessingFunction,
                    propagatorSettings_->getOutputSettings( ) );
        }
        else
        {
            std::shared_ptr< NonSequentialPropagationTerminationCondition > nonSequentialTerminations =
                    std::dynamic_pointer_cast< NonSequentialPropagationTerminationCondition >( propagationTerminationCondition_ );
            integrateEquations< SimulationResults, Eigen::Matrix< StateScalarType, Eigen::Dynamic, SimulationResults::number_of_columns >, TimeType >(
                    stateDerivativeFunction_,
                    processedInitialState ,
                    propagatorSettings_->getInitialTime( ),
                    integratorSettings_,
                    nonSequentialTerminations->getForwardPropagationTerminationCondition( ),
                    propagationResults,
                    dependentVariablesFunctions_,
                    statePostProcessingFunction,
                    propagatorSettings_->getOutputSettings( ) );

            integratorSettings_->initialTimeStep_ *= -1.0;
            integrateEquations< SimulationResults, Eigen::Matrix< StateScalarType, Eigen::Dynamic, SimulationResults::number_of_columns >, TimeType >(
                    stateDerivativeFunction_,
                    processedInitialState ,
                    propagatorSettings_->getInitialTime( ),
                    integratorSettings_,
                    nonSequentialTerminations->getBackwardPropagationTerminationCondition( ),
                    propagationResults,
                    dependentVariablesFunctions_,
                    statePostProcessingFunction,
                    propagatorSettings_->getOutputSettings( ) );
            integratorSettings_->initialTimeStep_ *= -1.0;
        }

        simulation_setup::setAreBodiesInPropagation( bodies_, false );
    }

    //! Function to perform steps necessary to reset all relevant models for the upcoming propagation
    /*
     *  Function to perform steps necessary to reset all relevant models for the upcoming propagation:
     *  - Whether to propagate dynamics and/or vatiational equations
     *  - Reset counter of function evaluations to zero
     *  - Reset termination conditions
     *  - Empty object holding the numerical simulation results of the previous run
     *  - Print messages to terminal, as requested by user settings
     */
    template< typename SimulationResults >
    void performPropagationPreProcessingSteps(
            const std::shared_ptr< SimulationResults > propagationResults )
    {
        // Reset functions
        dynamicsStateDerivative_->setPropagationSettings( std::vector< IntegratedStateType >( ), true, SimulationResults::is_variational );
        dynamicsStateDerivative_->resetFunctionEvaluationCounter( );
        dynamicsStateDerivative_->resetCumulativeFunctionEvaluationCounter( );
        resetPropagationTerminationConditions( );

        // Empty solution maps
        propagationResults->reset( );

        PropagationPrintingInterface< SimulationResults, StateScalarType, TimeType >::printSingleArcPrePropagationMessages(
                outputSettings_->getPrintSettings( ),
                outputSettings_->getPropagationStartHeader( ),
                propagationResults );
    }

    //! Function to perform steps necessary to finalize the propagation
    /*
     *  Function to perform steps necessary to finalize the propagation
     *  - Store number of function evaluations in the results object
     *  - Print messages to terminal, as requested by user settings
     *  - Update the environment (e.g. use numerical results to create tabulated ephemerides and similar for other dynamics)
     *    if requested by user
     */
    template< typename SimulationResults >
    void performPropagationPostProcessingSteps(
            const std::shared_ptr< SimulationResults > propagationResults )
    {
        // Retrieve number of cumulative function evaluations
        propagationResults->finalizePropagation( dynamicsStateDerivative_->getCumulativeNumberOfFunctionEvaluations( ) );
        PropagationPrintingInterface< SimulationResults, StateScalarType, TimeType >::printSingleArcPostPropagationMessages(
                outputSettings_->getPrintSettings( ),
                                               outputSettings_->getPropagationEndHeader( ),
                                               propagationResults );
        processNumericalEquationsOfMotionSolution( );
    }
};


//! Function to get a vector of initial states from a vector of propagator settings
/*!
 *  Function to get a vector of initial states from a vector of propagator settings.
 *  \param propagatorSettings List of propagator settings
 *  \return List of initial states, as retrieved from propagatorSettings list.
 */
template< typename StateScalarType = double >
std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1  > > getInitialStatesPerArc(
        const std::vector< std::shared_ptr< PropagatorSettings< StateScalarType > > > propagatorSettings )
{
    std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1  > > initialStatesList;
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
            throw std::runtime_error(
                        "Error when getting initial arc state from previous arc: no arc overlap" );
        }
        else
        {
            int currentIndex = 0;
            int initialTimeIndex = -1;

            std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > initialStateInterpolationMap;

            // Set sub-part of previous arc to interpolate for current arc
            for( typename std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >::
                 const_reverse_iterator previousArcIterator = previousArcDynamicsSolution.rbegin( );
                 previousArcIterator != previousArcDynamicsSolution.rend( ); previousArcIterator++ )
            {
                initialStateInterpolationMap[ previousArcIterator->first ] = previousArcIterator->second;
                if( initialTimeIndex < 0 )
                {
                    if( previousArcIterator->first <  currentArcInitialTime )
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
            currentArcInitialState =
                    std::make_shared< interpolators::LagrangeInterpolator<
                    TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >, long double > >(
                        initialStateInterpolationMap, 8 )->interpolate( currentArcInitialTime );

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
        const bool updateDependentVariableInterpolator = false  )
{
    std::shared_ptr< MultiArcPropagatorSettings< StateScalarType, TimeType > > multiArcPropagatorSettings =
            std::dynamic_pointer_cast< MultiArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings );
    if( multiArcPropagatorSettings == nullptr )
    {
        throw std::runtime_error( "Error in dynamics simulator (deprecated), input must be multi-arc." );
    }

    std::vector<std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > integratorSettingsList;
    if( integratorSettings.size( ) == 1 &&  multiArcPropagatorSettings->getSingleArcSettings( ).size( ) > 1 )
    {
        integratorSettingsList = std::vector<std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > >(
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
                std::cerr<<"Warning, multi-arc integrator settings, defined independently, and in propagator settings"<<std::endl;
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

    std::vector<std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > integratorSettingsList(
            propagationStartTimes.size( ), integratorSettings);

    std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > independentIntegratorSettingsList =
            utilities::cloneDuplicatePointers( integratorSettingsList );

    for( unsigned int i = 0; i < propagationStartTimes.size( ); i++ )
    {
        multiArcPropagatorSettings->getSingleArcSettings( ).at( i )->resetInitialTime(
                propagationStartTimes.at( i ) );
    }

    return validateDeprecatedMultiArcSettings(
                independentIntegratorSettingsList, propagatorSettings, clearNumericalSolutions, setIntegratedResult, updateDependentVariableInterpolator );
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
        const std::shared_ptr< MultiArcSimulationResults<SingleArcSimulationResults, StateScalarType, TimeType> > originalPropagationResults,
        const std::shared_ptr<SimulationResults> comparePropagationResults )
{

}

template< typename StateScalarType, typename TimeType >
void checkPropagationResultsObjectConsistency(
        const std::shared_ptr< MultiArcSimulationResults< SingleArcSimulationResults, StateScalarType, TimeType> > originalPropagationResults,
        const std::shared_ptr< MultiArcSimulationResults< SingleArcVariationalSimulationResults, StateScalarType, TimeType > > comparePropagationResults )
{
    if( originalPropagationResults->getSingleArcResults( ).size( ) != comparePropagationResults->getSingleArcResults( ).size( ) )
    {
        throw std::runtime_error( "Error when checking consistency of multi-arc dynamics results with variational input; results objects number of single arcs is inconsistent" );
    }

    for( unsigned int i = 0; i < originalPropagationResults->getSingleArcResults( ).size( ); i++ )
    {
        if( originalPropagationResults->getSingleArcResults( ) != comparePropagationResults->getSingleArcResults( )->getSingleArcResults( ) )
        {
            throw std::runtime_error( "Error when checking consistency of multi-arc dynamics results with variational input; results objects are incosistent" );
        }
    }
}
template< typename StateScalarType, typename TimeType >
void checkPropagationResultsObjectConsistency(
        const std::shared_ptr< MultiArcSimulationResults< SingleArcSimulationResults, StateScalarType, TimeType> > originalPropagationResults,
        const std::shared_ptr< MultiArcSimulationResults< SingleArcSimulationResults, StateScalarType, TimeType > > comparePropagationResults )
{
    if( originalPropagationResults != comparePropagationResults )
    {
        if( originalPropagationResults->getSingleArcResults( ).size( ) != comparePropagationResults->getSingleArcResults( ).size( ) )
        {
            throw std::runtime_error( "Error when checking consistency of multi-arc dynamics results with dynamics-only input; results objects number of single arcs is inconsistent" );
        }
        for( unsigned int i = 0; i < originalPropagationResults->getSingleArcResults( ).size( ); i++ )
        {
            if( originalPropagationResults->getSingleArcResults( ) != comparePropagationResults->getSingleArcResults( ) )
            {
                throw std::runtime_error( "Error when checking consistency of multi-arc dynamics results with dynamics-only input; results objects are incosistent" );
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
class MultiArcDynamicsSimulator: public DynamicsSimulator< StateScalarType, TimeType > {
public:

    typedef MultiArcSimulationResults<SingleArcSimulationResults, StateScalarType, TimeType> MultiArcResults;
    using DynamicsSimulator<StateScalarType, TimeType>::bodies_;

    MultiArcDynamicsSimulator(
            const simulation_setup::SystemOfBodies &bodies,
            const std::shared_ptr<MultiArcPropagatorSettings<StateScalarType, TimeType> > propagatorSettings,
            const bool areEquationsOfMotionToBeIntegrated = true ) :
            DynamicsSimulator<StateScalarType, TimeType>(
                    bodies, propagatorSettings ),
            multiArcPropagatorSettings_( propagatorSettings )
    {
        if ( multiArcPropagatorSettings_ == nullptr )
        {
            throw std::runtime_error( "Error when creating multi-arc dynamics simulator, input is not multi arc" );
        }
        else
        {

            std::vector<std::shared_ptr<SingleArcPropagatorSettings<StateScalarType, TimeType> > > singleArcSettings =
                    multiArcPropagatorSettings_->getSingleArcSettings( );

            // Create dynamics simulators
            std::vector<std::shared_ptr<SingleArcSimulationResults<StateScalarType, TimeType> > > singleArcResults;
            for ( unsigned int i = 0; i < singleArcSettings.size( ); i++ ) {
                singleArcDynamicsSimulators_.push_back(
                        std::make_shared<SingleArcDynamicsSimulator<StateScalarType, TimeType> >(
                                bodies, singleArcSettings.at( i ), false ));
                singleArcResults.push_back( singleArcDynamicsSimulators_.at( i )->getSingleArcPropagationResults( ));
                singleArcDynamicsSimulators_.at( i )->createAndSetIntegratedStateProcessors( );
            }

            std::vector< std::shared_ptr< SingleArcDependentVariablesInterface< TimeType > > > singleArcInterfaces;
            for ( unsigned int i = 0 ; i < singleArcSettings.size( ) ; i++ )
            {
                singleArcInterfaces.push_back( singleArcDynamicsSimulators_.at( i )->getSingleArcPropagationResults( )->getSingleArcDependentVariablesInterface( ) );
            }

            std::shared_ptr< MultiArcDependentVariablesInterface< TimeType > > dependentVariableInterface =
                std::make_shared< MultiArcDependentVariablesInterface< TimeType > >(
                    singleArcInterfaces, std::vector< double >( ), std::vector< double >( ) );

            propagationResults_ = std::make_shared<MultiArcResults>( singleArcResults, dependentVariableInterface );

            // Integrate equations of motion if required.
            if ( areEquationsOfMotionToBeIntegrated )
            {
                integrateEquationsOfMotion( multiArcPropagatorSettings_->getInitialStates( ));
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
    MultiArcDynamicsSimulator(
            const simulation_setup::SystemOfBodies &bodies,
            const std::shared_ptr<numerical_integrators::IntegratorSettings<TimeType> > integratorSettings,
            const std::shared_ptr<PropagatorSettings<StateScalarType> > propagatorSettings,
            const std::vector<double> arcStartTimes,
            const bool areEquationsOfMotionToBeIntegrated = true,
            const bool clearNumericalSolutions = true,
            const bool setIntegratedResult = true,
            const bool updateDependentVariableInterpolator = false ):
        MultiArcDynamicsSimulator( bodies, validateDeprecatedMultiArcSettings(
                                        integratorSettings, propagatorSettings, arcStartTimes,
                                        clearNumericalSolutions, setIntegratedResult, updateDependentVariableInterpolator ),
                                    areEquationsOfMotionToBeIntegrated ){ }


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
            const simulation_setup::SystemOfBodies &bodies,
            const std::vector<std::shared_ptr<numerical_integrators::IntegratorSettings<TimeType> > > integratorSettings,
            const std::shared_ptr<PropagatorSettings<StateScalarType> > propagatorSettings,
            const bool areEquationsOfMotionToBeIntegrated = true,
            const bool clearNumericalSolutions = true,
            const bool setIntegratedResult = true,
            const bool updateDependentVariableInterpolator = false ):
        MultiArcDynamicsSimulator( bodies, validateDeprecatedMultiArcSettings(
                                        integratorSettings, propagatorSettings,
                                        clearNumericalSolutions, setIntegratedResult, updateDependentVariableInterpolator ),
                                    areEquationsOfMotionToBeIntegrated ){ }

    //! Destructor
    ~MultiArcDynamicsSimulator( ) { }

    Eigen::Matrix<StateScalarType, Eigen::Dynamic, Eigen::Dynamic> getArcInitialState(
            const int arcIndex,
            const std::shared_ptr<MultiArcInitialStateProvider<StateScalarType> > stateProvider ) {
        bool initialStateFromPreviousArc = false;
        Eigen::Matrix<StateScalarType, Eigen::Dynamic, Eigen::Dynamic> currentArcInitialState =
                stateProvider->getArcInitialState( arcIndex, initialStateFromPreviousArc );
        if ( initialStateFromPreviousArc ) {
            currentArcInitialState.block( 0, currentArcInitialState.cols( ) - 1, currentArcInitialState.rows( ),
                                          1 ) = getArcInitialStateFromPreviousArcResult(
                    propagationResults_->getSingleArcResults( ).at(
                            arcIndex - 1 )->getEquationsOfMotionNumericalSolution( ),
                    singleArcDynamicsSimulators_.at( arcIndex )->getInitialPropagationTime( ));
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
    void integrateEquationsOfMotion(
            const Eigen::Matrix<StateScalarType, Eigen::Dynamic, Eigen::Dynamic> &concatenatedInitialStates ) {
        std::vector<Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> > splitInitialState;

        int currentIndex = 0;
        for ( unsigned int i = 0; i < singleArcDynamicsSimulators_.size( ); i++ ) {
            int currentSize = singleArcDynamicsSimulators_.at(
                    i )->getPropagatorSettings( )->getConventionalStateSize( );
            splitInitialState.push_back( concatenatedInitialStates.block( currentIndex, 0, currentSize, 1 ));
            currentIndex += currentSize;
        }

        if ( currentIndex != concatenatedInitialStates.rows( )) {
            throw std::runtime_error(
                    "Error when doing multi-arc integration, input state vector size is incompatible with settings" );
        }

        integrateEquationsOfMotion( splitInitialState );
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
    void integrateEquationsOfMotion(
            const std::vector<Eigen::Matrix<StateScalarType, Eigen::Dynamic, 1> > &initialStatesList )
    {
        integrateEquationsOfMotion < MultiArcResults > ( propagationResults_,
                std::make_shared<MultiArcInitialStateProvider<StateScalarType> >( initialStatesList ) );
    }


    template< typename MultiArcSimulationResults >
    void integrateEquationsOfMotion(
            const std::shared_ptr< MultiArcSimulationResults > propagationResults,
            const std::shared_ptr< MultiArcInitialStateProvider< StateScalarType > > initialStateProvider )
    {
        checkPropagationResultsObjectConsistency< StateScalarType, TimeType, MultiArcSimulationResults >(
                propagationResults_,
                propagationResults );

        Eigen::Matrix< StateScalarType, Eigen::Dynamic, MultiArcSimulationResults::single_arc_type::number_of_columns > currentArcInitialState;
        std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, MultiArcSimulationResults::single_arc_type::number_of_columns > > arcInitialStateList;

        initialStateProvider->restartPropagation( );
        propagationResults->restartPropagation( );
        propagationResults_->restartPropagation( );

        printPrePropagationMessages( );


        // Propagate dynamics for each arc
        for( unsigned int i = 0; i < singleArcDynamicsSimulators_.size( ); i++ )
        {
            currentArcInitialState = getArcInitialState( i, initialStateProvider );
            arcInitialStateList.push_back( currentArcInitialState );

            singleArcDynamicsSimulators_.at( i )->template integrateEquationsOfMotion<
                    typename MultiArcSimulationResults::single_arc_type >( currentArcInitialState, propagationResults->getSingleArcResults( ).at( i ) );
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

            multiArcPropagatorSettings_->resetInitialStatesList(
                        newInitialStates );
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
                std::map<IntegratedStateType, std::vector<std::shared_ptr<
                        SingleArcIntegratedStateProcessor<TimeType, StateScalarType> > > > singleArcIntegratedStatesProcessors;

                for ( unsigned int i = 0; i < singleArcDynamicsSimulators_.size( ); i++ ) {
                    std::map<IntegratedStateType, std::shared_ptr<
                            SingleArcIntegratedStateProcessor<TimeType, StateScalarType> > > currentArcStateProcessors =
                            singleArcDynamicsSimulators_.at( i )->getIntegratedStateProcessors( );

                    for ( auto itr: currentArcStateProcessors ) {
                        singleArcIntegratedStatesProcessors[ itr.first ].push_back( itr.second );
                    }
                }

                std::map<IntegratedStateType,
                        std::shared_ptr<MultiArcIntegratedStateProcessor<TimeType, StateScalarType> > > multiArcStateProcessors
                        = createMultiArcIntegratedStateProcessors( bodies_, propagationResults_->getArcStartTimes( ),
                                                                   singleArcIntegratedStatesProcessors );
                for ( auto itr: multiArcStateProcessors )
                {
                    itr.second->processIntegratedMultiArcStates(
                            propagationResults_->getConcatenatedEquationsOfMotionResults(
                                    multiArcPropagatorSettings_->getOutputSettings( )->getClearNumericalSolutions( )),
                            propagationResults_->getArcStartTimes( ));
                }
            }
            catch ( const std::exception &caughtException ) {
                std::cerr
                        << "Error occured when post-processing mulyi-arc integration results, and seting integrated states in environment, caught error is: "
                        << std::endl << std::endl;
                std::cerr << caughtException.what( ) << std::endl << std::endl;
                std::cerr
                        << "The problem may be that there is an insufficient number of data points (epochs) at which propagation results are produced for one or more arcs"
                        << std::endl;
                if ( multiArcPropagatorSettings_->getOutputSettings( )->getClearNumericalSolutions( )) {
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
            std::cout<<multiArcPropagatorSettings_->getOutputSettings( )->getPropagationStartHeader( )<<std::endl<<std::endl;
        }
    }

    void printPostPropagationMessages( )
    {
        if( multiArcPropagatorSettings_->getOutputSettings( )-> printAnyOutput( ) )
        {
            std::cout<<multiArcPropagatorSettings_->getOutputSettings( )->getPropagationEndHeader( )<<std::endl<<std::endl;
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
    std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
    getEquationsOfMotionNumericalSolution( )
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
    std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
    getEquationsOfMotionNumericalSolutionBase( )
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

    validateDeprecatedSingleArcSettings< StateScalarType, TimeType >(
                singleArcIntegratorSettings, hybridArcPropagatorSettings->getSingleArcPropagatorSettings( ) );
    validateDeprecatedMultiArcSettings< StateScalarType, TimeType >(
                multiArcIntegratorSettings, hybridArcPropagatorSettings->getMultiArcPropagatorSettings( ),
                clearNumericalSolutions, setIntegratedResult, updateDependentVariableInterpolator );

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
        const bool updateDependentVariableInterpolator = false  )
{
    std::shared_ptr< HybridArcPropagatorSettings< StateScalarType, TimeType > > hybridArcPropagatorSettings =
            std::dynamic_pointer_cast< HybridArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings );
    if( hybridArcPropagatorSettings == nullptr )
    {
        throw std::runtime_error( "Error in dynamics simulator (deprecated), input must be hybrid-arc." );
    }
    std::vector<std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > integratorSettingsList(
                arcStartTimes.size( ), multiArcIntegratorSettings);

    std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > independentIntegratorSettingsList =
            utilities::cloneDuplicatePointers( integratorSettingsList );
    for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
    {
        hybridArcPropagatorSettings->getMultiArcPropagatorSettings( )->getSingleArcSettings( ).at( i )->resetInitialTime(
                    arcStartTimes.at( i ) );
    }
    hybridArcPropagatorSettings->getSingleArcPropagatorSettings( )->resetInitialTime( singleArcIntegratorSettings->initialTimeDeprecated_ );
    return validateDeprecatedHybridArcSettings( singleArcIntegratorSettings, independentIntegratorSettingsList,
                                                propagatorSettings, clearNumericalSolutions, setIntegratedResult, updateDependentVariableInterpolator );
}

template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< HybridArcPropagatorSettings< StateScalarType, TimeType > > validateDeprecatedHybridArcSettings(
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
        const std::vector< double > arcStartTimes,
        const bool clearNumericalSolutions = true,
        const bool setIntegratedResult = true,
        const bool updateDependentVariableInterpolator = false )
{
    return validateDeprecatedHybridArcSettings(
                integratorSettings, utilities::deepcopyPointer( integratorSettings ), propagatorSettings,
                arcStartTimes, clearNumericalSolutions, setIntegratedResult, updateDependentVariableInterpolator );
}


//! Class for performing full numerical integration of a dynamical system, with a compbination of single and multi-arc propagations
/*!
 *  Class for performing full numerical integration of a dynamical system, with a compbination of single and multi-arc
 *  propagations. it is assumed that the single-arc propagations are not influence by the multi-arc propagations: first the
 *  single-arc propagation is done, followed by the multi-arc one (using the single-arc propagaton result for the environment)
 *  In this class, the governing equations are set once,  but can be re-integrated for different initial conditions using
 * the same instance of the class.
 */
template< typename StateScalarType = double, typename TimeType = double >
class HybridArcDynamicsSimulator: public DynamicsSimulator< StateScalarType, TimeType >
{
public:

    //! Using statemebts
    using DynamicsSimulator< StateScalarType, TimeType >::bodies_;

    typedef HybridArcSimulationResults< SingleArcSimulationResults, StateScalarType, TimeType > HybridArcResults;

    HybridArcDynamicsSimulator(
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< HybridArcPropagatorSettings< StateScalarType, TimeType > > hybridPropagatorSettings,
            const bool areEquationsOfMotionToBeIntegrated = true,
            const bool addSingleArcBodiesToMultiArcDynamics = false ):
        DynamicsSimulator< StateScalarType, TimeType >(
            bodies, hybridPropagatorSettings ),
        hybridPropagatorSettings_( hybridPropagatorSettings )
    {
        if( hybridPropagatorSettings_ == nullptr )
        {
            throw std::runtime_error( "Error when making HybridArcDynamicsSimulator, propagator settings are incompatible" );
        }

        singleArcDynamicsSize_ = hybridPropagatorSettings_->getSingleArcPropagatorSettings( )->getPropagatedStateSize( );
        multiArcDynamicsSize_ = hybridPropagatorSettings_->getMultiArcPropagatorSettings( )->getPropagatedStateSize( );

        if( !addSingleArcBodiesToMultiArcDynamics )
        {
            singleArcDynamicsSimulator_ = std::make_shared< SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                        bodies, hybridPropagatorSettings_->getSingleArcPropagatorSettings( ),
                        false, PredefinedSingleArcStateDerivativeModels< StateScalarType, TimeType >( ) );
            multiArcDynamicsSimulator_ = std::make_shared< MultiArcDynamicsSimulator< StateScalarType, TimeType > >(
                        bodies, hybridPropagatorSettings_->getMultiArcPropagatorSettings( ),
                        false );
            propagationResults_ = std::make_shared< HybridArcResults >(
                    singleArcDynamicsSimulator_->getSingleArcPropagationResults( ),
                    multiArcDynamicsSimulator_->getMultiArcPropagationResults( ) );
        }
        else
        {
            throw std::runtime_error( "Cannot yet add single-arc bodies to multi-arc propagation" );
        }

        if( areEquationsOfMotionToBeIntegrated )
        {
            integrateEquationsOfMotion( hybridPropagatorSettings_->getInitialStates( ) );
        }
    }

    //! Constructor of multi-arc simulator for same integration settings per arc.
    /*!
     *  Constructor of multi-arc simulator for same integration settings per arc.
     *  \param bodies Map of bodies (with names) of all bodies in integration.
     *  \param integratorSettings Integrator settings for numerical integrator, used for all arcs.
     *  \param propagatorSettings Propagator settings for dynamics (must be of multi arc type)
     *  \param arcStartTimes Times at which the separate arcs start, for the multi-arc case
     *  \param areEquationsOfMotionToBeIntegrated Boolean to denote whether equations of motion should be integrated at
     *  the end of the contructor or not.extern template class MultiArcDynamicsSimulator< double, double >;


     *  \param clearNumericalSolutions Boolean to determine whether to clear the raw numerical solution member variables
     *  after propagation and resetting ephemerides (default true).
     *  \param setIntegratedResult Boolean to determine whether to automatically use the integrated results to set
     *  ephemerides (default true).
     *  \param addSingleArcBodiesToMultiArcDynamics Boolean denoting whether to add single arc bodies to multi-arc
     *  dynamics (default false).
     */
    HybridArcDynamicsSimulator(
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const std::vector< double > arcStartTimes,
            const bool areEquationsOfMotionToBeIntegrated = true,
            const bool clearNumericalSolutions = true,
            const bool setIntegratedResult = true,
            const bool addSingleArcBodiesToMultiArcDynamics = false,
            const bool updateDependentVariableInterpolator = false ):
                HybridArcDynamicsSimulator( bodies, validateDeprecatedHybridArcSettings(
                                                integratorSettings, propagatorSettings, arcStartTimes,
                                                clearNumericalSolutions, setIntegratedResult, updateDependentVariableInterpolator ),
                                            areEquationsOfMotionToBeIntegrated, addSingleArcBodiesToMultiArcDynamics ){ }


//    HybridArcDynamicsSimulator(
//            const simulation_setup::SystemOfBodies& bodies,
//            const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > singleArcIntegratorSettings,
//            const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > multiArcIntegratorSettings,
//            const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
//            const std::vector< double > arcStartTimes,
//            const bool areEquationsOfMotionToBeIntegrated = true,
//            const bool clearNumericalSolutions = true,
//            const bool setIntegratedResult = true,
//            const bool addSingleArcBodiesToMultiArcDynamics = false,
//            const bool updateDependentVariableInterpolator = false ):
//        HybridArcDynamicsSimulator( bodies, validateDeprecatedHybridArcSettings(
//                                        singleArcIntegratorSettings, multiArcIntegratorSettings, propagatorSettings, arcStartTimes,
//                                        clearNumericalSolutions, setIntegratedResult ),
//                                    areEquationsOfMotionToBeIntegrated, addSingleArcBodiesToMultiArcDynamics, updateDependentVariableInterpolator ){ }
    //! Destructor
    ~HybridArcDynamicsSimulator( ){ }

    //! This function numerically (re-)integrates the equations of motion, using concatenated states for single and multi-arcs
    /*!
     *  This function numerically (re-)integrates the equations of motion, using the settings set through the constructor
     *  and a new initial state vector provided here. The raw results are set in the equationsOfMotionNumericalSolution_
     *  \param initialGlobalStates Initial state vector that is to be used for numerical integration. Note that this state
     *  should be in the correct frame (i.e. corresponding to centralBodies in propagatorSettings_), but not in the propagator-
     *  specific form (i.e Encke, Gauss, etc. for translational dynamics). The states for all arcs must be concatenated in
     *  order into a single Eigen Vector, starting with the single-arc states, followed by the mulit-arc states
     */
    void integrateEquationsOfMotion(
                const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& initialGlobalStates )
    {
        printPrePropagationMessages( );
        singleArcDynamicsSimulator_->integrateEquationsOfMotion(
                    initialGlobalStates.block( 0, 0, singleArcDynamicsSize_, 1 ) );
        multiArcDynamicsSimulator_->integrateEquationsOfMotion(
                    initialGlobalStates.block( singleArcDynamicsSize_, 0, multiArcDynamicsSize_, 1 ) );
        printPostPropagationMessages( );
    }

    //! This function updates the environment with the numerical solution of the propagation
    /*!
     *  This function updates the environment with the numerical solution of the propagation
     *  (no additional functionality in hybrid arc). Function may be used to process manually updtaed propagation results in
     *  single and/or multi-arc model
     */
    void processNumericalEquationsOfMotionSolution( )
    {
        singleArcDynamicsSimulator_->processNumericalEquationsOfMotionSolution( );
        multiArcDynamicsSimulator_->processNumericalEquationsOfMotionSolution( );
        if( hybridPropagatorSettings_->getOutputSettings( )->getUpdateDependentVariableInterpolator( ) )
        {
            propagationResults_->updateDependentVariableInterface( );
        }
    }

    //! Function to retrieve the single-arc dynamics simulator
    /*!
     * Function to retrieve the single-arc dynamics simulator
     * \return Single-arc dynamics simulator
     */
    std::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > getSingleArcDynamicsSimulator( )
    {
        return singleArcDynamicsSimulator_;
    }

    //! Function to retrieve the multi-arc dynamics simulator
    /*!
     * Function to retrieve the multi-arc dynamics simulator
     * \return Multi-arc dynamics simulator
     */
    std::shared_ptr< MultiArcDynamicsSimulator< StateScalarType, TimeType > > getMultiArcDynamicsSimulator( )
    {
        return multiArcDynamicsSimulator_;
    }

    //! Get whether the integration was completed successfully.
    /*!
     * Get whether the integration was completed successfully.
     * \return Whether the integration was completed successfully by reaching the termination condition.
     */
    virtual bool integrationCompletedSuccessfully( ) const
    {
        return propagationResults_->integrationCompletedSuccessfully( );
    }

    //! Function to return the numerical solution to the equations of motion (base class interface).
    /*!
     *  Function to return the numerical solution to the equations of motion for last numerical integration. First vector entry
     *  contains single-arc results. Each subsequent vector entry contains one of the multi-arcs. Key of map denotes time,
     *  values are full propagated state vectors.
     *  \return List of maps of history of numerically integrated states.
     */
    std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
    getEquationsOfMotionNumericalSolutionBase( )
    {
        return propagationResults_->getConcatenatedEquationsOfMotionResults( );
    }

    //! Function to return the numerical solution to the dependent variables (base class interface).
    /*!
     *  Function to return the numerical solution to the dependent variables for last numerical integration. First vector entry
     *  contains single-arc results. Each subsequent vector entry contains one of the multi-arcs. Key of map denotes time,
     *  values are dependent variables vectors
     *  \return List of maps of dependent variable history
     */
    std::vector< std::map< TimeType, Eigen::VectorXd > > getDependentVariableNumericalSolutionBase( )
    {
        return propagationResults_->getConcatenatedDependentVariableResults();
    }

    //! Function to return the map of cumulative computation time history that was saved during numerical propagation.
    /*!
     *  Function to return the map of cumulative computation time history that was saved during numerical propagation.  First vector
     *  entry contains single-arc results. Each subsequent vector entry contains one of the multi-arcs. Key of map denotes time,
     *  values are computation times.
     *  \return Vector is size 1, with entry: map of cumulative computation time history that was saved during numerical propagation.
     */
    std::vector< std::map< TimeType, double > > getCumulativeComputationTimeHistoryBase( )
    {
        return propagationResults_->getConcatenatedCumulativeComputationTimeHistory( );
    }

    void printPrePropagationMessages( )
    {
        if( hybridPropagatorSettings_->getOutputSettings( )->printAnyOutput( ) )
        {
            std::cout<<hybridPropagatorSettings_->getOutputSettings( )->getPropagationStartHeader( )<<std::endl<<std::endl;
        }
    }

    void printPostPropagationMessages( )
    {
        if( hybridPropagatorSettings_->getOutputSettings( )-> printAnyOutput( ) )
        {
            std::cout<<hybridPropagatorSettings_->getOutputSettings( )->getPropagationEndHeader( )<<std::endl<<std::endl;
        }
    }

    std::shared_ptr< SimulationResults< StateScalarType, TimeType > > getPropagationResults( )
    {
        return propagationResults_;
    }

    std::shared_ptr< HybridArcResults > getHybridArcPropagationResults( )
    {
        return propagationResults_;
    }

protected:

    std::shared_ptr< HybridArcPropagatorSettings< StateScalarType, TimeType > > hybridPropagatorSettings_;

    //! Object used to propagate single-arc dynamics
    std::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > singleArcDynamicsSimulator_;

    //! Object used to propagate multi-arc dynamics
    std::shared_ptr< MultiArcDynamicsSimulator< StateScalarType, TimeType > > multiArcDynamicsSimulator_;

    //! Size of single-arc (initial) state vector
    int singleArcDynamicsSize_;

    //! Size of multi-arc concatenated initial state vector
    int multiArcDynamicsSize_;

    std::shared_ptr< HybridArcResults > propagationResults_;

};




template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< PropagatorSettings< StateScalarType > > validateDeprecatePropagatorSettings(
        const std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > >& integratorSettings,
        const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings )
{
    if( std::dynamic_pointer_cast< propagators::SingleArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) != nullptr )
    {
        if( integratorSettings.size( ) == 0 )
        {
            throw std::runtime_error( "Error when validating deprecated propagator settings; did not find integrator settings for single-arc propagation" );;
        }
        return validateDeprecatedSingleArcSettings( integratorSettings.at( 0 ), propagatorSettings );
    }
    else if( std::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) != nullptr )
    {
        return validateDeprecatedMultiArcSettings( integratorSettings, propagatorSettings );
    }
    else if( std::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< StateScalarType > >( propagatorSettings ) != nullptr )
    {
        std::shared_ptr< propagators::HybridArcPropagatorSettings< StateScalarType > > hybridArcSettings =
                std::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< StateScalarType > >( propagatorSettings );

        if( integratorSettings.size( ) == 0 )
        {
            throw std::runtime_error( "Error when validating deprecated propagator settings; did not find integrator settings for hybrid-arc propagation" );;
        }
        validateDeprecatedSingleArcSettings< StateScalarType, TimeType >(
                    integratorSettings.at( 0 ), hybridArcSettings->getSingleArcPropagatorSettings( ) );
        validateDeprecatedMultiArcSettings< StateScalarType, TimeType >(
                    { integratorSettings.begin( ) + 1, integratorSettings.end( ) }, hybridArcSettings->getMultiArcPropagatorSettings( ) );
        return hybridArcSettings;
    }
    else
    {
        throw std::runtime_error( "Error when validating deprecated propagator settings" );
        return nullptr;
    }
}

extern template class SingleArcDynamicsSimulator< double, double >;
extern template class MultiArcDynamicsSimulator< double, double >;
extern template class HybridArcDynamicsSimulator< double, double >;



} // namespace propagators

} // namespace tudat


#endif // TUDAT_DYNAMICSSIMULATOR_H
