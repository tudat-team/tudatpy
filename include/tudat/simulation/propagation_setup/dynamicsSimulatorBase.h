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
#include "tudat/astro/propagators/singleArcPropagationController.h"
#include "tudat/simulation/propagation_setup/createStateDerivativeModel.h"
#include "tudat/simulation/propagation_setup/propagationResults.h"
#include "tudat/simulation/propagation_setup/createEnvironmentUpdater.h"
#include "tudat/simulation/propagation_setup/propagationTermination.h"
#include "tudat/astro/propagators/dynamicsStateDerivativeModel.h"
#include "tudat/math/interpolators/lagrangeInterpolator.h"
#include "tudat/simulation/propagation_setup/dependentVariablesInterface.h"
#include "tudat/simulation/propagation_setup/getInitialStateOfBody.h"

namespace tudat
{

namespace propagators
{
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
                    propagationResults->getStateTransitionMatrixSize( ) + propagationResults->getPropagatedSensitivityMatrixSize( ) + 1;
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

    bool getSetIntegratedResult( );

    //! This function updates the environment with the numerical solution of the propagation.
    /*!
     *  This function updates the environment with the numerical solution of the propagation. For instance, it sets
     *  the propagated translational dynamics solution as the new input for the Ephemeris object of the body that was
     *  propagated. This function is pure virtual and must be implemented in the derived class.
     */
    virtual void processNumericalEquationsOfMotionSolution( ) = 0;

    virtual std::shared_ptr< SimulationResults< StateScalarType, TimeType > > getPropagationResults( ) = 0;

    std::shared_ptr< DependentVariablesInterface< TimeType > > getDependentVariablesInterface( );

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

#if TUDAT_BUILD_TESTS_WITH_EXPLICIT_INSTANTIATIONS
extern template class DynamicsSimulator< double, double >;
#endif

}  // namespace propagators

}  // namespace tudat

#include "tudat/simulation/propagation_setup/dynamicsSimulatorBaseImplementation.h"

#endif  // TUDAT_DYNAMICSSIMULATOR_BASE_H
