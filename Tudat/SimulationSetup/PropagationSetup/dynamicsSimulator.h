/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include <boost/make_shared.hpp>
#include <boost/assign/list_of.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h"
#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"
#include "Tudat/Basics/utilities.h"
#include "Tudat/Astrodynamics/Propagators/nBodyStateDerivative.h"
#include "Tudat/Astrodynamics/Ephemerides/frameManager.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/compositeEphemeris.h"
#include "Tudat/Astrodynamics/Propagators/nBodyCowellStateDerivative.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationSettings.h"
#include "Tudat/SimulationSetup/PropagationSetup/setNumericallyIntegratedStates.h"
#include "Tudat/Astrodynamics/Propagators/integrateEquations.h"
#include "Tudat/SimulationSetup/PropagationSetup/createStateDerivativeModel.h"
#include "Tudat/SimulationSetup/PropagationSetup/createEnvironmentUpdater.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationTermination.h"
#include "Tudat/Astrodynamics/Propagators/dynamicsStateDerivativeModel.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"

namespace tudat
{

namespace propagators
{

//! Function to get the states of a set of bodies, w.r.t. some set of central bodies, at the requested time.
/*!
* Function to get the states of a set of bodies, w.r.t. some set of central bodies, at the requested time.
* \param bodiesToIntegrate List of bodies for which to retrieve state.
* \param centralBodies Origins w.r.t. which to retrieve states of bodiesToIntegrate.
* \param bodyMap List of bodies to use in simulations.
* \param initialTime Time at which to retrieve states.
* \param frameManager OBject with which to calculate frame origin translations.
* \return Initial state vector (with 6 Cartesian elements per body, in order of bodiesToIntegrate vector).
*/
template< typename TimeType = double, typename StateScalarType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialStatesOfBodies(
        const std::vector< std::string >& bodiesToIntegrate,
        const std::vector< std::string >& centralBodies,
        const simulation_setup::NamedBodyMap& bodyMap,
        const TimeType initialTime,
        const boost::shared_ptr< ephemerides::ReferenceFrameManager > frameManager )
{
    // Set initial states of bodies to integrate.
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > systemInitialState =
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( bodiesToIntegrate.size( ) * 6, 1 );
    boost::shared_ptr< ephemerides::Ephemeris > ephemerisOfCurrentBody;

    // Iterate over all bodies.
    for( unsigned int i = 0; i < bodiesToIntegrate.size( ) ; i++ )
    {
        ephemerisOfCurrentBody = bodyMap.at( bodiesToIntegrate.at( i ) )->getEphemeris( );

        if ( ! ephemerisOfCurrentBody )
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
            boost::shared_ptr< ephemerides::Ephemeris > correctionEphemeris =
                    frameManager->getEphemeris( ephemerisOfCurrentBody->getReferenceFrameOrigin( ), centralBodies.at( i ) );
            systemInitialState.segment( i * 6 , 6 ) -= correctionEphemeris->getTemplatedStateFromEphemeris<
                    StateScalarType, TimeType >( initialTime );
        }
    }
    return systemInitialState;
}


boost::shared_ptr< ephemerides::ReferenceFrameManager > createFrameManager(
        const simulation_setup::NamedBodyMap& bodyMap );

//! Function to get the states of a set of bodies, w.r.t. some set of central bodies, at the requested time.
/*!
* Function to get the states of a set of bodies, w.r.t. some set of central bodies, at the requested time, creates
* frameManager from input data.
* \param bodiesToIntegrate List of bodies for which to retrieve state.
* \param centralBodies Origins w.r.t. which to retrieve states of bodiesToIntegrate.
* \param bodyMap List of bodies to use in simulations.
* \param initialTime Time at which to retrieve states.
* \return Initial state vector (with 6 Cartesian elements per body, in order of bodiesToIntegrate vector).
*/
template< typename TimeType = double, typename StateScalarType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialStatesOfBodies(
        const std::vector< std::string >& bodiesToIntegrate,
        const std::vector< std::string >& centralBodies,
        const  simulation_setup::NamedBodyMap& bodyMap,
        const TimeType initialTime )
{
    // Create ReferenceFrameManager and call overloaded function.
    return getInitialStatesOfBodies< TimeType, StateScalarType >(
                bodiesToIntegrate, centralBodies, bodyMap, initialTime,
                createFrameManager( bodyMap ) );
}

//! Function to get the states of single body, w.r.t. some central body, at the requested time.
/*!
* Function to get the states of  single body, w.r.t. some central body, at the requested time. This function creates
* frameManager from input data to perform all required conversions.
* \param bodyToIntegrate Body for which to retrieve state
* \param centralBody Origin w.r.t. which to retrieve state of bodyToIntegrate.
* \param bodyMap List of bodies to use in simulations.
* \param initialTime Time at which to retrieve state.
* \return Initial state vector of bodyToIntegrate
*/
template< typename TimeType = double, typename StateScalarType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialStateOfBody(
        const std::string& bodyToIntegrate,
        const std::string& centralBody,
        const  simulation_setup::NamedBodyMap& bodyMap,
        const TimeType initialTime )
{
    return getInitialStatesOfBodies< TimeType, StateScalarType >(
                boost::assign::list_of( bodyToIntegrate ), boost::assign::list_of( centralBody ), bodyMap, initialTime );
}

//! Function to get the state of single body, w.r.t. some central body, at a set of requested times, concatanated into one vector.
/*!
* Function to get the states of  single body, w.r.t. some central body, at a set of requested times, concatanated into one vector.
* This function creates frameManager from input data to perform all required conversions.
* \param bodyToIntegrate Body for which to retrieve state
* \param centralBody Origin w.r.t. which to retrieve state of bodyToIntegrate.
* \param bodyMap List of bodies to use in simulations.
* \param arcStartTimes List of times at which to retrieve states.
* \return Initial state vectosr of bodyToIntegrate at requested times.
*/
template< typename TimeType = double, typename StateScalarType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialArcWiseStateOfBody(
        const std::string& bodyToIntegrate,
        const std::string& centralBody,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::vector< TimeType > arcStartTimes )
{
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialStates = Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero(
                6 * arcStartTimes.size( ), 1 );
    for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
    {
        initialStates.block( 6 * i, 0, 6, 1 ) = getInitialStateOfBody< double, StateScalarType >(
                    bodyToIntegrate, centralBody, bodyMap, arcStartTimes.at( i ) );
    }
    return initialStates;
}

//! Base class for performing full numerical integration of a dynamical system.
/*!
 *  Base class for performing full numerical integration of a dynamical system. Governing equations are set once,
 *  but can be re-integrated for different initial conditions using the same instance of the class.
 *  Derived classes define the specific kind of integration that is performed
 *  (single-arc/multi-arc/etc.)
 */
template< typename StateScalarType = double, typename TimeType = double >
class DynamicsSimulator
{
public:

    //! Constructor of simulator.
    /*!
     *  Constructor of simulator, constructs integrator and object for calculating each time step of integration.
     *  \param bodyMap Map of bodies (with names) of all bodies in integration.
     *  \param clearNumericalSolutions Boolean to determine whether to clear the raw numerical solution member variables
     *  after propagation and resetting ephemerides (default true).
     *  \param setIntegratedResult Boolean to determine whether to automatically use the integrated results to set
     *  ephemerides (default true).
     */
    DynamicsSimulator(
            const simulation_setup::NamedBodyMap& bodyMap,
            const bool clearNumericalSolutions = true,
            const bool setIntegratedResult = true ):
        bodyMap_( bodyMap ),
        clearNumericalSolutions_( clearNumericalSolutions ),
        setIntegratedResult_( setIntegratedResult ){ }

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

    virtual std::vector< std::map< TimeType, double > > getCummulativeComputationTimeHistoryBase( ) = 0;


    //! Function to get the map of named bodies involved in simulation.
    /*!
     *  Function to get the map of named bodies involved in simulation.
     *  \return Map of named bodies involved in simulation.
     */
    simulation_setup::NamedBodyMap getNamedBodyMap( )
    {
        return bodyMap_;
    }

    //! Function to reset the named body map.
    /*!
     *  Function to reset the named body map.
     *  \param bodyMap The new named body map.
     */
    void resetNamedBodyMap( const simulation_setup::NamedBodyMap& bodyMap )
    {
        bodyMap_ = bodyMap;
    }

    void resetSetIntegratedResult( const bool setIntegratedResult )
    {
        setIntegratedResult_ = setIntegratedResult;
    }

protected:

    //! This function updates the environment with the numerical solution of the propagation.
    /*!
     *  This function updates the environment with the numerical solution of the propagation. For instance, it sets
     *  the propagated translational dynamics solution as the new input for the Ephemeris object of the body that was
     *  propagated. This function is pure virtual and must be implemented in the derived class.
     */
    virtual void processNumericalEquationsOfMotionSolution( ) = 0;

    //!  Map of bodies (with names) of all bodies in integration.
    simulation_setup::NamedBodyMap bodyMap_;

    //! Boolean to determine whether to clear the raw numerical solution member variables after propagation and
    //! resetting ephemerides.
    bool clearNumericalSolutions_;

    //! Boolean to determine whether to automatically use the integrated results to set ephemerides.
    bool setIntegratedResult_;
};

//! Class for performing full numerical integration of a dynamical system in a single arc.
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

    using DynamicsSimulator< StateScalarType, TimeType >::bodyMap_;
    using DynamicsSimulator< StateScalarType, TimeType >::clearNumericalSolutions_;
    using DynamicsSimulator< StateScalarType, TimeType >::setIntegratedResult_;


    //! Constructor of simulator.
    /*!
     *  Constructor of simulator, constructs integrator and object for calculating each time step of integration.
     *  \param bodyMap Map of bodies (with names) of all bodies in integration.
     *  \param integratorSettings Settings for numerical integrator.
     *  \param propagatorSettings Settings for propagator.
     *  \param areEquationsOfMotionToBeIntegrated Boolean to denote whether equations of motion should be integrated
     *  immediately at the end of the contructor or not (default true).
     *  \param clearNumericalSolutions Boolean to determine whether to clear the raw numerical solution member variables
     *  after propagation and resetting ephemerides (default false).
     *  \param setIntegratedResult Boolean to determine whether to automatically use the integrated results to set
     *  ephemerides (default false).
     *  \param initialClockTime Initial clock time from which to determine cummulative computation time.
     *  By default now(), i.e. the moment at which this function is called.
     */
    SingleArcDynamicsSimulator(
            const simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const bool areEquationsOfMotionToBeIntegrated = true,
            const bool clearNumericalSolutions = false,
            const bool setIntegratedResult = false,
            const std::chrono::steady_clock::time_point initialClockTime = std::chrono::steady_clock::now( ) ):
        DynamicsSimulator< StateScalarType, TimeType >(
            bodyMap, clearNumericalSolutions, setIntegratedResult ),
        integratorSettings_( integratorSettings ),
        propagatorSettings_(
            boost::dynamic_pointer_cast< SingleArcPropagatorSettings< StateScalarType > >( propagatorSettings ) ),
        initialPropagationTime_( integratorSettings_->initialTime_ ), initialClockTime_( initialClockTime ),
        propagationTerminationReason_( boost::make_shared< PropagationTerminationDetails >( propagation_never_run ) )
    {
        if( propagatorSettings == NULL )
        {
            throw std::runtime_error( "Error in dynamics simulator, propagator settings not defined" );
        }
        else if( boost::dynamic_pointer_cast< SingleArcPropagatorSettings< StateScalarType > >( propagatorSettings ) == NULL )
        {
            throw std::runtime_error( "Error in dynamics simulator, input must be single-arc" );
        }

        if( integratorSettings == NULL )
        {
            throw std::runtime_error( "Error in dynamics simulator, integrator settings not defined" );
        }

        if( setIntegratedResult_ )
        {
            frameManager_ = createFrameManager( bodyMap );
            integratedStateProcessors_ = createIntegratedStateProcessors< TimeType, StateScalarType >(
                        propagatorSettings_, bodyMap_, frameManager_ );
        }

        environmentUpdater_ = createEnvironmentUpdaterForDynamicalEquations< StateScalarType, TimeType >(
                    propagatorSettings_, bodyMap_ );
        dynamicsStateDerivative_ = boost::make_shared< DynamicsStateDerivativeModel< TimeType, StateScalarType > >(
                    createStateDerivativeModels< StateScalarType, TimeType >(
                        propagatorSettings_, bodyMap_, initialPropagationTime_ ),
                    boost::bind( &EnvironmentUpdater< StateScalarType, TimeType >::updateEnvironment,
                                 environmentUpdater_, _1, _2, _3 ) );
        propagationTerminationCondition_ = createPropagationTerminationConditions(
                    propagatorSettings_->getTerminationSettings( ), bodyMap_, integratorSettings->initialTimeStep_ );

        if( propagatorSettings_->getDependentVariablesToSave( ) != NULL )
        {
            std::pair< boost::function< Eigen::VectorXd( ) >, std::map< int, std::string > > dependentVariableData =
                    createDependentVariableListFunction< TimeType, StateScalarType >(
                        propagatorSettings_->getDependentVariablesToSave( ), bodyMap_,
                        dynamicsStateDerivative_->getStateDerivativeModels( ) );
            dependentVariablesFunctions_ = dependentVariableData.first;
            dependentVariableIds_ = dependentVariableData.second;

            if( propagatorSettings_->getDependentVariablesToSave( )->printDependentVariableTypes_ )
            {
                std::cout << "Dependent variables being saved, output vectors contain: " << std::endl
                          << "Vector entry, Vector contents" << std::endl;
                utilities::printMapContents(
                            dependentVariableIds_ );
            }
        }

        stateDerivativeFunction_ =
                boost::bind( &DynamicsStateDerivativeModel< TimeType, StateScalarType >::computeStateDerivative,
                             dynamicsStateDerivative_, _1, _2 );
        doubleStateDerivativeFunction_ =
                boost::bind( &DynamicsStateDerivativeModel< TimeType, StateScalarType >::computeStateDoubleDerivative,
                             dynamicsStateDerivative_, _1, _2 );

        // Integrate equations of motion if required.
        if( areEquationsOfMotionToBeIntegrated )
        {
            integrateEquationsOfMotion( propagatorSettings_->getInitialStates( ) );
        }
    }

    //! Destructor
    ~SingleArcDynamicsSimulator( )
    { }

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

        equationsOfMotionNumericalSolution_.clear( );
        equationsOfMotionNumericalSolutionRaw_.clear( );

        dynamicsStateDerivative_->setPropagationSettings( std::vector< IntegratedStateType >( ), 1, 0 );
        dynamicsStateDerivative_->resetFunctionEvaluationCounter( );

        // Reset initial time to ensure consistency with multi-arc propagation.
        integratorSettings_->initialTime_ = this->initialPropagationTime_;

        // Integrate equations of motion numerically.
        propagationTerminationReason_ =
                EquationIntegrationInterface< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >, TimeType >::integrateEquations(
                    stateDerivativeFunction_, equationsOfMotionNumericalSolutionRaw_,
                    dynamicsStateDerivative_->convertFromOutputSolution(
                        initialStates, this->initialPropagationTime_ ), integratorSettings_,
                    propagationTerminationCondition_,
                    dependentVariableHistory_,
                    cummulativeComputationTimeHistory_,
                    dependentVariablesFunctions_,
                    propagatorSettings_->getPrintInterval( ),
                    initialClockTime_ );
        dynamicsStateDerivative_->convertNumericalStateSolutionsToOutputSolutions(
                    equationsOfMotionNumericalSolution_, equationsOfMotionNumericalSolutionRaw_ );

        if( this->setIntegratedResult_ )
        {
            processNumericalEquationsOfMotionSolution( );
        }
    }

    //! Function to return the map of state history of numerically integrated bodies.
    /*!
     * Function to return the map of state history of numerically integrated bodies.
     * \return Map of state history of numerically integrated bodies.
     */
    std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > getEquationsOfMotionNumericalSolution( )
    {
        return equationsOfMotionNumericalSolution_;
    }

    //! Function to return the map of dependent variable history that was saved during numerical propagation.
    /*!
     * Function to return the map of dependent variable history that was saved during numerical propagation.
     * \return Map of dependent variable history that was saved during numerical propagation.
     */
    std::map< TimeType, Eigen::VectorXd > getDependentVariableHistory( )
    {
        return dependentVariableHistory_;
    }

    //! Function to return the map of cummulative computation time history that was saved during numerical propagation.
    /*!
     * Function to return the map of cummulative computation time history that was saved during numerical propagation.
     * \return Map of cummulative computation time history that was saved during numerical propagation.
     */
    std::map< TimeType, double > getCummulativeComputationTimeHistory( )
    {
        return cummulativeComputationTimeHistory_;
    }

    //! Function to return the map of state history of numerically integrated bodies (base class interface).
    /*!
     * Function to return the map of state history of numerically integrated bodies (base class interface).
     * \return Vector is size 1, with entry: map of state history of numerically integrated bodies.
     */
    std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > getEquationsOfMotionNumericalSolutionBase( )
    {
        return std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >(
                    { getEquationsOfMotionNumericalSolution( ) } );
    }

    //! Function to return the map of dependent variable history that was saved during numerical propagation(base class interface)
    /*!
     * Function to return the map of dependent variable history that was saved during numerical propagation (base class interface)
     * \return Vector is size 1, with entry: map of dependent variable history that was saved during numerical propagation.
     */
    std::vector< std::map< TimeType, Eigen::VectorXd > > getDependentVariableNumericalSolutionBase( )
    {
        return std::vector< std::map< TimeType, Eigen::VectorXd > >(
                    { getDependentVariableHistory( ) } );
    }

    std::vector< std::map< TimeType, double > > getCummulativeComputationTimeHistoryBase( )
    {
        return std::vector< std::map< TimeType, double > >( { getCummulativeComputationTimeHistory( ) } );
    }


    //! Function to reset the environment from an externally generated state history.
    /*!
     * Function to reset the environment from an externally generated state history, the order of the entries in the
     * state vectors are proscribed by propagatorSettings
     * \param equationsOfMotionNumericalSolution Externally generated state history.
     * \param dependentVariableHistory Externally generated dependent variable history.
     */
    void manuallySetAndProcessRawNumericalEquationsOfMotionSolution(
            const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >&
            equationsOfMotionNumericalSolution,
            const std::map< TimeType, Eigen::VectorXd >& dependentVariableHistory)
    {
        equationsOfMotionNumericalSolution_ = equationsOfMotionNumericalSolution;
        dependentVariableHistory_ = dependentVariableHistory;
        processNumericalEquationsOfMotionSolution( );
    }

    //! Function to get the settings for the numerical integrator.
    /*!
     * Function to get the settings for the numerical integrator.
     * \return The settings for the numerical integrator.
     */
    boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > getIntegratorSettings( )
    {
        return integratorSettings_;
    }

    //! Function to get the function that performs a single state derivative function evaluation.
    /*!
     * Function to get the function that performs a single state derivative function evaluation.
     * \return Function that performs a single state derivative function evaluation.
     */
    boost::function< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >
    ( const TimeType, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >&) >
    getStateDerivativeFunction( )
    {
        return stateDerivativeFunction_;
    }

    //! Function to get the function that performs a single state derivative function evaluation with double precision.
    /*!
     * Function to get the function that performs a single state derivative function evaluation with double precision,
     * regardless of template arguments.
     * \return Function that performs a single state derivative function evaluation with double precision.
     */
    boost::function< Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >
    ( const double, const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& ) > getDoubleStateDerivativeFunction( )
    {
        return doubleStateDerivativeFunction_;
    }

    //! Function to get the settings for the propagator.
    /*!
     * Function to get the settings for the propagator.
     * \return The settings for the propagator.
     */
    boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > getPropagatorSettings( )
    {
        return propagatorSettings_;
    }

    //! Function to get the object that updates the environment.
    /*!
     * Function to get the object responsible for updating the environment based on the current state and time.
     * \return Object responsible for updating the environment based on the current state and time.
     */
    boost::shared_ptr< EnvironmentUpdater< StateScalarType, TimeType > > getEnvironmentUpdater( )
    {
        return environmentUpdater_;
    }

    //! Function to get the object that updates and returns state derivative
    /*!
     * Function to get the object that updates current environment and returns state derivative from single function call
     * \return Object that updates current environment and returns state derivative from single function call
     */
    boost::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > getDynamicsStateDerivative( )
    {
        return dynamicsStateDerivative_;
    }


    //! Function to retrieve the object defining when the propagation is to be terminated.
    /*!
     * Function to retrieve the object defining when the propagation is to be terminated.
     * \return Object defining when the propagation is to be terminated.
     */
    boost::shared_ptr< PropagationTerminationCondition > getPropagationTerminationCondition( )
    {
        return propagationTerminationCondition_;
    }

    //! Function to retrieve the list of object that process the integrated numerical solution by updating the environment
    /*!
     * Function to retrieve the List of object (per dynamics type) that process the integrated numerical solution by
     * updating the environment
     * \return List of object (per dynamics type) that process the integrated numerical solution by updating the environment
     */
    std::map< IntegratedStateType, std::vector< boost::shared_ptr<
    IntegratedStateProcessor< TimeType, StateScalarType > > > > getIntegratedStateProcessors( )
    {
        return integratedStateProcessors_;
    }


    //! Function to retrieve the event that triggered the termination of the last propagation
    /*!
     * Function to retrieve the event that triggered the termination of the last propagation
     * \return Event that triggered the termination of the last propagation
     */
    boost::shared_ptr< PropagationTerminationDetails > getPropagationTerminationReason( )
    {
        return propagationTerminationReason_;
    }

    //! Get whether the integration was completed successfully.
    /*!
     * @copybrief integrationCompletedSuccessfully
     * \return Whether the integration was completed successfully by reaching the termination condition.
     */
    virtual bool integrationCompletedSuccessfully( ) const
    {
        return ( propagationTerminationReason_->getPropagationTerminationReason( ) == termination_condition_reached );
    }


    //! Function to retrieve the dependent variables IDs
    /*!
     * Function to retrieve the dependent variables IDs
     * \return Map listing starting entry of dependent variables in output vector, along with associated ID
     */
    std::map< int, std::string > getDependentVariableIds( )
    {
        return dependentVariableIds_;
    }


    //! Function to retrieve initial time of propagation
    /*!
     * Function to retrieve initial time of propagation
     * \return Initial time of propagation
     */
    double getInitialPropagationTime( )
    {
        return this->initialPropagationTime_;
    }

    //! Function to retrieve the functions that compute the dependent variables at each time step
    /*!
     * Function to retrieve the functions that compute the dependent variables at each time step
     * \return Functions that compute the dependent variables at each time step
     */
    boost::function< Eigen::VectorXd( ) > getDependentVariablesFunctions( )
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
                    propagatorSettings_->getTerminationSettings(), bodyMap_,
                            integratorSettings_->initialTimeStep_ );
    }

protected:


    //! This function updates the environment with the numerical solution of the propagation.
    /*!
     *  This function updates the environment with the numerical solution of the propagation. It sets
     *  the propagated translational dynamics solution as the new input for the Ephemeris object of the body that was
     *  propagated.
     */
    void processNumericalEquationsOfMotionSolution( )
    {
        // Create and set interpolators for ephemerides
        resetIntegratedStates( equationsOfMotionNumericalSolution_, integratedStateProcessors_ );


        // Clear numerical solution if so required.
        if( clearNumericalSolutions_ )
        {
            equationsOfMotionNumericalSolution_.clear( );
        }

        for( simulation_setup::NamedBodyMap::const_iterator
             bodyIterator = bodyMap_.begin( );
             bodyIterator != bodyMap_.end( ); bodyIterator++ )
        {
            bodyIterator->second->updateConstantEphemerisDependentMemberQuantities( );
        }
    }


    //! List of object (per dynamics type) that process the integrated numerical solution by updating the environment
    std::map< IntegratedStateType, std::vector< boost::shared_ptr<
    IntegratedStateProcessor< TimeType, StateScalarType > > > > integratedStateProcessors_;

    //! Object responsible for updating the environment based on the current state and time.
    /*!
     *  Object responsible for updating the environment based on the current state and time. Calling the updateEnvironment
     * function automatically updates all dependent variables that are needed to calulate the state derivative.
     */
    boost::shared_ptr< EnvironmentUpdater< StateScalarType, TimeType > > environmentUpdater_;

    //! Interface object that updates current environment and returns state derivative from single function call.
    boost::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > dynamicsStateDerivative_;

    //! Function that performs a single state derivative function evaluation.
    /*!
     *  Function that performs a single state derivative function evaluation, will typically be set to
     *  DynamicsStateDerivativeModel< TimeType, StateScalarType >::computeStateDerivative function.
     *  Calling this function will first update the environment (using environmentUpdater_) and then calculate the
     *  full system state derivative.
     */
    boost::function< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >
    ( const TimeType, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& ) > stateDerivativeFunction_;

    //! Function that performs a single state derivative function evaluation with double precision.
    /*!
     *  Function that performs a single state derivative function evaluation with double precision
     *  \sa stateDerivativeFunction_
     */
    boost::function< Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >
    ( const double, const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& ) > doubleStateDerivativeFunction_;


    //! Settings for numerical integrator.
    boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings_;

    //! Settings for propagator.
    boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > propagatorSettings_;

    //! Object defining when the propagation is to be terminated.
    boost::shared_ptr< PropagationTerminationCondition > propagationTerminationCondition_;

    //! Function returning dependent variables (during numerical propagation)
    boost::function< Eigen::VectorXd( ) > dependentVariablesFunctions_;

    //! Map listing starting entry of dependent variables in output vector, along with associated ID.
    std::map< int, std::string > dependentVariableIds_;

    //! Object for retrieving ephemerides for transformation of reference frame (origins)
    boost::shared_ptr< ephemerides::ReferenceFrameManager > frameManager_;

    //! Map of state history of numerically integrated bodies.
    /*!
     *  Map of state history of numerically integrated bodies, i.e. the result of the numerical integration, transformed
     *  into the 'conventional form' (\sa SingleStateTypeDerivative::convertToOutputSolution). Key of map denotes time,
     *  values are concatenated vectors of integrated body states (order defined by propagatorSettings_).
     *  NOTE: this map is empty if clearNumericalSolutions_ is set to true.
     */
    std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > equationsOfMotionNumericalSolution_;

    std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > equationsOfMotionNumericalSolutionRaw_;

    //! Map of dependent variable history that was saved during numerical propagation.
    std::map< TimeType, Eigen::VectorXd > dependentVariableHistory_;

    //! Map of cummulative computation time history that was saved during numerical propagation.
    std::map< TimeType, double > cummulativeComputationTimeHistory_;

    //! Initial time of propagation
    double initialPropagationTime_;

    //!
    std::chrono::steady_clock::time_point initialClockTime_;

    //! Event that triggered the termination of the propagation
    boost::shared_ptr< PropagationTerminationDetails > propagationTerminationReason_;

};

//! Function to get a vector of initial states from a vector of propagator settings
/*!
 *  Function to get a vector of initial states from a vector of propagator settings.
 *  \param propagatorSettings List of propagator settings
 *  \return List of initial states, as retrieved from propagatorSettings list.
 */
template< typename StateScalarType = double >
std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1  > > getInitialStatesPerArc(
        const std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > propagatorSettings )
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
 *  \param previousArcDynamisSolution Numerical solution of previous arc
 *  \param currentArcInitialTime Start time of current arc
 *  \return Interpolated initial state of current arc
 */
template< typename StateScalarType = double, typename TimeType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getArcInitialStateFromPreviousArcResult(
        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& previousArcDynamisSolution,
        const double currentArcInitialTime )
{
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentArcInitialState;
    {
        // Check if overlap exists
        if( previousArcDynamisSolution.rbegin( )->first < currentArcInitialTime )
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
                 const_reverse_iterator previousArcIterator = previousArcDynamisSolution.rbegin( );
                 previousArcIterator != previousArcDynamisSolution.rend( ); previousArcIterator++ )
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
                    boost::make_shared< interpolators::LagrangeInterpolator<
                    TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >, long double > >(
                        initialStateInterpolationMap, 8 )->interpolate( currentArcInitialTime );

        }
    }
    return currentArcInitialState;
}

//! Class for performing full numerical integration of a dynamical system over multiple arcs.
/*!
 *  Class for performing full numerical integration of a dynamical system over multiple arcs, equations of motion are set up
 *  for each arc (and need not be equal for each arc). In this class, the governing equations are set once,
 *  but can be re-integrated for different initial conditions using the same instance of the class.
 */
template< typename StateScalarType = double, typename TimeType = double >
class MultiArcDynamicsSimulator: public DynamicsSimulator< StateScalarType, TimeType >
{
public:

    using DynamicsSimulator< StateScalarType, TimeType >::bodyMap_;
    using DynamicsSimulator< StateScalarType, TimeType >::clearNumericalSolutions_;

    //! Constructor of multi-arc simulator for same integration settings per arc.
    /*!
     *  Constructor of multi-arc simulator for same integration settings per arc.
     *  \param bodyMap Map of bodies (with names) of all bodies in integration.
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
            const simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const std::vector< double > arcStartTimes,
            const bool areEquationsOfMotionToBeIntegrated = true,
            const bool clearNumericalSolutions = true,
            const bool setIntegratedResult = true ):
        DynamicsSimulator< StateScalarType, TimeType >(
            bodyMap, clearNumericalSolutions, setIntegratedResult )
    {
        multiArcPropagatorSettings_ =
                boost::dynamic_pointer_cast< MultiArcPropagatorSettings< StateScalarType > >( propagatorSettings );
        if( multiArcPropagatorSettings_ == NULL )
        {
            throw std::runtime_error( "Error when creating multi-arc dynamics simulator, input is not multi arc" );
        }
        else
        {
            std::vector< boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > singleArcSettings =
                    multiArcPropagatorSettings_->getSingleArcSettings( );

            arcStartTimes_.resize( arcStartTimes.size( ) );

            if( singleArcSettings.size( ) != arcStartTimes.size( ) )
            {
                throw std::runtime_error( "Error when creating multi-arc dynamics simulator, input is inconsistent" );
            }
            // Create dynamics simulators
            for( unsigned int i = 0; i < singleArcSettings.size( ); i++ )
            {
                integratorSettings->initialTime_ = arcStartTimes.at( i );

                singleArcDynamicsSimulators_.push_back(
                            boost::make_shared< SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                                bodyMap, integratorSettings, singleArcSettings.at( i ), false, false, true ) );
                singleArcDynamicsSimulators_[ i ]->resetSetIntegratedResult( false );
            }

            equationsOfMotionNumericalSolution_.resize( arcStartTimes.size( ) );
            dependentVariableHistory_.resize( arcStartTimes.size( ) );
            cummulativeComputationTimeHistory_.resize( arcStartTimes.size( ) );
            propagationTerminationReasons_.resize( arcStartTimes.size( ) );

            // Integrate equations of motion if required.
            if( areEquationsOfMotionToBeIntegrated )
            {
                integrateEquationsOfMotion( multiArcPropagatorSettings_->getInitialStates( ) );
            }
        }
    }

    //! Constructor of multi-arc simulator for different integration settings per arc.
    /*!
         *  Constructor of multi-arc simulator for different integration settings per arc.
         *  \param bodyMap Map of bodies (with names) of all bodies in integration.
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
            const simulation_setup::NamedBodyMap& bodyMap,
            const std::vector< boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > > integratorSettings,
            const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const bool areEquationsOfMotionToBeIntegrated = true,
            const bool clearNumericalSolutions = true,
            const bool setIntegratedResult = true ):
        DynamicsSimulator< StateScalarType, TimeType >(
            bodyMap, clearNumericalSolutions, setIntegratedResult )
    {
        multiArcPropagatorSettings_ =
                boost::dynamic_pointer_cast< MultiArcPropagatorSettings< StateScalarType > >( propagatorSettings );
        if( multiArcPropagatorSettings_ == NULL )
        {
            throw std::runtime_error( "Error when creating multi-arc dynamics simulator, input is not multi arc" );
        }
        else
        {
            std::vector< boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > singleArcSettings =
                    multiArcPropagatorSettings_->getSingleArcSettings( );

            if( singleArcSettings.size( ) != integratorSettings.size( ) )
            {
                throw std::runtime_error( "Error when creating multi-arc dynamics simulator, input sizes are inconsistent" );
            }

            arcStartTimes_.resize( singleArcSettings.size( ) );

            // Create dynamics simulators
            for( unsigned int i = 0; i < singleArcSettings.size( ); i++ )
            {
                singleArcDynamicsSimulators_.push_back(
                            boost::make_shared< SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                                bodyMap, integratorSettings.at( i ), singleArcSettings.at( i ), false, false, true ) );
                singleArcDynamicsSimulators_[ i ]->resetSetIntegratedResult( false );
            }

            equationsOfMotionNumericalSolution_.resize( singleArcSettings.size( ) );
            dependentVariableHistory_.resize( singleArcSettings.size( ) );
            cummulativeComputationTimeHistory_.resize( singleArcSettings.size( ) );
            propagationTerminationReasons_.resize( singleArcSettings.size( ) );

            // Integrate equations of motion if required.
            if( areEquationsOfMotionToBeIntegrated )
            {
                integrateEquationsOfMotion( multiArcPropagatorSettings_->getInitialStates( ) );
            }
        }
    }

    //! Destructor
    ~MultiArcDynamicsSimulator( ) { }

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
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& concatenatedInitialStates )
    {
        std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > splitInitialState;

        int currentIndex = 0;
        for( unsigned int i = 0; i < singleArcDynamicsSimulators_.size( ); i++ )
        {
            int currentSize = singleArcDynamicsSimulators_.at( i )->getPropagatorSettings( )->getStateSize( );
            splitInitialState.push_back( concatenatedInitialStates.block( currentIndex, 0, currentSize, 1 ) );
            currentIndex += currentSize;
        }

        if( currentIndex != concatenatedInitialStates.rows( ) )
        {
            throw std::runtime_error( "Error when doing multi-arc integration, input state vector size is incompatible with settings" );
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
            const std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& initialStatesList )
    {
        // Clear existing solution (if any)
        for( unsigned int i = 0; i < equationsOfMotionNumericalSolution_.size( ); i++ )
        {
            equationsOfMotionNumericalSolution_.at( i ).clear( );
        }

        for( unsigned int i = 0; i < dependentVariableHistory_.size( ); i++ )
        {
            dependentVariableHistory_.at( i ).clear( );
        }

        for( unsigned int i = 0; i < cummulativeComputationTimeHistory_.size( ); i++ )
        {
            cummulativeComputationTimeHistory_.at( i ).clear( );
        }


        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentArcInitialState;
        std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > arcInitialStateList;
        bool updateInitialStates = false;

        // Propagate dynamics for each arc
        for( unsigned int i = 0; i < singleArcDynamicsSimulators_.size( ); i++ )
        {
            // Get arc initial state. If initial state is NaN, this signals that the initial state is to be taken from previous
            // arc
            if( ( i == 0 ) || ( !linear_algebra::doesMatrixHaveNanEntries( initialStatesList.at( i ) ) ) )
            {
                currentArcInitialState = initialStatesList.at( i );
            }
            else
            {
                currentArcInitialState = getArcInitialStateFromPreviousArcResult(
                            equationsOfMotionNumericalSolution_.at( i - 1 ),
                            singleArcDynamicsSimulators_.at( i )->getInitialPropagationTime( ) );

                // If arc initial state is taken from previous arc, this indicates that the initial states in propagator settings
                // need to be updated.
                updateInitialStates = true;
            }
            arcInitialStateList.push_back( currentArcInitialState );

            singleArcDynamicsSimulators_.at( i )->integrateEquationsOfMotion( currentArcInitialState );
            equationsOfMotionNumericalSolution_[ i ] =
                    singleArcDynamicsSimulators_.at( i )->getEquationsOfMotionNumericalSolution( );
            dependentVariableHistory_[ i ] =
                    singleArcDynamicsSimulators_.at( i )->getDependentVariableHistory( );
            cummulativeComputationTimeHistory_[ i ] =
                    singleArcDynamicsSimulators_.at( i )->getCummulativeComputationTimeHistory( );
            propagationTerminationReasons_[ i ] = singleArcDynamicsSimulators_.at( i )->getPropagationTerminationReason( );
            arcStartTimes_[ i ] = equationsOfMotionNumericalSolution_[ i ].begin( )->first;
        }


        if( updateInitialStates )
        {
            multiArcPropagatorSettings_->resetInitialStatesList(
                        arcInitialStateList );
        }

        if( this->setIntegratedResult_ )
        {
            processNumericalEquationsOfMotionSolution( );
        }
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
        return equationsOfMotionNumericalSolution_;
    }

    //! Function to return the numerical solution of the dependent variables
    /*!
     *  Function to return the numerical solution of the dependent variables for last numerical integration. Each vector entry
     *  denotes one arc. Key of map denotes time, values are dependent variable vectors
     *  \return List of maps of dependent variable history
     */
    std::vector< std::map< TimeType, Eigen::VectorXd > > getDependentVariableHistory( )
    {
        return dependentVariableHistory_;
    }

    std::vector< std::map< TimeType, double > > getCummulativeComputationTimeHistory( )
    {
        return cummulativeComputationTimeHistory_;
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

    std::vector< std::map< TimeType, double > > getCummulativeComputationTimeHistoryBase( )
    {
        return getCummulativeComputationTimeHistory( );
    }


    //! Function to reset the environment using an externally provided list of (numerically integrated) states
    /*!
     *  Function to reset the environment using an externally provided list of (numerically integrated) states, for instance
     *  provided by a variational equations solver.
     *  \param equationsOfMotionNumericalSolution Vector of state histories
     *  (externally provided equationsOfMotionNumericalSolution_)
     *  \param dependentVariableHistory Vector of dependent variable histories
     *  (externally provided dependentVariableHistory_)
     */
    void manuallySetAndProcessRawNumericalEquationsOfMotionSolution(
            std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >&
            equationsOfMotionNumericalSolution,
            std::vector< std::map< TimeType, Eigen::VectorXd > >&
            dependentVariableHistory)
    {
        // Set equationsOfMotionNumericalSolution_
        equationsOfMotionNumericalSolution_.resize( equationsOfMotionNumericalSolution.size( ) );

        for( unsigned int i = 0; i < equationsOfMotionNumericalSolution.size( ); i++ )
        {
            equationsOfMotionNumericalSolution_[ i ].clear( );
            equationsOfMotionNumericalSolution_[ i ] = equationsOfMotionNumericalSolution[ i ];
            arcStartTimes_[ i ] = equationsOfMotionNumericalSolution_[ i ].begin( )->first;

        }

        // Reset environment with new states.
        processNumericalEquationsOfMotionSolution( );

        dependentVariableHistory_.resize( dependentVariableHistory.size( ) );

        for( unsigned int i = 0; i < dependentVariableHistory.size( ); i++ )
        {
            dependentVariableHistory_[ i ].clear( );
            dependentVariableHistory_[ i ] = dependentVariableHistory[ i ];
        }
    }

    //! Function to get the list of DynamicsStateDerivativeModel objects used for each arc
    /*!
     * Function to get the list of DynamicsStateDerivativeModel objects used for each arc
     * \return List of DynamicsStateDerivativeModel objects used for each arc
     */
    std::vector< boost::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > > getDynamicsStateDerivative( )
    {
        std::vector< boost::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > > dynamicsStateDerivatives;
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
    std::vector< boost::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > > getSingleArcDynamicsSimulators( )
    {
        return singleArcDynamicsSimulators_;
    }

    //! Function to retrieve the current state and end times of the arcs
    /*!
     * Function to retrieve the current state and end times of the arcs
     * \return The current state and end times of the arcs
     */
    std::vector< double > getArcStartTimes( )
    {
        return arcStartTimes_;
    }

    //! Get whether the integration was completed successfully.
    /*!
     * @copybrief integrationCompletedSuccessfully
     * \return Whether the integration was completed successfully by reaching the termination condition.
     */
    virtual bool integrationCompletedSuccessfully( ) const
    {
        for ( const boost::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > >
              singleArcDynamicsSimulator : singleArcDynamicsSimulators_ )
        {
            if ( ! singleArcDynamicsSimulator->integrationCompletedSuccessfully( ) )
            {
                return false;
            }
        }
        return true;
    }


protected:

    //! This function updates the environment with the numerical solution of the propagation.
    /*!
     *  This function updates the environment with the numerical solution of the propagation. It sets
     *  the propagated dynamics solution as the new input for e.g., the ephemeris object of the boies that were
     *  propagated (for translational states).
     */
    void processNumericalEquationsOfMotionSolution( )
    {
        resetIntegratedMultiArcStatesWithEqualArcDynamics(
                    equationsOfMotionNumericalSolution_,
                    singleArcDynamicsSimulators_.at( 0 )->getIntegratedStateProcessors( ), arcStartTimes_ );

        if( clearNumericalSolutions_ )
        {
            for( unsigned int i = 0; i < equationsOfMotionNumericalSolution_.size( ); i++ )
            {
                equationsOfMotionNumericalSolution_.at( i ).clear( );
            }
            equationsOfMotionNumericalSolution_.clear( );
        }
    }

    //! List of maps of state history of numerically integrated states.
    /*!
     *  List of maps of state history of numerically integrated states. Each entry in the list contains data on a single arc.
     *  Key of map denotes time, values are concatenated vectors of body states in order of bodiesToIntegrate
     */
    std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > equationsOfMotionNumericalSolution_;

    //! List of maps of dependent variable history that was saved during numerical propagation.
    std::vector< std::map< TimeType, Eigen::VectorXd > > dependentVariableHistory_;

    std::vector< std::map< TimeType, double > > cummulativeComputationTimeHistory_;

    //! Objects used to compute the dynamics of the sepatrate arcs
    std::vector< boost::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > > singleArcDynamicsSimulators_;

    //! List of start times of each arc. NOTE: This list is updated after every propagation.
    std::vector< double > arcStartTimes_;

    //! Event that triggered the termination of the propagation
    std::vector< boost::shared_ptr< PropagationTerminationDetails > > propagationTerminationReasons_;

    //! Propagator settings used by this objec
    boost::shared_ptr< MultiArcPropagatorSettings< StateScalarType > > multiArcPropagatorSettings_;


};

} // namespace propagators

} // namespace tudat


#endif // TUDAT_DYNAMICSSIMULATOR_H
