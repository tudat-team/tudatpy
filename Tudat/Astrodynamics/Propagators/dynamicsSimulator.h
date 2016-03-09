#ifndef TUDAT_DYNAMICSSIMULATOR_H
#define TUDAT_DYNAMICSSIMULATOR_H

#include <boost/make_shared.hpp>
#include <boost/assign/list_of.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h"
#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"

#include "Tudat/Astrodynamics/Propagators/nBodyStateDerivative.h"
#include "Tudat/Astrodynamics/Ephemerides/frameManager.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/compositeEphemeris.h"
#include "Tudat/Astrodynamics/Propagators/nBodyCowellStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/propagationSettings.h"
#include "Tudat/Astrodynamics/Propagators/setNumericallyIntegratedStates.h"
#include "Tudat/Astrodynamics/Propagators/integrateEquations.h"
#include "Tudat/Astrodynamics/Propagators/createStateDerivativeModel.h"
#include "Tudat/Astrodynamics/Propagators/createEnvironmentUpdater.h"
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
        // Get body initial state from ephemeris
        systemInitialState.segment( i * 6 , 6 ) = ephemerisOfCurrentBody->getTemplatedStateFromEphemeris<
                StateScalarType, TimeType >( initialTime );

        // Correct initial state if integration origin and ephemeris origin are not equal.
        ephemerisOfCurrentBody = bodyMap.at( bodiesToIntegrate.at( i ) )->getEphemeris( );
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
    return getInitialStatesOfBodies( bodiesToIntegrate, centralBodies, bodyMap, initialTime,
                                     boost::make_shared< ephemerides::ReferenceFrameManager >( bodyMap ) );
}

//! Function to get the states of single, w.r.t. some set of central body, at the requested time.
/*!
* Function to get the states of single, w.r.t. some set of central body, at the requested time., creates
* frameManager from input data.
* \param bodyToIntegrate Bodies for which to retrieve state
* \param centralBody Origins w.r.t. which to retrieve state of bodyToIntegrate.
* \param bodyMap List of bodies to use in simulations.
* \param initialTime Time at which to retrieve states.
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

//! Base class for performing full numerical integration of a dynamical system.
/*!
 *  Base class for performing full numerical integration of a dynamical system. Governing equations are set once,
 *  but can be re-integrated for different initial conditions using the same instance of the class.
 *  Derived classes define the specific kind of integration that is performed
 *  (single-arc/multi-arc; dynamics/variational equations, etc.)
 */
template< typename StateScalarType = double, typename TimeType = double >
class DynamicsSimulator
{
public:

    //! Constructor of simulator.
    /*!
     *  Constructor of simulator, constructs integrator and object for calculating each time step of integration.
     *  \param bodyMap Map of bodies (with names) of all bodies in integration.
     *  \param integratorSettings Settings for numerical integrator.
     *  \param propagatorSettings Settings for propagator.
     *  \param clearNumericalSolutions Boolean to determine whether to clear the raw numerical solution member variables
     *  after propagation and resetting ephemerides.
     *  \param setIntegratedResult Boolean to determine whether to automatically use the integrated results to set
     *  ephemerides.
     */

    DynamicsSimulator(
            const  simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const bool clearNumericalSolutions = true,
            const bool setIntegratedResult = true ):
        bodyMap_( bodyMap ), integratorSettings_( integratorSettings ),
        propagatorSettings_( propagatorSettings ), clearNumericalSolutions_( clearNumericalSolutions ),
        setIntegratedResult_( setIntegratedResult )
    {
        frameManager_ = boost::make_shared< ephemerides::ReferenceFrameManager >( bodyMap );
        integratedStateProcessors_ = createIntegratedStateProcessors< TimeType, StateScalarType >(
                    propagatorSettings_, bodyMap_, frameManager_ );
        environmentUpdater_ = createEnvironmentUpdaterForDynamicalEquations< StateScalarType, TimeType >(
                    propagatorSettings_, bodyMap_ );
        dynamicsStateDerivative_ = boost::make_shared< DynamicsStateDerivativeModel< TimeType, StateScalarType > >(
                    createStateDerivativeModels< StateScalarType, TimeType >(
                        propagatorSettings_, bodyMap_ ), environmentUpdater_ );
        stateDerivativeFunction_ =
                boost::bind( &DynamicsStateDerivativeModel< TimeType, StateScalarType >::computeStateDerivative,
                             dynamicsStateDerivative_, _1, _2 );
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

    //! Function to get the settings for the propagator.
    /*!
     * Function to get the settings for the propagator.
     * \return The settings for the propagator.
     */
    boost::shared_ptr< PropagatorSettings< StateScalarType > > getPropagatorSettings( )
    {
        return propagatorSettings_;
    }

    //! Function to get the object resposible for updating the environment based on the current state and time.
    /*!
     * Function to get the object resposible for updating the environment based on the current state and time.
     * \return Object resposible for updating the environment based on the current state and time.
     */
    boost::shared_ptr< EnvironmentUpdater< StateScalarType, TimeType > > getEnvironmentUpdater( )
    {
        return environmentUpdater_;
    }

    //! Function to get the object that updates current environment and returns state derivative from single function call.
    /*!
     * Function to get the object that updates current environment and returns state derivative from single function call
     * \return Object that updates current environment and returns state derivative from single function call
     */
    boost::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > getDynamicsStateDerivative( )
    {
        return dynamicsStateDerivative_;
    }

    //! Function to get the map of named bodies involved in simulation.
    /*!
     *  Function to get the map of named bodies involved in simulation.
     *  \return Map of named bodies involved in simulation.
     */
    simulation_setup::NamedBodyMap getNamedBodyMap( )
    {
        return bodyMap_;
    }

protected:

    //! This function updates the environment with the numerical solution of the propagation.
    /*!
     *  This function updates the environment with the numerical solution of the propagation. For instance, it sets
     *  the propagated translational dynamics solution as the new input for the Ephemeris object of the body that was
     *  propagated. This function is pure virtual and must be implemented in the derived class.
     */
    virtual void processNumericalEquationsOfMotionSolution( ) = 0;

    //! List of object (per dynamics type) that process the integrated numerical solution by updating the environment
    std::map< IntegratedStateType, std::vector< boost::shared_ptr<
    IntegratedStateProcessor< TimeType, StateScalarType > > > > integratedStateProcessors_;

    //! Object resposible for updating the environment based on the current state and time.
    /*!
     *  Object resposible for updating the environment based on the current state and time. Calling the updateEnvironment
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


    //!  Map of bodies (with names) of all bodies in integration.
    simulation_setup::NamedBodyMap bodyMap_;

    //! Settings for numerical integrator.
    boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings_;

    //! Settings for propagator.
    boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings_;

    //! Object for retrieving ephemerides for transformation of reference frame (origins)
    boost::shared_ptr< ephemerides::ReferenceFrameManager > frameManager_;

    //! Boolean to determine whether to clear the raw numerical solution member variables after propagation and
    //! resetting ephemerides.
    bool clearNumericalSolutions_;

    //! Boolean to determine whether to automatically use the integrated results to set ephemerides.
    bool setIntegratedResult_;


};

//! Class for performing full numerical integration of a dynamical system in a single arc..
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
    using DynamicsSimulator< StateScalarType, TimeType >::environmentUpdater_;
    using DynamicsSimulator< StateScalarType, TimeType >::dynamicsStateDerivative_;
    using DynamicsSimulator< StateScalarType, TimeType >::clearNumericalSolutions_;
    using DynamicsSimulator< StateScalarType, TimeType >::stateDerivativeFunction_;
    using DynamicsSimulator< StateScalarType, TimeType >::integratorSettings_;
    using DynamicsSimulator< StateScalarType, TimeType >::propagatorSettings_;
    using DynamicsSimulator< StateScalarType, TimeType >::integratedStateProcessors_;


    //! Constructor of simulator.
    /*!
     *  Constructor of simulator, constructs integrator and object for calculating each time step of integration.
     *  \param bodyMap Map of bodies (with names) of all bodies in integration.
     *  \param integratorSettings Settings for numerical integrator.
     *  \param propagatorSettings Settings for propagator.
     *  \param areEquationsOfMotionToBeIntegrated Boolean to denote whether equations of motion should be integrated
     *  immediately at the end of the contructor or not.
     *  \param clearNumericalSolutions Boolean to determine whether to clear the raw numerical solution member variables
     *  after propagation and resetting ephemerides.
     *  \param setIntegratedResult Boolean to determine whether to automatically use the integrated results to set
     *  ephemerides.
     */
    SingleArcDynamicsSimulator(
            const  simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const bool areEquationsOfMotionToBeIntegrated = true,
            const bool clearNumericalSolutions = true,
            const bool setIntegratedResult = true ):
        DynamicsSimulator< StateScalarType, TimeType >(
            bodyMap, integratorSettings, propagatorSettings, clearNumericalSolutions, setIntegratedResult )
    {
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
        // Update body properties that are independent of integrated states.
        for( std::map< std::string, boost::shared_ptr< simulation_setup::Body > >::iterator bodyIterator = bodyMap_.begin( );
             bodyIterator != bodyMap_.end( ); bodyIterator++ )
        {
            std::cerr<<( "Error, updateConstantEphemerisIndependentMemberQuantities not set" )<<std::endl;
            //bodyIterator->second->updateConstantEphemerisIndependentMemberQuantities( );
        }

        equationsOfMotionNumericalSolution_.clear( );
        dynamicsStateDerivative_->setPropagationSettings( std::vector< IntegratedStateType >( ), 1 );


        // Integrate equations of motion numerically.
        equationsOfMotionNumericalSolution_ =
                integrateEquations< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >, TimeType >(
                    stateDerivativeFunction_, dynamicsStateDerivative_->convertFromOutputSolution(
                        initialStates, integratorSettings_->initialTime_ ), integratorSettings_ );
        equationsOfMotionNumericalSolution_ = dynamicsStateDerivative_->
                convertNumericalStateSolutionsToOutputSolutions( equationsOfMotionNumericalSolution_ );

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

    //! Function to reset the environment from an externally generated state history.
    /*!
     * Function to reset the environment from an externally generated state history, the order of the entries in the
     * state vectors are proscribed by propagatorSettings
     * \param equationsOfMotionNumericalSolution Externally generated state history.
     */
    void manuallySetAndProcessRawNumericalEquationsOfMotionSolution(
            const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >&
            equationsOfMotionNumericalSolution )
    {
        equationsOfMotionNumericalSolution_ = equationsOfMotionNumericalSolution;
        processNumericalEquationsOfMotionSolution( );
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
    }

    //! Map of state history of numerically integrated bodies.
    /*!
     *  Map of state history of numerically integrated bodies, i.e. the result of the numerical integration, transformed
     *  into the 'conventional form' (\sa SingleStateTypeDerivative::convertToOutputSolution). Key of map denotes time,
     *  values are concatenated vectors of integrated body states (order defined by propagatorSettings_).
     *  NOTEL this map is empty if clearNumericalSolutions_ is set to true.
     */
    std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > equationsOfMotionNumericalSolution_;
};

}

}


#endif // TUDAT_DYNAMICSSIMULATOR_H
