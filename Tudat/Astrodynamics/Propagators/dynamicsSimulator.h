#ifndef DYNAMICSSIMULATOR2_H
#define DYNAMICSSIMULATOR2_H

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
#include "Tudat/Astrodynamics/Propagators/hybridStateDerivativeModel.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"

namespace tudat
{

namespace propagators
{

template< typename TimeType = double, typename StateScalarType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > setInitialStatesOfBodies(
        const std::vector< std::string >& bodiesToIntegrate,
        const std::vector< std::string >& centralBodies,
        const  simulation_setup::NamedBodyMap& bodyMap,
        const TimeType initialTime )
{
    // Set initial states of bodies to integrate.
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > systemInitialState =
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( bodiesToIntegrate.size( ) * 6, 1 );
    boost::shared_ptr< ephemerides::Ephemeris > ephemerisOfCurrentBody;

    ephemerides::ReferenceFrameManager frameManager = ephemerides::ReferenceFrameManager(
                bodyMap );

    for( unsigned int i = 0; i < bodiesToIntegrate.size( ) ; i++ )
    {
        ephemerisOfCurrentBody = bodyMap.at( bodiesToIntegrate.at( i ) )->getEphemeris( );
        systemInitialState.segment( i * 6 , 6 ) = ephemerisOfCurrentBody->getTemplatedStateFromEphemeris<
                StateScalarType, TimeType >( initialTime );

        if( centralBodies.at( i ) != ephemerisOfCurrentBody->getReferenceFrameOrigin( ) )
        {
            boost::shared_ptr< ephemerides::Ephemeris > correctionEphemeris =
                    frameManager.getEphemeris( ephemerisOfCurrentBody->getReferenceFrameOrigin( ), centralBodies.at( i ) );
            systemInitialState.segment( i * 6 , 6 ) -= correctionEphemeris->getTemplatedStateFromEphemeris<
                    StateScalarType, TimeType >( initialTime );
        }
    }
    return systemInitialState;
}

template< typename TimeType = double, typename StateScalarType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > setInitialStateOfBody(
        const std::string& bodyToIntegrate,
        const std::string& centralBody,
        const  simulation_setup::NamedBodyMap& bodyMap,
        const TimeType initialTime )
{
    return setInitialStatesOfBodies< TimeType, StateScalarType >(
                boost::assign::list_of( bodyToIntegrate ), boost::assign::list_of( centralBody ), bodyMap, initialTime );
}

template< typename StateScalarType = double, typename TimeType = double >
class DynamicsSimulatorBase
{
public:
    DynamicsSimulatorBase( const  simulation_setup::NamedBodyMap& bodyMap ): bodyMap_( bodyMap )
    {
        frameManager_ = boost::make_shared< ephemerides::ReferenceFrameManager >( bodyMap );
    }

    virtual ~DynamicsSimulatorBase( ){ }

    virtual void integrateEquationsOfMotion(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& initialGlobalStates ) = 0;

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

    //! Object for retrieving ephemerides for transformation of reference frame (origins)
    /*!
     *  Object for retrieving ephemerides for transformation of reference frame (origins)
     */
    boost::shared_ptr< ephemerides::ReferenceFrameManager > frameManager_;

    //!  Map of bodies (with names) of all bodies in integration.
    /*!
     *   Map of bodies (with names) of all bodies in integration.
     */
     simulation_setup::NamedBodyMap bodyMap_;

};

//! Class to manage integration of equations of motion of bodies undergoing accelerations.
/*!
 *  Class to manage integration of equations of motion of bodies undergoing accelerations. Equations of motion are set once, can be re-integrated
 *  for different initial conditions.
 */
template< typename StateScalarType = double, typename TimeType = double >
class DynamicsSimulator: public DynamicsSimulatorBase< StateScalarType, TimeType >
{
public:

    using DynamicsSimulatorBase< StateScalarType, TimeType >::frameManager_;
    using DynamicsSimulatorBase< StateScalarType, TimeType >::bodyMap_;

    //! Constructor of simulator.
    /*!
     *  Constructor of simulator, constructs integrator and object for calculating each time step of integration.
     *  \param bodyMap Map of bodies (with names) of all bodies in integration.
     *  \param integratorSettings Settings for numerical integrator.
     *  \param propagatorSettings Settings for propagator.
     *  \param areEquationsOfMotionToBeIntegrated Boolean to denote whether equations of motion should be integrated at the end of the contructor or not.
     */
    DynamicsSimulator(
            const  simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const bool clearNumericalSolutions = true,
            const bool setIntegratedResult = true ):
        DynamicsSimulatorBase< StateScalarType, TimeType >( bodyMap ),  integratorSettings_( integratorSettings ),
        propagatorSettings_( propagatorSettings ), clearNumericalSolutions_( clearNumericalSolutions ), setIntegratedResult_( setIntegratedResult )
    {
        using namespace state_derivative_models;

        integratedStateProcessors_ = createIntegratedStateProcessors< TimeType, StateScalarType >(  propagatorSettings_, bodyMap_, frameManager_ );
        environmentUpdater_ = createEnvironmentUpdaterForDynamicalEquations< StateScalarType, TimeType >( propagatorSettings_, bodyMap_ );
        dynamicsStateDerivative_ = boost::make_shared< HybridStateDerivativeModel< TimeType, StateScalarType > >(
                    createStateDerivativeModels< StateScalarType, TimeType >( propagatorSettings_, bodyMap_, integratorSettings_->initialTime_ ),
                    environmentUpdater_ );

        // Create function evaluation for differential equations.
        stateDerivativeFunction_ =
                boost::bind( &HybridStateDerivativeModel< TimeType, StateScalarType >::computeStateDerivative,
                             dynamicsStateDerivative_, _1, _2 );
        doubleStateDerivativeFunction_ =
                boost::bind( &HybridStateDerivativeModel< TimeType, StateScalarType >::computeStateDoubleDerivative,
                             dynamicsStateDerivative_, _1, _2 );
    }

    //! Virtual destructor
    /*!
     *  Virtual destructor
     */
    virtual ~DynamicsSimulator( ) { }

    //! This function numerically (re-)integrates the equations of motion.
    /*!
     *  This function numerically (re-)integrates the equations of motion, using the settings set through the constructor
     *  and a new initial state vector provided here. The raw results are set in the equationsOfMotionNumericalSolution_
     *  \param initialState Initial state vector that is to be used for numerical integration.
     */
    virtual void integrateEquationsOfMotion(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& initialGlobalStates ) = 0;


    boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > getIntegratorSettings( )
    {
        return integratorSettings_;
    }

    boost::function< Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >
    ( const double, const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& ) > getDoubleStateDerivativeFunction( )
    {
        return doubleStateDerivativeFunction_;
    }

    boost::function< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >
    ( const TimeType, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >&) > getStateDerivativeFunction( )
    {
        return stateDerivativeFunction_;
    }

    boost::shared_ptr< PropagatorSettings< StateScalarType > > getPropagatorSettings( )
    {
        return propagatorSettings_;
    }

    boost::shared_ptr< EnvironmentUpdater< StateScalarType, TimeType > > getEnvironmentUpdater( )
    {
        return environmentUpdater_;
    }

    boost::shared_ptr< HybridStateDerivativeModel< TimeType, StateScalarType > > getDynamicsStateDerivative( )
    {
        return dynamicsStateDerivative_;
    }

protected:

    virtual void processRawNumericalEquationsOfMotionSolution( ) = 0;

    std::map< IntegratedStateType, std::vector< boost::shared_ptr< IntegratedStateProcessor< TimeType, StateScalarType > > > > integratedStateProcessors_;

    boost::shared_ptr< EnvironmentUpdater< StateScalarType, TimeType > > environmentUpdater_;

    //! Interface object that updates current environment and return state derivative from single function call
    /*!
     *  Interface object that updates current environment and return state derivative from single function call
     */
    boost::shared_ptr< HybridStateDerivativeModel< TimeType, StateScalarType > > dynamicsStateDerivative_;

    bool clearNumericalSolutions_;

    bool setIntegratedResult_;


    boost::function< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >
    ( const TimeType, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& ) > stateDerivativeFunction_;

    boost::function< Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >
    ( const double, const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& ) > doubleStateDerivativeFunction_;




    //! Settings for numerical integrator.
    /*!
     *  Settings for numerical integrator.
     */
    boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings_;

    //! Settings for propagator.
    /*!
     *  Settings for propagator.
     */
    boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings_;


};

template< typename StateScalarType = double, typename TimeType = double >
class SingleArcDynamicsSimulator: public DynamicsSimulator< StateScalarType, TimeType >
{

public:

    using DynamicsSimulatorBase< StateScalarType, TimeType >::bodyMap_;

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
     *  \param areEquationsOfMotionToBeIntegrated Boolean to denote whether equations of motion should be integrated at the end of the contructor or not.
     */
    SingleArcDynamicsSimulator(
            const  simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const bool areEquationsOfMotionToBeIntegrated = true,
            const bool clearNumericalSolutions = true,
            const bool setIntegratedResult = true ):DynamicsSimulator< StateScalarType, TimeType >(
            bodyMap, integratorSettings, propagatorSettings, clearNumericalSolutions, setIntegratedResult )
    {
        if( propagatorSettings->isPropagatorMultiArc_ )
        {
            std::cerr<<"Error when making single arc dynamics simulator, cannot process propagator settings with multiple integration arcs"<<std::endl;
        }

        if( areEquationsOfMotionToBeIntegrated )
        {
            // Integrate equations of motion.
            integrateEquationsOfMotion( propagatorSettings_->getInitialStates( ) );
        }
    }

    //! Virtual destructor
    /*!
     *  Virtual destructor
     */
    ~SingleArcDynamicsSimulator( )
    {
        std::cerr<<"Deleting SingleArcDynamicsSimulator"<<std::endl;
    }

    //! This function numerically (re-)integrates the equations of motion.
    /*!
     *  This function numerically (re-)integrates the equations of motion, using the settings set through the constructor
     *  and a new initial state vector provided here. The raw results are set in the equationsOfMotionNumericalSolution_
     *  \param initialState Initial state vector that is to be used for numerical integration.
     */
    void integrateEquationsOfMotion(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& initialGlobalStates )
    {
        using namespace state_derivative_models;

        // Update body properties that are independent of integrated states.
        for( std::map< std::string, boost::shared_ptr< simulation_setup::Body > >::iterator bodyIterator = bodyMap_.begin( );
             bodyIterator != bodyMap_.end( ); bodyIterator++ )
        {
            //throw std::runtime_error( "Error, updateConstantEphemerisIndependentMemberQuantities not set" );
            //bodyIterator->second->updateConstantEphemerisIndependentMemberQuantities( );
        }

        equationsOfMotionNumericalSolution_.clear( );
        dynamicsStateDerivative_->setPropagationSettings( std::vector< IntegratedStateType >( ), 1 );


        // Integrate equations of motion numerically.
        equationsOfMotionNumericalSolution_ =
                integrateEquations< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >, TimeType >(
                    stateDerivativeFunction_, dynamicsStateDerivative_->convertFromOutputSolution(
                        initialGlobalStates, integratorSettings_->initialTime_ ), integratorSettings_ );
        equationsOfMotionNumericalSolution_ = convertNumericalStateSolutionsToOutputSolutions(
                    equationsOfMotionNumericalSolution_, dynamicsStateDerivative_ );
        if( this->setIntegratedResult_ )
        {
            processRawNumericalEquationsOfMotionSolution( );
        }
    }

    std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > getEquationsOfMotionNumericalSolution( )
    {
        return equationsOfMotionNumericalSolution_;
    }

    void manuallySetAndProcessRawNumericalEquationsOfMotionSolution(
            const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& equationsOfMotionNumericalSolution )
    {
        equationsOfMotionNumericalSolution_ = equationsOfMotionNumericalSolution;
        processRawNumericalEquationsOfMotionSolution( );
    }

protected:
    void processRawNumericalEquationsOfMotionSolution( )
    {
        // Create and set interpolators for ephemerides
        resetIntegratedStates( equationsOfMotionNumericalSolution_, integratedStateProcessors_ );

        if( clearNumericalSolutions_ )
        {
            equationsOfMotionNumericalSolution_.clear( );
        }
    }

    //! Map of state history of numerically integrated bodies.
    /*!
     *  Map of state history of numerically integrated bodies. Key of map denotes time, values are concatenated vectors of
     *  body states in order of bodiesToIntegrate
     */
    std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > equationsOfMotionNumericalSolution_;
};

template< typename StateScalarType = double, typename TimeType = double >
class HybridDynamicsSimulator: public DynamicsSimulatorBase< StateScalarType, TimeType >
{
public:

    HybridDynamicsSimulator(
            const  simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const std::vector< boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > > propagatorSettings,
            const bool areEquationsOfMotionToBeIntegrated = true,
            const bool clearNumericalSolutions = true ):
        DynamicsSimulatorBase< StateScalarType, TimeType >( bodyMap )
    {

        for( unsigned int i = 0; i < propagatorSettings.size( ); i++ )
        {
            dynamicsSimulators_.push_back(
                        boost::make_shared< DynamicsSimulator< StateScalarType, TimeType > >(
                            integratorSettings, propagatorSettings.at( i ), areEquationsOfMotionToBeIntegrated,
                            clearNumericalSolutions ) );
        }
    }

    ~HybridDynamicsSimulator( ){ }

    void integrateEquationsOfMotion(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& initialGlobalStates )
    {
        for( unsigned int i = 0; i < dynamicsSimulators_.size( ); i++ )
        {
            dynamicsSimulators_.at( i )->integrateEquationsOfMotion(
                        initialGlobalStates.block( initialStateSegments_.at( i ).first, 0,
                                                   initialStateSegments_.at( i ).second, 1 ) );
        }
    }

private:

    std::vector< DynamicsSimulator< StateScalarType, TimeType > > dynamicsSimulators_;

    std::vector< std::pair< int, int > > initialStateSegments_;
};

template< typename StateScalarType = double, typename TimeType = double >
std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >
getForwardAndBackwardIntegrationDifference(
        boost::shared_ptr< DynamicsSimulator< StateScalarType, TimeType > > dynamicsSimulator,
        const double interpolationTimeStep )
{
    return getForwardAndBackwardIntegrationDifference< StateScalarType, TimeType >(
                dynamicsSimulator->getStateDerivativeFunction( ), dynamicsSimulator->getPropagatorSettings( )->getInitialStates( ),
                dynamicsSimulator->getIntegratorSettings( ), interpolationTimeStep );
}

}

}


#endif // DYNAMICSMANAGER_H
