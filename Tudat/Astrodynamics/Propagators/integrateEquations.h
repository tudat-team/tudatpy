#ifndef INTEGRATEEQUATIONS_H
#define INTEGRATEEQUATIONS_H

#include <Eigen/Core>

#include <map>

#include "Tudat/Mathematics/NumericalIntegrators/numericalIntegrator.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"
#include "Tudat/Astrodynamics/Propagators/environmentUpdater.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"

namespace tudat
{

namespace propagators
{

//! Function to numerically integrate a given equation
/*!
 *  Function to numerically integrate a given equation of a single independent variable and the current state
 *  \param stateDerivativeFunction Function returning the state derivative from current time and state.
 *  \param initialState Initial state
 *  \param integratorSettings Settings for numerical integrator.
 *  \return History of numerical states (first of pair) and derivatives of states (second of pair) given as maps with time as key.
 */
template< typename StateType = Eigen::MatrixXd, typename TimeType = double >
std::map< TimeType, StateType > integrateEquations(
        boost::function< StateType( const TimeType, const StateType&) > stateDerivativeFunction,
        const StateType initialState,
        boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings )
{
    using namespace tudat::numerical_integrators;


    // Create numerical integrator.
    boost::shared_ptr< NumericalIntegrator< TimeType, StateType, StateType > > integrator =
            createIntegrator< TimeType, StateType >( stateDerivativeFunction, initialState, integratorSettings );

    // Get Initial state and time.
    TimeType currentTime = integratorSettings->initialTime_;
    StateType newState = initialState;

    // Initialization of numerical solutions for variational equations.
    std::map< TimeType, StateType > solutionHistory;
    solutionHistory[ currentTime ] = newState;

    // Check if numerical integration is forward or backwrd.
    TimeType timeStepSign = 1.0L;
    if( integratorSettings->initialTimeStep_ < 0.0 )
    {
        timeStepSign = -1.0L;
    }

    // Set initial time step and total integration time.
    TimeType timeStep = integratorSettings->initialTimeStep_;
    TimeType endTime = integratorSettings->endTime_;
    TimeType previousTime = currentTime;

    // Perform first integration step.
    newState = integrator->performIntegrationStep( timeStep );

    currentTime = integrator->getCurrentIndependentVariable( );

    timeStep = timeStepSign * integrator->getNextStepSize( );
    solutionHistory[ currentTime ] = newState;

    int printIndex = 0;
    int printFrequency = integratorSettings->printFrequency_;
    // Perform numerical integration steps until end time reached.
    while( timeStepSign * static_cast< TimeType >( currentTime ) < timeStepSign * static_cast< TimeType >( endTime ) )
    {
        previousTime = currentTime;

        // Perform integration step.
        newState = integrator->performIntegrationStep( timeStep );
        currentTime = integrator->getCurrentIndependentVariable( );
        timeStep = timeStepSign * integrator->getNextStepSize( );

        printIndex++;
        printIndex = printIndex % printFrequency;

        // Save integration result in map.
        if( printIndex == 0 )
        {
            solutionHistory[ currentTime ] = newState;
        }


        // Print solutions
//        if( ( static_cast<int>( std::fabs( currentTime - integratorSettings->initialTime_ ) ) % static_cast< int >( 1.0E6 ) ) <
//                ( static_cast< int >( std::fabs( previousTime - integratorSettings->initialTime_ ) ) % static_cast<int>( 1.0E6 ) )  )
//        {
//            std::cout<<std::setprecision( 10 )<<timeStep<<" "<<currentTime<<std::endl;
//        }
    }

    return solutionHistory;
}

template< typename StateScalarType = double, typename TimeType = double >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > updateEnvironmentAndCalculateStateDerivative(
        const boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > stateDerivativeModel,
        const boost::shared_ptr< EnvironmentUpdater< StateScalarType, TimeType > > environmentUpdater,
        const TimeType currentTime, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > currentState,
        const boost::function< double( const double ) > timeConverterToUpdateTime = &basic_astrodynamics::doDummyTimeConversion< TimeType > )
{
    std::map< IntegratedStateType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > updateStateInput;
    updateStateInput[ stateDerivativeModel->getIntegratedStateType( ) ] = stateDerivativeModel->convertCurrentStateToGlobalRepresentation(
                currentState, currentTime );
    environmentUpdater->updateEnvironment( timeConverterToUpdateTime( currentTime ), updateStateInput );
    return stateDerivativeModel->calculateSystemStateDerivative( currentTime, currentState );
}

template< typename StateScalarType = double, typename TimeType = double >
std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > integrateEquations(
        const boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > stateDerivativeModel,
        const boost::shared_ptr< EnvironmentUpdater< StateScalarType, TimeType > > environmentUpdater,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialState,
        boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const boost::function< double( const double ) > timeConverterToUpdateTime = &basic_astrodynamics::doDummyTimeConversion< TimeType >  )
{
    boost::function< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >(
        const TimeType, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >&) > stateDerivativeFunction =
            boost::bind( &updateEnvironmentAndCalculateStateDerivative< StateScalarType, TimeType >, stateDerivativeModel, environmentUpdater, _1, _2,
                         timeConverterToUpdateTime );
    return integrateEquations< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >, TimeType >(
                stateDerivativeFunction, initialState, integratorSettings );
}


template< typename StateScalarType = double, typename TimeType = double >
std::pair< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >,
std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > >
integrateEquationsForwardAndBackward( boost::function< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >(
                                          const TimeType, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >&) > stateDerivativeFunction,
                                      const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > initialState,
                                      boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings )
{
    using namespace tudat::numerical_integrators;

    std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > forwardIntegrated =
            integrateEquations< StateScalarType, TimeType >(
                stateDerivativeFunction, initialState, integratorSettings );

    TimeType originalStartTime = integratorSettings->initialTime_;
    TimeType originalEndTime = integratorSettings->endTime_;
    double originalTimeStep = integratorSettings->initialTimeStep_;

    typename std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >::iterator outputMapIterator =
            forwardIntegrated.end( );
    --outputMapIterator;

    integratorSettings->initialTime_ = outputMapIterator->first;
    integratorSettings->endTime_ = originalStartTime;
    integratorSettings->initialTimeStep_ = -originalTimeStep * 1.2;

    std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > backwardIntegrated =
            integrateEquations< StateScalarType, TimeType >(
                stateDerivativeFunction, outputMapIterator->second, integratorSettings );

    integratorSettings->initialTime_ = originalStartTime;
    integratorSettings->endTime_ = originalEndTime;
    integratorSettings->initialTimeStep_ = originalTimeStep;

    return std::make_pair( forwardIntegrated, backwardIntegrated );
}

template< typename StateScalarType = double, typename TimeType = double >
std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >
getForwardAndBackwardIntegrationDifference(
        boost::function< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >(
            const TimeType, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& ) > stateDerivativeFunction,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > initialState,
        boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const double interpolationTimeStep )
{
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > StateType;

    typedef std::map< TimeType, StateType > StateMap;
    std::pair< StateMap, StateMap > forwardBackwardIntegrationResult =
            integrateEquationsForwardAndBackward( stateDerivativeFunction, initialState, integratorSettings );

    boost::shared_ptr< interpolators::LagrangeInterpolator< TimeType, StateType > > forwardInterpolator =
            boost::make_shared< interpolators::LagrangeInterpolator< TimeType, StateType > >( forwardBackwardIntegrationResult.first, 8 );
    boost::shared_ptr< interpolators::LagrangeInterpolator< TimeType, StateType > > backwardInterpolator =
            boost::make_shared< interpolators::LagrangeInterpolator< TimeType, StateType > >( forwardBackwardIntegrationResult.second, 8 );

    StateMap stateDifference;

    TimeType currentTime = integratorSettings->initialTime_;
    while( currentTime < integratorSettings->endTime_ )
    {
        forwardInterpolator->interpolate( currentTime );
        backwardInterpolator->interpolate( currentTime );
        stateDifference[ currentTime ] = forwardInterpolator->interpolate( currentTime ) - backwardInterpolator->interpolate( currentTime );
        currentTime += interpolationTimeStep;
    }

    return stateDifference;
}

}

}
#endif // INTEGRATEEQUATIONS_H
