#include "Tudat/Astrodynamics/Propagators/integrateEquations.h"

namespace tudat
{

namespace propagators
{

template< typename StateType = Eigen::MatrixXd, typename TimeType = double >
void integrateEquations(
        boost::function< StateType( const TimeType, const StateType& ) > stateDerivativeFunction,
        std::map< TimeType, StateType >& solutionHistory,
        const StateType initialState,
        const boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const boost::function< bool( const double ) > stopPropagationFunction,
        std::map< TimeType, Eigen::VectorXd >& dependentVariableHistory,
        const boost::function< Eigen::VectorXd( ) > dependentVariableFunction =
        boost::function< Eigen::VectorXd( ) >( ),
        const TimeType printInterval = TUDAT_NAN )
{
    // Create numerical integrator.
    boost::shared_ptr< numerical_integrators::NumericalIntegrator< TimeType, StateType, StateType > > integrator =
            numerical_integrators::createIntegrator< TimeType, StateType >(
                stateDerivativeFunction, initialState, integratorSettings );

    integrateEquations< StateType, TimeType >(
                integrator, integratorSettings->initialTimeStep_, stopPropagationFunction, solutionHistory,
                dependentVariableHistory,
                dependentVariableFunction,
                integratorSettings->saveFrequency_, printInterval );

}

}

}
