#include "Tudat/Astrodynamics/Propagators/integrateEquations.h"

namespace tudat
{

namespace propagators
{

template std::shared_ptr< PropagationTerminationDetails > integrateEquationsFromIntegrator<
Eigen::MatrixXd, double, double >(
        const std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::MatrixXd, Eigen::MatrixXd, double > > integrator,
        const double initialTimeStep,
        const std::shared_ptr< PropagationTerminationCondition > propagationTerminationCondition,
        std::map< double, Eigen::MatrixXd >& solutionHistory,
        std::map< double, Eigen::VectorXd >& dependentVariableHistory,
        std::map< double, double >& cumulativeComputationTimeHistory,
        const std::function< Eigen::VectorXd( ) > dependentVariableFunction,
        const std::function< void( Eigen::MatrixXd& ) > statePostProcessingFunction,
        const int saveFrequency,
        const double printInterval,
        const std::chrono::steady_clock::time_point initialClockTime );

template std::shared_ptr< PropagationTerminationDetails > integrateEquationsFromIntegrator<
Eigen::VectorXd, double, double >(
        const std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::VectorXd, Eigen::VectorXd, double > > integrator,
        const double initialTimeStep,
        const std::shared_ptr< PropagationTerminationCondition > propagationTerminationCondition,
        std::map< double, Eigen::VectorXd >& solutionHistory,
        std::map< double, Eigen::VectorXd >& dependentVariableHistory,
        std::map< double, double >& cumulativeComputationTimeHistory,
        const std::function< Eigen::VectorXd( ) > dependentVariableFunction,
        const std::function< void( Eigen::VectorXd& ) > statePostProcessingFunction,
        const int saveFrequency,
        const double printInterval,
        const std::chrono::steady_clock::time_point initialClockTime );

} // namespace propagators

} // namespace tudat
