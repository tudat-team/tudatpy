#include "Tudat/SimulationSetup/EstimationSetup/createNumericalSimulator.h"

namespace tudat
{

namespace simulation_setup
{

template std::shared_ptr< propagators::SingleArcVariationalEquationsSolver< double, double > >
createSingleArcVariationalEquationsSolver< double, double >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

template std::shared_ptr< propagators::SingleArcDynamicsSimulator< double, double > >
createSingleArcDynamicsSimulator< double, double >(
        const  simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const bool areEquationsOfMotionToBeIntegrated,
        const bool clearNumericalSolutions,
        const bool setIntegratedResult );

template std::shared_ptr< propagators::MultiArcVariationalEquationsSolver< double, double > >
createMultiArcVariationalEquationsSolver< double, double >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

template std::shared_ptr< propagators::HybridArcVariationalEquationsSolver< double, double > >
createHybridArcVariationalEquationsSolver< double, double >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );


template std::shared_ptr< propagators::VariationalEquationsSolver< double, double > >
createVariationalEquationsSolver< double, double >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )


template std::shared_ptr< propagators::SingleArcDynamicsSimulator< double, Time > >
createSingleArcDynamicsSimulator< double, Time >(
        const  simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< Time > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const bool areEquationsOfMotionToBeIntegrated,
        const bool clearNumericalSolutions,
        const bool setIntegratedResult );

template std::shared_ptr< propagators::SingleArcDynamicsSimulator< long double, double > >
createSingleArcDynamicsSimulator< long double, double >(
        const  simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< long double > > propagatorSettings,
        const bool areEquationsOfMotionToBeIntegrated,
        const bool clearNumericalSolutions,
        const bool setIntegratedResult );

template std::shared_ptr< propagators::SingleArcDynamicsSimulator< long double, Time > >
createSingleArcDynamicsSimulator< long double, Time >(
        const  simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< Time > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< long double > > propagatorSettings,
        const bool areEquationsOfMotionToBeIntegrated,
        const bool clearNumericalSolutions,
        const bool setIntegratedResult );

template std::shared_ptr< propagators::SingleArcVariationalEquationsSolver< double, Time > >
createSingleArcVariationalEquationsSolver< double, Time >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< Time > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

template std::shared_ptr< propagators::SingleArcVariationalEquationsSolver< long double, Time > >
createSingleArcVariationalEquationsSolver< long double, Time >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< Time > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< long double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  long double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

template std::shared_ptr< propagators::SingleArcVariationalEquationsSolver< long double, double > >
createSingleArcVariationalEquationsSolver< long double, double >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< long double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  long double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );



template std::shared_ptr< propagators::HybridArcVariationalEquationsSolver< double, Time > >
createHybridArcVariationalEquationsSolver< double, Time >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< Time > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

template std::shared_ptr< propagators::HybridArcVariationalEquationsSolver< long double, Time > >
createHybridArcVariationalEquationsSolver< long double, Time >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< Time > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< long double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  long double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

template std::shared_ptr< propagators::HybridArcVariationalEquationsSolver< long double, double > >
createHybridArcVariationalEquationsSolver< long double, double >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< long double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  long double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );



template std::shared_ptr< propagators::MultiArcVariationalEquationsSolver< double, Time > >
createMultiArcVariationalEquationsSolver< double, Time >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< Time > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

template std::shared_ptr< propagators::MultiArcVariationalEquationsSolver< long double, Time > >
createMultiArcVariationalEquationsSolver< long double, Time >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< Time > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< long double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  long double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

template std::shared_ptr< propagators::MultiArcVariationalEquationsSolver< long double, double > >
createMultiArcVariationalEquationsSolver< long double, double >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< long double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  long double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );




template std::shared_ptr< propagators::VariationalEquationsSolver< double, Time > >
createVariationalEquationsSolver< double, Time >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< Time > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

template std::shared_ptr< propagators::VariationalEquationsSolver< long double, Time > >
createVariationalEquationsSolver< long double, Time >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< Time > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< long double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  long double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

template std::shared_ptr< propagators::VariationalEquationsSolver< long double, double > >
createVariationalEquationsSolver< long double, double >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< long double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  long double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

#endif

} // namespace simulation_setup

} // namespace tudat

