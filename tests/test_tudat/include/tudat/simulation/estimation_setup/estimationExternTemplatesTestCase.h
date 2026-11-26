/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_ESTIMATIONEXTERNTEMPLATESTESTCASE_H
#define TUDAT_ESTIMATIONEXTERNTEMPLATESTESTCASE_H

#include "tudat/basics/timeType.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/estimation_setup/createObservationManager.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"

// Explicit instantiations of the propagation and estimation stack for the state scalar and time types that the unit tests use
// next to < double, double > (which the Tudat libraries instantiate). Each group is instantiated in its own file of
// tudat_estimation_test_case_support (only built with the estimation tools), so that neither a test nor a single
// instantiation file compiles the whole stack.
#define TUDAT_TEST_DYNAMICS_SIMULATOR_INSTANTIATIONS( PREFIX, StateScalarType, TimeType )              \
    PREFIX template class tudat::propagators::DynamicsSimulator< StateScalarType, TimeType >;          \
    PREFIX template class tudat::propagators::SingleArcDynamicsSimulator< StateScalarType, TimeType >; \
    PREFIX template class tudat::propagators::MultiArcDynamicsSimulator< StateScalarType, TimeType >;  \
    PREFIX template class tudat::propagators::HybridArcDynamicsSimulator< StateScalarType, TimeType >;

#define TUDAT_TEST_VARIATIONAL_EQUATIONS_SOLVER_INSTANTIATIONS( PREFIX, StateScalarType, TimeType )                                  \
    PREFIX template class tudat::propagators::VariationalEquationsSolver< StateScalarType, TimeType >;                               \
    PREFIX template class tudat::propagators::SingleArcVariationalEquationsSolver< StateScalarType, TimeType >;                      \
    PREFIX template class tudat::propagators::MultiArcVariationalEquationsSolver< StateScalarType, TimeType >;                       \
    PREFIX template class tudat::propagators::HybridArcVariationalEquationsSolver< StateScalarType, TimeType >;                      \
    PREFIX template std::shared_ptr< tudat::propagators::VariationalEquationsSolver< StateScalarType, TimeType > >                   \
    tudat::simulation_setup::createVariationalEquationsSolver< StateScalarType, TimeType >(                                          \
            const tudat::simulation_setup::SystemOfBodies& bodies,                                                                   \
            const std::shared_ptr< tudat::propagators::PropagatorSettings< StateScalarType > > propagatorSettings,                   \
            const std::shared_ptr< tudat::estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate, \
            const bool integrateEquationsOnCreation );

#define TUDAT_TEST_ORBIT_DETERMINATION_INSTANTIATIONS( PREFIX, StateScalarType, TimeType )                 \
    PREFIX template class tudat::simulation_setup::OrbitDeterminationManager< StateScalarType, TimeType >; \
    PREFIX template struct tudat::simulation_setup::CovarianceAnalysisOutput< StateScalarType, TimeType >; \
    PREFIX template struct tudat::simulation_setup::EstimationOutput< StateScalarType, TimeType >;

#define TUDAT_TEST_OBSERVATION_MANAGER_INSTANTIATIONS( PREFIX, StateScalarType, TimeType )                                            \
    PREFIX template std::shared_ptr< tudat::observation_models::ObservationManagerBase< StateScalarType, TimeType > >                 \
    tudat::observation_models::createObservationManagerBase< StateScalarType, TimeType >(                                             \
            const tudat::observation_models::ObservableType observableType,                                                           \
            const std::vector< std::shared_ptr< tudat::observation_models::ObservationModelSettings > > observationModelSettingsList, \
            const tudat::simulation_setup::SystemOfBodies& bodies,                                                                    \
            const std::shared_ptr< tudat::estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,  \
            const std::shared_ptr< tudat::propagators::CombinedStateTransitionAndSensitivityMatrixInterface >                         \
                    stateTransitionMatrixInterface,                                                                                   \
            const std::shared_ptr< tudat::propagators::DependentVariablesInterface< TimeType > > dependentVariablesInterface );       \
    PREFIX template std::map< tudat::observation_models::ObservableType,                                                              \
                              std::shared_ptr< tudat::observation_models::ObservationManagerBase< StateScalarType, TimeType > > >     \
    tudat::observation_models::createObservationManagersBase< StateScalarType, TimeType >(                                            \
            const std::vector< std::shared_ptr< tudat::observation_models::ObservationModelSettings > >& observationSettingsList,     \
            const tudat::simulation_setup::SystemOfBodies& bodies,                                                                    \
            const std::shared_ptr< tudat::estimatable_parameters::EstimatableParameterSet< StateScalarType > > fullParameters,        \
            const std::shared_ptr< tudat::propagators::CombinedStateTransitionAndSensitivityMatrixInterface >                         \
                    stateTransitionMatrixInterface,                                                                                   \
            const std::shared_ptr< tudat::propagators::DependentVariablesInterface< TimeType > > dependentVariablesInterface );

#define TUDAT_TEST_ESTIMATION_INSTANTIATIONS( PREFIX, StateScalarType, TimeType )               \
    TUDAT_TEST_DYNAMICS_SIMULATOR_INSTANTIATIONS( PREFIX, StateScalarType, TimeType )           \
    TUDAT_TEST_VARIATIONAL_EQUATIONS_SOLVER_INSTANTIATIONS( PREFIX, StateScalarType, TimeType ) \
    TUDAT_TEST_ORBIT_DETERMINATION_INSTANTIATIONS( PREFIX, StateScalarType, TimeType )          \
    TUDAT_TEST_OBSERVATION_MANAGER_INSTANTIATIONS( PREFIX, StateScalarType, TimeType )

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS && TUDAT_BUILD_WITH_ESTIMATION_TOOLS
TUDAT_TEST_ESTIMATION_INSTANTIATIONS( extern, long double, tudat::Time )
TUDAT_TEST_ESTIMATION_INSTANTIATIONS( extern, long double, double )
TUDAT_TEST_ESTIMATION_INSTANTIATIONS( extern, double, tudat::Time )
#endif

#endif  // TUDAT_ESTIMATIONEXTERNTEMPLATESTESTCASE_H
