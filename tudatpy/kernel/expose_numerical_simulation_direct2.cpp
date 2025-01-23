/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "docstrings.h"
#include "scalarTypes.h"

#include "expose_numerical_simulation.h"

#include "expose_numerical_simulation/expose_environment_setup.h"
#include "expose_numerical_simulation/expose_estimation_setup.h"
#include "expose_numerical_simulation/expose_propagation_setup.h"

#include "expose_numerical_simulation/expose_environment.h"
#include "expose_numerical_simulation/expose_estimation.h"
#include "expose_numerical_simulation/expose_propagation.h"

#include "tudat/basics/timeType.h"
#include "tudat/astro/basic_astro/dateTime.h"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace tp = tudat::propagators;
namespace tss = tudat::simulation_setup;
namespace tni = tudat::numerical_integrators;
namespace tep = tudat::estimatable_parameters;
namespace tom = tudat::observation_models;
namespace tba = tudat::basic_astrodynamics;

namespace tudatpy {


namespace numerical_simulation {

void expose_numerical_simulation_direct2(py::module &m) {



    //TODO: Remove variationalOnlyIntegratorSettings
    py::class_<
            tp::SingleArcVariationalEquationsSolver<STATE_SCALAR_TYPE, TIME_TYPE>,
            std::shared_ptr<tp::SingleArcVariationalEquationsSolver<STATE_SCALAR_TYPE, TIME_TYPE>>>(m, "SingleArcVariationalSimulator",
                                                                                      get_docstring("SingleArcVariationalSimulator").c_str() )
            .def(py::init<
                 const tudat::simulation_setup::SystemOfBodies&,
                 const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings<TIME_TYPE>>,
                 const std::shared_ptr< tp::PropagatorSettings<STATE_SCALAR_TYPE>>,
                 const std::shared_ptr< tep::EstimatableParameterSet< STATE_SCALAR_TYPE > >,
                 const bool,
                 const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > >,
                 const bool,
                 const bool,
                 const bool,
                 const bool>(),
                 py::arg("bodies"),
                 py::arg("integrator_settings"),
                 py::arg("propagator_settings"),
                 py::arg("estimated_parameters"),
                 py::arg("integrate_equations_concurrently") = true,
                 py::arg("variational_only_integrator_settings") = std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< TIME_TYPE > >( ),
                 py::arg("clear_numerical_solutions") = false,
                 py::arg("integrate_on_creation") = true,
                 py::arg("set_integrated_result") = false,
                 py::arg("print_dependent_variable_data") = true,
                 get_docstring("SingleArcVariationalSimulator.ctor").c_str() )
            .def("integrate_equations_of_motion_only",
                 &tp::SingleArcVariationalEquationsSolver<STATE_SCALAR_TYPE, TIME_TYPE>::integrateDynamicalEquationsOfMotionOnly,
                 py::arg("initial_states"),
                 get_docstring("SingleArcVariationalSimulator.integrate_equations_of_motion_only").c_str() )
            .def("integrate_full_equations",
                 &tp::SingleArcVariationalEquationsSolver<STATE_SCALAR_TYPE, TIME_TYPE>::integrateVariationalAndDynamicalEquations,
                 py::arg("initial_states"),
                 py::arg("integrate_equations_concurrently") = true,
                 get_docstring("SingleArcVariationalSimulator.integrate_full_equations").c_str() )
            .def_property("parameter_vector",
                          &tp::SingleArcVariationalEquationsSolver<STATE_SCALAR_TYPE, TIME_TYPE>::getParametersToEstimate,
                          &tp::SingleArcVariationalEquationsSolver<STATE_SCALAR_TYPE, TIME_TYPE>::resetParameterEstimate,
                          get_docstring("SingleArcVariationalSimulator.parameter_vector").c_str() )
            .def_property_readonly("variational_equations_history",
                                   &tp::SingleArcVariationalEquationsSolver<STATE_SCALAR_TYPE, TIME_TYPE>::getNumericalVariationalEquationsSolution,
                                   get_docstring("SingleArcVariationalSimulator.variational_equations_history").c_str() )
            .def_property_readonly("state_transition_matrix_history",
                                   &tp::SingleArcVariationalEquationsSolver<STATE_SCALAR_TYPE, TIME_TYPE>::getStateTransitionMatrixSolution,
                                   get_docstring("SingleArcVariationalSimulator.state_transition_matrix_history").c_str() )
            .def_property_readonly("sensitivity_matrix_history",
                                   &tp::SingleArcVariationalEquationsSolver<STATE_SCALAR_TYPE, TIME_TYPE>::getSensitivityMatrixSolution,
                                   get_docstring("SingleArcVariationalSimulator.sensitivity_matrix_history").c_str() )
            .def_property_readonly("state_history",
                                   &tp::SingleArcVariationalEquationsSolver<STATE_SCALAR_TYPE, TIME_TYPE>::getEquationsOfMotionSolution,
                                   get_docstring("SingleArcVariationalSimulator.state_history").c_str() )
            .def_property_readonly("dynamics_simulator",
                                   &tp::SingleArcVariationalEquationsSolver<STATE_SCALAR_TYPE, TIME_TYPE>::getDynamicsSimulator,
                                   get_docstring("SingleArcVariationalSimulator.dynamics_simulator").c_str() );


};

}// namespace numerical_simulation
}// namespace tudatpy
