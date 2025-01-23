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

void expose_numerical_simulation_direct3(py::module &m) {





    py::class_<
            tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE>,
            std::shared_ptr<tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE>>>(m, "Estimator",
                                                                             get_docstring("Estimator").c_str() )
            .def(py::init<
                 const tss::SystemOfBodies&,
                 const std::shared_ptr< tep::EstimatableParameterSet< STATE_SCALAR_TYPE > >,
                 const std::vector< std::shared_ptr< tom::ObservationModelSettings > >&,
                 const std::shared_ptr< tp::PropagatorSettings< STATE_SCALAR_TYPE > >,
                 const bool >( ),
                 py::arg("bodies"),
                 py::arg("estimated_parameters"),
                 py::arg("observation_settings"),
                 py::arg("propagator_settings"),
                 py::arg("integrate_on_creation") = true,
                 get_docstring("Estimator.ctor").c_str() )
            .def_property_readonly("observation_simulators",
                                   &tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE>::getObservationSimulators,
                                   get_docstring("Estimator.observation_simulators").c_str() )
            .def_property_readonly("observation_managers",
                                   &tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE>::getObservationManagers,
                                   get_docstring("Estimator.observation_managers").c_str() )
            .def_property_readonly("state_transition_interface",
                                   &tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE>::getStateTransitionAndSensitivityMatrixInterface,
                                   get_docstring("Estimator.state_transition_interface").c_str() )
            .def("perform_estimation",
                 &tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE>::estimateParameters,
                 py::arg( "estimation_input" ),
                 get_docstring("Estimator.perform_estimation").c_str() )
            .def("compute_covariance",
                 &tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE>::computeCovariance,
                 py::arg( "covariance_analysis_input" ),
                 get_docstring("Estimator.compute_covariance").c_str() )
            .def_property_readonly("variational_solver",
                                   &tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE>::getVariationalEquationsSolver,
                                   get_docstring("Estimator.variational_solver").c_str() );
};

}// namespace numerical_simulation
}// namespace tudatpy
