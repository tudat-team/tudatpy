/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_estimation_setup.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

namespace py = pybind11;
namespace tep = tudat::estimatable_parameters;
namespace tp = tudat::propagators;
namespace tss = tudat::simulation_setup;
namespace tni = tudat::numerical_integrators;


namespace tudatpy {

void expose_estimation_setup(py::module &m) {

    py::enum_<tep::EstimatebleParametersEnum >(m, "EstimatebleParameterTypes")
            .value("arc_wise_initial_body_state", tep::EstimatebleParametersEnum::arc_wise_initial_body_state)
            .value("initial_body_state", tep::EstimatebleParametersEnum::initial_body_state)
            .value("initial_rotational_body_state", tep::EstimatebleParametersEnum::initial_rotational_body_state)
            .value("gravitational_parameter", tep::EstimatebleParametersEnum::gravitational_parameter)
            .value("constant_drag_coefficient", tep::EstimatebleParametersEnum::constant_drag_coefficient)
            .value("radiation_pressure_coefficient", tep::EstimatebleParametersEnum::radiation_pressure_coefficient)
            .value("arc_wise_radiation_pressure_coefficient", tep::EstimatebleParametersEnum::arc_wise_radiation_pressure_coefficient)
            .value("spherical_harmonics_cosine_coefficient_block", tep::EstimatebleParametersEnum::arc_wise_initial_body_state)
            .value("spherical_harmonics_sine_coefficient_block", tep::EstimatebleParametersEnum::spherical_harmonics_sine_coefficient_block)
            .value("constant_rotation_rate", tep::EstimatebleParametersEnum::constant_rotation_rate)
            .value("rotation_pole_position", tep::EstimatebleParametersEnum::rotation_pole_position)
            .value("constant_additive_observation_bias", tep::EstimatebleParametersEnum::constant_additive_observation_bias)
            .value("arcwise_constant_additive_observation_bias", tep::EstimatebleParametersEnum::arcwise_constant_additive_observation_bias)
            .value("constant_relative_observation_bias", tep::EstimatebleParametersEnum::constant_relative_observation_bias)
            .value("arcwise_constant_relative_observation_bias", tep::EstimatebleParametersEnum::arcwise_constant_relative_observation_bias)
            .value("ppn_parameter_gamma", tep::EstimatebleParametersEnum::ppn_parameter_gamma)
            .value("ppn_parameter_beta", tep::EstimatebleParametersEnum::ppn_parameter_beta)
            .value("ground_station_position", tep::EstimatebleParametersEnum::ground_station_position)
            .value("equivalence_principle_lpi_violation_parameter", tep::EstimatebleParametersEnum::equivalence_principle_lpi_violation_parameter)
            .value("empirical_acceleration_coefficients", tep::EstimatebleParametersEnum::empirical_acceleration_coefficients)
            .value("arc_wise_empirical_acceleration_coefficients", tep::EstimatebleParametersEnum::arc_wise_empirical_acceleration_coefficients)
            .value("full_degree_tidal_love_number", tep::EstimatebleParametersEnum::full_degree_tidal_love_number)
            .value("single_degree_variable_tidal_love_number", tep::EstimatebleParametersEnum::single_degree_variable_tidal_love_number)
            .value("direct_dissipation_tidal_time_lag", tep::EstimatebleParametersEnum::direct_dissipation_tidal_time_lag)
            .value("mean_moment_of_inertia", tep::EstimatebleParametersEnum::mean_moment_of_inertia)
            .value("arc_wise_constant_drag_coefficient", tep::EstimatebleParametersEnum::arc_wise_constant_drag_coefficient)
            .value("periodic_spin_variation", tep::EstimatebleParametersEnum::periodic_spin_variation)
            .value("polar_motion_amplitude", tep::EstimatebleParametersEnum::polar_motion_amplitude)
            .value("core_factor", tep::EstimatebleParametersEnum::core_factor)
            .value("free_core_nutation_rate", tep::EstimatebleParametersEnum::free_core_nutation_rate)
            .value("desaturation_delta_v_values", tep::EstimatebleParametersEnum::desaturation_delta_v_values)
            .export_values();

    py::class_<tep::EstimatableParameterSettings,
            std::shared_ptr<tep::EstimatableParameterSettings>>(m, "EstimatableParameterSettings")
            .def(py::init<
                 const std::string,
                 const tep::EstimatebleParametersEnum,
                 const std::string>(),
                 py::arg("associated_body"),
                 py::arg("parameter_type"),
                 py::arg("point_on_body_id") = "");

    //TODO: Remove variationalOnlyIntegratorSettings
    py::class_<
            tp::SingleArcVariationalEquationsSolver<double, double>,
            std::shared_ptr<tp::SingleArcVariationalEquationsSolver<double, double>>>(m, "SingleArcVariationalEquationsSolver")
            .def(py::init<
                 const tudat::simulation_setup::NamedBodyMap&,
                 const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings<double>>,
                 const std::shared_ptr< tp::PropagatorSettings<double>>,
                 const std::shared_ptr< tep::EstimatableParameterSet< double > >,
                 const bool,
                 const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > >,
                 const bool,
                 const bool,
                 const bool >(),
                 py::arg("body_map"),
                 py::arg("integrator_settings"),
                 py::arg("propagator_settings"),
                 py::arg("estimated_parameters"),
                 py::arg("integrate_equations_concurrently") = true,
                 py::arg("variational_only_integrator_settings") = std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > >( ),
                 py::arg("clear_numerical_solutions") = false,
                 py::arg("integrate_on_creation") = false,
                 py::arg("set_integrated_result") = false )
            .def("integrate_equations_of_motion_only",
                 &tp::SingleArcVariationalEquationsSolver<double, double>::integrateDynamicalEquationsOfMotionOnly,
                 py::arg("initial_states"))
            .def("integrate_full_equations",
                 &tp::SingleArcVariationalEquationsSolver<double, double>::integrateVariationalAndDynamicalEquations,
                 py::arg("initial_states"),
                 py::arg("integrate_equations_concurrently"))
            .def("get_numerical_variational_equations_solution",
                 &tp::SingleArcVariationalEquationsSolver<double, double>::getNumericalVariationalEquationsSolution)
            .def("get_dynamics_simulator",
                 &tp::SingleArcVariationalEquationsSolver<double, double>::getDynamicsSimulator)
            .def("reset_parameter_estimate",
                 &tp::SingleArcVariationalEquationsSolver<double, double>::resetParameterEstimate);
}

}// namespace tudatpy
