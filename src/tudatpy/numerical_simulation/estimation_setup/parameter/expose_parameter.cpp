/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/simulation/estimation_setup/createEstimatableParameters.h>


namespace py = pybind11;

namespace tep = tudat::estimatable_parameters;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tba = tudat::basic_astrodynamics;


namespace tudatpy {
    namespace numerical_simulation {
        namespace estimation_setup {
            namespace parameter {

                PYBIND11_MODULE(expose_parameter, m) {
                    py::enum_<tep::EstimatebleParametersEnum>(
                        m, "EstimatableParameterTypes",
                        R"doc(Enumeration of model parameters that are available for estimation.
In order to establish a parameter estimation settings for a parameter of a certain type, use the factory function dedicated to this parameter type.
Note that not all of the listed types might be accessible via factory functions in the python interface yet.



	:member arc_wise_initial_body_state_type:
	:member initial_body_state_type:
	:member initial_rotational_body_state_type:
	:member constant_drag_coefficient_type:
	:member arc_wise_constant_drag_coefficient_type:
	:member radiation_pressure_coefficient_type:
	:member arc_wise_radiation_pressure_coefficient_type:
	:member empirical_acceleration_coefficients_type:
	:member arc_wise_empirical_acceleration_coefficients_type:
	:member desaturation_delta_v_values_type:
	:member gravitational_parameter_type:
	:member spherical_harmonics_cosine_coefficient_block_type:
	:member spherical_harmonics_sine_coefficient_block_type:
	:member mean_moment_of_inertia_type:
	:member constant_rotation_rate_type:
	:member rotation_pole_position_type:
	:member polar_motion_amplitude_type:
	:member core_factor_type:
	:member free_core_nutation_rate_type:
	:member periodic_spin_variation_type:
	:member constant_additive_observation_bias_type:
	:member arcwise_constant_additive_observation_bias_type:
	:member constant_relative_observation_bias_type:
	:member arcwise_constant_relative_observation_bias_type:
	:member ground_station_position_type:
	:member full_degree_tidal_love_number_type:
	:member single_degree_variable_tidal_love_number_type:
	:member direct_dissipation_tidal_time_lag_type:
	:member ppn_parameter_gamma_type:
	:member ppn_parameter_beta_type:
	:member equivalence_principle_lpi_violation_parameter_type:
)doc")
                        .value("arc_wise_initial_body_state_type",
                               tep::EstimatebleParametersEnum::
                                   arc_wise_initial_body_state)
                        .value(
                            "initial_body_state_type",
                            tep::EstimatebleParametersEnum::initial_body_state)
                        .value("initial_rotational_body_state_type",
                               tep::EstimatebleParametersEnum::
                                   initial_rotational_body_state)
                        .value("gravitational_parameter_type",
                               tep::EstimatebleParametersEnum::
                                   gravitational_parameter)
                        .value("constant_drag_coefficient_type",
                               tep::EstimatebleParametersEnum::
                                   constant_drag_coefficient)
                        .value("radiation_pressure_coefficient_type",
                               tep::EstimatebleParametersEnum::
                                   radiation_pressure_coefficient)
                        .value("arc_wise_radiation_pressure_coefficient_type",
                               tep::EstimatebleParametersEnum::
                                   arc_wise_radiation_pressure_coefficient)
                        .value(
                            "spherical_harmonics_cosine_coefficient_block_type",
                            tep::EstimatebleParametersEnum::
                                arc_wise_initial_body_state)
                        .value(
                            "spherical_harmonics_sine_coefficient_block_type",
                            tep::EstimatebleParametersEnum::
                                spherical_harmonics_sine_coefficient_block)
                        .value("constant_rotation_rate_type",
                               tep::EstimatebleParametersEnum::
                                   constant_rotation_rate)
                        .value("rotation_pole_position_type",
                               tep::EstimatebleParametersEnum::
                                   rotation_pole_position)
                        .value("constant_additive_observation_bias_type",
                               tep::EstimatebleParametersEnum::
                                   constant_additive_observation_bias)
                        .value(
                            "arcwise_constant_additive_observation_bias_type",
                            tep::EstimatebleParametersEnum::
                                arcwise_constant_additive_observation_bias)
                        .value("constant_relative_observation_bias_type",
                               tep::EstimatebleParametersEnum::
                                   constant_relative_observation_bias)
                        .value(
                            "arcwise_constant_relative_observation_bias_type",
                            tep::EstimatebleParametersEnum::
                                arcwise_constant_relative_observation_bias)
                        .value(
                            "ppn_parameter_gamma_type",
                            tep::EstimatebleParametersEnum::ppn_parameter_gamma)
                        .value(
                            "ppn_parameter_beta_type",
                            tep::EstimatebleParametersEnum::ppn_parameter_beta)
                        .value("ground_station_position_type",
                               tep::EstimatebleParametersEnum::
                                   ground_station_position)
                        .value(
                            "equivalence_principle_lpi_violation_parameter_"
                            "type",
                            tep::EstimatebleParametersEnum::
                                equivalence_principle_lpi_violation_parameter)
                        .value(
                            "empirical_acceleration_coefficients_type",
                            tep::EstimatebleParametersEnum::
                                empirical_acceleration_coefficients)  // TO
                                                                      // EXPOSE
                        .value(
                            "arc_wise_empirical_acceleration_coefficients_type",
                            tep::EstimatebleParametersEnum::
                                arc_wise_empirical_acceleration_coefficients)  // TO EXPOSE
                        .value("full_degree_tidal_love_number_type",
                               tep::EstimatebleParametersEnum::
                                   full_degree_tidal_love_number)  // TO EXPOSE
                        .value("single_degree_variable_tidal_love_number_type",
                               tep::EstimatebleParametersEnum::
                                   single_degree_variable_tidal_love_number)
                        .value("direct_dissipation_tidal_time_lag_type",
                               tep::EstimatebleParametersEnum::
                                   direct_dissipation_tidal_time_lag)
                        .value("mean_moment_of_inertia_type",
                               tep::EstimatebleParametersEnum::
                                   mean_moment_of_inertia)
                        .value("arc_wise_constant_drag_coefficient_type",
                               tep::EstimatebleParametersEnum::
                                   arc_wise_constant_drag_coefficient)
                        .value("periodic_spin_variation_type",
                               tep::EstimatebleParametersEnum::
                                   periodic_spin_variation)
                        .value("polar_motion_amplitude_type",
                               tep::EstimatebleParametersEnum::
                                   polar_motion_amplitude)
                        .value("core_factor_type",
                               tep::EstimatebleParametersEnum::core_factor)
                        .value("free_core_nutation_rate_type",
                               tep::EstimatebleParametersEnum::
                                   free_core_nutation_rate)
                        .value("desaturation_delta_v_values_type",
                               tep::EstimatebleParametersEnum::
                                   desaturation_delta_v_values)
                        .value("constant_time_drift_observation_bias_type",
                               tep::EstimatebleParametersEnum::
                                   constant_time_drift_observation_bias)
                        .value("arc_wise_time_drift_observation_bias_type",
                               tep::EstimatebleParametersEnum::
                                   arc_wise_time_drift_observation_bias)
                        .value("inverse_tidal_quality_factor_type",
                               tep::EstimatebleParametersEnum::
                                   inverse_tidal_quality_factor)
                        .export_values();

                    py::enum_<tba::EmpiricalAccelerationComponents>(
                        m, "EmpiricalAccelerationComponents", "")
                        .value("radial_empirical_acceleration_component",
                               tba::EmpiricalAccelerationComponents::
                                   radial_empirical_acceleration_component)
                        .value("along_track_empirical_acceleration_component",
                               tba::EmpiricalAccelerationComponents::
                                   along_track_empirical_acceleration_component)
                        .value(
                            "across_track_empirical_acceleration_component",
                            tba::EmpiricalAccelerationComponents::
                                across_track_empirical_acceleration_component)
                        .export_values();

                    py::enum_<tba::EmpiricalAccelerationFunctionalShapes>(
                        m, "EmpiricalAccelerationFunctionalShapes", "")
                        .value("constant_empirical",
                               tba::EmpiricalAccelerationFunctionalShapes::
                                   constant_empirical)
                        .value("sine_empirical",
                               tba::EmpiricalAccelerationFunctionalShapes::
                                   sine_empirical)
                        .value("cosine_empirical",
                               tba::EmpiricalAccelerationFunctionalShapes::
                                   cosine_empirical)
                        .export_values();

                    py::class_<tep::CustomAccelerationPartialSettings,
                               std::shared_ptr<
                                   tep::CustomAccelerationPartialSettings>>(
                        m, "CustomAccelerationPartialSettings", "");


                    m.def("custom_analytical_partial",
                          &tep::analyticalAccelerationPartialSettings,
                          py::arg("analytical_partial_function"),
                          py::arg("body_undergoing_acceleration"),
                          py::arg("body_exerting_acceleration"),
                          py::arg("acceleration_type"), "");

                    m.def("custom_numerical_partial",
                          &tep::numericalAccelerationPartialSettings,
                          py::arg("parameter_perturbation"),
                          py::arg("body_undergoing_acceleration"),
                          py::arg("body_exerting_acceleration"),
                          py::arg("acceleration_type"),
                          py::arg("environment_updates") =
                              std::map<tp::EnvironmentModelsToUpdate,
                                       std::vector<std::string>>(),
                          "");

                    py::class_<
                        tep::EstimatableParameterSettings,
                        std::shared_ptr<tep::EstimatableParameterSettings>>(
                        m, "EstimatableParameterSettings",
                        R"doc(Base class to defining settings of parameter to be estimated.

	Functional (base) class for settings of model parameter to be estimated.
	Settings of simple parameters types are managed via this class, more complex parameter types are handled by specialised derivates of this class.
	Instances of either base or derived class can be created via dedicated factory functions.

)doc")
                        .def_readwrite("custom_partial_settings",
                                       &tep::EstimatableParameterSettings::
                                           customPartialSettings_);


                    // ###############    Initial States
                    // ################################

                    m.def(
                        "initial_states",
                        &tss::getInitialStateParameterSettings<double>,
                        py::arg("propagator_settings"), py::arg("bodies"),
                        py::arg("arc_initial_times") = std::vector<double>(),
                        R"doc(Function for defining parameter settings for initial state parameters.

	Factory function for creating a (linear sensitivity) parameter settings object for initial state parameters.
	The factory function uses the propagator settings to determine which type of initial state parameter (single/multi/hybrid-arc; translational/rotational/... dynamics) is to be estimated,
	e.g. if a single-arc translational state propagator is defined, the function will automatically create the parameters for the associated initial state parameter

	Note: These function return lists of parameter settings objects.
	This means that which the return of this function cannot simply be added to the parameter settings objects of single parameters in a list creation statement.
	Instead, list concatenation is recommended. Please see the following example:

	.. code-block:: python

	   # define single estimatable parameters
	   single_parameter_1 = ...
	   single_parameter_2 = ...
	   ...

	   # bad: list creation statement --> will result in nested list, undesired!
	   list_of_all_parameters = [estimation_setup.parameter.initial_states(...), single_parameter_1, single_parameter_2, ...]

	   # better: list concatenation --> will result in simple list, desired!
	   list_of_all_parameters = estimation_setup.parameter.initial_states(...) + [single_parameter_1, single_parameter_2, ...]


	:param propagator_settings:
		Object containing the consolidated propagation settings of the simulation in the context of which the given model parameters are to be estimated.

	:param bodies:
		Object consolidating all bodies and environment models that constitute the physical environment.

	:param arc_initial_times:
		Initial times of arcs, only required if arc-wise propagation settings are passed via the `propagator_settings` object.

	:return:
		List of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` objects, one per component of each initial state in the simulation.

)doc");


                    // ###############    Vehicle Model Parameters
                    // ################################

                    m.def(
                        "constant_drag_coefficient",
                        &tep::constantDragCoefficient, py::arg("body"),
                        R"doc(Function for defining parameter settings for constant drag coefficients.

	Factory function for creating a (linear sensitivity) parameter settings object for constant drag coefficients.
	Using the constant drag coefficient as an estimatable parameter requires:
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.constant` aerodynamic interface to be defined for the body specified by the ``body`` parameter
	  * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.aerodynamic` acceleration


	:param body:
		Name of the body, with whose drag acceleration model the estimatable parameter is associated.

	:return:
		:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's constant drag coefficient.

)doc");

                    m.def(
                        "arcwise_constant_drag_coefficient",
                        &tep::arcwiseDragCoefficient, py::arg("body"),
                        py::arg("arc_initial_times"),
                        R"doc(Function for defining parameter settings for arc-wise constant drag coefficients.

	Factory function for creating (linear sensitivity) parameter settings object for arc-wise constant drag coefficients (arc-wise version of :func:``~tudatpy.numerical_simulation.estimation_setup.parameter.constant_drag_coefficient`).
	Using the arc-wise constant drag coefficient as an estimatable parameter requires:
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.constant` aerodynamic interface to be defined for the body specified by the ``body`` parameter
	  * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.aerodynamic` acceleration

	Note: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the drag coefficient may, but need not, correspond to the arcs used for a multi-arc propagation.


	:param body:
		Name of the body, with whose drag acceleration model the estimatable parameter is associated.

	:param arc_initial_times:
		List of times at which the arcs over which the drag coefficient is to be estimated will start.

	:return:
		Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.ArcWiseDragCoefficientEstimatableParameterSettings` class
		for arc-wise treatment of the specified body's constant drag coefficient.

)doc");

                    m.def(
                        "radiation_pressure_coefficient",
                        &tep::radiationPressureCoefficient, py::arg("body"),
                        R"doc(Function for defining parameter settings for radiation pressure coefficients.

	Factory function for creating a (linear sensitivity) parameter settings object for a radiation pressure coefficient.
	Using the radiation pressure coefficient as an estimatable parameter requires
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.radiation_pressure.cannonball` radiation pressure interface to be defined for the body specified by the ``body`` parameter
	  * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.cannonball_radiation_pressure` acceleration


	:param body:
		Name of the body, with whose radiation pressure acceleration model the estimatable parameter is associated.

	:return:
		:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's radiation pressure coefficient.

)doc");

                    m.def(
                        "arcwise_radiation_pressure_coefficient",
                        &tep::arcwiseRadiationPressureCoefficient,
                        py::arg("body"), py::arg("arc_initial_times"),
                        R"doc(Function for defining parameter settings for arc-wise radiation pressure coefficients.

	Factory function for creating a (linear sensitivity) parameter settings object for arc-wise radiation pressure coefficients (arc-wise version of :func:``~tudatpy.numerical_simulation.estimation_setup.parameter.radiation_pressure_coefficient`).
	Using the radiation pressure coefficient as an estimatable parameter requires
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.radiation_pressure.cannonball` radiation pressure interface to be defined for the body specified by the ``body`` parameter
	  * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.cannonball_radiation_pressure` acceleration
	The radiation pressure coefficient is defined according to the universal convention for a cannonball model and thus no further model information is given.

	Note: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the radiation pressure coefficient may, but need not, correspond to the arcs used for a multi-arc propagation.


	:param body:
		Name of the body, with whose radiation pressure acceleration model the estimatable parameter is associated.

	:param arc_initial_times:
		List of times at which the arcs over which the radiation pressure coefficient is to be estimated will start.

	:return:
		Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.ArcWiseRadiationPressureCoefficientEstimatableParameterSettings` class
		for arc-wise treatment of the specified body's radiation pressure coefficient.

)doc");


                    m.def("radiation_pressure_target_direction_scaling",
                          &tep::radiationPressureTargetDirectionScaling,
                          py::arg("target_body"), py::arg("source_body"), "");

                    m.def(
                        "radiation_pressure_target_perpendicular_direction_"
                        "scaling",
                        &tep::
                            radiationPressureTargetPerpendicularDirectionScaling,
                        py::arg("target_body"), py::arg("source_body"), "");

                    m.def(
                        "constant_empirical_acceleration_terms",
                        &tep::constantEmpiricalAccelerationMagnitudes,
                        py::arg("body"), py::arg("centralBody"),
                        R"doc(Function for defining parameter settings for constant empirical acceleration terms.

	As :func:`~tudatpy.numerical_simulation.estimation_setup.parameter.empirical_accelerations`, but only using the constant R, S and W components (no sine or cosine term estimation). This function is added as a function of convenience


	:param body:
		Name of the body, with whose empirical acceleration model the estimatable parameter is associated.

	:param centralBody:
		Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration

	:return:
		Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EmpiricalAccelerationEstimatableParameterSettings` class
		for the specified body's empirical acceleration terms.

)doc");

                    m.def(
                        "full_empirical_acceleration_terms",
                        &tep::empiricalAccelerationMagnitudesFull,
                        py::arg("body"), py::arg("centralBody"),
                        R"doc(Function for defining parameter settings for empirical acceleration magnitudes.

	Factory function for creating a (linear sensitivity) parameter settings object for empirical acceleration magnitudes.
	Using the empirical acceleration terms as estimatable parameters requires
	  * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.empirical` acceleration, which include constant (in RSW frame) terms


	:param body:
		Name of the body, with whose empirical acceleration model the estimatable parameter is associated.

	:param centralBody:
		Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration

	:param acceleration_components:
		Dictionary of components of the empirical acceleration which are to be estimated. There are two 'degrees of freedom' in these components: the direction of the acceleration (e.g. R, S or W direction) and the temporal signature (constant, sine of true anomaly or cosine of true anomaly). With this input, any subset may be selected. This parameter is a dictionary, with the key denoting the direction of the acceleration, and the value a list of the temporal signatures to estimate for this empirical acceleration direction.

	:return:
		Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EmpiricalAccelerationEstimatableParameterSettings` class
		for the specified body's empirical acceleration terms.

)doc");

                    m.def(
                        "empirical_accelerations",
                        &tep::empiricalAccelerationMagnitudes, py::arg("body"),
                        py::arg("centralBody"),
                        py::arg("acceleration_components"),
                        R"doc(Function for defining parameter settings for empirical acceleration magnitudes.

	Factory function for creating a (linear sensitivity) parameter settings object for empirical acceleration magnitudes.
	Using the empirical acceleration terms as estimatable parameters requires
	  * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.empirical` acceleration, which include constant (in RSW frame) terms


	:param body:
		Name of the body, with whose empirical acceleration model the estimatable parameter is associated.

	:param centralBody:
		Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration

	:param acceleration_components:
		Dictionary of components of the empirical acceleration which are to be estimated. There are two 'degrees of freedom' in these components: the direction of the acceleration (e.g. R, S or W direction) and the temporal signature (constant, sine of true anomaly or cosine of true anomaly). With this input, any subset may be selected. This parameter is a dictionary, with the key denoting the direction of the acceleration, and the value a list of the temporal signatures to estimate for this empirical acceleration direction.

	:return:
		Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EmpiricalAccelerationEstimatableParameterSettings` class
		for the specified body's empirical acceleration terms.

)doc");

                    m.def(
                        "arcwise_empirical_accelerations",
                        &tep::arcWiseEmpiricalAccelerationMagnitudes,
                        py::arg("body"), py::arg("centralBody"),
                        py::arg("acceleration_components"),
                        py::arg("arc_start_times"),
                        R"doc(Function for defining parameter settings for arc-wise empirical acceleration magnitudes.

	Factory function for creating a (linear sensitivity) parameter settings object for arc-wise empirical acceleration magnitudes (arc-wise version of :func:``~tudatpy.numerical_simulation.estimation_setup.parameter.empirical_accelerations`).
	Using the empirical acceleration terms as estimatable parameters requires
	  * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.empirical` acceleration, which include constant (in RSW frame) terms

	Note: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the radiation pressure coefficient may, but need not, correspond to the arcs used for a multi-arc propagation.


	:param body:
		Name of the body, with whose empirical acceleration model the estimatable parameter is associated.

	:param centralBody:
		Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration

	:param acceleration_components:
		Dictionary of components of the empirical acceleration which are to be estimated. There are two 'degrees of freedom' in these components: the direction of the acceleration (e.g. R, S or W direction) and the temporal signature (constant, sine of true anomaly or cosine of true anomaly). With this input, any subset may be selected. This parameter is a dictionary, with the key denoting the direction of the acceleration, and the value a list of the temporal signatures to estimate for this empirical acceleration direction.

	:param arc_initial_times:
		List of times at which the arcs over which the empirical accelerations are to be estimated will start.

	:return:
		Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EmpiricalAccelerationEstimatableParameterSettings` class
		for the specified body's arc-wise empirical acceleration terms.

)doc");

                    m.def(
                        "arcwise_constant_empirical_acceleration_terms",
                        &tep::constantArcWiseEmpiricalAccelerationMagnitudes,
                        py::arg("body"), py::arg("centralBody"),
                        py::arg("arc_start_times"),
                        R"doc(Function for defining parameter settings for arc-wise constant empirical acceleration terms.

	As :func:`~tudatpy.numerical_simulation.estimation_setup.parameter.arcwise_empirical_accelerations`, but only using the constant R, S and W components (no sine or cosine term estimation). This function is added as a function of convenience


	:param body:
		Name of the body, with whose empirical acceleration model the estimatable parameter is associated.

	:param centralBody:
		Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration

	:param arc_initial_times:
		List of times at which the arcs over which the empirical accelerations are to be estimated will start.

	:return:
		Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EmpiricalAccelerationEstimatableParameterSettings` class
		for the specified body's arc-wise constant empirical acceleration terms.

)doc");

                    m.def(
                        "quasi_impulsive_shots", &tep::quasiImpulsiveShots,
                        py::arg("body"),
                        R"doc(Function for defining parameter settings for quasi-impulsive shots.

	Factory function for creating a (linear sensitivity) parameter settings object for so-called 'quasi-impulsive shots', such as desaturation maneuvers. With this parameter, the total :math:`\Delta V` vector of a set of such shots/maneuvers can be estimated.
	Using the quasi-impulsive shots as an estimatable parameter requires
	  * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.quasi_impulsive_shots_acceleration` acceleration
	Note: this parameter considers *all* shots/maneuvers used in the above acceleration model, and estimates the value of the 'delta_v_values' input of this acceleration.


	:param body:
		Name of the body, with which the quasi-impulsive shot estimatable parameter is associated.

	:return:
		:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's quasi-impulsive shots

)doc");


                    // ###############    Gravity Model Parameters
                    // ################################

                    m.def(
                        "gravitational_parameter", &tep::gravitationalParameter,
                        py::arg("body"),
                        R"doc(Function for defining parameter settings for a massive body's gravitational parameter.

	Factory function for creating a (linear sensitivity) parameter settings object for the gravitational parameter of massive bodies.
	Using the gravitational parameter as estimatable parameter requires
	  * The body specified by the ``body`` parameter to be endowed with a gravity field (see :ref:`\`\`gravity_field\`\`` module for options)
	  * Any dynamical or observational model to depend on the gravitational parameter of the body specified by the ``body`` parameter


	:param body:
		Name of the body, with whose gravitational model the estimatable parameter is associated.

	:return:
		:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's gravitational parameter.

)doc");

                    m.def(
                        "spherical_harmonics_c_coefficients",
                        py::overload_cast<const std::string, const int,
                                          const int, const int, const int>(
                            &tep::sphericalHarmonicsCosineBlock),
                        py::arg("body"), py::arg("minimum_degree"),
                        py::arg("minimum_order"), py::arg("maximum_degree"),
                        py::arg("maximum_order"),
                        R"doc(Function for defining parameter settings for the cosine coefficients of body's spherical harmonics gravitational model.

	Factory function for creating a (linear sensitivity) parameter settings object for the spherical harmonics cosine-coefficients (:math:`\bar{C}_{lm}`) of a hody with a spherical harmonic gravity field. Using this function, a 'full' set of spherical harmonic coefficients between an minimum/maximum degree/order are estimated. For instance, for minimum degree/order of 2/0, and maximum degree/order 4/4, all spherical harmonic cosine coefficients of degrees 2, 3 and 4 are estimated. If the maximum degree/order is set to 4/2, only coefficients with an order of 0, 1 and 2 are included. The entries in the parameter are sorted first by degree, and then by order (both in ascending order)
	Using the spherical harmonics cosine coefficients as estimatable parameter requires:
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field.spherical_harmonic` (or derived) gravity model to be defined for the body specified by the ``body`` parameter
	  * Any dynamical or observational model to depend on the estimated cosine coefficients of the body specified by the ``body`` parameter. Typically, this dependency will be a :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.spherical_harmonic` acceleration


	:param body:
		Name of the body, with whose gravitational model the estimatable parameters are associated.

	:param minimum_degree:
		Minimum degree of c-coefficients to be included.
	:param minimum_order:
		Minimum order of c-coefficients to be included.
	:param maximum_degree:
		Maximum degree of c-coefficients to be included.
	:param maximum_order:
		Maximum order of c-coefficients to be included.
	:return:
		Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.SphericalHarmonicEstimatableParameterSettings` class
		for the applicable spherical harmonics c-coefficients of the specified body's gravitational model.

)doc");

                    m.def(
                        "spherical_harmonics_c_coefficients_block",
                        py::overload_cast<const std::string,
                                          std::vector<std::pair<int, int>>>(
                            &tep::sphericalHarmonicsCosineBlock),
                        py::arg("body"), py::arg("block_indices"),
                        R"doc(Function for defining parameter settings for the cosine coefficients of body's spherical harmonics gravitational model.

	As :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.spherical_harmonics_c_coefficients`, but with a manually defined set of coefficients.


	:param body:
		Name of the body, with whose gravitational model the estimatable parameters are associated.

	:param block_indices:
		List of block indices. The length of this list can be arbitrary, as long as the pairs are unique.
		For each pair, the first value is the degree and the second the order of the coefficient to be included.

	:return:
		Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.SphericalHarmonicEstimatableParameterSettings` class
		for the applicable spherical harmonics c-coefficients of the specified body's gravitational model.

)doc");

                    m.def(
                        "spherical_harmonics_s_coefficients",
                        py::overload_cast<const std::string, const int,
                                          const int, const int, const int>(
                            &tep::sphericalHarmonicsSineBlock),
                        py::arg("body"), py::arg("minimum_degree"),
                        py::arg("minimum_order"), py::arg("maximum_degree"),
                        py::arg("maximum_order"),
                        R"doc(Function for defining parameter settings for the sine coefficients of body's spherical harmonics gravitational model.

	Factory function for creating a (linear sensitivity) parameter settings object for the spherical harmonics sine-coefficients (:math:`\bar{S}_{lm}`) of a hody with a spherical harmonic gravity field. Using this function, a 'full' set of spherical harmonic coefficients between an minimum/maximum degree/order are estimated. For instance, for minimum degree/order of 2/1 (there is no order 0 sine coefficient), and maximum degree/order 4/4, all spherical harmonic sine coefficients of degrees 2, 3 and 4 are estimated. If the maximum degree/order is set to 4/2, only coefficients with an order of 1 and 2 are included. The entries in the parameter are sorted first by degree, and then by order (both in ascending order)
	Using the spherical harmonics cosine coefficients as estimatable parameter requires:
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field.spherical_harmonic` (or derived) gravity model to be defined for the body specified by the ``body`` parameter
	  * Any dynamical or observational model to depend on the estimated cosine coefficients of the body specified by the ``body`` parameter. Typically, this dependency will be a :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.spherical_harmonic` acceleration


	:param body:
		Name of the body, with whose gravitational model the estimatable parameters are associated.

	:param minimum_degree:
		Minimum degree of s-coefficients to be included.
	:param minimum_order:
		Minimum order of s-coefficients to be included.
	:param maximum_degree:
		Maximum degree of s-coefficients to be included.
	:param maximum_order:
		Maximum order of s-coefficients to be included.
	:return:
		Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.SphericalHarmonicEstimatableParameterSettings` class
		for the applicable spherical harmonics s-coefficients of the specified body's gravitational model.

)doc");

                    m.def(
                        "spherical_harmonics_s_coefficients_block",
                        py::overload_cast<const std::string,
                                          std::vector<std::pair<int, int>>>(
                            &tep::sphericalHarmonicsSineBlock),
                        py::arg("body"), py::arg("block_indices"),
                        R"doc(Function for defining parameter settings for the sine coefficients of body's spherical harmonics gravitational model.

	As :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.spherical_harmonics_s_coefficients`, but with a manually defined set of coefficients.


	:param body:
		Name of the body, with whose gravitational model the estimatable parameters are associated.

	:param block_indices:
		List of block indices. The length of this list can be arbitrary, as long as the pairs are unique.
		For each pair, the first value is the degree and the second the order of the coefficient to be included.

	:return:
		Instance of :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.SphericalHarmonicEstimatableParameterSettings` class
		for the applicable spherical harmonics s-coefficients of the specified body's gravitational model.

)doc");


                    // ###############    Rotation Model Parameters
                    // ################################

                    m.def(
                        "mean_moment_of_inertia", &tep::meanMomentOfInertia,
                        py::arg("body"),
                        R"doc(Function for defining parameter settings for a body's mean moment of inertia.

	Factory function for creating a (linear sensitivity) parameter settings object for a body's mean moment of inertia.
	In most cases, the mean moment of inertia will not influence the dynamics/observation directly and sensitivity to this parameter will not be included. The dynamics/observation will be sensitive to this parameter if the rotational dynamics of a relevant body is estimated.
	Using the mean moment of inertia as estimatable parameter requires
	  * The estimation of an initial rotational state of the body specified by the ``body`` parameter


	:param body:
		Name of the body, with whose body model the estimatable parameter is associated.

	:return:
		:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's mean moment of inertia.

)doc");

                    m.def(
                        "constant_rotation_rate", &tep::constantRotationRate,
                        py::arg("body"),
                        R"doc(Function for defining parameter settings for a body's constant rotation rate.

	Factory function for creating a (linear sensitivity) parameter settings object for a body's constant rotation rate parameter.
	Using the constant rotation rate as estimatable parameter requires:
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.simple` or :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.simple_from_spice` rotation model specified by the ``body`` parameter
	  * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter


	:param body:
		Name of the body, with whose rotation model the estimatable parameter is associated.

	:return:
		:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's constant spin rate.

)doc");

                    m.def(
                        "rotation_pole_position", &tep::rotationPolePosition,
                        py::arg("body"),
                        R"doc(Function for defining parameter settings for a body's rotation pole position.

	Factory function for creating a (linear sensitivity) parameter settings object for a body's rotation pole position, parameterized by the constant pole rotation angles (:math:`alpha` and :math:`\delta`).
	Using the rotation pole position as estimatable parameter requires:
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.simple` or :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.simple_from_spice` rotation model specified by the ``body`` parameter
	  * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter


	:param body:
		Name of the body, with whose rotation model the estimatable parameter is associated.

	:return:
		:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's rotation pole position.

)doc");

                    m.def(
                        "core_factor", &tep::coreFactor, py::arg("body"),
                        R"doc(Function for defining parameter settings for a body's core factor.

	Factory function for creating a (linear sensitivity) parameter settings object for a body's core factor.
	Using the core factor as estimatable parameter requires
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.mars_high_accuracy` rotation model specified by the ``body`` parameter
	  * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter


	:param body:
		Name of the body, with whose rotation model the estimatable parameter is associated.

	:return:
		:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's core factor.

)doc");

                    m.def(
                        "free_core_nutation_rate", &tep::freeCoreNutationRate,
                        py::arg("body"),
                        R"doc(Function for defining parameter settings for a body's free core nutation rate.

	Factory function for creating a (linear sensitivity) parameter settings object for a body's free core nutation rate.
	Using the free core nutation rate as estimatable parameter requires
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.mars_high_accuracy` rotation model specified by the ``body`` parameter
	  * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter


	:param body:
		Name of the body, with whose rotation model the estimatable parameter is associated.

	:return:
		:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's free core nutation rate.

)doc");

                    m.def(
                        "periodic_spin_variations",
                        &tep::periodicSpinVariations, py::arg("body"),
                        R"doc(Function for defining parameter settings for a body's periodic spin variations.

	Factory function for creating a (linear sensitivity) parameter settings object for a body's periodic spin variation parameters.
	Using the mean moment of inertia as estimatable parameter requires
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.mars_high_accuracy` rotation model specified by the ``body`` parameter
	  * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter


	:param body:
		Name of the body, with whose rotation model the estimatable parameter is associated.

	:return:
		:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's periodic spin variations.

)doc");

                    m.def(
                        "polar_motion_amplitudes", &tep::polarMotionAmplitudes,
                        py::arg("body"),
                        R"doc(Function for defining parameter settings for a body's polar motion amplitudes.

	Factory function for creating a (linear sensitivity) parameter settings object for a body's polar motion amplitudes.
	Using the polar motion amplitudes as estimatable parameter requires
	  * A :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.mars_high_accuracy` rotation model specified by the ``body`` parameter
	  * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter


	:param body:
		Name of the body, with whose rotation model the estimatable parameter is associated.

	:return:
		:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified body's polar motion amplitudes.

)doc");


                    // ###############   Observation Model Parameters
                    // ################################

                    m.def(
                        "absolute_observation_bias", &tep::observationBias,
                        py::arg("link_ends"), py::arg("observable_type"),
                        R"doc(Function for defining parameter settings for an absolute observation bias.

	Factory function for creating a (linear sensitivity) parameter settings object for an observation's absolute bias parameter.
	Using the absolute observation bias as estimatable parameter requires:
	  * The observation model (corresponding to the `link_ends` and `observable_type`) to include an absolute bias (:func:`~tudatpy.numerical_simulation.estimation_setup.observation.absolute_bias`)


	:param link_ends:
		Set of link ends that define the geometry of the biased observations.

	:param observable_type:
		Observable type of the biased observations.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.ConstantObservationBiasEstimatableParameterSettings`
		for the specified observation's arc-wise absolute bias.

)doc");

                    m.def(
                        "arcwise_absolute_observation_bias",
                        &tep::arcwiseObservationBias, py::arg("link_ends"),
                        py::arg("observable_type"), py::arg("arc_start_times"),
                        py::arg("time_link_end"),
                        R"doc(Function for defining parameter settings for arc-wise absolute observation bias.

	Factory function for creating a (linear sensitivity) parameter settings object for the arc-wise treatment of an observation's absolute bias parameter.
	Using the arc-wise absolute observation bias as estimatable parameter requires
	  * The observation model (corresponding to the `link_ends` and `observable_type`) to include an arc-wise absolute bias (:func:`~tudatpy.numerical_simulation.estimation_setup.observation.arcwise_absolute_bias`)


	:param link_ends:
		Set of link ends that define the geometry of the biased observations.

	:param observable_type:
		Observable type of the biased observations.
	:param arc_start_times:
		List of times at which the arcs over which the bias is to be estimated will start.
	:param time_link_end:
		The link end type (transmitter, receiver, etc.) at which the arc_start_times is evaluated.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.ArcWiseConstantObservationBiasEstimatableParameterSettings`
		for the specified observation's arc-wise absolute bias.

)doc");

                    m.def(
                        "relative_observation_bias",
                        &tep::relativeObservationBias, py::arg("link_ends"),
                        py::arg("observable_type"),
                        R"doc(Function for defining parameter settings for an relative observation bias.

	Factory function for creating a (linear sensitivity) parameter settings object for an observation's relative bias parameter.
	Using the relative observation bias as estimatable parameter requires
	  * The observation model (corresponding to the `link_ends` and `observable_type`) to include a relative bias (:func:`~tudatpy.numerical_simulation.estimation_setup.observation.relative_bias`)
	Note: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the bias may, but need not, correspond to the arcs used for a multi-arc propagation.


	:param link_ends:
		Set of link ends that define the geometry of the biased observations.

	:param observable_type:
		Observable type of the biased observations.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.ConstantObservationBiasEstimatableParameterSettings`
		for the specified observation's arc-wise relative bias.

)doc");

                    m.def(
                        "arcwise_relative_observation_bias",
                        &tep::arcwiseRelativeObservationBias,
                        py::arg("link_ends"), py::arg("observable_type"),
                        py::arg("arc_start_times"), py::arg("time_link_end"),
                        R"doc(Function for defining parameter settings for arc-wise absolute observation bias.

	Factory function for creating a (linear sensitivity) parameter settings object for the arc-wise treatment of an observation's relative bias parameter.
	Using the arc-wise relative observation bias as estimatable parameter requires
	  * The observation model (corresponding to the `link_ends` and `observable_type`) to include an arc-wise relative bias (:func:`~tudatpy.numerical_simulation.estimation_setup.observation.arcwise_relative_bias`)
	Note: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the bias may, but need not, correspond to the arcs used for a multi-arc propagation.


	:param link_ends:
		Set of link ends that define the geometry of the biased observations.

	:param observable_type:
		Observable type of the biased observations.
	:param arc_start_times:
		List of times at which the arcs over which the bias is to be estimated will start.
	:param time_link_end:
		The link end type (transmitter, receiver, etc.) at which the arc_start_times is evaluated.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.ArcWiseConstantObservationBiasEstimatableParameterSettings`
		for the specified observation's arc-wise relative bias.

)doc");

                    m.def("time_drift_observation_bias",
                          &tep::timeDriftObservationBias, py::arg("link_ends"),
                          py::arg("observable_type"), py::arg("ref_epoch"),
                          py::arg("time_link_end"));

                    m.def("arcwise_time_drift_observation_bias",
                          &tep::arcwiseTimeDriftObservationBias,
                          py::arg("link_ends"), py::arg("observable_type"),
                          py::arg("arc_start_times"), py::arg("ref_epochs"),
                          py::arg("time_link_end"));

                    m.def(
                        "ground_station_position", &tep::groundStationPosition,
                        py::arg("body"), py::arg("ground_station_name"),
                        R"doc(Function for defining parameter settings for ground station position bias.

	Factory function for creating a (linear sensitivity) parameter settings object for a ground station's body-fixed Cartesian position.
	Using the ground station position bias as estimatable parameter requires:
	  * At least one observation model to rely on the specified ground station


	:param body:
		Body name identifying the body, with which the ground station is associated.
	:param ground_station_name:
		Name which identifies the position-biased ground station.
	:return:
		:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for the specified ground station's position bias.

)doc");


                    // ###############  Tidal Model Parameters
                    // ################################

                    m.def("direct_tidal_dissipation_time_lag",
                          py::overload_cast<const std::string&,
                                            const std::string&>(
                              &tep::directTidalDissipationLagTime),
                          py::arg("body"), py::arg("deforming_body"), "");

                    m.def("direct_tidal_dissipation_time_lag",
                          py::overload_cast<const std::string&,
                                            const std::vector<std::string>&>(
                              &tep::directTidalDissipationLagTime),
                          py::arg("body"), py::arg("deforming_body"), "");

                    m.def("inverse_tidal_quality_factor",
                          py::overload_cast<const std::string&,
                                            const std::string&>(
                              &tep::inverseTidalQualityFactor),
                          py::arg("body"), py::arg("deforming_body"), "");

                    m.def("inverse_tidal_quality_factor",
                          py::overload_cast<const std::string&,
                                            const std::vector<std::string>&>(
                              &tep::inverseTidalQualityFactor),
                          py::arg("body"), py::arg("deforming_body"), "");

                    m.def("order_invariant_k_love_number",
                          py::overload_cast<const std::string&, const int,
                                            const std::string, const bool>(
                              &tep::orderInvariantKLoveNumber),
                          py::arg("deformed_body"), py::arg("degree"),
                          py::arg("deforming_body"),
                          py::arg("use_complex_love_number") = 0, "");

                    m.def("order_invariant_k_love_number",
                          py::overload_cast<const std::string&, const int,
                                            const std::vector<std::string>&,
                                            const bool>(
                              &tep::orderInvariantKLoveNumber),
                          py::arg("deformed_body"), py::arg("degree"),
                          py::arg("deforming_bodies"),
                          py::arg("use_complex_love_number") = 0, "");

                    m.def("order_invariant_k_love_number",
                          py::overload_cast<const std::string&, const int,
                                            const bool>(
                              &tep::orderInvariantKLoveNumber),
                          py::arg("deformed_body"), py::arg("degree"),
                          py::arg("use_complex_love_number") = 0, "");

                    m.def("order_varying_k_love_number",
                          py::overload_cast<const std::string&, const int,
                                            const std::vector<int>&,
                                            const std::string, const bool>(
                              &tep::orderVaryingKLoveNumber),
                          py::arg("deformed_body"), py::arg("degree"),
                          py::arg("orders"), py::arg("deforming_body"),
                          py::arg("use_complex_love_number") = 0, "");

                    m.def("order_varying_k_love_number",
                          py::overload_cast<const std::string&, const int,
                                            const std::vector<int>&,
                                            const std::vector<std::string>&,
                                            const bool>(
                              &tep::orderVaryingKLoveNumber),
                          py::arg("deformed_body"), py::arg("degree"),
                          py::arg("orders"), py::arg("deforming_bodies"),
                          py::arg("use_complex_love_number") = 0, "");

                    m.def(
                        "order_varying_k_love_number",
                        py::overload_cast<const std::string&, const int,
                                          const std::vector<int>&, const bool>(
                            &tep::orderVaryingKLoveNumber),
                        py::arg("deformed_body"), py::arg("degree"),
                        py::arg("orders"),
                        py::arg("use_complex_love_number") = 0, "");

                    m.def("polynomial_gravity_field_variation_amplitudes",
                          &tep::polynomialGravityFieldVariationParameter,
                          py::arg("body_name"),
                          py::arg("cosine_indices_per_power"),
                          py::arg("sine_indices_per_power"), "");

                    m.def(
                        "monomial_gravity_field_variation_amplitudes",
                        &tep::
                            polynomialSinglePowerGravityFieldVariationParameter,
                        py::arg("body_name"), py::arg("power"),
                        py::arg("cosine_indices"), py::arg("sine_indices"), "");

                    m.def(
                        "monomial_full_block_gravity_field_variation_"
                        "amplitudes",
                        &tep::
                            polynomialSinglePowerFullBlockGravityFieldVariationParameter,
                        py::arg("body_name"), py::arg("power"),
                        py::arg("minimum_degree"), py::arg("minimum_order"),
                        py::arg("maximum_degree"), py::arg("maximum_order"),
                        "");

                    m.def("scaled_longitude_libration_amplitude",
                          &tep::scaledLongitudeLibrationAmplitude,
                          py::arg("body_name"), "");

                    m.def("yarkovsky_parameter", &tep::yarkovskyParameter,
                          py::arg("body_name"),
                          py::arg("central_body_name") = "Sun", "");

                    m.def("custom_parameter", &tep::customParameterSettings,
                          py::arg("custom_id"), py::arg("parameter_size"),
                          py::arg("get_parameter_function"),
                          py::arg("set_parameter_function"), "");


                    // ###############  Global (GR) Model Parameters
                    // ################################

                    m.def(
                        "ppn_parameter_gamma", &tep::ppnParameterGamma,
                        R"doc(Function for defining parameter settings for post-newtonian gamma parameter.

	Factory function for creating a (linear sensitivity) parameter settings object for a global PPN :math:`\gamma` parameter.
	Using the post-newtonian gamma parameter as estimatable parameter requires at least one of the following:
	 * An acceleration model depending on this parameter, such as :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.relativistic_correction`
	 * An observation model with a light-time correction depending on this parameter, such as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.first_order_relativistic_light_time_correction`

	:return:
		:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for a global post-newtonian :math:`\gamma` parameter.

)doc");

                    m.def(
                        "ppn_parameter_beta", &tep::ppnParameterBeta,
                        R"doc(Function for defining parameter settings for post-newtonian beta parameter.

	Factory function for creating a (linear sensitivity) parameter settings object for a global PPN :math:`\beta` parameter.
	Using the post-newtonian gamma parameter as estimatable parameter requires at least one of the following:
	 * An acceleration model depending on this parameter, such as :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.relativistic_correction`
	 * An observation model with a light-time correction depending on this parameter (none yet implemented)

	:return:
		:class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` object for a global post-newtonian :math:`\beta` parameter.

)doc");
                }

            }  // namespace parameter
        }  // namespace estimation_setup
    }  // namespace numerical_simulation
}  // namespace tudatpy
