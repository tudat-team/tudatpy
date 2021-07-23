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

#include "tudat/astro/propagators/propagateCovariance.h"

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
namespace tom = tudat::observation_models;

namespace tudat{
namespace simulation_setup{

std::shared_ptr< PodInput< > > createPodInput(
        const std::vector< std::tuple< tom::ObservableType, tom::LinkEnds, Eigen::VectorXd,
        std::vector< double >, tom::LinkEndType > >& tudatpyObservationsAndTimes,
          const int numberOfEstimatedParameters,
          const Eigen::MatrixXd inverseOfAprioriCovariance = Eigen::MatrixXd::Zero( 0, 0 ),
          const Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterDeviationEstimate =
        Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( 0, 0 ) )
{
    auto tudatObservationsAndTimes = getTudatCompatibleObservationsAndTimes( tudatpyObservationsAndTimes );
    return std::make_shared< PodInput< > >( tudatObservationsAndTimes, numberOfEstimatedParameters, inverseOfAprioriCovariance,
                                         initialParameterDeviationEstimate );
}

OrbitDeterminationManager< double, double > createOrbitDeterminationManager(
        const tss::SystemOfBodies& bodies,
        const std::shared_ptr< tep::EstimatableParameterSet< double > > parameterSet,
        const std::vector< std::pair< tom::LinkEnds, std::shared_ptr< tom::ObservationSettings > > >&
        observationSettings,
        const std::shared_ptr< tni::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< tp::PropagatorSettings< double > > propagatorSettings,
        const bool propagateOnCreation )
{
    tom::SortedObservationSettingsMap observationSettingsList;

        for( int i = 0; i < observationSettings.size( ); i++ )
        {
            observationSettingsList[ observationSettings.at( i ).second->observableType_ ][
                    observationSettings.at( i ).first ] = observationSettings.at( i ).second;
        }

    return OrbitDeterminationManager< double, double >(
                bodies, parameterSet, observationSettingsList, integratorSettings, propagatorSettings,
                propagateOnCreation );
}

} // namespace estimation_setup
} // namespace tudat



namespace tudatpy {
namespace simulation {
namespace estimation_setup {

    void expose_observation_setup(py::module &m) {

        py::enum_<tom::LinkEndType>(m, "LinkEndType")
                .value("unidentified_link_end", tom::LinkEndType::unidentified_link_end)
                .value("transmitter", tom::LinkEndType::transmitter)
                .value("reflector1", tom::LinkEndType::reflector1)
                .value("retransmitter", tom::LinkEndType::retransmitter)
                .value("reflector2", tom::LinkEndType::reflector2)
                .value("reflector3", tom::LinkEndType::reflector3)
                .value("reflector4", tom::LinkEndType::reflector4)
                .value("receiver", tom::LinkEndType::receiver)
                .value("observed_body", tom::LinkEndType::observed_body)
                .export_values();

        py::enum_<tom::ObservableType>(m, "ObservableType")
                .value("one_way_range_type", tom::ObservableType::one_way_range)
                .value("angular_position_type", tom::ObservableType::angular_position)
                .value("position_observable_type", tom::ObservableType::position_observable)
                .value("one_way_doppler_type", tom::ObservableType::one_way_doppler)
                .value("one_way_differenced_range_type", tom::ObservableType::one_way_differenced_range)
                .value("n_way_range_type", tom::ObservableType::n_way_range)
                .value("two_way_doppler_type", tom::ObservableType::two_way_doppler)
                .value("euler_angle_313_observable_type", tom::ObservableType::euler_angle_313_observable)
                .value("velocity_observable_type", tom::ObservableType::velocity_observable)
                .export_values();

        py::class_<tom::ObservationSettings,
                std::shared_ptr<tom::ObservationSettings>>(
                m, "ObservationSettings");

        m.def("one_way_range",
              &tom::oneWayRangeSettings,
              py::arg("light_time_correction_settings") = std::vector<std::shared_ptr<tom::LightTimeCorrectionSettings> >(),
              py::arg("bias_settings") = nullptr);

        m.def("angular",
              &tom::angularPositionSettings,
              py::arg("light_time_correction_settings") = std::vector<std::shared_ptr<tom::LightTimeCorrectionSettings> >(),
              py::arg("bias_settings") = nullptr);


        m.def("position",
              &tom::positionObservableSettings,
              py::arg("bias_settings") = nullptr);

        m.def("velocity",
              &tom::velocityObservableSettings,
              py::arg("bias_settings") = nullptr);

        m.def("313_euler_angle",
              &tom::eulerAngle313ObservableSettings,
              py::arg("bias_settings") = nullptr);

        m.def("one_way_open_loop_doppler",
              &tom::oneWayOpenLoopDoppler,
              py::arg("light_time_correction_settings") = std::vector<std::shared_ptr<tom::LightTimeCorrectionSettings> >(),
              py::arg("bias_settings") = nullptr,
              py::arg("transmitter_proper_time_rate_settings") = nullptr,
              py::arg("receiver_proper_time_rate_settings") = nullptr);

        m.def("two_way_open_loop_doppler",
              &tom::twoWayOpenLoopDoppler,
              py::arg("uplink_doppler_settings"),
              py::arg("downlink_doppler_settings"),
              py::arg("bias_settings") = nullptr);

        m.def("one_way_closed_loop_doppler",
              py::overload_cast<
                      const double,
                      const std::vector<std::shared_ptr<tom::LightTimeCorrectionSettings> >,
                      const std::shared_ptr<tom::ObservationBiasSettings> >(&tom::oneWayClosedLoopDoppler),
              py::arg("integration_time"),
              py::arg("light_time_correction_settings") = std::vector<std::shared_ptr<tom::LightTimeCorrectionSettings> >(),
              py::arg("bias_settings") = nullptr);

        m.def("one_way_closed_loop_doppler",
              py::overload_cast<
                      const std::function<double(const double)>,
                      const std::vector<std::shared_ptr<tom::LightTimeCorrectionSettings> >,
                      const std::shared_ptr<tom::ObservationBiasSettings> >(&tom::oneWayClosedLoopDoppler),
              py::arg("integration_time_function"),
              py::arg("light_time_correction_settings") = std::vector<std::shared_ptr<tom::LightTimeCorrectionSettings> >(),
              py::arg("bias_settings") = nullptr);

        m.def("n_way_range",
              &tom::nWayRange,
              py::arg("one_way_range_settings"),
              py::arg("bias_settings") = nullptr,
              py::arg("retransmission_times_function") = nullptr);

        py::class_<tom::LightTimeCorrectionSettings,
                std::shared_ptr<tom::LightTimeCorrectionSettings>>(
                m, "LightTimeCorrectionSettings");

        m.def("first_order_relativistic_light_time_correction",
              &tom::firstOrderRelativisticLightTimeCorrectionSettings,
              py::arg("perturbing_bodies"));

        py::class_<tom::ObservationBiasSettings,
                std::shared_ptr<tom::ObservationBiasSettings>>(
                m, "ObservationBiasSettings");

        m.def("bias",
              &tom::constantAbsoluteBias,
              py::arg("bias_value"));

        m.def("relative_bias",
              &tom::constantRelativeBias,
              py::arg("bias_value"));

        m.def("arcwise_bias",
              py::overload_cast<
                      const std::vector<double> &,
                      const std::vector<Eigen::VectorXd> &,
                      const tom::LinkEndType>(&tom::arcWiseAbsoluteBias),
              py::arg("arc_start_times"),
              py::arg("bias_values"),
              py::arg("time_link_end"));

        m.def("arcwise_bias",
              py::overload_cast<
                      const std::map<double, Eigen::VectorXd> &,
                      const tom::LinkEndType>(&tom::arcWiseAbsoluteBias),
              py::arg("bias_values_per_start_time"),
              py::arg("time_link_end"));

        m.def("arcwise_relative_bias",
              py::overload_cast<
                      const std::vector<double> &,
                      const std::vector<Eigen::VectorXd> &,
                      const tom::LinkEndType>(&tom::arcWiseRelativeBias),
              py::arg("arc_start_times"),
              py::arg("bias_values"),
              py::arg("time_link_end"));

        m.def("arcwise_relative_bias",
              py::overload_cast<
                      const std::map<double, Eigen::VectorXd> &,
                      const tom::LinkEndType>(&tom::arcWiseRelativeBias),
              py::arg("bias_values_per_start_time"),
              py::arg("time_link_end"));

        m.def("combined_bias",
              &tom::multipleObservationBiasSettings,
              py::arg("bias_list"));


        py::enum_<tom::ObservationViabilityType>(m, "ObservationViabilityType")
                .value("minimum_elevation_angle", tom::ObservationViabilityType::minimum_elevation_angle)
                .value("body_avoidance_angle", tom::ObservationViabilityType::body_avoidance_angle)
                .value("body_occultation", tom::ObservationViabilityType::body_occultation)
                .export_values();


        py::class_<tom::ObservationViabilitySettings,
                std::shared_ptr<tom::ObservationViabilitySettings>>(
                m, "ObservationViabilitySettings")
                .def(py::init<const tom::ObservationViabilityType,
                             const std::pair<std::string, std::string>,
                             const std::string,
                             const double>(),
                     py::arg("viability_type"),
                     py::arg("associated_link_end"),
                     py::arg("string_input"),
                     py::arg("double_input"));

        m.def("elevation_angle_viability",
              &tom::elevationAngleViabilitySettings,
              py::arg("link_end"),
              py::arg("elevationAngle"));


        m.def("body_avoidance_viability",
              &tom::bodyAvoidanceAngleViabilitySettings,
              py::arg("link_end"),
              py::arg("body_to_avoid"),
              py::arg("avoidance_angle"));

        m.def("body_occultation_viability",
              &tom::bodyOccultationViabilitySettings,
              py::arg("link_end"),
              py::arg("occulting_body"));

        m.def("create_observation_simulators",
              py::overload_cast<const tom::ObservationSettingsVector &, const tss::SystemOfBodies &>(
                      &tom::createObservationSimulators<double, double>),
              py::arg("observation_settings"),
              py::arg("bodies"));

    }

    void expose_estimated_parameter_setup(py::module &m) {

        py::enum_<tep::EstimatebleParametersEnum>(m, "EstimatebleParameterTypes")
                .value("arc_wise_initial_body_state_type", tep::EstimatebleParametersEnum::arc_wise_initial_body_state)
                .value("initial_body_state_type", tep::EstimatebleParametersEnum::initial_body_state)
                .value("initial_rotational_body_state_type",
                       tep::EstimatebleParametersEnum::initial_rotational_body_state)
                .value("gravitational_parameter_type", tep::EstimatebleParametersEnum::gravitational_parameter)
                .value("constant_drag_coefficient_type", tep::EstimatebleParametersEnum::constant_drag_coefficient)
                .value("radiation_pressure_coefficient_type",
                       tep::EstimatebleParametersEnum::radiation_pressure_coefficient)
                .value("arc_wise_radiation_pressure_coefficient_type",
                       tep::EstimatebleParametersEnum::arc_wise_radiation_pressure_coefficient)
                .value("spherical_harmonics_cosine_coefficient_block_type",
                       tep::EstimatebleParametersEnum::arc_wise_initial_body_state)
                .value("spherical_harmonics_sine_coefficient_block_type",
                       tep::EstimatebleParametersEnum::spherical_harmonics_sine_coefficient_block)
                .value("constant_rotation_rate_type", tep::EstimatebleParametersEnum::constant_rotation_rate)
                .value("rotation_pole_position_type", tep::EstimatebleParametersEnum::rotation_pole_position)
                .value("constant_additive_observation_bias_type",
                       tep::EstimatebleParametersEnum::constant_additive_observation_bias)
                .value("arcwise_constant_additive_observation_bias_type",
                       tep::EstimatebleParametersEnum::arcwise_constant_additive_observation_bias)
                .value("constant_relative_observation_bias_type",
                       tep::EstimatebleParametersEnum::constant_relative_observation_bias)
                .value("arcwise_constant_relative_observation_bias_type",
                       tep::EstimatebleParametersEnum::arcwise_constant_relative_observation_bias)
                .value("ppn_parameter_gamma_type", tep::EstimatebleParametersEnum::ppn_parameter_gamma)
                .value("ppn_parameter_beta_type", tep::EstimatebleParametersEnum::ppn_parameter_beta)
                .value("ground_station_position_type", tep::EstimatebleParametersEnum::ground_station_position)
                .value("equivalence_principle_lpi_violation_parameter_type",
                       tep::EstimatebleParametersEnum::equivalence_principle_lpi_violation_parameter)
                .value("empirical_acceleration_coefficients_type",
                       tep::EstimatebleParametersEnum::empirical_acceleration_coefficients)// TO EXPOSE
                .value("arc_wise_empirical_acceleration_coefficients_type",
                       tep::EstimatebleParametersEnum::arc_wise_empirical_acceleration_coefficients)// TO EXPOSE
                .value("full_degree_tidal_love_number_type",
                       tep::EstimatebleParametersEnum::full_degree_tidal_love_number)// TO EXPOSE
                .value("single_degree_variable_tidal_love_number_type",
                       tep::EstimatebleParametersEnum::single_degree_variable_tidal_love_number)
                .value("direct_dissipation_tidal_time_lag_type",
                       tep::EstimatebleParametersEnum::direct_dissipation_tidal_time_lag)
                .value("mean_moment_of_inertia_type", tep::EstimatebleParametersEnum::mean_moment_of_inertia)
                .value("arc_wise_constant_drag_coefficient_type",
                       tep::EstimatebleParametersEnum::arc_wise_constant_drag_coefficient)
                .value("periodic_spin_variation_type", tep::EstimatebleParametersEnum::periodic_spin_variation)
                .value("polar_motion_amplitude_type", tep::EstimatebleParametersEnum::polar_motion_amplitude)
                .value("core_factor_type", tep::EstimatebleParametersEnum::core_factor)
                .value("free_core_nutation_rate_type", tep::EstimatebleParametersEnum::free_core_nutation_rate)
                .value("desaturation_delta_v_values_type", tep::EstimatebleParametersEnum::desaturation_delta_v_values)
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

        m.def("initial_states",
              &tss::getInitialStateParameterSettings<double>,
              py::arg("propagator_settings"), py::arg("bodies"));

        m.def("gravitational_parameter",
              &tep::gravitationalParameter,
              py::arg("body_name"));


        m.def("spherical_harmonics_c_coefficients",
              py::overload_cast<const std::string,
                      const int,
                      const int,
                      const int,
                      const int>(&tep::sphericalHarmonicsCosineBlock),
              py::arg("body_name"),
              py::arg("minimum_degree"),
              py::arg("minimum_order"),
              py::arg("maximum_degree"),
              py::arg("maximum_order"));

        m.def("spherical_harmonics_c_coefficients",
              py::overload_cast<const std::string,
                      std::vector<std::pair<int, int> > >(&tep::sphericalHarmonicsCosineBlock),
              py::arg("body_name"),
              py::arg("block_indices"));

        m.def("spherical_harmonics_s_coefficients",
              py::overload_cast<const std::string,
                      const int,
                      const int,
                      const int,
                      const int>(&tep::sphericalHarmonicsSineBlock),
              py::arg("body_name"),
              py::arg("minimum_degree"),
              py::arg("minimum_order"),
              py::arg("maximum_degree"),
              py::arg("maximum_order"));

        m.def("spherical_harmonics_s_coefficients",
              py::overload_cast<const std::string,
                      std::vector<std::pair<int, int> > >(&tep::sphericalHarmonicsSineBlock),
              py::arg("body_name"),
              py::arg("block_indices"));


        m.def("constant_drag_coefficient",
              &tep::constantDragCoefficient,
              py::arg("body_name"));


        m.def("radiation_pressure_coefficient",
              &tep::radiationPressureCoefficient,
              py::arg("body_name"));

        m.def("arcwise_radiation_pressure_coefficient",
              &tep::arcwiseRadiationPressureCoefficient,
              py::arg("body_name"),
              py::arg("arc_initial_times"));

        m.def("arcwise_drag_coefficient",
              &tep::arcwiseDragCoefficient,
              py::arg("body_name"),
              py::arg("arc_initial_times"));

        m.def("constant_rotation_rate",
              &tep::constantRotationRate,
              py::arg("body_name"));

        m.def("rotation_pole_position",
              &tep::rotationPolePosition,
              py::arg("body_name"));

        m.def("observation_bias",
              &tep::observationBias,
              py::arg("link_ends"),
              py::arg("observable_type"));

        m.def("relative_observation_bias",
              &tep::relativeObservationBias,
              py::arg("link_ends"),
              py::arg("observable_type"));

        m.def("arcwise_observation_bias",
              &tep::arcwiseObservationBias,
              py::arg("link_ends"),
              py::arg("observable_type"),
              py::arg("arc_start_times"),
              py::arg("time_link_end"));

        m.def("arcwise_relative_observation_bias",
              &tep::arcwiseRelativeObservationBias,
              py::arg("link_ends"),
              py::arg("observable_type"),
              py::arg("arc_start_times"),
              py::arg("time_link_end"));


        m.def("ppn_parameter_gamma",
              &tep::ppnParameterGamma);

        m.def("ppn_parameter_beta",
              &tep::ppnParameterBeta);

        m.def("ground_station_position",
              &tep::groundStationPosition,
              py::arg("body_name"),
              py::arg("ground_station_name"));

        m.def("direct_tidal_dissipation_time_lag",
              py::overload_cast<const std::string &, const std::string &>(
                      &tep::directTidalDissipationLagTime),
              py::arg("body_name"),
              py::arg("deforming_body"));

        m.def("direct_tidal_dissipation_time_lag",
              py::overload_cast<const std::string &, const std::vector<std::string> &>(
                      &tep::directTidalDissipationLagTime),
              py::arg("body_name"),
              py::arg("deforming_body"));

        m.def("mean_moment_of_inertia",
              &tep::meanMomentOfInertia,
              py::arg("body_name"));

        m.def("order_invariant_k_love_number",
              py::overload_cast<const std::string &,
                      const int,
                      const std::string,
                      const bool>(&tep::orderInvariantKLoveNumber),
              py::arg("deformed_body"),
              py::arg("degree"),
              py::arg("deforming_body"),
              py::arg("use_complex_love_number") = 0);

        m.def("order_invariant_k_love_number",
              py::overload_cast<const std::string &,
                      const int,
                      const std::vector<std::string> &,
                      const bool>(&tep::orderInvariantKLoveNumber),
              py::arg("deformed_body"),
              py::arg("degree"),
              py::arg("deforming_bodies"),
              py::arg("use_complex_love_number") = 0);

        m.def("order_invariant_k_love_number",
              py::overload_cast<const std::string &,
                      const int,
                      const bool>(&tep::orderInvariantKLoveNumber),
              py::arg("deformed_body"),
              py::arg("degree"),
              py::arg("use_complex_love_number") = 0);


        m.def("order_varying_k_love_number",
              py::overload_cast<const std::string &,
                      const int,
                      const std::vector<int> &,
                      const std::string,
                      const bool>(&tep::orderVaryingKLoveNumber),
              py::arg("deformed_body"),
              py::arg("degree"),
              py::arg("orders"),
              py::arg("deforming_body"),
              py::arg("use_complex_love_number") = 0);

        m.def("order_varying_k_love_number",
              py::overload_cast<const std::string &,
                      const int,
                      const std::vector<int> &,
                      const std::vector<std::string> &,
                      const bool>(&tep::orderVaryingKLoveNumber),
              py::arg("deformed_body"),
              py::arg("degree"),
              py::arg("orders"),
              py::arg("deforming_bodies"),
              py::arg("use_complex_love_number") = 0);

        m.def("order_varying_k_love_number",
              py::overload_cast<const std::string &,
                      const int,
                      const std::vector<int> &,
                      const bool>(&tep::orderVaryingKLoveNumber),
              py::arg("deformed_body"),
              py::arg("degree"),
              py::arg("orders"),
              py::arg("use_complex_love_number") = 0);

        m.def("constant_empirical_acceleration_terms",
              &tep::constantEmpiricalAccelerationMagnitudes,
              py::arg("body"),
              py::arg("centralBody"));
    }

    void expose_estimation_setup(py::module &m) {

        //TODO: Remove variationalOnlyIntegratorSettings
        py::class_<
                tp::SingleArcVariationalEquationsSolver<double, double>,
                std::shared_ptr<tp::SingleArcVariationalEquationsSolver<double, double>>>(m,
                                                                                          "SingleArcVariationalEquationsSolver")
                .def(py::init<
                             const tudat::simulation_setup::SystemOfBodies &,
                             const std::shared_ptr<tudat::numerical_integrators::IntegratorSettings<double>>,
                             const std::shared_ptr<tp::PropagatorSettings<double>>,
                             const std::shared_ptr<tep::EstimatableParameterSet<double> >,
                             const bool,
                             const std::shared_ptr<tudat::numerical_integrators::IntegratorSettings<double> >,
                             const bool,
                             const bool,
                             const bool>(),
                     py::arg("bodies"),
                     py::arg("integrator_settings"),
                     py::arg("propagator_settings"),
                     py::arg("estimated_parameters"),
                     py::arg("integrate_equations_concurrently") = true,
                     py::arg("variational_only_integrator_settings") = std::shared_ptr<tudat::numerical_integrators::IntegratorSettings<double> >(),
                     py::arg("clear_numerical_solutions") = false,
                     py::arg("integrate_on_creation") = true,
                     py::arg("set_integrated_result") = false)
                .def("integrate_equations_of_motion_only",
                     &tp::SingleArcVariationalEquationsSolver<double, double>::integrateDynamicalEquationsOfMotionOnly,
                     py::arg("initial_states"))
                .def("integrate_full_equations",
                     &tp::SingleArcVariationalEquationsSolver<double, double>::integrateVariationalAndDynamicalEquations,
                     py::arg("initial_states"),
                     py::arg("integrate_equations_concurrently"))
                .def("get_variational_equations_solution",
                     &tp::SingleArcVariationalEquationsSolver<double, double>::getNumericalVariationalEquationsSolution)
                .def("get_state_transition_matrix_solution",
                     &tp::SingleArcVariationalEquationsSolver<double, double>::getStateTransitionMatrixSolution)
                .def_property_readonly("state_transition_matrix_history",
                                       &tp::SingleArcVariationalEquationsSolver<double, double>::getStateTransitionMatrixSolution)
                .def("get_sensitivity_matrix_solution",
                     &tp::SingleArcVariationalEquationsSolver<double, double>::getSensitivityMatrixSolution)
                .def_property_readonly("sensitivity_matrix_history",
                                       &tp::SingleArcVariationalEquationsSolver<double, double>::getSensitivityMatrixSolution)
                .def("get_equations_of_motion_solution",
                     &tp::SingleArcVariationalEquationsSolver<double, double>::getEquationsOfMotionSolution)
                .def_property_readonly("state_history",
                                       &tp::SingleArcVariationalEquationsSolver<double, double>::getEquationsOfMotionSolution)
                .def("get_dynamics_simulator",
                     &tp::SingleArcVariationalEquationsSolver<double, double>::getDynamicsSimulator)
                .def("reset_parameter_estimate",
                     &tp::SingleArcVariationalEquationsSolver<double, double>::resetParameterEstimate);

        py::class_<
                tp::CombinedStateTransitionAndSensitivityMatrixInterface,
                std::shared_ptr<tp::CombinedStateTransitionAndSensitivityMatrixInterface>>(
                m, "CombinedStateTransitionAndSensitivityMatrixInterface")
                .def("state_transition_sensitivity_at_epoch",
                     &tp::CombinedStateTransitionAndSensitivityMatrixInterface::
                     getCombinedStateTransitionAndSensitivityMatrix,
                     py::arg("time"))
                .def("full_state_transition_sensitivity_at_epoch",
                     &tp::CombinedStateTransitionAndSensitivityMatrixInterface::
                     getCombinedStateTransitionAndSensitivityMatrix,
                     py::arg("time"))
                .def_property_readonly(
                        "state_transition_size",
                        &tp::CombinedStateTransitionAndSensitivityMatrixInterface::getStateTransitionMatrixSize)
                .def_property_readonly(
                        "sensitivity_size",
                        &tp::CombinedStateTransitionAndSensitivityMatrixInterface::getSensitivityMatrixSize)
                .def_property_readonly(
                        "full_parameter_size",
                        &tp::CombinedStateTransitionAndSensitivityMatrixInterface::getFullParameterVectorSize);;

        m.def("propagate_covariance",
              py::overload_cast<const Eigen::MatrixXd &,
                      const std::shared_ptr<tp::CombinedStateTransitionAndSensitivityMatrixInterface>,
                      const std::vector<double> >(&tp::propagateCovariance),
              py::arg("initial_covariance"),
              py::arg("state_transition_interface"),
              py::arg("output_times"));

        py::class_<tep::EstimatableParameterSet<double>,
                std::shared_ptr<tep::EstimatableParameterSet<double>>>(m, "EstimatableParameterSet")
                .def_property_readonly("parameter_set_size",
                                       &tep::EstimatableParameterSet<double>::getEstimatedParameterSetSize)
                .def_property_readonly("initial_states_size",
                                       &tep::EstimatableParameterSet<double>::getInitialDynamicalStateParameterSize)
                .def_property_readonly("initial_multi_arc_states_size",
                                       &tep::EstimatableParameterSet<double>::getInitialDynamicalSingleArcStateParameterSize)
                .def_property_readonly("initial_multi_arc_states_size",
                                       &tep::EstimatableParameterSet<double>::getInitialDynamicalMultiArcStateParameterSize)
                .def_property_readonly("constraints_size",
                                       &tep::EstimatableParameterSet<double>::getConstraintSize)
                .def_property_readonly("values",
                                       &tep::EstimatableParameterSet<double>::getFullParameterValues<double>)
                .def("reset_values",
                     &tep::EstimatableParameterSet<double>::resetParameterValues<double>,
                     py::arg("new_parameter_values"))
                .def("indices_for_parameter_type",
                     &tep::EstimatableParameterSet<double>::getIndicesForParameterType,
                     py::arg("parameter_type"));
        py::class_<
                tss::EstimationConvergenceChecker,
                std::shared_ptr<tss::EstimationConvergenceChecker>>(m, "EstimationConvergenceChecker")
                .def(py::init<const unsigned int,
                             const double,
                             const double,
                             const int>(),
                     py::arg("maximum_iterations") = 5,
                     py::arg("minimum_residual_change") = 0.0,
                     py::arg("minimum_residual") = 0.0,
                     py::arg("number_of_iterations_without_improvement") = 2);

        py::class_<
                tss::PodInput<double, double>,
                std::shared_ptr<tss::PodInput<double, double>>>(m, "PodInput")
                .def(py::init<
                             const tss::PodInput<double, double>::PodInputDataType &,
                             const int,
                             const Eigen::MatrixXd,
                             const Eigen::VectorXd>(),
                     py::arg("observations_and_times"),
                     py::arg("parameter_size"),
                     py::arg("inverse_apriori_covariance") = Eigen::MatrixXd::Zero(0, 0),
                     py::arg("apriori_parameter_correction") = Eigen::VectorXd(0))
                .def("set_constant_weight",
                     &tss::PodInput<double, double>::setConstantWeightsMatrix,
                     py::arg("weight"))
                .def("set_constant_weight_per_observable",
                     &tss::PodInput<double, double>::setConstantPerObservableWeightsMatrix,
                     py::arg("weight_per_observable"))
                .def("set_constant_weight_per_observable_and_link_end",
                     &tss::PodInput<double, double>::setConstantPerObservableAndLinkEndsWeights,
                     py::arg("weight_per_observable_and_link"))
                .def("define_estimation_settings",
                     &tss::PodInput<double, double>::defineEstimationSettings,
                     py::arg("reintegrate_equations_on_first_iteration") = true,
                     py::arg("reintegrate_variational_equations") = true,
                     py::arg("save_design_matrix") = true,
                     py::arg("print_output_to_terminal") = true,
                     py::arg("save_residuals_and_parameters_per_iteration") = true,
                     py::arg("save_state_history_per_iteration") = false);

        py::class_<
                tss::PodOutput<double, double>,
                std::shared_ptr<tss::PodOutput<double, double>>>(m, "PodOutput")
                .def_property_readonly("inverse_covariance",
                                       &tss::PodOutput<double, double>::getUnnormalizedInverseCovarianceMatrix)
                .def_property_readonly("covariance",
                                       &tss::PodOutput<double, double>::getUnnormalizedCovarianceMatrix)
                .def_property_readonly("formal_errors",
                                       &tss::PodOutput<double, double>::getFormalErrorVector)
                .def_property_readonly("correlations",
                                       &tss::PodOutput<double, double>::getCorrelationMatrix)
                .def_property_readonly("residual_history",
                                       &tss::PodOutput<double, double>::getResidualHistoryMatrix)
                .def_property_readonly("parameter_history",
                                       &tss::PodOutput<double, double>::getParameterHistoryMatrix)
                .def_property_readonly("design_matrix",
                                       &tss::PodOutput<double, double>::getUnnormalizedPartialDerivatives);

        py::class_<
                tss::OrbitDeterminationManager<double, double>,
                std::shared_ptr<tss::OrbitDeterminationManager<double, double>>>(m, "OrbitDeterminationManager")
                .def(py::init<const tss::SystemOfBodies &,
                             const std::shared_ptr<tep::EstimatableParameterSet<double> >,
                             const tom::SortedObservationSettingsMap &,
                             const std::shared_ptr<tni::IntegratorSettings<double> >,
                             const std::shared_ptr<tp::PropagatorSettings<double> >,
                             const bool>(),
                     py::arg("bodies"),
                     py::arg("estimated_parameters"),
                     py::arg("observation_settings"),
                     py::arg("integrator_settings"),
                     py::arg("propagator_settings"),
                     py::arg("integrate_on_creation") = true)
                .def(py::init<const tss::SystemOfBodies &,
                             const std::shared_ptr<tep::EstimatableParameterSet<double> >,
                             const tom::ObservationSettingsMap &,
                             const std::shared_ptr<tni::IntegratorSettings<double> >,
                             const std::shared_ptr<tp::PropagatorSettings<double> >,
                             const bool>(),
                     py::arg("bodies"),
                     py::arg("estimated_parameters"),
                     py::arg("observation_settings"),
                     py::arg("integrator_settings"),
                     py::arg("propagator_settings"),
                     py::arg("integrate_on_creation") = true)
                .def_property_readonly("observation_simulators",
                                       &tss::OrbitDeterminationManager<double, double>::getObservationSimulators)
                .def_property_readonly("observation_managers",
                                       &tss::OrbitDeterminationManager<double, double>::getObservationManagers)
                .def("perform_estimation",
                     &tss::OrbitDeterminationManager<double, double>::estimateParameters,
                     py::arg("estimation_input"),
                     py::arg("convergence_checker") = std::make_shared<tss::EstimationConvergenceChecker>());

        m.def("create_orbit_determination_manager",
              &tss::createOrbitDeterminationManager,
              py::arg("bodies"),
              py::arg("estimated_parameters"),
              py::arg("observation_settings"),
              py::arg("integrator_settings"),
              py::arg("propagator_settings"),
              py::arg("integrate_on_creation") = true);

        m.def("pod_input",
              &tss::createPodInput,
              py::arg("observations_and_times"),
              py::arg("number_of_parameters"),
              py::arg("inverse_a_priori_covariance") = Eigen::MatrixXd(0, 0),
              py::arg("initial_parameter_deviation") = Eigen::VectorXd(0));


        m.def("create_parameters_to_estimate",
              &tss::createParametersToEstimate<double>,
              py::arg("parameter_settings"),
              py::arg("bodies"),
              py::arg("propagator_settings") =
                      std::shared_ptr<tp::PropagatorSettings<double> >());

        auto parameter_setup = m.def_submodule("parameter");
        expose_estimated_parameter_setup(parameter_setup);

        auto observation_setup = m.def_submodule("observation_setup");
        expose_observation_setup(observation_setup);
    }


}// namespace estimation_setup
}// namespace simulation
}// namespace tudatpy
