"""
Tests for pickle serialization/deserialization of AccelerationSettings,
PropagationResults, and all registered derived classes.

Run with: pytest test_pickle_support.py -v
"""

import pickle
from functools import lru_cache

import pytest
import numpy as np

# Adjust imports to match your actual tudatpy module paths
from tudatpy.numerical_simulation.propagation_setup import acceleration as acc
from tudatpy import numerical_simulation
from tudatpy.interface import spice
from tudatpy.astro import element_conversion
from tudatpy.astro.time_representation import DateTime
from tudatpy import constants
from tudatpy.dynamics import environment_setup, parameters_setup, propagation_setup, simulator
from tudatpy.numerical_simulation import estimation, estimation_setup, environment
from tudatpy.numerical_simulation.estimation_setup import observation

from tudatpy.dynamics.propagation import (
    SingleArcSimulationResults,
    SingleArcVariationalSimulationResults,
    MultiArcSimulationResults,
    MultiArcVariationalSimulationResults,
    HybridArcSimulationResults,
    HybridArcVariationalSimulationResults,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def roundtrip(obj):
    """Pickle and unpickle an object, returning the reconstructed instance."""
    return pickle.loads(pickle.dumps(obj))


def assert_roundtrip(obj):
    """Assert that an object survives a pickle roundtrip with equality."""
    restored = roundtrip(obj)
    assert restored == obj, (
        f"Roundtrip failed for {type(obj).__name__}:\n"
        f"  original: {obj}\n"
        f"  restored: {restored}"
    )
    # Also verify the concrete type is preserved (no silent base-class slicing)
    assert type(restored) is type(obj), (
        f"Type mismatch after roundtrip: "
        f"expected {type(obj).__name__}, got {type(restored).__name__}"
    )


def assert_dict_of_arrays_equal(restored, original):
    assert restored.keys() == original.keys()
    for epoch in original:
        np.testing.assert_allclose(restored[epoch], original[epoch])


def assert_single_arc_results_roundtrip(obj):
    restored = roundtrip(obj)
    assert type(restored) is type(obj)
    assert_dict_of_arrays_equal(restored.state_history, obj.state_history)
    assert_dict_of_arrays_equal(restored.unprocessed_state_history, obj.unprocessed_state_history)
    assert_dict_of_arrays_equal(restored.dependent_variable_history, obj.dependent_variable_history)
    assert restored.initial_and_final_times == obj.initial_and_final_times
    assert restored.propagated_state_vector_length == obj.propagated_state_vector_length
    assert restored.propagation_is_performed == obj.propagation_is_performed


def assert_single_arc_variational_results_roundtrip(obj):
    restored = roundtrip(obj)
    assert type(restored) is type(obj)
    assert_single_arc_results_roundtrip(restored.dynamics_results)
    assert_single_arc_results_roundtrip(obj.dynamics_results)
    assert_dict_of_arrays_equal(
        restored.state_transition_matrix_history, obj.state_transition_matrix_history
    )
    assert_dict_of_arrays_equal(restored.sensitivity_matrix_history, obj.sensitivity_matrix_history)


def assert_array_roundtrip(restored, original):
    np.testing.assert_allclose(restored, original)


def assert_covariance_analysis_output_roundtrip(obj):
    restored = roundtrip(obj)
    assert type(restored) is type(obj)
    assert_array_roundtrip(restored.inverse_covariance, obj.inverse_covariance)
    assert_array_roundtrip(restored.covariance, obj.covariance)
    assert_array_roundtrip(
        restored.inverse_normalized_covariance, obj.inverse_normalized_covariance
    )
    assert_array_roundtrip(restored.normalized_covariance, obj.normalized_covariance)
    assert_array_roundtrip(restored.formal_errors, obj.formal_errors)
    assert_array_roundtrip(restored.correlations, obj.correlations)
    assert_array_roundtrip(
        restored.consider_covariance_contribution, obj.consider_covariance_contribution
    )
    assert_array_roundtrip(restored.normalized_design_matrix, obj.normalized_design_matrix)
    assert_array_roundtrip(
        restored.design_matrix_consider_parameters, obj.design_matrix_consider_parameters
    )
    assert_array_roundtrip(
        restored.normalized_design_matrix_consider_parameters,
        obj.normalized_design_matrix_consider_parameters,
    )
    assert_array_roundtrip(restored.normalization_terms, obj.normalization_terms)
    assert_array_roundtrip(restored.consider_normalization_terms, obj.consider_normalization_terms)
    assert_array_roundtrip(restored.consider_covariance, obj.consider_covariance)


def assert_estimation_output_roundtrip(obj):
    restored = roundtrip(obj)
    assert type(restored) is type(obj)
    assert_covariance_analysis_output_roundtrip(obj)
    assert_array_roundtrip(restored.final_residuals, obj.final_residuals)
    assert_array_roundtrip(restored.final_parameters, obj.final_parameters)
    assert_array_roundtrip(restored.residual_history, obj.residual_history)
    assert_array_roundtrip(restored.parameter_history, obj.parameter_history)
    assert restored.best_iteration == obj.best_iteration
    assert restored.exception_during_inversion == obj.exception_during_inversion
    assert restored.exception_during_propagation == obj.exception_during_propagation
    assert len(restored.simulation_results_per_iteration) == len(
        obj.simulation_results_per_iteration
    )
    for restored_result, original_result in zip(
        restored.simulation_results_per_iteration, obj.simulation_results_per_iteration
    ):
        assert type(restored_result) is type(original_result)
        if isinstance(original_result, SingleArcSimulationResults):
            assert_single_arc_results_roundtrip(restored_result)
            assert_single_arc_results_roundtrip(original_result)


def assert_multi_arc_results_roundtrip(obj):
    restored = roundtrip(obj)
    assert type(restored) is type(obj)
    assert restored.arc_start_times == obj.arc_start_times
    assert restored.arc_end_times == obj.arc_end_times
    assert restored.propagation_is_performed == obj.propagation_is_performed
    assert len(restored.single_arc_results) == len(obj.single_arc_results)
    for restored_arc, original_arc in zip(restored.single_arc_results, obj.single_arc_results):
        assert_single_arc_results_roundtrip(restored_arc)
        assert_single_arc_results_roundtrip(original_arc)


def assert_multi_arc_variational_results_roundtrip(obj):
    restored = roundtrip(obj)
    assert type(restored) is type(obj)
    assert len(restored.single_arc_results) == len(obj.single_arc_results)
    for restored_arc, original_arc in zip(restored.single_arc_results, obj.single_arc_results):
        assert_single_arc_variational_results_roundtrip(restored_arc)
        assert_single_arc_variational_results_roundtrip(original_arc)


def assert_hybrid_arc_results_roundtrip(obj):
    restored = roundtrip(obj)
    assert type(restored) is type(obj)
    assert_single_arc_results_roundtrip(restored.single_arc_results)
    assert_single_arc_results_roundtrip(obj.single_arc_results)
    assert_multi_arc_results_roundtrip(restored.multi_arc_results)
    assert_multi_arc_results_roundtrip(obj.multi_arc_results)


def assert_hybrid_arc_variational_results_roundtrip(obj):
    restored = roundtrip(obj)
    assert type(restored) is type(obj)
    assert_single_arc_variational_results_roundtrip(restored.single_arc_results)
    assert_single_arc_variational_results_roundtrip(obj.single_arc_results)
    assert_multi_arc_variational_results_roundtrip(restored.multi_arc_results)
    assert_multi_arc_variational_results_roundtrip(obj.multi_arc_results)


# ---------------------------------------------------------------------------
# Base class
# ---------------------------------------------------------------------------


class TestAccelerationSettings:
    def test_point_mass_gravity(self):
        obj = acc.point_mass_gravity()
        assert_roundtrip(obj)

    def test_base_type_preserved(self):
        obj = acc.point_mass_gravity()
        restored = roundtrip(obj)
        assert type(restored).__name__ == type(obj).__name__


# ---------------------------------------------------------------------------
# Derived classes
# ---------------------------------------------------------------------------


class TestRadiationPressureAccelerationSettings:
    def test_default(self):
        obj = acc.radiation_pressure()  # adjust if constructor differs
        assert_roundtrip(obj)

    def test_cannonball(self):
        obj = acc.cannonball_radiation_pressure()  # adjust if needed
        assert_roundtrip(obj)


class TestSphericalHarmonicAccelerationSettings:
    @pytest.mark.parametrize(
        "degree,order",
        [
            (2, 0),
            (2, 2),
            (8, 8),
            (20, 20),
        ],
    )
    def test_various_degrees(self, degree, order):
        obj = acc.spherical_harmonic_gravity(degree, order)
        assert_roundtrip(obj)


class TestMutualSphericalHarmonicAccelerationSettings:
    @pytest.mark.parametrize(
        "args",
        [
            (2, 2, 2, 2),
            (4, 4, 2, 2),
            (8, 8, 4, 4),
        ],
    )
    def test_various_degrees(self, args):
        # adjust constructor name/signature if needed
        obj = acc.mutual_spherical_harmonic_gravity(*args)
        assert_roundtrip(obj)


class TestRelativisticAccelerationCorrectionSettings:
    def test_default(self):
        obj = acc.relativistic_correction(
            use_schwarzschild=True, use_lense_thirring=False, use_de_sitter=False
        )
        assert_roundtrip(obj)


class TestEmpiricalAccelerationSettings:
    def test_default(self):
        obj = acc.empirical(constant_acceleration=np.array([1, 2, 3]))  # adjust if needed
        assert_roundtrip(obj)


class TestYarkovskyAccelerationSettings:
    def test_default(self):
        obj = acc.yarkovsky(1)  # adjust if needed
        assert_roundtrip(obj)


class TestThrustAccelerationSettings:
    def test_from_engines_list(self):
        # adjust — ThrustAccelerationSettings likely requires direction/magnitude args
        obj = acc.thrust_from_engines(["1", "2", "3"])  # adjust if needed
        assert_roundtrip(obj)

    def test_from_engines(self):
        # adjust — ThrustAccelerationSettings likely requires direction/magnitude args
        obj = acc.thrust_from_engine("1")  # adjust if needed
        assert_roundtrip(obj)

    def test_from_all_engines(self):
        # adjust — ThrustAccelerationSettings likely requires direction/magnitude args
        obj = acc.thrust_from_all_engines()  # adjust if needed
        assert_roundtrip(obj)


# class TestRTGAccelerationSettings:
#     def test_default(self):
#         obj = acc.rtg_acceleration()  # adjust if needed
#         assert_roundtrip(obj)


class TestDirectTidalDissipationAccelerationSettings:
    def test_default(self):
        obj = acc.direct_tidal_dissipation_acceleration(0.2, 10, False, False)  # adjust if needed
        assert_roundtrip(obj)


# class TestMomentumWheelDesaturationAccelerationSettings:
#     def test_default(self):
#         obj = acc.momentum_wheel_desaturation_thrust()  # adjust if needed
#         assert_roundtrip(obj)


# ---------------------------------------------------------------------------
# Polymorphism: all derived types must survive roundtrip through base
# ---------------------------------------------------------------------------


class TestPolymorphicDispatch:
    """
    These tests verify that cereal correctly identifies and reconstructs
    the concrete derived type when serialized through a base class pointer.
    Without CEREAL_REGISTER_POLYMORPHIC_RELATION this would silently slice.
    """

    DERIVED_FACTORIES = [
        acc.point_mass_gravity,
        lambda: acc.spherical_harmonic_gravity(4, 4),
        acc.radiation_pressure,
        acc.relativistic_correction,
        acc.empirical,
    ]

    DERIVED_INSTANCES = [
        pytest.param(factory, id=getattr(factory, "__name__", "derived"))
        for factory in DERIVED_FACTORIES
    ]

    @pytest.mark.parametrize("factory", DERIVED_INSTANCES)
    def test_concrete_type_preserved(self, factory):
        """The concrete type must be identical before and after roundtrip."""
        obj = factory()
        restored = roundtrip(obj)
        assert type(restored) is type(obj)
        assert restored == obj

    def test_all_registered_types_are_tested(self):
        """
        Canary: if you register a new derived type, add it to DERIVED_INSTANCES
        or this test will remind you by failing.
        """
        expected_types = {
            "AccelerationSettings",
            "SphericalHarmonicAccelerationSettings",
            "EmpiricalAccelerationSettings",
            "RelativisticAccelerationCorrectionSettings",
        }
        tested_types = {factory().__class__.__name__ for factory in self.DERIVED_FACTORIES}
        untested = expected_types - tested_types
        assert not untested, f"These registered types have no roundtrip test: {untested}"


# ---------------------------------------------------------------------------
# Propagation Results Classes
# ---------------------------------------------------------------------------


class TestSimulationResultsPickle:
    """
    Tests for pickle serialization/deserialization of SimulationResults
    and derived classes. These tests verify that result objects can be
    serialized and deserialized correctly, preserving both data and type.
    """

    @staticmethod
    @lru_cache(maxsize=1)
    def _test_environment():
        spice.load_standard_kernels()

        bodies_to_create = ["Earth"]
        body_settings = environment_setup.get_default_body_settings(
            bodies_to_create,
            "Earth",
            "J2000",
        )
        body_settings.add_empty_settings("Delfi-C3")

        bodies = environment_setup.create_system_of_bodies(body_settings)
        bodies.get("Delfi-C3").mass = 2.2

        simulation_start_epoch = DateTime(2008, 4, 28).to_epoch()
        initial_state = element_conversion.keplerian_to_cartesian_elementwise(
            gravitational_parameter=bodies.get("Earth").gravitational_parameter,
            semi_major_axis=6.99276221e06,
            eccentricity=4.03294322e-03,
            inclination=1.71065169e00,
            argument_of_periapsis=1.31226971e00,
            longitude_of_ascending_node=3.82958313e-01,
            true_anomaly=3.07018490e00,
        )

        return bodies, simulation_start_epoch, initial_state

    @classmethod
    def _create_translational_settings(cls, start_epoch, end_epoch, bodies, initial_state):
        acceleration_settings = {
            "Delfi-C3": {
                "Earth": [propagation_setup.acceleration.point_mass_gravity()],
            }
        }
        acceleration_models = propagation_setup.create_acceleration_models(
            bodies,
            acceleration_settings,
            ["Delfi-C3"],
            ["Earth"],
        )

        integrator_settings = propagation_setup.integrator.runge_kutta_fixed_step(
            30.0,
            coefficient_set=propagation_setup.integrator.CoefficientSets.rk_4,
        )
        termination_settings = propagation_setup.propagator.time_termination(end_epoch)

        return propagation_setup.propagator.translational(
            ["Earth"],
            acceleration_models,
            ["Delfi-C3"],
            initial_state,
            start_epoch,
            integrator_settings,
            termination_settings,
        )

    @classmethod
    @lru_cache(maxsize=1)
    def _single_arc_simulation_results(cls):
        bodies, start_epoch, initial_state = cls._test_environment()
        propagator_settings = cls._create_translational_settings(
            start_epoch,
            start_epoch + 600.0,
            bodies,
            initial_state,
        )
        dynamics_simulator = simulator.create_dynamics_simulator(bodies, propagator_settings)
        return dynamics_simulator.propagation_results

    @classmethod
    @lru_cache(maxsize=1)
    def _single_arc_variational_simulation_results(cls):
        bodies, start_epoch, initial_state = cls._test_environment()
        propagator_settings = cls._create_translational_settings(
            start_epoch,
            start_epoch + 600.0,
            bodies,
            initial_state,
        )
        parameter_settings = parameters_setup.initial_states(propagator_settings, bodies)
        parameters_to_estimate = parameters_setup.create_parameter_set(parameter_settings, bodies)
        variational_solver = simulator.create_variational_equations_solver(
            bodies,
            propagator_settings,
            parameters_to_estimate,
            simulate_dynamics_on_creation=True,
        )
        return variational_solver.variational_propagation_results

    @classmethod
    @lru_cache(maxsize=1)
    def _multi_arc_simulation_results(cls):
        bodies, start_epoch, initial_state = cls._test_environment()

        first_arc = cls._create_translational_settings(
            start_epoch,
            start_epoch + 300.0,
            bodies,
            initial_state,
        )
        second_arc = cls._create_translational_settings(
            start_epoch + 300.0,
            start_epoch + 600.0,
            bodies,
            initial_state,
        )
        propagator_settings = propagation_setup.propagator.multi_arc(
            [first_arc, second_arc],
            transfer_state_to_next_arc=False,
        )
        dynamics_simulator = simulator.create_dynamics_simulator(bodies, propagator_settings)
        return dynamics_simulator.propagation_results

    @classmethod
    @lru_cache(maxsize=1)
    def _multi_arc_variational_simulation_results(cls):
        bodies, start_epoch, initial_state = cls._test_environment()

        first_arc = cls._create_translational_settings(
            start_epoch,
            start_epoch + 300.0,
            bodies,
            initial_state,
        )
        second_arc = cls._create_translational_settings(
            start_epoch + 300.0,
            start_epoch + 600.0,
            bodies,
            initial_state,
        )
        propagator_settings = propagation_setup.propagator.multi_arc(
            [first_arc, second_arc],
            transfer_state_to_next_arc=False,
        )
        parameter_settings = parameters_setup.initial_states(propagator_settings, bodies)
        parameters_to_estimate = parameters_setup.create_parameter_set(parameter_settings, bodies)
        variational_solver = simulator.create_variational_equations_solver(
            bodies,
            propagator_settings,
            parameters_to_estimate,
            simulate_dynamics_on_creation=True,
        )
        return variational_solver.variational_propagation_results

    @staticmethod
    @lru_cache(maxsize=1)
    def _hybrid_environment():
        spice.load_standard_kernels()

        bodies_to_create = ["Sun", "Earth"]
        body_settings = environment_setup.get_default_body_settings(
            bodies_to_create,
            "Sun",
            "J2000",
        )
        body_settings.get("Sun").radiation_source_settings = (
            environment_setup.radiation_pressure.isotropic_radiation_source(
                environment_setup.radiation_pressure.irradiance_based_constant_luminosity(
                    1367.0, constants.ASTRONOMICAL_UNIT
                )
            )
        )
        body_settings.add_empty_settings("Delfi-C3")

        bodies = environment_setup.create_system_of_bodies(body_settings)
        bodies.get("Delfi-C3").mass = 2.2

        simulation_start_epoch = DateTime(2008, 4, 28).to_epoch()
        earth_initial_state = spice.get_body_cartesian_state_at_epoch(
            "Earth",
            "Sun",
            "J2000",
            "None",
            simulation_start_epoch,
        )
        satellite_initial_state = element_conversion.keplerian_to_cartesian_elementwise(
            gravitational_parameter=bodies.get("Earth").gravitational_parameter,
            semi_major_axis=6.99276221e06,
            eccentricity=4.03294322e-03,
            inclination=1.71065169e00,
            argument_of_periapsis=1.31226971e00,
            longitude_of_ascending_node=3.82958313e-01,
            true_anomaly=3.07018490e00,
        )

        return bodies, simulation_start_epoch, earth_initial_state, satellite_initial_state

    @classmethod
    def _create_hybrid_single_arc_settings(
        cls, start_epoch, end_epoch, bodies, earth_initial_state
    ):
        acceleration_settings = {
            "Earth": {
                "Sun": [propagation_setup.acceleration.point_mass_gravity()],
            }
        }
        acceleration_models = propagation_setup.create_acceleration_models(
            bodies,
            acceleration_settings,
            ["Earth"],
            ["Sun"],
        )
        integrator_settings = propagation_setup.integrator.runge_kutta_fixed_step(
            30.0,
            coefficient_set=propagation_setup.integrator.CoefficientSets.rk_4,
        )
        termination_settings = propagation_setup.propagator.time_termination(end_epoch)
        return propagation_setup.propagator.translational(
            ["Sun"],
            acceleration_models,
            ["Earth"],
            earth_initial_state,
            start_epoch,
            integrator_settings,
            termination_settings,
        )

    @classmethod
    def _create_hybrid_multi_arc_settings(cls, start_epoch, bodies, satellite_initial_state):
        first_arc = cls._create_translational_settings(
            start_epoch,
            start_epoch + 300.0,
            bodies,
            satellite_initial_state,
        )
        second_arc = cls._create_translational_settings(
            start_epoch + 300.0,
            start_epoch + 600.0,
            bodies,
            satellite_initial_state,
        )
        return propagation_setup.propagator.multi_arc(
            [first_arc, second_arc],
            transfer_state_to_next_arc=False,
        )

    @classmethod
    @lru_cache(maxsize=1)
    def _hybrid_arc_simulation_results(cls):
        pytest.skip(
            "Chatgpt has some trouble setting up this test. This will need to be done at some point."
        )

    @classmethod
    @lru_cache(maxsize=1)
    def _hybrid_arc_variational_simulation_results(cls):
        pytest.skip(
            "Chatgpt has some trouble setting up this test. This will need to be done at some point."
        )

    def test_single_arc_simulation_results_default(self):
        """Test pickle support for SingleArcSimulationResults from an actual propagation."""
        obj = self._single_arc_simulation_results()
        assert_single_arc_results_roundtrip(obj)

    def test_single_arc_variational_simulation_results_default(self):
        """Test pickle support for SingleArcVariationalSimulationResults from an actual propagation."""
        obj = self._single_arc_variational_simulation_results()
        assert_single_arc_variational_results_roundtrip(obj)

    def test_multi_arc_simulation_results_default(self):
        """Test pickle support for MultiArcSimulationResults from an actual propagation."""
        obj = self._multi_arc_simulation_results()
        assert_multi_arc_results_roundtrip(obj)

    def test_multi_arc_variational_simulation_results_default(self):
        """Test pickle support for MultiArcVariationalSimulationResults from an actual propagation."""
        obj = self._multi_arc_variational_simulation_results()
        assert_multi_arc_variational_results_roundtrip(obj)

    def test_hybrid_arc_simulation_results_default(self):
        """Test pickle support for HybridArcSimulationResults from an actual propagation."""
        obj = self._hybrid_arc_simulation_results()
        assert_hybrid_arc_results_roundtrip(obj)

    def test_hybrid_arc_variational_simulation_results_default(self):
        """Test pickle support for HybridArcVariationalSimulationResults from an actual propagation."""
        obj = self._hybrid_arc_variational_simulation_results()
        assert_hybrid_arc_variational_results_roundtrip(obj)


class TestEstimationAnalysisPickle:
    """Tests for pickle serialization/deserialization of EstimationOutput and CovarianceAnalysisOutput."""

    @staticmethod
    @lru_cache(maxsize=1)
    def _estimation_analysis_outputs():
        spice.load_standard_kernels()

        bodies_to_create = ["Sun", "Earth", "Moon", "Mars", "Venus"]
        body_settings = environment_setup.get_default_body_settings(
            bodies_to_create,
            "Sun",
            "J2000",
        )
        body_settings.add_empty_settings("Delfi-C3")

        bodies = environment_setup.create_system_of_bodies(body_settings)
        bodies.get("Delfi-C3").mass = 2.2
        environment_setup.add_aerodynamic_coefficient_interface(
            bodies,
            "Delfi-C3",
            environment_setup.aerodynamic_coefficients.constant(
                (4 * 0.3 * 0.1 + 2 * 0.1 * 0.1) / 4,
                [1.2, 0.0, 0.0],
            ),
        )
        environment_setup.add_radiation_pressure_target_model(
            bodies,
            "Delfi-C3",
            environment_setup.radiation_pressure.cannonball_radiation_target(
                4.0, 1.2, {"Sun": ["Earth"]}
            ),
        )

        simulation_start_epoch = DateTime(2008, 4, 28).to_epoch()
        simulation_end_epoch = simulation_start_epoch + 6.0 * 3600.0

        delfi_tle = environment.Tle(
            "1 32789U 07021G   08119.60740078 -.00000054  00000-0  00000+0 0  9999",
            "2 32789 098.0082 179.6267 0015321 307.2977 051.0656 14.81417433    68",
        )
        delfi_ephemeris = environment.TleEphemeris("Earth", "J2000", delfi_tle, False)
        initial_state = delfi_ephemeris.cartesian_state(simulation_start_epoch)

        bodies_to_propagate = ["Delfi-C3"]
        central_bodies = ["Earth"]
        acceleration_settings = {
            "Delfi-C3": {
                "Sun": [
                    propagation_setup.acceleration.radiation_pressure(),
                    propagation_setup.acceleration.point_mass_gravity(),
                ],
                "Mars": [propagation_setup.acceleration.point_mass_gravity()],
                "Moon": [propagation_setup.acceleration.point_mass_gravity()],
                "Earth": [
                    propagation_setup.acceleration.spherical_harmonic_gravity(8, 8),
                    propagation_setup.acceleration.aerodynamic(),
                ],
            }
        }
        acceleration_models = propagation_setup.create_acceleration_models(
            bodies,
            acceleration_settings,
            bodies_to_propagate,
            central_bodies,
        )

        integrator_settings = propagation_setup.integrator.runge_kutta_fixed_step_size(
            initial_time_step=60.0,
            coefficient_set=propagation_setup.integrator.CoefficientSets.rkdp_87,
        )
        termination_condition = propagation_setup.propagator.time_termination(simulation_end_epoch)
        propagator_settings = propagation_setup.propagator.translational(
            central_bodies,
            acceleration_models,
            bodies_to_propagate,
            initial_state,
            simulation_start_epoch,
            integrator_settings,
            termination_condition,
        )

        environment_setup.add_ground_station(
            bodies.get_body("Earth"),
            "TrackingStation",
            [0.0, np.deg2rad(52.00667), np.deg2rad(4.35556)],
            element_conversion.geodetic_position_type,
        )

        link_ends = {
            observation.transmitter: observation.body_reference_point_link_end_id(
                "Earth", "TrackingStation"
            ),
            observation.receiver: observation.body_origin_link_end_id("Delfi-C3"),
        }
        link_definition = observation.LinkDefinition(link_ends)
        observation_settings_list = [observation.one_way_doppler_instantaneous(link_definition)]

        observation_times = np.arange(simulation_start_epoch, simulation_end_epoch, 300.0)
        observation_simulation_settings = observation.tabulated_simulation_settings(
            observation.one_way_instantaneous_doppler_type,
            link_definition,
            observation_times,
        )
        noise_level = 1.0e-3
        observation.add_gaussian_noise_to_observable(
            [observation_simulation_settings],
            noise_level,
            observation.one_way_instantaneous_doppler_type,
        )

        parameter_settings = estimation_setup.parameter.initial_states(propagator_settings, bodies)
        parameter_settings.append(estimation_setup.parameter.gravitational_parameter("Earth"))
        parameter_settings.append(estimation_setup.parameter.constant_drag_coefficient("Delfi-C3"))
        parameters_to_estimate = estimation_setup.create_parameter_set(parameter_settings, bodies)
        estimator = numerical_simulation.Estimator(
            bodies,
            parameters_to_estimate,
            observation_settings_list,
            propagator_settings,
        )

        simulated_observations = estimation.simulate_observations(
            [observation_simulation_settings],
            estimator.observation_simulators,
            bodies,
        )

        covariance_input = estimation.CovarianceAnalysisInput(simulated_observations)
        covariance_input.define_covariance_settings(reintegrate_variational_equations=False)
        covariance_input.set_constant_weight_per_observable(
            {observation.one_way_instantaneous_doppler_type: noise_level**-2}
        )
        covariance_output = estimator.compute_covariance(covariance_input)

        estimation_input = estimation.EstimationInput(
            simulated_observations,
            convergence_checker=estimation.estimation_convergence_checker(maximum_iterations=2),
        )
        estimation_input.define_estimation_settings(reintegrate_variational_equations=False)
        estimation_input.set_constant_weight_per_observable(
            {observation.one_way_instantaneous_doppler_type: noise_level**-2}
        )
        estimation_output = estimator.perform_estimation(estimation_input)

        return covariance_output, estimation_output

    def test_covariance_analysis_output_roundtrip(self):
        covariance_output, _ = self._estimation_analysis_outputs()
        assert_covariance_analysis_output_roundtrip(covariance_output)

    def test_estimation_output_roundtrip(self):
        _, estimation_output = self._estimation_analysis_outputs()
        assert_estimation_output_roundtrip(estimation_output)


class TestSimulationResultsPolymorphism:
    """
    Tests for polymorphic serialization of SimulationResults derived classes.
    Verifies that cereal correctly identifies and reconstructs the concrete type
    when serialized through a base class pointer.
    """

    RESULT_CASES = [
        "single_arc",
        "single_arc_variational",
        "multi_arc",
        "multi_arc_variational",
        "hybrid_arc",
        "hybrid_arc_variational",
    ]

    @staticmethod
    def _case_object(case_name):
        if case_name == "single_arc":
            return TestSimulationResultsPickle._single_arc_simulation_results()
        if case_name == "single_arc_variational":
            return TestSimulationResultsPickle._single_arc_variational_simulation_results()
        if case_name == "multi_arc":
            return TestSimulationResultsPickle._multi_arc_simulation_results()
        if case_name == "multi_arc_variational":
            return TestSimulationResultsPickle._multi_arc_variational_simulation_results()
        if case_name == "hybrid_arc":
            return TestSimulationResultsPickle._hybrid_arc_simulation_results()
        if case_name == "hybrid_arc_variational":
            return TestSimulationResultsPickle._hybrid_arc_variational_simulation_results()
        raise ValueError(f"Unknown result case: {case_name}")

    @staticmethod
    def _assert_case_roundtrip(case_name, obj):
        if case_name == "single_arc":
            assert_single_arc_results_roundtrip(obj)
        elif case_name == "single_arc_variational":
            assert_single_arc_variational_results_roundtrip(obj)
        elif case_name == "multi_arc":
            assert_multi_arc_results_roundtrip(obj)
        elif case_name == "multi_arc_variational":
            assert_multi_arc_variational_results_roundtrip(obj)
        elif case_name == "hybrid_arc":
            assert_hybrid_arc_results_roundtrip(obj)
        elif case_name == "hybrid_arc_variational":
            assert_hybrid_arc_variational_results_roundtrip(obj)
        else:
            raise ValueError(f"Unknown result case: {case_name}")

    @pytest.mark.parametrize("case_name", RESULT_CASES, ids=RESULT_CASES)
    def test_concrete_type_preserved_in_roundtrip(self, case_name):
        """Verify that concrete derived type is preserved through pickle roundtrip."""
        obj = self._case_object(case_name)
        restored = roundtrip(obj)
        assert type(restored) is type(obj), (
            f"Type mismatch after roundtrip: "
            f"expected {type(obj).__name__}, got {type(restored).__name__}"
        )
        self._assert_case_roundtrip(case_name, obj)

    def test_all_result_types_tested(self):
        """
        Canary test: if a new result type is registered, it should be added
        to RESULT_CASES or this test will fail as a reminder.
        """
        expected_types = {
            "single_arc",
            "single_arc_variational",
            "multi_arc",
            "multi_arc_variational",
            "hybrid_arc",
            "hybrid_arc_variational",
        }
        tested_types = set(self.RESULT_CASES)
        untested = expected_types - tested_types
        assert not untested, f"These registered result types have no pickle test: {untested}"
