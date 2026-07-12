"""
Comprehensive pickle roundtrip tests for all serializable classes.

Tests that every Python-exposed serializable class survives a pickle
roundtrip with equality preserved and concrete type preserved.

Run with: pytest test_serialization_pickle.py -v
"""

import pickle
import tempfile
import os

import pytest
import numpy as np

from tudatpy.astro.time_representation import Time
from tudatpy.estimation.observable_models_setup import links
from tudatpy.estimation.observations_setup import ancillary_settings
from tudatpy.estimation.observations_setup import (
    observations_dependent_variables as obs_dep_var,
)
from tudatpy.estimation.observations import SingleObservationSet, ObservationCollection
from tudatpy.math import root_finders
from tudatpy.dynamics.propagation_setup import dependent_variable as dep_var
from tudatpy.dynamics.propagation_setup import torque
from tudatpy.dynamics.propagation_setup import propagator
from tudatpy.dynamics.propagation_setup import acceleration as acc
from tudatpy.dynamics.environment_setup.aerodynamic_coefficients import (
    AerodynamicsReferenceFrames,
)
from tudatpy.dynamics.environment_setup.gravity_field_variation import (
    BodyDeformationTypes,
)
from tudatpy.dynamics.environment_setup import gravity_field_variation
from tudatpy.estimation.observable_models_setup import model_settings

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def roundtrip(obj):
    """Pickle and unpickle an object, returning the reconstructed instance."""
    return pickle.loads(pickle.dumps(obj))


def assert_roundtrip(obj):
    """Assert that an object survives a pickle roundtrip with equality and type preservation."""
    restored = roundtrip(obj)
    assert restored == obj, (
        f"Roundtrip failed for {type(obj).__name__}:\n"
        f"  original: {obj}\n"
        f"  restored: {restored}"
    )
    assert type(restored) is type(obj), (
        f"Type mismatch after roundtrip: "
        f"expected {type(obj).__name__}, got {type(restored).__name__}"
    )


def assert_roundtrip_fails(obj):
    """Assert that pickle roundtrip is expected to fail (e.g. for custom/lambda types)."""
    with pytest.raises(Exception):
        pickle.loads(pickle.dumps(obj))


# ===========================================================================
# Acceleration Settings
# ===========================================================================


class TestAccelerationSettingsPickle:
    """Pickle roundtrip tests for AccelerationSettings and all derived types."""

    def test_point_mass_gravity(self):
        assert_roundtrip(acc.point_mass_gravity())

    def test_radiation_pressure(self):
        assert_roundtrip(acc.radiation_pressure())

    def test_cannonball_radiation_pressure(self):
        assert_roundtrip(acc.cannonball_radiation_pressure())

    def test_spherical_harmonic_gravity(self):
        assert_roundtrip(acc.spherical_harmonic_gravity(4, 4))

    def test_mutual_spherical_harmonic_gravity(self):
        assert_roundtrip(acc.mutual_spherical_harmonic_gravity(2, 2, 2, 2))

    def test_relativistic_correction(self):
        assert_roundtrip(
            acc.relativistic_correction(
                use_schwarzschild=True, use_lense_thirring=False, use_de_sitter=False
            )
        )

    def test_empirical(self):
        assert_roundtrip(acc.empirical(constant_acceleration=np.array([1.0, 2.0, 3.0])))

    def test_yarkovsky(self):
        assert_roundtrip(acc.yarkovsky(0.1))

    def test_thrust_from_engines(self):
        assert_roundtrip(acc.thrust_from_engines(["engine1", "engine2"]))

    def test_thrust_from_engine(self):
        assert_roundtrip(acc.thrust_from_engine("engine1"))

    def test_thrust_from_all_engines(self):
        assert_roundtrip(acc.thrust_from_all_engines())

    def test_rtg(self):
        assert_roundtrip(acc.rtg(np.array([0.1, 0.0, 0.0]), 0.1, 0.0))

    def test_direct_tidal_dissipation(self):
        assert_roundtrip(acc.direct_tidal_dissipation_acceleration(0.2, 10, False, False))

    def test_momentum_wheel_desaturation_expected_failure(self):
        """MomentumWheelDesaturationAccelerationSettings has no Python factory or constructor."""
        pass

    def test_custom_acceleration_expected_failure(self):
        """CustomAccelerationSettings contains a lambda and is expected to fail roundtrip."""
        obj = acc.custom_acceleration(lambda t: np.array([1.0, 2.0, 3.0]))
        assert_roundtrip_fails(obj)


# ===========================================================================
# Torque Settings
# ===========================================================================


class TestTorqueSettingsPickle:
    """Pickle roundtrip tests for TorqueSettings and derived types."""

    def test_aerodynamic_torque(self):
        assert_roundtrip(torque.aerodynamic())

    def test_radiation_pressure_torque(self):
        assert_roundtrip(torque.radiation_pressure_torque())

    def test_second_degree_gravitational_torque(self):
        assert_roundtrip(torque.second_degree_gravitational())

    def test_spherical_harmonic_gravitational_torque(self):
        assert_roundtrip(torque.spherical_harmonic_gravitational(2, 2))

    def test_custom_torque_expected_failure(self):
        """CustomTorqueSettings contains a lambda and is expected to fail roundtrip."""
        obj = torque.custom_torque(lambda t: np.array([1.0, 0.0, 0.0]))
        assert_roundtrip_fails(obj)


# ===========================================================================
# Propagation Termination Settings
# ===========================================================================


class TestPropagationTerminationSettingsPickle:
    """Pickle roundtrip tests for PropagationTerminationSettings and derived types."""

    def test_time_termination(self):
        assert_roundtrip(propagator.time_termination(1000.0))

    def test_cpu_time_termination(self):
        assert_roundtrip(propagator.cpu_time_termination(60.0))

    def test_dependent_variable_termination(self):
        assert_roundtrip(
            propagator.dependent_variable_termination(
                dep_var.altitude("Earth", "Earth"), 100.0, False
            )
        )

    def test_hybrid_termination(self):
        assert_roundtrip(
            propagator.hybrid_termination(
                [
                    propagator.time_termination(500.0),
                    propagator.cpu_time_termination(30.0),
                ],
                fulfill_single_condition=True,
            )
        )

    def test_non_sequential_termination(self):
        assert_roundtrip(
            propagator.non_sequential_termination(
                propagator.time_termination(1000.0),
                propagator.time_termination(-1000.0),
            )
        )

    def test_custom_termination_expected_failure(self):
        """PropagationCustomTerminationSettings contains a lambda and is expected to fail."""
        obj = propagator.custom_termination(lambda t: True)
        assert_roundtrip_fails(obj)


# ===========================================================================
# Dependent Variable Settings (Propagation)
# ===========================================================================


class TestDependentVariableSettingsPickle:
    """Pickle roundtrip tests for VariableSettings and all derived save-settings types."""

    def test_variable_settings_base(self):
        assert_roundtrip(dep_var.altitude("Earth", "Earth"))

    def test_single_dependent_variable(self):
        assert_roundtrip(dep_var.altitude("Earth", "Earth"))

    def test_single_acceleration(self):
        assert_roundtrip(dep_var.single_acceleration(acc.point_mass_gravity_type, "Earth", "Moon"))

    def test_spherical_harmonic_terms_acceleration(self):
        assert_roundtrip(dep_var.spherical_harmonic_terms_acceleration("Earth", "Moon", [(2, 2)]))

    def test_single_torque(self):
        assert_roundtrip(dep_var.single_torque(torque.aerodynamic_type, "Earth", "Moon"))

    def test_intermediate_aerodynamic_rotation_matrix(self):
        assert_roundtrip(
            dep_var.intermediate_aerodynamic_rotation_matrix_variable(
                "Earth",
                AerodynamicsReferenceFrames.inertial_frame,
                AerodynamicsReferenceFrames.body_frame,
            )
        )

    def test_body_aerodynamic_orientation_angle(self):
        assert_roundtrip(dep_var.angle_of_attack("Earth", "Sun"))

    def test_local_wind_velocity(self):
        assert_roundtrip(dep_var.local_wind_velocity("Earth", "Earth"))

    def test_control_surface_deflection(self):
        assert_roundtrip(dep_var.control_surface_deflection("Earth", "aileron"))

    def test_single_gravity_field_variation_acceleration(self):
        assert_roundtrip(
            dep_var.single_gravity_field_variation_acceleration(
                "Earth", "Moon", BodyDeformationTypes.basic_solid_body
            )
        )

    def test_single_per_term_gravity_field_variation_acceleration(self):
        assert_roundtrip(
            dep_var.single_per_term_gravity_field_variation_acceleration(
                "Earth", "Moon", [(2, 2)], BodyDeformationTypes.basic_solid_body
            )
        )

    def test_acceleration_partial_wrt_state(self):
        assert_roundtrip(
            dep_var.acceleration_partial_wrt_body_translational_state(
                acc.point_mass_gravity_type, "Earth", "Moon", "Earth"
            )
        )

    def test_total_acceleration_partial_wrt_state(self):
        assert_roundtrip(
            dep_var.total_acceleration_partial_wrt_body_translational_state("Earth", "Earth")
        )

    def test_minimum_constellation_distance(self):
        assert_roundtrip(dep_var.minimum_body_distance("Earth", ["Moon", "Mars"]))

    def test_minimum_constellation_station_distance(self):
        assert_roundtrip(
            dep_var.minimum_visible_station_body_distances(
                "Earth", "station1", ["Moon", "Mars"], 10.0
            )
        )

    def test_cross_section(self):
        assert_roundtrip(
            dep_var.actual_cross_section("Earth", "Sun", "cannon_ball_radiation_pressure")
        )

    def test_illuminated_panel_fraction(self):
        assert_roundtrip(dep_var.illuminated_panel_fraction("Earth", "Sun"))

    def test_total_gravity_field_variation_acceleration(self):
        assert_roundtrip(dep_var.total_gravity_field_variation_acceleration("Earth", "Moon"))

    def test_custom_dependent_variable_expected_failure(self):
        """CustomDependentVariableSaveSettings contains a lambda and is expected to fail."""
        obj = dep_var.custom_dependent_variable(lambda: np.zeros(3), 3)
        assert_roundtrip_fails(obj)


# ===========================================================================
# Observation Dependent Variable Settings
# ===========================================================================


class TestObservationDependentVariableSettingsPickle:
    """Pickle roundtrip tests for ObservationDependentVariableSettings and derived types."""

    def test_base_observation_dependent_variable(self):
        assert_roundtrip(obs_dep_var.elevation_angle_dependent_variable(links.receiver))

    def test_station_angle_observation_dependent_variable(self):
        assert_roundtrip(obs_dep_var.azimuth_angle_dependent_variable(links.receiver))

    def test_interlink_observation_dependent_variable(self):
        assert_roundtrip(obs_dep_var.integration_time_dependent_variable())

    def test_ancillary_observation_dependent_variable(self):
        assert_roundtrip(obs_dep_var.retransmission_delays_dependent_variable())

    def test_light_time_correction_components(self):
        assert_roundtrip(obs_dep_var.link_end_epochs_dependent_variable())


# ===========================================================================
# Observation Ancillary Simulation Settings
# ===========================================================================


class TestObservationAncillarySimulationSettingsPickle:
    """Pickle roundtrip tests for ObservationAncillarySimulationSettings."""

    def test_doppler_ancillary_settings(self):
        assert_roundtrip(ancillary_settings.doppler_ancillary_settings())


# ===========================================================================
# Link Types
# ===========================================================================


class TestLinkTypesPickle:
    """Pickle roundtrip tests for LinkEndId and LinkDefinition."""

    def test_link_end_id(self):
        assert_roundtrip(links.body_origin_link_end_id("Earth"))

    def test_link_definition(self):
        link_ends = {
            links.transmitter: links.body_origin_link_end_id("Earth"),
            links.receiver: links.body_origin_link_end_id("Delfi-C3"),
        }
        assert_roundtrip(links.LinkDefinition(link_ends))


# ===========================================================================
# Observation Collection / SingleObservationSet
# ===========================================================================

# Disabled until Obscollection refactor is done
# class TestObservationCollectionPickle:
#     """Pickle roundtrip tests for SingleObservationSet and ObservationCollection."""

#     def test_single_observation_set(self):
#         link_ends = {
#             links.transmitter: links.body_origin_link_end_id("Earth"),
#             links.receiver: links.body_origin_link_end_id("Delfi-C3"),
#         }
#         link_def = links.LinkDefinition(link_ends)
#         obs = [np.array([1.0])]
#         times = [Time(0, 0.0)]
#         obj = SingleObservationSet(
#             model_settings.one_way_range_type,
#             link_def,
#             obs,
#             times,
#             links.receiver,
#         )
#         assert_roundtrip(obj)

#     def test_observation_collection(self):
#         link_ends = {
#             links.transmitter: links.body_origin_link_end_id("Earth"),
#             links.receiver: links.body_origin_link_end_id("Delfi-C3"),
#         }
#         link_def = links.LinkDefinition(link_ends)
#         obs = [np.array([1.0])]
#         times = [Time(0, 0.0)]
#         single_set = SingleObservationSet(
#             model_settings.one_way_range_type,
#             link_def,
#             obs,
#             times,
#             links.receiver,
#         )
#         obj = ObservationCollection([single_set])
#         assert_roundtrip(obj)


# ===========================================================================
# Root Finder Settings
# ===========================================================================


class TestRootFinderSettingsPickle:
    """Pickle roundtrip tests for RootFinderSettings."""

    def test_bisection(self):
        assert_roundtrip(root_finders.bisection())

    def test_newton_raphson(self):
        assert_roundtrip(root_finders.newton_raphson())

    def test_halley(self):
        assert_roundtrip(root_finders.halley())

    def test_secant(self):
        assert_roundtrip(root_finders.secant())


# ===========================================================================
# Time
# ===========================================================================


class TestTimePickle:
    """Pickle roundtrip tests for Time."""

    def test_time_zero(self):
        assert_roundtrip(Time(0, 0.0))

    def test_time_nonzero(self):
        assert_roundtrip(Time(1, 0.5))

    def test_time_negative(self):
        assert_roundtrip(Time(-1, 0.0))


# ===========================================================================
# Gravity Field Variation Settings
# ===========================================================================


class TestGravityFieldVariationSettingsPickle:
    """Pickle roundtrip tests for GravityFieldVariationSettings."""

    def test_solid_body_tide(self):
        assert_roundtrip(gravity_field_variation.solid_body_tide("Moon", 0.3, 2))

    def test_solid_body_tide_degree_variable_k(self):
        assert_roundtrip(
            gravity_field_variation.solid_body_tide_degree_variable_k("Moon", {2: 0.3, 3: 0.2})
        )

    def test_solid_body_tide_degree_order_variable_k(self):
        assert_roundtrip(
            gravity_field_variation.solid_body_tide_degree_order_variable_k("Moon", {2: [0.3, 0.3]})
        )

    def test_single_period_periodic(self):
        assert_roundtrip(
            gravity_field_variation.single_period_periodic(
                cosine_coefficient_amplitude_cosine_time=np.zeros((1, 1)),
                cosine_coefficient_amplitude_sine_time=np.zeros((1, 1)),
                sine_coefficient_amplitude_cosine_time=np.zeros((1, 1)),
                sine_coefficient_amplitude_sine_time=np.zeros((1, 1)),
                angular_frequency=1.0e-4,
                reference_epoch=Time(0, 0.0),
                minimum_degree=2,
                minimum_order=0,
            )
        )

    def test_single_power_polynomial(self):
        assert_roundtrip(
            gravity_field_variation.single_power_polynomial(
                cosine_amplitudes=np.zeros((1, 1)),
                sine_amplitudes=np.zeros((1, 1)),
                polynomial_power=1,
                reference_epoch=Time(0, 0.0),
                minimum_degree=2,
                minimum_order=0,
            )
        )

    @pytest.mark.xfail(
        reason="Love numbers nested dict with tuple keys does not survive pickle roundtrip"
    )
    def test_mode_coupled_solid_body_tide(self):
        """ModeCoupledSolidBodyGravityFieldVariationSettings."""
        love_numbers = {(2, 0): {(2, 0): 0.3}}
        assert_roundtrip(
            gravity_field_variation.mode_coupled_solid_body_tide(["Moon"], love_numbers)
        )


# ===========================================================================
# Polymorphic Dispatch Tests
# ===========================================================================


class TestPolymorphicDispatch:
    """
    Verify that cereal correctly identifies and reconstructs the concrete
    derived type when serialized through a base class pointer.
    """

    # --- Acceleration polymorphic ---
    ACCELERATION_FACTORIES = [
        ("point_mass_gravity", lambda: acc.point_mass_gravity()),
        ("spherical_harmonic_gravity", lambda: acc.spherical_harmonic_gravity(4, 4)),
        ("radiation_pressure", lambda: acc.radiation_pressure()),
        (
            "relativistic_correction",
            lambda: acc.relativistic_correction(
                use_schwarzschild=True, use_lense_thirring=False, use_de_sitter=False
            ),
        ),
        ("empirical", lambda: acc.empirical(constant_acceleration=np.array([1, 2, 3]))),
        ("yarkovsky", lambda: acc.yarkovsky(0.1)),
        ("thrust_from_engines", lambda: acc.thrust_from_engines(["engine1"])),
        ("rtg", lambda: acc.rtg(np.array([0.1, 0.0, 0.0]), 0.1, 0.0)),
        (
            "direct_tidal_dissipation",
            lambda: acc.direct_tidal_dissipation_acceleration(0.2, 10, False, False),
        ),
    ]

    @pytest.mark.parametrize(
        "name,factory",
        ACCELERATION_FACTORIES,
        ids=lambda x: x[0] if isinstance(x, tuple) else x,
    )
    def test_acceleration_concrete_type_preserved(self, name, factory):
        obj = factory()
        assert_roundtrip(obj)

    # --- Torque polymorphic ---
    TORQUE_FACTORIES = [
        ("aerodynamic", lambda: torque.aerodynamic()),
        ("radiation_pressure_torque", lambda: torque.radiation_pressure_torque()),
        ("second_degree_gravitational", lambda: torque.second_degree_gravitational()),
        (
            "spherical_harmonic_gravitational",
            lambda: torque.spherical_harmonic_gravitational(2, 2),
        ),
    ]

    @pytest.mark.parametrize(
        "name,factory",
        TORQUE_FACTORIES,
        ids=lambda x: x[0] if isinstance(x, tuple) else x,
    )
    def test_torque_concrete_type_preserved(self, name, factory):
        obj = factory()
        assert_roundtrip(obj)

    # --- Termination polymorphic ---
    TERMINATION_FACTORIES = [
        ("time_termination", lambda: propagator.time_termination(1000.0)),
        ("cpu_time_termination", lambda: propagator.cpu_time_termination(60.0)),
        (
            "dependent_variable_termination",
            lambda: propagator.dependent_variable_termination(
                dep_var.altitude("Earth", "Earth"), 100.0, False
            ),
        ),
        (
            "hybrid_termination",
            lambda: propagator.hybrid_termination(
                [propagator.time_termination(500.0)], fulfill_single_condition=True
            ),
        ),
        (
            "non_sequential_termination",
            lambda: propagator.non_sequential_termination(
                propagator.time_termination(1000.0),
                propagator.time_termination(-1000.0),
            ),
        ),
    ]

    @pytest.mark.parametrize(
        "name,factory",
        TERMINATION_FACTORIES,
        ids=lambda x: x[0] if isinstance(x, tuple) else x,
    )
    def test_termination_concrete_type_preserved(self, name, factory):
        obj = factory()
        assert_roundtrip(obj)

    # --- Dependent variable polymorphic ---
    DEPVAR_FACTORIES = [
        ("altitude", lambda: dep_var.altitude("Earth", "Earth")),
        (
            "single_acceleration",
            lambda: dep_var.single_acceleration(acc.point_mass_gravity_type, "Earth", "Moon"),
        ),
        (
            "spherical_harmonic_terms",
            lambda: dep_var.spherical_harmonic_terms_acceleration("Earth", "Moon", [(2, 2)]),
        ),
        (
            "single_torque",
            lambda: dep_var.single_torque(torque.aerodynamic_type, "Earth", "Moon"),
        ),
        (
            "intermediate_aero_rotation",
            lambda: dep_var.intermediate_aerodynamic_rotation_matrix_variable(
                "Earth",
                AerodynamicsReferenceFrames.inertial_frame,
                AerodynamicsReferenceFrames.body_frame,
            ),
        ),
        ("body_aero_angle", lambda: dep_var.angle_of_attack("Earth", "Sun")),
        ("local_wind", lambda: dep_var.local_wind_velocity("Earth", "Earth")),
        (
            "control_surface",
            lambda: dep_var.control_surface_deflection("Earth", "aileron"),
        ),
        (
            "single_gravity_field_var",
            lambda: dep_var.single_gravity_field_variation_acceleration(
                "Earth", "Moon", BodyDeformationTypes.basic_solid_body
            ),
        ),
        (
            "single_per_term_gravity_field_var",
            lambda: dep_var.single_per_term_gravity_field_variation_acceleration(
                "Earth", "Moon", [(2, 2)], BodyDeformationTypes.basic_solid_body
            ),
        ),
        (
            "acceleration_partial",
            lambda: dep_var.acceleration_partial_wrt_body_translational_state(
                acc.point_mass_gravity_type, "Earth", "Moon", "Earth"
            ),
        ),
        (
            "total_acceleration_partial",
            lambda: dep_var.total_acceleration_partial_wrt_body_translational_state(
                "Earth", "Earth"
            ),
        ),
        (
            "min_body_distance",
            lambda: dep_var.minimum_body_distance("Earth", ["Moon", "Mars"]),
        ),
        (
            "min_station_distance",
            lambda: dep_var.minimum_visible_station_body_distances(
                "Earth", "station1", ["Moon", "Mars"], 10.0
            ),
        ),
        (
            "cross_section",
            lambda: dep_var.actual_cross_section("Earth", "Sun", "cannon_ball_radiation_pressure"),
        ),
        (
            "illuminated_panel_fraction",
            lambda: dep_var.illuminated_panel_fraction("Earth", "Sun"),
        ),
        (
            "total_gravity_field_var",
            lambda: dep_var.total_gravity_field_variation_acceleration("Earth", "Moon"),
        ),
    ]

    @pytest.mark.parametrize(
        "name,factory",
        DEPVAR_FACTORIES,
        ids=lambda x: x[0] if isinstance(x, tuple) else x,
    )
    def test_depvar_concrete_type_preserved(self, name, factory):
        obj = factory()
        assert_roundtrip(obj)

    # --- Root finder polymorphic ---
    ROOT_FINDER_FACTORIES = [
        ("bisection", lambda: root_finders.bisection()),
        ("newton_raphson", lambda: root_finders.newton_raphson()),
        ("halley", lambda: root_finders.halley()),
        ("secant", lambda: root_finders.secant()),
    ]

    @pytest.mark.parametrize(
        "name,factory",
        ROOT_FINDER_FACTORIES,
        ids=lambda x: x[0] if isinstance(x, tuple) else x,
    )
    def test_root_finder_concrete_type_preserved(self, name, factory):
        obj = factory()
        assert_roundtrip(obj)


# ===========================================================================
# Simulation-Based Pickle Tests
# ===========================================================================


class TestSimulationBasedPickle:
    """
    Pickle roundtrip tests that require running a real propagation.
    Tests PropagationTerminationDetails, SimulationResults, and related types.
    """

    @staticmethod
    def _run_propagation():
        """Run a simple propagation and return the dynamics simulator."""
        from tudatpy.interface import spice
        from tudatpy.astro import element_conversion
        from tudatpy.astro.time_representation import DateTime
        from tudatpy.dynamics import environment_setup, propagation_setup, simulator

        spice.load_standard_kernels()
        bodies = environment_setup.create_system_of_bodies(
            environment_setup.get_default_body_settings(["Earth"], "Earth", "J2000")
        )
        bodies.add_empty_settings("Delfi-C3")
        bodies.get("Delfi-C3").mass = 2.2

        start = DateTime(2008, 4, 28).to_epoch()
        initial_state = element_conversion.keplerian_to_cartesian_elementwise(
            gravitational_parameter=bodies.get("Earth").gravitational_parameter,
            semi_major_axis=6.99276221e06,
            eccentricity=4.03294322e-03,
            inclination=1.71065169e00,
            argument_of_periapsis=1.31226971e00,
            longitude_of_ascending_node=3.82958313e-01,
            true_anomaly=3.07018490e00,
        )
        propagator_settings = propagation_setup.propagator.translational(
            ["Earth"],
            propagation_setup.create_acceleration_models(
                bodies,
                {"Delfi-C3": {"Earth": [propagation_setup.acceleration.point_mass_gravity()]}},
                ["Delfi-C3"],
                ["Earth"],
            ),
            ["Delfi-C3"],
            initial_state,
            start,
            propagation_setup.integrator.runge_kutta_fixed_step(
                30.0, coefficient_set=propagation_setup.integrator.CoefficientSets.rk_4
            ),
            propagation_setup.propagator.time_termination(start + 600.0),
        )
        return simulator.create_dynamics_simulator(bodies, propagator_settings)

    def test_propagation_termination_details(self):
        """Pickle roundtrip for PropagationTerminationDetails."""
        pytest.importorskip("spice")
        sim = self._run_propagation()
        obj = sim.propagation_results.termination_details
        assert_roundtrip(obj)

    def test_single_arc_simulation_results(self):
        """Pickle roundtrip for SingleArcSimulationResults."""
        pytest.importorskip("spice")
        sim = self._run_propagation()
        obj = sim.propagation_results
        assert_roundtrip(obj)

    def test_propagation_termination_details_from_hybrid(self):
        """Pickle roundtrip for PropagationTerminationDetailsFromHybridCondition."""
        pytest.importorskip("spice")
        from tudatpy.interface import spice
        from tudatpy.astro import element_conversion
        from tudatpy.astro.time_representation import DateTime
        from tudatpy.dynamics import environment_setup, propagation_setup, simulator

        spice.load_standard_kernels()
        bodies = environment_setup.create_system_of_bodies(
            environment_setup.get_default_body_settings(["Earth"], "Earth", "J2000")
        )
        bodies.add_empty_settings("Delfi-C3")
        bodies.get("Delfi-C3").mass = 2.2

        start = DateTime(2008, 4, 28).to_epoch()
        initial_state = element_conversion.keplerian_to_cartesian_elementwise(
            gravitational_parameter=bodies.get("Earth").gravitational_parameter,
            semi_major_axis=6.99276221e06,
            eccentricity=4.03294322e-03,
            inclination=1.71065169e00,
            argument_of_periapsis=1.31226971e00,
            longitude_of_ascending_node=3.82958313e-01,
            true_anomaly=3.07018490e00,
        )
        propagator_settings = propagation_setup.propagator.translational(
            ["Earth"],
            propagation_setup.create_acceleration_models(
                bodies,
                {"Delfi-C3": {"Earth": [propagation_setup.acceleration.point_mass_gravity()]}},
                ["Delfi-C3"],
                ["Earth"],
            ),
            ["Delfi-C3"],
            initial_state,
            start,
            propagation_setup.propagator.hybrid_termination(
                [
                    propagation_setup.propagator.time_termination(start + 600.0),
                    propagation_setup.propagator.cpu_time_termination(600.0),
                ],
                fulfill_single_condition=True,
            ),
        )
        sim = simulator.create_dynamics_simulator(bodies, propagator_settings)
        obj = sim.propagation_results.termination_details
        assert_roundtrip(obj)

    def test_single_arc_variational_simulation_results(self):
        """Pickle roundtrip for SingleArcVariationalSimulationResults."""
        pytest.importorskip("spice")
        from tudatpy.interface import spice
        from tudatpy.astro import element_conversion
        from tudatpy.astro.time_representation import DateTime
        from tudatpy.dynamics import environment_setup, propagation_setup, simulator

        spice.load_standard_kernels()
        bodies = environment_setup.create_system_of_bodies(
            environment_setup.get_default_body_settings(["Earth"], "Earth", "J2000")
        )
        bodies.add_empty_settings("Delfi-C3")
        bodies.get("Delfi-C3").mass = 2.2

        start = DateTime(2008, 4, 28).to_epoch()
        initial_state = element_conversion.keplerian_to_cartesian_elementwise(
            gravitational_parameter=bodies.get("Earth").gravitational_parameter,
            semi_major_axis=6.99276221e06,
            eccentricity=4.03294322e-03,
            inclination=1.71065169e00,
            argument_of_periapsis=1.31226971e00,
            longitude_of_ascending_node=3.82958313e-01,
            true_anomaly=3.07018490e00,
        )
        propagator_settings = propagation_setup.propagator.translational(
            ["Earth"],
            propagation_setup.create_acceleration_models(
                bodies,
                {"Delfi-C3": {"Earth": [propagation_setup.acceleration.point_mass_gravity()]}},
                ["Delfi-C3"],
                ["Earth"],
            ),
            ["Delfi-C3"],
            initial_state,
            start,
            propagation_setup.integrator.runge_kutta_fixed_step(
                30.0, coefficient_set=propagation_setup.integrator.CoefficientSets.rk_4
            ),
            propagation_setup.propagator.time_termination(start + 600.0),
        )
        variational_equations = propagation_setup.variational_equations(
            bodies,
            ["Delfi-C3"],
            [propagation_setup.variation.estimated_bodies_settings(["Delfi-C3"])],
            propagate_on_creation=True,
        )
        sim = simulator.create_dynamics_simulator(
            bodies, propagator_settings, variational_equations
        )
        obj = sim.propagation_results
        assert_roundtrip(obj)


# ===========================================================================
# Untestable Classes (explanatory notes)
# ===========================================================================
#
# The following classes from the serialization inventory cannot be tested:
#
# 1. Not exposed as pybind11 classes (C++ internal only):
#    - FixedTimeHodographicShapingOptimisationProblem
#    - HodographicShapingOptimisationProblem
#    - ObservationDependentVariableBookkeeping
#    - ModelInterpolationSettings
#    - PropagatorType (exists as enum TranslationalPropagatorType instead)
#
# 2. Exposed but lack __eq__ or __reduce__:
#    - InterpolatorSettings / LagrangeInterpolatorSettings
#      (value types with direct constructors, no serialization macros)
#
# 3. Require full multi-arc / hybrid-arc propagation setup:
#    - MultiArcSimulationResults
#    - HybridArcSimulationResults
#      (need dedicated multi-arc/hybrid test infrastructure)
#
# 4. Require full estimation pipeline:
#    - CovarianceAnalysisOutput
#    - EstimationOutput
#      (need Estimator with observations, too heavy for unit tests)
#
# 5. Missing Python factory:
#    - MomentumWheelDesaturationAccelerationSettings
#      (registered in bindings but no public factory function)
