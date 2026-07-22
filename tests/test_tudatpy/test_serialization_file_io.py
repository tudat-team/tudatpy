"""
File-format roundtrip tests for serializable classes that expose
save_to_json / load_from_json and/or save_to_binary / load_from_binary.

For each class that implements file-format methods, we:
  1. Create an object
  2. Save it to a temporary file via save_to_json / save_to_binary
  3. Load it back via load_from_json / load_from_binary
  4. Verify equality and type preservation

File extensions (appended automatically by the file writers):
  - Binary: .tudat
  - JSON:   .json

Run with: pytest test_serialization_file_io.py -v
"""

import pickle
import tempfile
import os

import pytest
import numpy as np

from tudatpy.astro.time_representation import Time
from tudatpy.estimation.observable_models_setup import links
from tudatpy.estimation.observations_setup import ancillary_settings
from tudatpy.estimation.observations_setup import observations_dependent_variables as obs_dep_var
from tudatpy.estimation.observations import (
    SingleObservationSet,
    ObservationCollection,
)  # noqa: F401 — skipped during refactoring
from tudatpy.math import root_finders
from tudatpy.dynamics.propagation_setup import dependent_variable as dep_var
from tudatpy.dynamics.propagation_setup import torque
from tudatpy.dynamics.propagation_setup import propagator
from tudatpy.dynamics.propagation_setup import acceleration as acc
from tudatpy.dynamics import environment
from tudatpy.dynamics.environment_setup.aerodynamic_coefficients import AerodynamicsReferenceFrames
from tudatpy.dynamics.environment_setup.gravity_field_variation import BodyDeformationTypes
from tudatpy.estimation.observable_models_setup import model_settings
from tudatpy.dynamics.propagation import (
    PropagationTerminationDetails,
    SimulationResults,
    SingleArcSimulationResults,
    SingleArcVariationalSimulationResults,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _temp_path(suffix=""):
    """Return a temporary file path (stem, without extension)."""
    fd, path = tempfile.mkstemp(suffix=suffix)
    os.close(fd)
    os.unlink(path)
    return path


def assert_json_roundtrip(obj, **save_kwargs):
    """Save object to JSON, load it back, and verify equality + type."""
    path = _temp_path()
    try:
        obj.save_to_json(path, **save_kwargs)
        loaded = type(obj).load_from_json(path, **save_kwargs)
        assert loaded == obj, f"JSON roundtrip failed for {type(obj).__name__}"
        assert type(loaded) is type(obj), (
            f"Type mismatch after JSON roundtrip: "
            f"expected {type(obj).__name__}, got {type(loaded).__name__}"
        )
    finally:
        for ext in (".json", ".tudat"):
            full = path + ext
            if os.path.exists(full):
                os.unlink(full)


def assert_json_poly_roundtrip(obj, base_class):
    """
    Polymorphic JSON roundtrip: save object through base class,
    load through base class, verify equality + concrete type.
    Uses base_class.__eq__ because derived __eq__ may be type-strict.
    Type preservation is verified implicitly by C++ dynamic_cast in equals().
    """
    path = _temp_path()
    try:
        obj.save_to_json(path)
        loaded = base_class.load_from_json(path)
        assert base_class.__eq__(
            loaded, obj
        ), f"JSON polymorphic roundtrip failed for {type(obj).__name__}"
    finally:
        for ext in (".json", ".tudat"):
            full = path + ext
            if os.path.exists(full):
                os.unlink(full)


def assert_binary_roundtrip(obj, **save_kwargs):
    """Save object to binary, load it back, and verify equality + type."""
    path = _temp_path()
    try:
        obj.save_to_binary(path, **save_kwargs)
        loaded = type(obj).load_from_binary(path, **save_kwargs)
        assert loaded == obj, f"Binary roundtrip failed for {type(obj).__name__}"
        assert type(loaded) is type(obj), (
            f"Type mismatch after binary roundtrip: "
            f"expected {type(obj).__name__}, got {type(loaded).__name__}"
        )
    finally:
        for ext in (".json", ".tudat"):
            full = path + ext
            if os.path.exists(full):
                os.unlink(full)


def assert_binary_poly_roundtrip(obj, base_class):
    """
    Polymorphic binary roundtrip: save object through base class,
    load through base class, verify equality + concrete type.
    Uses base_class.__eq__ because derived __eq__ may be type-strict.
    Type preservation is verified implicitly by C++ dynamic_cast in equals().
    """
    path = _temp_path()
    try:
        obj.save_to_binary(path)
        loaded = base_class.load_from_binary(path)
        assert base_class.__eq__(
            loaded, obj
        ), f"Binary polymorphic roundtrip failed for {type(obj).__name__}"
    finally:
        for ext in (".json", ".tudat"):
            full = path + ext
            if os.path.exists(full):
                os.unlink(full)


# ===========================================================================
# 1. ACCELERATION SETTINGS  (polymorphic, has FILE_IO)
# ===========================================================================


class TestAccelerationSettingsFileIO:
    """JSON and binary file roundtrip for AccelerationSettings (+ derived types)."""

    @pytest.mark.parametrize(
        "obj",
        [
            acc.point_mass_gravity(),
            acc.spherical_harmonic_gravity(4, 4),
            acc.radiation_pressure(),
            acc.cannonball_radiation_pressure(),
            acc.relativistic_correction(
                use_schwarzschild=True, use_lense_thirring=False, use_de_sitter=False
            ),
            pytest.param(
                acc.empirical(constant_acceleration=np.array([1.0, 2.0, 3.0])),
                id="empirical",
            ),
            acc.yarkovsky(0.1),
            acc.thrust_from_engines(["engine1"]),
            acc.rtg(np.array([0.1, 0.0, 0.0]), 0.1, 0.0),
            acc.direct_tidal_dissipation_acceleration(0.2, 10, False, False),
            acc.quasi_impulsive_shots_acceleration(
                [0.0], [np.array([0.1, 0.0, 0.0])], Time(0, 0.0), Time(0, 0.0)
            ),
        ],
        ids=lambda obj: type(obj).__name__,
    )
    def test_json_roundtrip(self, obj):
        assert_json_roundtrip(obj)

    @pytest.mark.parametrize(
        "obj",
        [
            acc.point_mass_gravity(),
            acc.spherical_harmonic_gravity(4, 4),
            acc.radiation_pressure(),
            acc.cannonball_radiation_pressure(),
            acc.relativistic_correction(
                use_schwarzschild=True, use_lense_thirring=False, use_de_sitter=False
            ),
            pytest.param(
                acc.empirical(constant_acceleration=np.array([1.0, 2.0, 3.0])),
                id="empirical",
            ),
            acc.yarkovsky(0.1),
            acc.thrust_from_engines(["engine1"]),
            acc.rtg(np.array([0.1, 0.0, 0.0]), 0.1, 0.0),
            acc.direct_tidal_dissipation_acceleration(0.2, 10, False, False),
            acc.quasi_impulsive_shots_acceleration(
                [0.0], [np.array([0.1, 0.0, 0.0])], Time(0, 0.0), Time(0, 0.0)
            ),
        ],
        ids=lambda obj: type(obj).__name__,
    )
    def test_binary_roundtrip(self, obj):
        assert_binary_roundtrip(obj)

    def test_polymorphic_json_through_base(self):
        """Poly dispatch: save derived, load through AccelerationSettings base."""
        obj = acc.spherical_harmonic_gravity(4, 4)
        assert_json_poly_roundtrip(obj, acc.AccelerationSettings)

    def test_polymorphic_binary_through_base(self):
        """Poly dispatch: save derived, load through AccelerationSettings base."""
        obj = acc.empirical(constant_acceleration=np.array([1.0, 2.0, 3.0]))
        assert_binary_poly_roundtrip(obj, acc.AccelerationSettings)


# ===========================================================================
# 2. TORQUE SETTINGS  (polymorphic, has FILE_IO)
# ===========================================================================


class TestTorqueSettingsFileIO:
    """JSON and binary file roundtrip for TorqueSettings (+ derived types)."""

    TORQUE_OBJS = [
        torque.aerodynamic(),
        torque.radiation_pressure_torque(),
        torque.second_degree_gravitational(),
        torque.spherical_harmonic_gravitational(2, 2),
    ]

    @pytest.mark.parametrize("obj", TORQUE_OBJS, ids=lambda obj: type(obj).__name__)
    def test_json_roundtrip(self, obj):
        assert_json_roundtrip(obj)

    @pytest.mark.parametrize("obj", TORQUE_OBJS, ids=lambda obj: type(obj).__name__)
    def test_binary_roundtrip(self, obj):
        assert_binary_roundtrip(obj)

    def test_polymorphic_json_through_base(self):
        obj = torque.spherical_harmonic_gravitational(2, 2)
        assert_json_poly_roundtrip(obj, torque.TorqueSettings)

    def test_polymorphic_binary_through_base(self):
        obj = torque.spherical_harmonic_gravitational(2, 2)
        assert_binary_poly_roundtrip(obj, torque.TorqueSettings)


# ===========================================================================
# 3. PROPAGATION TERMINATION SETTINGS  (polymorphic, has FILE_IO)
# ===========================================================================


class TestPropagationTerminationSettingsFileIO:
    """JSON and binary file roundtrip for PropagationTerminationSettings (+ derived types)."""

    TERMINATION_OBJS = [
        propagator.time_termination(1000.0),
        propagator.cpu_time_termination(60.0),
        propagator.dependent_variable_termination(dep_var.altitude("Earth", "Earth"), 100.0, False),
        propagator.hybrid_termination(
            [propagator.time_termination(500.0), propagator.cpu_time_termination(30.0)],
            fulfill_single_condition=True,
        ),
        propagator.non_sequential_termination(
            propagator.time_termination(1000.0), propagator.time_termination(-1000.0)
        ),
    ]

    @pytest.mark.parametrize("obj", TERMINATION_OBJS, ids=lambda obj: type(obj).__name__)
    def test_json_roundtrip(self, obj):
        assert_json_roundtrip(obj)

    @pytest.mark.parametrize("obj", TERMINATION_OBJS, ids=lambda obj: type(obj).__name__)
    def test_binary_roundtrip(self, obj):
        assert_binary_roundtrip(obj)

    def test_polymorphic_json_through_base(self):
        obj = propagator.hybrid_termination(
            [propagator.time_termination(500.0)], fulfill_single_condition=True
        )
        assert_json_poly_roundtrip(obj, propagator.PropagationTerminationSettings)

    def test_polymorphic_binary_through_base(self):
        obj = propagator.non_sequential_termination(
            propagator.time_termination(1000.0), propagator.time_termination(-1000.0)
        )
        assert_binary_poly_roundtrip(obj, propagator.PropagationTerminationSettings)


# ===========================================================================
# 4. DEPENDENT VARIABLE SETTINGS  (polymorphic, has FILE_IO)
# ===========================================================================


class TestVariableSettingsFileIO:
    """JSON and binary file roundtrip for VariableSettings (+ derived types)."""

    DEPVAR_OBJS = [
        dep_var.altitude("Earth", "Earth"),
        dep_var.single_acceleration(acc.point_mass_gravity_type, "Earth", "Moon"),
        dep_var.spherical_harmonic_terms_acceleration("Earth", "Moon", [(2, 2)]),
        dep_var.single_torque(torque.aerodynamic_type, "Earth", "Moon"),
        dep_var.intermediate_aerodynamic_rotation_matrix_variable(
            "Earth",
            AerodynamicsReferenceFrames.inertial_frame,
            AerodynamicsReferenceFrames.body_frame,
        ),
        dep_var.angle_of_attack("Earth", "Sun"),
        dep_var.local_wind_velocity("Earth", "Earth"),
        dep_var.control_surface_deflection("Earth", "aileron"),
        dep_var.single_gravity_field_variation_acceleration(
            "Earth", "Moon", BodyDeformationTypes.basic_solid_body
        ),
        dep_var.single_per_term_gravity_field_variation_acceleration(
            "Earth", "Moon", [(2, 2)], BodyDeformationTypes.basic_solid_body
        ),
        dep_var.acceleration_partial_wrt_body_translational_state(
            acc.point_mass_gravity_type, "Earth", "Moon", "Earth"
        ),
        dep_var.total_acceleration_partial_wrt_body_translational_state("Earth", "Earth"),
        dep_var.minimum_body_distance("Earth", ["Moon", "Mars"]),
        dep_var.minimum_visible_station_body_distances("Earth", "station1", ["Moon", "Mars"], 10.0),
        dep_var.actual_cross_section("Earth", "Sun", "cannon_ball_radiation_pressure"),
        dep_var.illuminated_panel_fraction("Earth", "Sun"),
        dep_var.total_gravity_field_variation_acceleration("Earth", "Moon"),
    ]

    @pytest.mark.parametrize("obj", DEPVAR_OBJS, ids=lambda obj: type(obj).__name__)
    def test_json_roundtrip(self, obj):
        assert_json_poly_roundtrip(obj, dep_var.VariableSettings)

    @pytest.mark.parametrize("obj", DEPVAR_OBJS, ids=lambda obj: type(obj).__name__)
    def test_binary_roundtrip(self, obj):
        assert_binary_poly_roundtrip(obj, dep_var.VariableSettings)

    def test_polymorphic_json_through_base(self):
        obj = dep_var.single_acceleration(acc.point_mass_gravity_type, "Earth", "Moon")
        assert_json_poly_roundtrip(obj, dep_var.VariableSettings)

    def test_polymorphic_binary_through_base(self):
        obj = dep_var.total_gravity_field_variation_acceleration("Earth", "Moon")
        assert_binary_poly_roundtrip(obj, dep_var.VariableSettings)


# ===========================================================================
# 5. OBSERVATION DEPENDENT VARIABLE SETTINGS  (has FILE_IO)
# ===========================================================================


class TestObservationDependentVariableSettingsFileIO:
    """JSON and binary file roundtrip for ObservationDependentVariableSettings."""

    OBS_DEPVAR_OBJS = [
        obs_dep_var.elevation_angle_dependent_variable(links.receiver),
        obs_dep_var.azimuth_angle_dependent_variable(links.receiver),
        obs_dep_var.integration_time_dependent_variable(),
        obs_dep_var.retransmission_delays_dependent_variable(),
        obs_dep_var.link_end_epochs_dependent_variable(),
    ]

    @pytest.mark.parametrize("obj", OBS_DEPVAR_OBJS, ids=lambda obj: type(obj).__name__)
    def test_json_roundtrip(self, obj):
        assert_json_roundtrip(obj)

    @pytest.mark.parametrize("obj", OBS_DEPVAR_OBJS, ids=lambda obj: type(obj).__name__)
    def test_binary_roundtrip(self, obj):
        assert_binary_roundtrip(obj)


# ===========================================================================
# 6. OBSERVATION ANCILLARY SIMULATION SETTINGS  (has FILE_IO)
# ===========================================================================


class TestObservationAncillarySettingsFileIO:
    """JSON and binary file roundtrip for ObservationAncillarySimulationSettings."""

    def test_json_roundtrip(self):
        obj = ancillary_settings.doppler_ancillary_settings()
        assert_json_roundtrip(obj)

    def test_binary_roundtrip(self):
        obj = ancillary_settings.doppler_ancillary_settings()
        assert_binary_roundtrip(obj)


# ===========================================================================
# 7. ROOT FINDER SETTINGS  (has FILE_IO)
# ===========================================================================


class TestRootFinderSettingsFileIO:
    """JSON and binary file roundtrip for RootFinderSettings."""

    ROOT_FINDER_OBJS = [
        root_finders.bisection(),
        root_finders.newton_raphson(),
        root_finders.halley(),
        root_finders.secant(),
    ]

    @pytest.mark.parametrize("obj", ROOT_FINDER_OBJS, ids=lambda obj: type(obj).__name__)
    def test_json_roundtrip(self, obj):
        assert_json_roundtrip(obj)

    @pytest.mark.parametrize("obj", ROOT_FINDER_OBJS, ids=lambda obj: type(obj).__name__)
    def test_binary_roundtrip(self, obj):
        assert_binary_roundtrip(obj)


# ===========================================================================
# 8. LINK TYPES  (has FILE_IO)
# ===========================================================================


class TestLinkTypesFileIO:
    """JSON and binary file roundtrip for LinkEndId and LinkDefinition."""

    def test_link_end_id_json(self):
        assert_json_roundtrip(links.body_origin_link_end_id("Earth"))

    def test_link_end_id_binary(self):
        assert_binary_roundtrip(links.body_origin_link_end_id("Earth"))

    def test_link_definition_json(self):
        link_ends = {
            links.transmitter: links.body_origin_link_end_id("Earth"),
            links.receiver: links.body_origin_link_end_id("Delfi-C3"),
        }
        assert_json_roundtrip(links.LinkDefinition(link_ends))

    def test_link_definition_binary(self):
        link_ends = {
            links.transmitter: links.body_origin_link_end_id("Earth"),
            links.receiver: links.body_origin_link_end_id("Delfi-C3"),
        }
        assert_binary_roundtrip(links.LinkDefinition(link_ends))


# ===========================================================================
# 9. OBSERVATION COLLECTION / SINGLEOBSERVATIONSET  (has BINARY_IO only)
# ===========================================================================


class TestObservationCollectionFileIO:
    """Binary file roundtrip for SingleObservationSet and ObservationCollection.

    .. note::
        Skipped — SingleObservationSet and ObservationCollection are undergoing
        major refactoring and the constructor API is not stable.
    """

    @staticmethod
    def _make_single_set():
        pytest.skip("SingleObservationSet/ObservationCollection under refactoring")

    def test_single_observation_set_binary(self):
        pytest.skip("SingleObservationSet under refactoring")

    def test_observation_collection_binary(self):
        pytest.skip("ObservationCollection under refactoring")


# ===========================================================================
# 10. TIME  (has BINARY_IO only)
# ===========================================================================


class TestTimeFileIO:
    """Binary file roundtrip for Time."""

    def test_time_zero(self):
        assert_binary_roundtrip(Time(0, 0.0))

    def test_time_nonzero(self):
        assert_binary_roundtrip(Time(1, 0.5))

    def test_time_negative(self):
        assert_binary_roundtrip(Time(-1, 0.0))


# ===========================================================================
# 11. PROPAGATION TERMINATION DETAILS  (custom save_binary/load_binary)
# ===========================================================================


class TestPropagationTerminationDetailsFileIO:
    """
    Binary file roundtrip for PropagationTerminationDetails.
    These use custom save_binary/load_binary lambdas (not the macros).
    """

    def test_termination_details_binary(self):
        """Requires a real propagation to create a valid instance."""
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
        dynamics_simulator = simulator.create_dynamics_simulator(bodies, propagator_settings)
        obj = dynamics_simulator.propagation_results.termination_details

        path = _temp_path()
        try:
            obj.save_binary(path)
            loaded = PropagationTerminationDetails.load_binary(path)
            assert loaded == obj
            assert type(loaded) is type(obj)
        finally:
            for ext in (".json", ".tudat"):
                full = path + ext
                if os.path.exists(full):
                    os.unlink(full)


# ===========================================================================
# 12. SIMULATION RESULTS  (custom save_binary/load_binary)
# ===========================================================================


class TestSimulationResultsFileIO:
    """
    Binary file roundtrip for SimulationResults and derived classes.
    These use custom save_binary/load_binary lambdas through the base class.
    """

    @staticmethod
    def _single_arc_results():
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
        dynamics_simulator = simulator.create_dynamics_simulator(bodies, propagator_settings)
        return dynamics_simulator.propagation_results

    def test_single_arc_simulation_results_binary(self):
        """Save/load SingleArcSimulationResults through binary file."""
        obj = self._single_arc_results()
        path = _temp_path()
        try:
            obj.save_binary(path)
            loaded = SimulationResults.load_binary(path)
            assert type(loaded) is type(
                obj
            ), f"Expected {type(obj).__name__}, got {type(loaded).__name__}"
            # Compare state dictionaries
            np.testing.assert_allclose(
                list(loaded.state_history.values()),
                list(obj.state_history.values()),
            )
            assert loaded.propagated_state_vector_length == obj.propagated_state_vector_length
            assert loaded.propagation_is_performed == obj.propagation_is_performed
        finally:
            for ext in (".json", ".tudat"):
                full = path + ext
                if os.path.exists(full):
                    os.unlink(full)

    def test_single_arc_variational_simulation_results_binary(self):
        """Save/load SingleArcVariationalSimulationResults through binary file."""
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
        assert type(obj) is SingleArcVariationalSimulationResults

        path = _temp_path()
        try:
            obj.save_to_binary(path)
            loaded = SimulationResults.load_from_binary(path)
            assert type(loaded) is type(obj)
            assert loaded.propagated_state_vector_length == obj.propagated_state_vector_length
            assert loaded.propagation_is_performed == obj.propagation_is_performed
        finally:
            for ext in (".json", ".tudat"):
                full = path + ext
                if os.path.exists(full):
                    os.unlink(full)


# ===========================================================================
# 13. HYBRID TERMINATION DETAILS  (custom save_binary/load_binary)
# ===========================================================================


class TestHybridTerminationDetailsFileIO:
    """
    Binary file roundtrip for PropagationTerminationDetailsFromHybridCondition.
    Uses hybrid termination to create the derived type, saves/loads through base.
    """

    def test_hybrid_termination_details_binary(self):
        """Save/load PropagationTerminationDetailsFromHybridCondition through binary."""
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
        dynamics_simulator = simulator.create_dynamics_simulator(bodies, propagator_settings)
        obj = dynamics_simulator.propagation_results.termination_details

        path = _temp_path()
        try:
            obj.save_binary(path)
            loaded = PropagationTerminationDetails.load_binary(path)
            assert loaded == obj
            assert type(loaded) is type(obj)
        finally:
            for ext in (".json", ".tudat"):
                full = path + ext
                if os.path.exists(full):
                    os.unlink(full)


# ===========================================================================
# Notes on untestable classes
# ===========================================================================
#
# The following classes from the serialization inventory cannot be tested
# for file-format roundtrip:
#
# 1. Have only C++ binary IO (no Python save_to_json / save_to_binary):
#    - All concrete AccelerationSettings derived types
#      (inherit file_io from AccelerationSettings base, tested via polymorphic)
#    - All concrete TorqueSettings derived types
#      (inherit file_io from TorqueSettings base, tested via polymorphic)
#    - All concrete PropagationTerminationSettings derived types
#      (inherit file_io from base, tested via polymorphic)
#    - All concrete VariableSettings derived types
#      (inherit file_io from base, tested via polymorphic)
#    - All concrete ObservationDependentVariableSettings derived types
#      (inherit file_io from base, tested via polymorphic)
#
# 2. Require full estimation pipeline (save_binary/load_binary):
#    - CovarianceAnalysisOutput
#    - EstimationOutput
#      (need Estimator with observations, too heavy for unit tests)
#
# 3. Require full multi-arc / hybrid-arc propagation setup:
#    - MultiArcSimulationResults
#    - HybridArcSimulationResults
#      (need dedicated multi-arc/hybrid test infrastructure)
#
# 4. Not exposed as pybind11 classes:
#    - ObservationDependentVariableBookkeeping
#    - ModelInterpolationSettings
#    - PropagatorType (exists as enum TranslationalPropagatorType instead)
#
# 5. No file_io methods exposed:
#    - InterpolatorSettings / LagrangeInterpolatorSettings
#    - GravityFieldVariationSettings (base and all derived)
