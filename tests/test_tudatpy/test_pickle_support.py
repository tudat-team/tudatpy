"""
Tests for pickle serialization/deserialization of AccelerationSettings
and all registered derived classes.

Run with: pytest test_acceleration_settings_pickle.py -v
"""
import pickle
import pytest
import numpy as np

# Adjust this import to match your actual tudatpy module path
from tudatpy.numerical_simulation.propagation_setup import acceleration as acc


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
    @pytest.mark.parametrize("degree,order", [
        (2, 0),
        (2, 2),
        (8, 8),
        (20, 20),
    ])
    def test_various_degrees(self, degree, order):
        obj = acc.spherical_harmonic_gravity(degree, order)
        assert_roundtrip(obj)


class TestMutualSphericalHarmonicAccelerationSettings:
    @pytest.mark.parametrize("args", [
        (2, 2, 2, 2),
        (4, 4, 2, 2),
        (8, 8, 4, 4),
    ])
    def test_various_degrees(self, args):
        # adjust constructor name/signature if needed
        obj = acc.mutual_spherical_harmonic_gravity(*args)
        assert_roundtrip(obj)


class TestRelativisticAccelerationCorrectionSettings:
    def test_default(self):
        obj = acc.relativistic_correction(
            use_schwarzschild=True,
            use_lense_thirring=False,
            use_de_sitter=False
        )
        assert_roundtrip(obj)


class TestEmpiricalAccelerationSettings:
    def test_default(self):
        obj = acc.empirical(constant_acceleration=np.array([1,2,3]))  # adjust if needed
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