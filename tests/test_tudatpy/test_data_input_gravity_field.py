"""Gravity records, consuming settings and assembled-model binding contracts."""

from pathlib import Path
import pickle

import numpy as np
import pytest

from tudatpy.data_input.environment_data import earth_orientation
from tudatpy.data_input.environment_data import gravity_field as gravity_data
from tudatpy.dynamics import environment_setup
from tudatpy.dynamics.environment_setup import ephemeris, gravity_field, rotation_model
from tudatpy.dynamics.environment_setup import gravity_field_variation as variation
from tudatpy.math import interpolators

# Existing pinned provider extracts also supply the native scientific tests.
TEST_DATA = Path(__file__).parents[1] / "test_tudat"


# Provider conventions must survive adaptation without silently changing coefficients.
def test_icgem_records_settings_and_assembled_tide_system():
    """Adapt a provider field and retain an explicit tide declaration through assembly."""
    data = gravity_data.read_icgem_gravity_field(
        str(TEST_DATA / "reference_files/icgem_gravity/gfz_rl061_ddk4_202104_degree2.gfc")
    )
    data.validate()
    assert data.coefficients[(0, 0)].cosine == 1.0
    assert data.coefficients[(2, 0)].formal_errors is not None
    assert data.tide_system == gravity_data.TideSystem.tide_free

    settings = gravity_field.from_icgem(data, "ITRS")

    original = settings.normalized_cosine_coefficients.copy()

    assert settings.tide_system == data.tide_system
    settings.tide_system = gravity_data.TideSystem.zero_tide
    np.testing.assert_array_equal(settings.normalized_cosine_coefficients, original)

    bodies_settings = environment_setup.BodyListSettings("SSB", "GCRS")

    bodies_settings.add_empty_settings("Earth")
    bodies_settings.get("Earth").ephemeris_settings = ephemeris.constant(np.zeros(6), "SSB", "GCRS")
    bodies_settings.get("Earth").gravity_field_settings = settings
    bodies = environment_setup.create_system_of_bodies(bodies_settings)

    assert bodies.get("Earth").gravity_field_model.tide_system == settings.tide_system

    # Temporal records retain printed dates and cannot silently pass through a static adapter.
    temporal = gravity_data.read_icgem_gravity_field(
        str(TEST_DATA / "reference_files/icgem_gravity/eigen6s4_v2_c20_1950_2006.gfc")
    )

    assert isinstance(temporal.temporal_coefficients[0], gravity_data.IcgemTemporalCoefficient)
    assert temporal.temporal_coefficients[0].reference_date
    assert temporal.temporal_coefficients[0].validity_end_date
    with pytest.raises(RuntimeError):
        gravity_field.from_icgem(temporal, "ITRS")


# One construction exercises all native input types and the existing deferred sampler.
def gravity_settings():
    """Supply physical source constants and a shared interior sampling interval."""
    time = (45701.75 - 51544.5) * 86400.0 + 54.184

    sampling = variation.ModelInterpolationSettings(
        time - 20.0, time + 20.0, 10.0, interpolators.linear_interpolation()
    )

    # The ocean-pole table contains integrals; its forcing constants are separate inputs.
    ocean = gravity_data.read_ocean_pole_gravity_data(
        str(TEST_DATA / "reference_files/ocean_pole_gravity/desaiscopolecoef_degree180.txt")
    )

    constants = variation.OceanPoleGravityConstants()

    constants.angular_velocity = 7.292115e-5
    constants.gravitational_constant = 6.673e-11
    constants.ocean_density = 1025.0
    constants.equatorial_gravity = 9.7803278
    constants.gravitational_parameter = 3.986004418e14
    constants.reference_radius = 6378136.6
    constants.gamma2 = complex(0.687, 0.0036)

    potential = gravity_data.read_harmonic_potential_data(
        str(TEST_DATA / "reference_files/harmonic_potential/fes2004-degree3.potential"),
        gravity_data.HarmonicPotentialProfile.fes2004,
    )

    potential.validate()
    assert len(potential.constituents[0].multipliers) == 6
    metadata = variation.HarmonicGravityMetadata()

    metadata.component = variation.HarmonicGravityComponent.ocean_mass
    metadata.gravitational_parameter = 3.986004415e14
    metadata.reference_radius = 6378136.3
    metadata.fully_normalized = True
    metadata.revision = "pinned FES2004 binding extract"
    metadata.loading = variation.HarmonicGravityLoading.included

    settings = [
        variation.iers_2010_solid_body_tide(sampling),
        variation.solid_pole_tide(sampling),
        variation.ocean_pole_tide_from_data(ocean, [0.0, 0.0, -0.3075], constants, 2, 2, sampling),
        variation.harmonic_potential_from_data(
            potential, metadata, 2, 2, ["255.555"], None, sampling
        ),
    ]

    return time, settings


# The Moon's Earth-relative ephemeris also checks deferred origin conversion in the public factory.
def assembled_corrections(time, variations):
    """Retrieve correction blocks from the public factory with deferred ephemeris sampling."""
    settings = environment_setup.BodyListSettings("SSB", "GCRS")

    for name in ("Earth", "Moon", "Sun"):
        settings.add_empty_settings(name)
    settings.get("Earth").ephemeris_settings = ephemeris.constant(
        [1e8, 0.0, 0.0, 0.0, 0.0, 0.0], "SSB", "GCRS"
    )
    settings.get("Moon").ephemeris_settings = ephemeris.constant(
        [4e8, 2e8, -1e8, 0.0, 0.0, 0.0], "Earth", "GCRS"
    )
    settings.get("Sun").ephemeris_settings = ephemeris.constant(
        [1e11, -1.3e11, 2e10, 0.0, 0.0, 0.0], "SSB", "GCRS"
    )
    settings.get("Moon").gravity_field_settings = gravity_field.central(4.902800066e12)
    settings.get("Sun").gravity_field_settings = gravity_field.central(1.32712440018e20)
    cosine = np.zeros((5, 5))

    cosine[0, 0] = 1.0
    field = gravity_field.spherical_harmonic(
        3.986004415e14, 6378136.3, cosine, np.zeros((5, 5)), "ITRS"
    )

    field.tide_system = gravity_data.TideSystem.tide_free
    settings.get("Earth").gravity_field_settings = field
    eop = earth_orientation.read_earth_orientation_data(
        str(TEST_DATA / "src/io/earth_orientation_fixtures/modern_c04.txt")
    )

    settings.get("Earth").rotation_model_settings = rotation_model.gcrs_to_itrs(eop_data=eop)
    settings.get("Earth").gravity_field_variation_settings = variations
    bodies = environment_setup.create_system_of_bodies(settings)

    models = bodies.get("Earth").gravity_field_model.gravity_field_variation_models

    assert len(models) == len(variations)
    corrections = []

    for model in models:
        c, s = model.calculate_spherical_harmonics_corrections(time)

        assert c.shape == (
            model.maximum_degree - model.minimum_degree + 1,
            model.maximum_order - model.minimum_order + 1,
        )
        assert np.isfinite(c).all() and np.isfinite(s).all()
        assert np.linalg.norm(c) + np.linalg.norm(s) > 0.0
        corrections.append((c, s))
    return corrections


def test_conventional_and_harmonic_settings_assemble_corrections():
    """Reader records and consuming settings must cross module boundaries without manual runtime models."""
    time, settings = gravity_settings()

    assembled_corrections(time, settings)


@pytest.mark.skipif(
    not hasattr(variation.GravityFieldVariationSettings, "__setstate__"),
    reason="Tudat serialization is disabled in this build",
)
def test_gravity_settings_pickle_preserves_assembled_corrections():
    """Existing cereal support preserves real provider inputs, final inventory and computed blocks exactly."""
    time, settings = gravity_settings()

    restored = pickle.loads(pickle.dumps(settings))

    assert all(a == b and type(a) is type(b) for a, b in zip(settings, restored))
    before = assembled_corrections(time, settings)

    after = assembled_corrections(time, restored)

    for original, reconstructed in zip(before, after):
        np.testing.assert_array_equal(original[0], reconstructed[0])
        np.testing.assert_array_equal(original[1], reconstructed[1])


def test_aod1b_records_and_time22_final_inventory_settings():
    """Specialized provider loaders must return the existing consuming settings types."""
    data = gravity_data.read_aod1b_gravity_field(
        str(TEST_DATA / "reference_files/aod1b_gravity/AOD1B_2026-01-01_degree2_extract_07.asc")
    )

    data.validate()
    assert data.release == 7
    assert len(data.samples[0].coefficients[(2, 0)].values) == 2
    aod = variation.aod1b_from_data([data], gravity_data.Aod1bComponent.atmosphere, 2, 2)

    assert isinstance(aod, variation.Aod1bGravityFieldVariationSettings)
    assert aod.component == gravity_data.Aod1bComponent.atmosphere
    time22 = variation.time22(
        str(TEST_DATA / "reference_files/time22-v1_1"),
        variation.HarmonicGravityComponent.direct_atmospheric_mass,
        2,
        2,
    )

    assert isinstance(time22, variation.HarmonicGravityFieldVariationSettings)
    assert len(time22.constituents) == 16
    assert time22.metadata.component == variation.HarmonicGravityComponent.direct_atmospheric_mass
    if hasattr(variation.GravityFieldVariationSettings, "__setstate__"):
        assert pickle.loads(pickle.dumps(aod)) == aod
        assert pickle.loads(pickle.dumps(time22)) == time22
