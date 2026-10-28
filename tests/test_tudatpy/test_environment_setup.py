import numpy as np
import pytest

from tudatpy.dynamics import environment, environment_setup
from tudatpy.interface import spice


def test_deprecated_mass_properties_alias_preserves_behavior(capfd):
    body_settings = environment_setup.BodyListSettings("SSB", "J2000")
    body_settings.add_empty_settings("Canonical")
    body_settings.add_empty_settings("Legacy")
    bodies = environment_setup.create_system_of_bodies(body_settings)
    settings = environment_setup.rigid_body.constant_rigid_body_properties(125.0)

    environment_setup.add_rigid_body_properties(
        bodies=bodies,
        body_name="Canonical",
        rigid_body_property_settings=settings,
    )
    assert capfd.readouterr().err == ""

    environment_setup.add_mass_properties_model(
        bodies=bodies,
        body_name="Legacy",
        mass_property_settings=settings,
    )
    warning = capfd.readouterr().err
    assert "add_mass_properties_model is deprecated" in warning
    assert "tudatpy.dynamics.environment_setup.add_rigid_body_properties" in warning
    assert bodies.get("Legacy").mass == pytest.approx(bodies.get("Canonical").mass)
    assert bodies.get("Legacy").mass == pytest.approx(125.0)


def test_sgp4_from_tle_elements_matches_tle_lines():
    spice.load_standard_kernels()
    tle_from_lines = environment.Tle(
        "1 32789U 07021G   08119.60740078 -.00000054  00000-0  00000+0 0  9999",
        "2 32789 098.0082 179.6267 0015321 307.2977 051.0656 14.81417433    68",
    )
    tle_from_elements = environment.Tle(
        reference_epoch=tle_from_lines.reference_epoch,
        b_star=tle_from_lines.b_star,
        inclination=tle_from_lines.inclination,
        right_ascension=tle_from_lines.right_ascension,
        eccentricity=tle_from_lines.eccentricity,
        argument_of_perigee=tle_from_lines.argument_of_perigee,
        mean_anomaly=tle_from_lines.mean_anomaly,
        mean_motion=tle_from_lines.mean_motion,
    )
    ephemeris_settings = environment_setup.ephemeris.sgp4(tle_from_elements)

    epoch = tle_from_lines.reference_epoch + 3600.0
    state_from_lines = environment.TleEphemeris(tle=tle_from_lines).cartesian_state(epoch)
    state_from_elements = environment.TleEphemeris(tle=ephemeris_settings.tle).cartesian_state(
        epoch
    )
    np.testing.assert_allclose(state_from_elements, state_from_lines, rtol=0.0, atol=1.0e-6)
