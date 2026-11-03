import pytest

from tudatpy.dynamics import environment_setup


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
