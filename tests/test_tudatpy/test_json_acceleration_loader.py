import json

import pytest

from tudatpy.dynamics import propagation_setup
from tudatpy.json_interface import (
    ACCELERATION_CONTRACT_PATH,
    JSONSettingsValidationError,
    load_acceleration_settings,
    validate_contract_against_module,
)

PERTURBED_SATELLITE_ACCELERATIONS = {
    "Delfi-C3": {
        "Sun": [
            {"radiation_pressure": {}},
            {"point_mass_gravity": {}},
        ],
        "Earth": [
            {
                "spherical_harmonic_gravity": {
                    "maximum_degree": 5,
                    "maximum_order": 5,
                }
            },
            {"aerodynamic": {}},
        ],
        "Moon": [{"point_mass_gravity": {}}],
        "Mars": [{"point_mass_gravity": {}}],
        "Venus": [{"point_mass_gravity": {}}],
    }
}


def _write_json(path, value):
    path.write_text(json.dumps(value), encoding="utf-8")
    return path


def test_acceleration_contract_matches_tudat_api():
    assert (
        validate_contract_against_module(
            ACCELERATION_CONTRACT_PATH,
            "tudatpy.dynamics.propagation_setup.acceleration",
            (
                "tudatpy.dynamics.environment_setup",
                "tudatpy.dynamics.propagation_setup",
            ),
        )
        == 21
    )


def test_loads_perturbed_satellite_accelerations_with_tudat_api(tmp_path):
    settings_path = _write_json(tmp_path / "accelerations.json", PERTURBED_SATELLITE_ACCELERATIONS)

    result = load_acceleration_settings(settings_path)

    assert result.keys() == {"Delfi-C3"}
    assert result["Delfi-C3"].keys() == {
        "Sun",
        "Earth",
        "Moon",
        "Mars",
        "Venus",
    }
    assert all(
        isinstance(setting, propagation_setup.acceleration.AccelerationSettings)
        for settings in result["Delfi-C3"].values()
        for setting in settings
    )
    assert isinstance(
        result["Delfi-C3"]["Earth"][0],
        propagation_setup.acceleration.SphericalHarmonicAccelerationSettings,
    )


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        ({"maximum_degree": 5}, "missing required properties: maximum_order"),
        (
            {"maximum_degree": 5, "maximum_order": 5, "extra": 1},
            "contains undeclared properties: extra",
        ),
        (
            {"maximum_degree": 5.0, "maximum_order": 5},
            "maximum_degree must be an int",
        ),
    ],
)
def test_rejects_invalid_acceleration_settings(tmp_path, arguments, message):
    settings_path = _write_json(
        tmp_path / "accelerations.json",
        {"Delfi-C3": {"Earth": [{"spherical_harmonic_gravity": arguments}]}},
    )

    with pytest.raises(JSONSettingsValidationError, match=message):
        load_acceleration_settings(settings_path)
