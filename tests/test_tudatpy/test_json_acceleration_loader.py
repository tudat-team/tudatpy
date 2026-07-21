import json
import re
from pathlib import Path

import pytest

from tudatpy.dynamics import environment_setup, propagation_setup
from tudatpy.json_interface import (
    ACCELERATION_CONTRACT_PATH,
    JSONSettingsValidationError,
    load_acceleration_settings,
    load_body_list_settings,
    load_settings,
    load_translational_propagator_settings,
    validate_all_contracts,
    validate_contract_against_module,
)
from tudatpy.json_interface._contract import module_name_from_contract

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
    expected = len(json.loads(ACCELERATION_CONTRACT_PATH.read_text(encoding="utf-8")))
    assert validate_contract_against_module(ACCELERATION_CONTRACT_PATH) == expected


def test_all_contracts_match_tudat_api():
    contract_root = ACCELERATION_CONTRACT_PATH.parents[1]
    expected = sum(
        len(json.loads(path.read_text(encoding="utf-8"))) for path in contract_root.rglob("*.json")
    )
    assert validate_all_contracts() == expected


def test_nested_contracts_are_limited_to_compatible_settings_families():
    contract_root = ACCELERATION_CONTRACT_PATH.parents[1]
    integrators = json.loads(
        (contract_root / "propagation_setup/integrator/integrator_settings.json").read_text(
            encoding="utf-8"
        )
    )
    terminations = json.loads(
        (contract_root / "propagation_setup/propagator/termination_settings.json").read_text(
            encoding="utf-8"
        )
    )

    assert "step_size_validation" not in integrators
    assert set(terminations) == {
        "time_termination",
        "cpu_time_termination",
        "dependent_variable_termination",
        "hybrid_termination",
        "non_sequential_termination",
    }


def test_all_contract_entries_are_in_api_docs():
    repository_root = Path(__file__).parents[2]
    contract_root = repository_root / "src/tudatpy/json_interface/contracts"
    docs_root = repository_root / "docs/tudatpy/source/dynamics"
    pattern = re.compile(
        r"^\.\. autofunction:: tudatpy\.dynamics\."
        r"((?:environment_setup|propagation_setup)\.[\w.]+)$",
        re.MULTILINE,
    )
    documented = {
        name
        for group in ("environment_setup", "propagation_setup")
        for path in [docs_root / f"{group}.rst", *(docs_root / group).glob("*.rst")]
        for name in pattern.findall(path.read_text(encoding="utf-8"))
    }
    contracted = {
        f"{module_name_from_contract(path).removeprefix('tudatpy.dynamics.')}.{factory}"
        for path in contract_root.rglob("*.json")
        for factory in json.loads(path.read_text(encoding="utf-8"))
    }
    assert contracted <= documented


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


def test_factory_with_optional_callable_uses_tudat_default(tmp_path):
    contract_path = (
        ACCELERATION_CONTRACT_PATH.parents[1] / "environment_setup" / "rotation_model.json"
    )
    settings_path = _write_json(
        tmp_path / "rotation.json",
        {
            "orbital_state_direction_based": {
                "central_body": "Earth",
                "is_colinear_with_velocity": True,
                "direction_is_opposite_to_vector": False,
                "base_frame": "J2000",
            }
        },
    )

    result = load_settings(settings_path, contract_path)

    assert isinstance(result, environment_setup.rotation_model.RotationModelSettings)


def test_body_list_settings_support_relative_json_references(tmp_path):
    _write_json(
        tmp_path / "vehicle.json",
        {
            "body_settings": {
                "body_name": "Vehicle",
                "use_default_settings": False,
                "rigid_body_settings": {"constant_rigid_body_properties": {"mass": 2.2}},
            }
        },
    )
    settings_path = _write_json(
        tmp_path / "body_list.json",
        {
            "body_list_settings": {
                "body_settings": {"Vehicle": {"$ref": "vehicle.json"}},
                "global_frame_origin": "SSB",
                "global_frame_orientation": "ECLIPJ2000",
            }
        },
    )

    result = load_body_list_settings(settings_path)

    assert isinstance(result, environment_setup.BodyListSettings)
    assert result.get("Vehicle").rigid_body_settings is not None


def test_translational_settings_create_acceleration_models_from_json(tmp_path):
    body_settings_path = _write_json(
        tmp_path / "bodies.json",
        {
            "body_list_settings": {
                "body_settings": {
                    "Earth": {
                        "body_settings": {
                            "body_name": "Earth",
                            "use_default_settings": False,
                            "gravity_field_settings": {
                                "central": {"gravitational_parameter": 3.986004418e14}
                            },
                        }
                    },
                    "Vehicle": {
                        "body_settings": {
                            "body_name": "Vehicle",
                            "use_default_settings": False,
                            "rigid_body_settings": {
                                "constant_rigid_body_properties": {"mass": 2.2}
                            },
                        }
                    },
                },
                "global_frame_origin": "SSB",
                "global_frame_orientation": "ECLIPJ2000",
            }
        },
    )
    bodies = environment_setup.create_system_of_bodies(load_body_list_settings(body_settings_path))
    propagator_path = _write_json(
        tmp_path / "propagator.json",
        {
            "translational_from_acceleration_settings": {
                "central_bodies": ["Earth"],
                "acceleration_settings": {"Vehicle": {"Earth": [{"point_mass_gravity": {}}]}},
                "bodies_to_integrate": ["Vehicle"],
                "initial_states": [7000000, 0, 0, 0, 7500, 0],
                "initial_time": 0,
                "integrator_settings": {
                    "runge_kutta_fixed_step": {
                        "time_step": 10,
                        "coefficient_set": "rk_4",
                    }
                },
                "termination_settings": {"time_termination": {"termination_time": 100}},
            }
        },
    )

    result = load_translational_propagator_settings(propagator_path, bodies)

    assert isinstance(
        result,
        propagation_setup.propagator.TranslationalStatePropagatorSettings,
    )


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        ({"unknown_acceleration": {}}, "unknown settings factory"),
        ({"maximum_degree": 5}, "missing required properties: maximum_order"),
        (
            {"maximum_degree": 5, "maximum_order": 5, "extra": 1},
            "contains undeclared properties: extra",
        ),
        (
            {"maximum_degree": 5.0, "maximum_order": 5},
            "maximum_degree must be an int",
        ),
        (
            {"maximum_degree": None, "maximum_order": 5},
            "maximum_degree must not be null",
        ),
    ],
)
def test_rejects_invalid_acceleration_settings(tmp_path, arguments, message):
    definition = (
        arguments
        if "unknown_acceleration" in arguments
        else {"spherical_harmonic_gravity": arguments}
    )
    settings_path = _write_json(
        tmp_path / "accelerations.json",
        {"Delfi-C3": {"Earth": [definition]}},
    )

    with pytest.raises(JSONSettingsValidationError, match=message):
        load_acceleration_settings(settings_path)
