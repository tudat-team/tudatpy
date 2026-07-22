"""Contract, loader, and public-API regression tests for the JSON interface."""

import json
import re
import sys
from pathlib import Path
from types import ModuleType

import pytest

from tudatpy import astro, math as tudat_math
from tudatpy.dynamics import environment_setup, propagation_setup
from tudatpy.kernel.dynamics import environment_setup as kernel_environment_setup
from tudatpy.kernel.dynamics.propagation_setup import propagator as kernel_propagator
from tudatpy.json_interface import (
    ACCELERATION_CONTRACT_PATH,
    CONTRACT_ROOT,
    JSONSettingsValidationError,
    load_acceleration_settings,
    load_body_list_settings,
    load_single_arc_propagator_settings,
    load_settings,
    load_translational_propagator_settings,
    validate_all_contracts,
    validate_contract_against_module,
)
from tudatpy.json_interface._contract import module_name_from_contract

PERTURBED_SATELLITE_ACCELERATIONS = {
    "Delfi-C3": {
        "Sun": ["radiation_pressure", "point_mass_gravity"],
        "Earth": [
            {
                "spherical_harmonic_gravity": {
                    "maximum_degree": 5,
                    "maximum_order": 5,
                }
            },
            "aerodynamic",
        ],
        "Moon": "point_mass_gravity",
        "Mars": {"point_mass_gravity": {}},
        "Venus": {"point_mass_gravity": None},
    }
}


def _write_json(path, value):
    """Write a value as a temporary JSON document.

    Parameters
    ----------
    path : pathlib.Path
        Destination file.
    value : Any
        JSON-serializable value to write.

    Returns
    -------
    pathlib.Path
        The supplied destination path, allowing creation and use in one expression.
    """

    path.write_text(json.dumps(value), encoding="utf-8")
    return path


def test_acceleration_contract_matches_tudat_api():
    """Keep every acceleration contract entry aligned with the exposed API.

    Returns
    -------
    None
        The test asserts that every acceleration entry validates successfully.
    """

    expected = len(json.loads(ACCELERATION_CONTRACT_PATH.read_text(encoding="utf-8")))
    assert validate_contract_against_module(ACCELERATION_CONTRACT_PATH) == expected


def test_all_contracts_match_tudat_api():
    """Keep all packaged contracts aligned with their corresponding modules.

    Returns
    -------
    None
        The test asserts that the validated count equals the packaged entry count.
    """

    contract_root = ACCELERATION_CONTRACT_PATH.parents[1]
    expected = sum(
        len(json.loads(path.read_text(encoding="utf-8"))) for path in contract_root.rglob("*.json")
    )
    assert validate_all_contracts() == expected


def test_nested_contracts_are_limited_to_compatible_settings_families():
    """Restrict nested inputs to compatible settings families.

    Returns
    -------
    None
        The test asserts the exact permitted termination family and excludes a
        step-size validation factory from the integrator-settings family.
    """

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
    """Require every contracted factory to appear in public API documentation.

    Returns
    -------
    None
        The test asserts that contracted qualified names are a subset of the
        documented factory names.
    """

    repository_root = Path(__file__).parents[2]
    contract_root = repository_root / "src/tudatpy/json_interface/contracts"
    docs_root = repository_root / "docs/tudatpy/source"
    pattern = re.compile(
        r"^\.\. autofunction:: (tudatpy\.(?:dynamics|math|astro)\.[\w.]+)$",
        re.MULTILINE,
    )
    documented = {
        name
        for path in docs_root.rglob("*.rst")
        for name in pattern.findall(path.read_text(encoding="utf-8"))
    }
    contracted = {
        f"{module_name_from_contract(path)}.{factory}"
        for path in contract_root.rglob("*.json")
        for factory in json.loads(path.read_text(encoding="utf-8"))
    }
    assert contracted <= documented


def test_json_composition_factories_are_public_tudat_api():
    """Ensure JSON composition uses ordinary kernel-exposed Tudat factories.

    Returns
    -------
    None
        The test asserts object identity between public imports and kernel
        bindings for every added composition factory.
    """

    assert environment_setup.body_settings is kernel_environment_setup.body_settings
    assert environment_setup.body_list_settings is kernel_environment_setup.body_list_settings
    assert (
        propagation_setup.propagator.translational_from_acceleration_settings
        is kernel_propagator.translational_from_acceleration_settings
    )
    assert (
        propagation_setup.propagator.rotational_from_torque_settings
        is kernel_propagator.rotational_from_torque_settings
    )
    assert (
        propagation_setup.propagator.mass_from_mass_rate_settings
        is kernel_propagator.mass_from_mass_rate_settings
    )
    for factory_name in (
        "propagation_print_settings",
        "single_arc_processing_settings",
        "multi_arc_processing_settings",
        "hybrid_arc_processing_settings",
    ):
        assert getattr(propagation_setup.propagator, factory_name) is getattr(
            kernel_propagator, factory_name
        )


def test_external_settings_types_load_from_their_public_modules(tmp_path):
    """Create interpolation, root-finder, and position-element settings from JSON.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Pytest-managed directory containing the temporary settings documents.

    Returns
    -------
    None
        The test asserts that each external contract creates the corresponding
        public Tudat type and that position-element enum conversion is accepted.
    """

    interpolator = load_settings(
        _write_json(
            tmp_path / "interpolator.json", {"lagrange_interpolation": {"number_of_points": 6}}
        ),
        CONTRACT_ROOT / "math/interpolators/interpolator_settings.json",
    )
    assert isinstance(interpolator, tudat_math.interpolators.InterpolatorSettings)

    root_finder = load_settings(
        _write_json(
            tmp_path / "root_finder.json",
            {"bisection": {"absolute_variable_tolerance": 1.0e-12}},
        ),
        CONTRACT_ROOT / "math/root_finders/root_finder_settings.json",
    )
    assert isinstance(root_finder, tudat_math.root_finders.RootFinderSettings)

    station = load_settings(
        _write_json(
            tmp_path / "station.json",
            {
                "basic_station": {
                    "station_name": "Test station",
                    "station_nominal_position": [6378137.0, 0.0, 0.0],
                    "station_position_element_type": "cartesian_position_type",
                }
            },
        ),
        CONTRACT_ROOT / "environment_setup/ground_station/ground_station_settings.json",
    )
    assert isinstance(station, environment_setup.ground_station.GroundStationSettings)
    assert astro.element_conversion.PositionElementTypes.cartesian_position_type is not None


@pytest.mark.parametrize(
    "definition",
    [
        "linear_interpolation",
        {"linear_interpolation": {}},
        {"linear_interpolation": None},
    ],
)
def test_no_argument_factory_shorthand_is_generic(tmp_path, definition):
    """Accept every no-argument spelling through the generic settings loader.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Pytest-managed directory containing the temporary settings document.
    definition : str or dict[str, object]
        One of the equivalent bare-string, empty-object, or null factory forms.

    Returns
    -------
    None
        The test asserts that each generic spelling creates an interpolator
        setting without using acceleration-specific loading code.
    """

    result = load_settings(
        _write_json(tmp_path / "interpolator.json", definition),
        CONTRACT_ROOT / "math/interpolators/interpolator_settings.json",
    )

    assert isinstance(result, tudat_math.interpolators.InterpolatorSettings)


def test_unknown_bare_factory_name_is_rejected(tmp_path):
    """Reject an unknown root string instead of returning unconverted JSON.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Pytest-managed directory containing the temporary settings document.

    Returns
    -------
    None
        The test succeeds only when the generic loader reports the unknown
        factory name against the selected contract.
    """

    settings_path = _write_json(tmp_path / "unknown.json", "unknown_factory")

    with pytest.raises(JSONSettingsValidationError, match="unknown settings factory"):
        load_settings(
            settings_path,
            CONTRACT_ROOT / "math/interpolators/interpolator_settings.json",
        )


def test_singleton_array_shorthand_is_generic(tmp_path):
    """Accept one factory object where a contracted settings array is expected.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Pytest-managed directory containing the temporary termination document.

    Returns
    -------
    None
        The test asserts that a single termination setting is converted to the
        array required by ``hybrid_termination`` without acceleration-specific
        handling.
    """

    result = load_settings(
        _write_json(
            tmp_path / "termination.json",
            {
                "hybrid_termination": {
                    "termination_settings": {"time_termination": {"termination_time": 100.0}},
                    "fulfill_single_condition": True,
                }
            },
        ),
        CONTRACT_ROOT / "propagation_setup/propagator/termination_settings.json",
    )

    assert isinstance(
        result,
        propagation_setup.propagator.PropagationTerminationSettings,
    )


def test_processing_settings_load_from_json(tmp_path):
    """Create all propagation processing-settings variants from JSON.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Pytest-managed directory containing temporary settings documents.

    Returns
    -------
    None
        The test asserts concrete public types for single-, multi-, and
        hybrid-arc processing settings, including nested print settings.
    """

    definitions = (
        (
            "single_arc_processing_settings.json",
            {
                "single_arc_processing_settings": {
                    "results_save_frequency_in_steps": 4,
                    "print_settings": {
                        "propagation_print_settings": {"print_termination_reason": True}
                    },
                }
            },
            propagation_setup.propagator.SingleArcPropagatorProcessingSettings,
        ),
        (
            "multi_arc_processing_settings.json",
            {"multi_arc_processing_settings": {"print_output_on_first_arc_only": True}},
            propagation_setup.propagator.MultiArcPropagatorProcessingSettings,
        ),
        (
            "hybrid_arc_processing_settings.json",
            {"hybrid_arc_processing_settings": {"print_state_type_start": True}},
            propagation_setup.propagator.HybridArcPropagatorProcessingSettings,
        ),
    )

    for contract_name, definition, expected_type in definitions:
        result = load_settings(
            _write_json(tmp_path / contract_name, definition),
            CONTRACT_ROOT / "propagation_setup/propagator" / contract_name,
        )
        assert isinstance(result, expected_type)


def test_body_settings_contract_uses_rigid_body_settings_for_mass():
    """Represent body mass through rigid-body settings instead of a special input.

    Returns
    -------
    None
        The test asserts that ``constant_mass`` is absent and the rigid-body
        property points to its nested contract.
    """

    contract = json.loads(
        (
            ACCELERATION_CONTRACT_PATH.parents[1] / "environment_setup" / "body_settings.json"
        ).read_text(encoding="utf-8")
    )["body_settings"]

    assert "constant_mass" not in contract["properties"]
    assert contract["properties"]["rigid_body_settings"] == "environment_setup.rigid_body"


def test_contract_defaults_control_factory_arguments(tmp_path):
    """Verify that finite contract defaults are passed to the selected factory.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Pytest-managed directory for the temporary contract and settings files.

    Returns
    -------
    None
        The test asserts that the contract default overrides the test factory's
        Python default.
    """

    module_name = "_tudatpy_json_contract_default_test"
    module = ModuleType(module_name)
    module.settings = lambda value=1: value
    sys.modules[module_name] = module
    try:
        contract_path = _write_json(
            tmp_path / "contract.json",
            {
                "settings": {
                    "properties": {"value": "int"},
                    "optional": {"value": 7},
                }
            },
        )
        settings_path = _write_json(tmp_path / "settings.json", {"settings": {}})

        assert load_settings(settings_path, contract_path, module_name=module_name) == 7
    finally:
        del sys.modules[module_name]


def test_unsupported_contract_input_is_rejected(tmp_path):
    """Reject a required input that the JSON interface cannot construct.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Pytest-managed directory for the temporary contract and settings files.

    Returns
    -------
    None
        The test succeeds only when loading raises the expected validation error.
    """

    module_name = "_tudatpy_json_unsupported_input_test"
    module = ModuleType(module_name)
    module.settings = lambda external_input: external_input
    sys.modules[module_name] = module
    try:
        contract_path = _write_json(
            tmp_path / "contract.json",
            {
                "settings": {
                    "properties": {"external_input": "any"},
                    "unsupported": ["external_input"],
                }
            },
        )
        settings_path = _write_json(
            tmp_path / "settings.json",
            {"settings": {"external_input": {"value": 1}}},
        )

        with pytest.raises(JSONSettingsValidationError, match="not yet supported"):
            load_settings(settings_path, contract_path, module_name=module_name)
    finally:
        del sys.modules[module_name]


def test_optional_unsupported_input_uses_factory_default(tmp_path):
    """Confirm an omitted unsupported optional input uses the factory default.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Pytest-managed directory for the temporary contract and settings files.

    Returns
    -------
    None
        The test asserts that no contract placeholder is passed to the factory.
    """

    module_name = "_tudatpy_json_optional_unsupported_input_test"
    module = ModuleType(module_name)
    module.settings = lambda external_input="factory default": external_input
    sys.modules[module_name] = module
    try:
        contract_path = _write_json(
            tmp_path / "contract.json",
            {
                "settings": {
                    "properties": {"external_input": "ExternalType"},
                    "optional": {"external_input": "contract default"},
                    "unsupported": ["external_input"],
                }
            },
        )
        settings_path = _write_json(tmp_path / "settings.json", {"settings": {}})

        assert (
            load_settings(settings_path, contract_path, module_name=module_name)
            == "factory default"
        )
    finally:
        del sys.modules[module_name]


def test_untyped_inputs_are_flagged_unsupported_and_coma_is_excluded():
    """Require unsupported markers for untyped inputs and exclude Coma factories.

    Returns
    -------
    None
        The test asserts that every ``any`` property is explicitly unsupported
        and that Coma factories are absent from the JSON contracts.
    """

    contract_root = ACCELERATION_CONTRACT_PATH.parents[1]
    for path in contract_root.rglob("*.json"):
        for entry in json.loads(path.read_text(encoding="utf-8")).values():
            unsupported = set(entry.get("unsupported", []))
            any_properties = {
                name
                for name, type_name in entry["properties"].items()
                if type_name.startswith("any")
            }
            assert any_properties <= unsupported

    atmosphere = json.loads(
        (contract_root / "environment_setup/atmosphere/atmosphere_settings.json").read_text(
            encoding="utf-8"
        )
    )
    wind = json.loads(
        (contract_root / "environment_setup/atmosphere/wind_model_settings.json").read_text(
            encoding="utf-8"
        )
    )
    assert "coma_model_from_poly_data" not in atmosphere
    assert "coma_model_from_stokes_data" not in atmosphere
    assert "coma_wind_model" not in wind


def test_loads_perturbed_satellite_accelerations_with_tudat_api(tmp_path):
    """Create the perturbed-satellite acceleration map with real Tudat factories.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Pytest-managed directory containing the generated acceleration document.

    Returns
    -------
    None
        The test asserts the nested mapping shape and concrete Tudat settings
        types returned by the loader.
    """

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
    assert all(isinstance(settings, list) for settings in result["Delfi-C3"].values())
    assert len(result["Delfi-C3"]["Sun"]) == 2
    assert all(len(result["Delfi-C3"][body]) == 1 for body in ("Moon", "Mars", "Venus"))
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
    """Omit an optional callable so Tudat retains its own implementation.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Pytest-managed directory containing the rotation-model document.

    Returns
    -------
    None
        The test asserts successful construction without a JSON callable value.
    """

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
    """Resolve a linked body definition relative to its enclosing JSON file.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Pytest-managed directory containing the root and referenced JSON files.

    Returns
    -------
    None
        The test asserts that the link produces body-list and rigid-body settings.
    """

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
    """Build acceleration models while loading translational settings.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Pytest-managed directory containing temporary environment and propagator
        JSON documents.

    Returns
    -------
    None
        The test asserts the configured vehicle mass and the returned Tudat
        translational propagator-settings type.
    """

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
    assert bodies.get("Vehicle").mass == pytest.approx(2.2)
    propagator_path = _write_json(
        tmp_path / "propagator.json",
        {
            "translational_from_acceleration_settings": {
                "central_bodies": ["Earth"],
                "acceleration_settings": {"Vehicle": {"Earth": "point_mass_gravity"}},
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


def test_rotational_and_mass_settings_create_models_from_json(tmp_path):
    """Build rotational and mass propagator settings from model settings.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Pytest-managed directory containing temporary environment and propagator
        JSON documents.

    Returns
    -------
    None
        The test asserts that both from-settings composition factories return
        their concrete public propagator-settings types.
    """

    body_settings_path = _write_json(
        tmp_path / "bodies.json",
        {
            "body_list_settings": {
                "body_settings": {
                    "Vehicle": {
                        "body_settings": {
                            "body_name": "Vehicle",
                            "use_default_settings": False,
                            "rigid_body_settings": {
                                "constant_rigid_body_properties": {"mass": 2.2}
                            },
                        }
                    }
                },
                "global_frame_origin": "SSB",
                "global_frame_orientation": "ECLIPJ2000",
            }
        },
    )
    bodies = environment_setup.create_system_of_bodies(load_body_list_settings(body_settings_path))
    common = {
        "initial_time": 0.0,
        "integrator_settings": {
            "runge_kutta_fixed_step": {"time_step": 10.0, "coefficient_set": "rk_4"}
        },
        "termination_settings": {"time_termination": {"termination_time": 100.0}},
    }

    rotational_path = _write_json(
        tmp_path / "rotational.json",
        {
            "rotational_from_torque_settings": {
                "torque_settings": {},
                "bodies_to_integrate": ["Vehicle"],
                "initial_states": [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                **common,
            }
        },
    )
    rotational = load_single_arc_propagator_settings(rotational_path, bodies)
    assert isinstance(rotational, propagation_setup.propagator.RotationalStatePropagatorSettings)

    mass_path = _write_json(
        tmp_path / "mass.json",
        {
            "mass_from_mass_rate_settings": {
                "bodies_with_mass_to_propagate": ["Vehicle"],
                "mass_rate_settings": {},
                "initial_body_masses": [2.2],
                **common,
            }
        },
    )
    mass = load_single_arc_propagator_settings(mass_path, bodies)
    assert isinstance(mass, propagation_setup.propagator.MassPropagatorSettings)


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
    """Reject malformed acceleration factory input.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Pytest-managed directory containing the generated settings document.
    arguments : dict[str, Any]
        Parameterized invalid factory definition or argument mapping.
    message : str
        Expected fragment of the resulting validation error.

    Returns
    -------
    None
        The test succeeds only when loading raises the matching validation error.
    """

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
