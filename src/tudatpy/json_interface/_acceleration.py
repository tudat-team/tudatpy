"""Load acceleration settings from contract-driven JSON input."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ._contract import (
    JSONSettingsValidationError,
    create_settings_object,
    expect_object,
    load_contract,
    read_json_object,
)

ACCELERATION_CONTRACT_PATH = (
    Path(__file__).parent / "contracts" / "propagation_setup" / "acceleration.json"
)


def load_acceleration_settings(
    settings_path: str | Path,
    contract_path: str | Path = ACCELERATION_CONTRACT_PATH,
) -> dict[str, dict[str, list[Any]]]:
    """Validate JSON input and create Tudat acceleration settings objects."""

    from tudatpy.dynamics import environment_setup, propagation_setup

    acceleration_module = propagation_setup.acceleration
    type_modules = (environment_setup, propagation_setup)
    contract = load_contract(contract_path, acceleration_module, type_modules)
    document = read_json_object(settings_path, "accelerations")

    result: dict[str, dict[str, list[Any]]] = {}
    for body_undergoing, exerting_bodies_value in document.items():
        exerting_bodies = expect_object(exerting_bodies_value, f"accelerations.{body_undergoing}")
        result[body_undergoing] = {}

        for body_exerting, definitions in exerting_bodies.items():
            path = f"accelerations.{body_undergoing}.{body_exerting}"
            if not isinstance(definitions, list):
                raise JSONSettingsValidationError(
                    f"{path} must be an array, got {type(definitions).__name__}"
                )

            result[body_undergoing][body_exerting] = [
                create_settings_object(
                    definition,
                    contract,
                    acceleration_module,
                    type_modules,
                    f"{path}[{index}]",
                )
                for index, definition in enumerate(definitions)
            ]

    return result
