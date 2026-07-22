"""Load acceleration settings from contract-driven JSON input."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ._contract import JSONSettingsValidationError, load_settings

ACCELERATION_CONTRACT_PATH = (
    Path(__file__).parent / "contracts" / "propagation_setup" / "acceleration.json"
)


def load_acceleration_settings(
    settings_path: str | Path,
    contract_path: str | Path = ACCELERATION_CONTRACT_PATH,
) -> dict[str, dict[str, list[Any]]]:
    """Create a Tudat acceleration-settings mapping from a JSON document.

    Parameters
    ----------
    settings_path : str or pathlib.Path
        Path to the JSON document. Its outer keys name the bodies undergoing
        acceleration, the next level names the bodies exerting acceleration,
        Each value is either a list of acceleration factory definitions or one
        definition using the singleton-list shorthand. A zero-input factory may
        be written as its bare name, as ``{"factory": null}``, or as
        ``{"factory": {}}``.
    contract_path : str or pathlib.Path, optional
        Contract used to validate and convert every acceleration definition.
        The packaged acceleration contract is used by default.

    Returns
    -------
    dict[str, dict[str, list[Any]]]
        Nested mapping accepted by
        ``propagation_setup.create_acceleration_models``. Every list item is a
        Tudat ``AccelerationSettings`` instance created by the contracted
        acceleration factory.

    Raises
    ------
    JSONSettingsValidationError
        If either JSON file cannot be read, a reference is invalid, the input
        does not satisfy the contract, or Tudat rejects a factory invocation.
    """

    # Acceleration documents have a domain-specific outer mapping, but their
    # leaves use the same one-key factory representation as every other
    # contract. The generic recursive loader can therefore construct the full
    # tree; this wrapper only fixes the appropriate contract and public return
    # shape for callers.
    converted = load_settings(settings_path, contract_path)

    # The generic converter preserves the user's singleton shorthand. Tudat's
    # selected-acceleration map always requires lists, so normalize only this
    # domain-specific outer structure after all factory calls have succeeded.
    from tudatpy.dynamics.propagation_setup.acceleration import AccelerationSettings

    normalized: dict[str, dict[str, list[Any]]] = {}
    for body_undergoing, exerting_bodies in converted.items():
        if not isinstance(exerting_bodies, dict):
            raise JSONSettingsValidationError(
                f"settings.{body_undergoing} must map exerting bodies to accelerations"
            )
        normalized[body_undergoing] = {}
        for body_exerting, settings in exerting_bodies.items():
            settings_list = settings if isinstance(settings, list) else [settings]
            if any(not isinstance(setting, AccelerationSettings) for setting in settings_list):
                raise JSONSettingsValidationError(
                    f"settings.{body_undergoing}.{body_exerting} must contain "
                    "acceleration settings factories"
                )
            normalized[body_undergoing][body_exerting] = settings_list
    return normalized
