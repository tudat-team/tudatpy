"""Load acceleration settings from contract-driven JSON input."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ._contract import load_settings

ACCELERATION_CONTRACT_PATH = (
    Path(__file__).parent / "contracts" / "propagation_setup" / "acceleration.json"
)


def load_acceleration_settings(
    settings_path: str | Path,
    contract_path: str | Path = ACCELERATION_CONTRACT_PATH,
) -> dict[str, dict[str, list[Any]]]:
    """Validate JSON input and create Tudat acceleration settings objects."""

    return load_settings(settings_path, contract_path)
