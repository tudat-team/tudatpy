"""Load propagator settings from contract-driven JSON input."""

from pathlib import Path

from ._contract import load_settings

SINGLE_ARC_PROPAGATOR_CONTRACT_PATH = (
    Path(__file__).parent
    / "contracts"
    / "propagation_setup"
    / "propagator"
    / "single_arc_propagator_settings.json"
)


def load_translational_propagator_settings(
    settings_path,
    bodies,
    contract_path=SINGLE_ARC_PROPAGATOR_CONTRACT_PATH,
):
    """Create translational settings, injecting bodies for acceleration-model creation."""

    return load_settings(
        settings_path,
        contract_path,
        factory_context={"translational_from_acceleration_settings": {"bodies": bodies}},
    )
