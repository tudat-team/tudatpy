"""Load body-list settings from contract-driven JSON input."""

from pathlib import Path

from ._contract import load_settings, read_json_object

BODY_LIST_SETTINGS_CONTRACT_PATH = (
    Path(__file__).parent / "contracts" / "environment_setup" / "body_list_settings.json"
)


def load_body_list_settings(settings_path, contract_path=BODY_LIST_SETTINGS_CONTRACT_PATH):
    """Create Tudat BodyListSettings from a JSON file."""

    document = read_json_object(settings_path, "settings")
    arguments = document.get("body_list_settings", {})
    context = {}
    if isinstance(arguments, dict):
        origin = arguments.get("global_frame_origin")
        orientation = arguments.get("global_frame_orientation")
        if isinstance(origin, str) and isinstance(orientation, str):
            context = {
                "body_settings": {
                    "base_frame_origin": origin,
                    "base_frame_orientation": orientation,
                }
            }
    return load_settings(settings_path, contract_path, factory_context=context)
