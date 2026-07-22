"""Load body-list settings from contract-driven JSON input."""

from pathlib import Path

from ._contract import load_settings, read_json_object

BODY_LIST_SETTINGS_CONTRACT_PATH = (
    Path(__file__).parent / "contracts" / "environment_setup" / "body_list_settings.json"
)


def load_body_list_settings(settings_path, contract_path=BODY_LIST_SETTINGS_CONTRACT_PATH):
    """Create Tudat body-list settings from linked JSON environment models.

    Parameters
    ----------
    settings_path : str or pathlib.Path
        Path to a JSON document containing one ``body_list_settings`` factory.
        Per-body definitions and their model settings may be included through
        relative ``$ref`` objects.
    contract_path : str or pathlib.Path, optional
        Contract governing the body-list document. The packaged body-list
        contract is used by default.

    Returns
    -------
    tudatpy.dynamics.environment_setup.BodyListSettings
        Settings ready for
        ``environment_setup.create_system_of_bodies``.

    Raises
    ------
    JSONSettingsValidationError
        If the document or a linked file is invalid, does not satisfy its
        contract, or a contracted Tudat factory rejects its arguments.
    """

    # A default body is constructed before the enclosing BodyListSettings
    # object exists. Its default ephemeris and rotation settings must already
    # use the global frame declared at the body-list level. The loader therefore
    # extracts those two values and injects them into every nested body_settings
    # call as context instead of duplicating them in each body JSON document.
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

    # The second read is intentional: load_settings owns the complete generic
    # validation/construction path, while the first read above only determines
    # wrapper-specific runtime context.
    return load_settings(settings_path, contract_path, factory_context=context)
