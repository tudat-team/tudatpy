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


def load_single_arc_propagator_settings(
    settings_path,
    bodies,
    contract_path=SINGLE_ARC_PROPAGATOR_CONTRACT_PATH,
):
    """Create single-arc propagator settings from a JSON document.

    Parameters
    ----------
    settings_path : str or pathlib.Path
        Path to a JSON definition using a supported single-arc propagator
        factory. Translational, rotational, and mass propagation can supply
        settings rather than already-created runtime models.
    bodies : tudatpy.dynamics.environment.SystemOfBodies
        Runtime body system used by Tudat to create acceleration, torque, or
        mass-rate models from the contracted settings.
    contract_path : str or pathlib.Path, optional
        Contract governing single-arc propagator settings. The packaged
        single-arc contract is used by default.

    Returns
    -------
    tudatpy.dynamics.propagation_setup.propagator.SingleArcPropagatorSettings
        Fully composed translational, rotational, mass, or multitype propagator
        settings ready for ``simulator.create_dynamics_simulator``.

    Raises
    ------
    JSONSettingsValidationError
        If JSON validation or nested factory construction fails.
    RuntimeError
        If Tudat cannot create the requested models from the supplied body
        system and contracted settings.
    """

    # SystemOfBodies contains live ephemerides, frames, and environment models,
    # so serializing it would violate the settings-only role of this interface.
    # It is injected only into public composition factories, which convert the
    # contracted settings into runtime models and delegate to the corresponding
    # ordinary propagator factories.
    return load_settings(
        settings_path,
        contract_path,
        factory_context={
            factory_name: {"bodies": bodies}
            for factory_name in (
                "translational_from_acceleration_settings",
                "rotational_from_torque_settings",
                "mass_from_mass_rate_settings",
            )
        },
    )


def load_translational_propagator_settings(
    settings_path,
    bodies,
    contract_path=SINGLE_ARC_PROPAGATOR_CONTRACT_PATH,
):
    """Create translational propagator settings from a JSON document.

    Parameters
    ----------
    settings_path : str or pathlib.Path
        Path to a JSON ``translational_from_acceleration_settings`` definition.
    bodies : tudatpy.dynamics.environment.SystemOfBodies
        Runtime body system used to create acceleration models.
    contract_path : str or pathlib.Path, optional
        Contract governing single-arc propagator settings.

    Returns
    -------
    tudatpy.dynamics.propagation_setup.propagator.TranslationalStatePropagatorSettings
        Fully composed translational propagator settings.

    Raises
    ------
    JSONSettingsValidationError
        If JSON validation or nested factory construction fails.
    RuntimeError
        If Tudat cannot create acceleration models from the supplied settings.
    """

    return load_single_arc_propagator_settings(settings_path, bodies, contract_path)
