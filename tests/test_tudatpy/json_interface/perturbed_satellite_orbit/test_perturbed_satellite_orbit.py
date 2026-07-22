"""End-to-end propagation using only contract-governed JSON settings."""

from pathlib import Path

from tudatpy.dynamics import environment_setup, simulator
from tudatpy.interface import spice
from tudatpy.json_interface import (
    load_body_list_settings,
    load_translational_propagator_settings,
)

settings_directory = Path(__file__).parent


def test_perturbed_satellite_orbit_from_json():
    """Reproduce the perturbed Delfi-C3 propagation from top-level JSON files.

    Returns
    -------
    None
        The test succeeds after the simulator produces all 8,641 expected state
        samples and otherwise fails through an exception or assertion.
    """

    # Kernel loading initializes runtime ephemeris data; it does not construct a
    # Tudat settings object and therefore remains an explicit runtime operation.
    spice.load_standard_kernels()

    # This file links the individual default bodies and the single file holding
    # every Delfi-C3 environment setting.
    body_list_settings = load_body_list_settings(settings_directory / "body_list_settings.json")
    bodies = environment_setup.create_system_of_bodies(body_list_settings)

    # The bodies provide runtime context for converting the linked acceleration
    # settings into acceleration models. All other propagation settings are
    # included from JSON by propagator_settings.json.
    propagator_settings = load_translational_propagator_settings(
        settings_directory / "propagator_settings.json", bodies
    )

    # Simulator creation is the execution step; every required settings object
    # has already been produced through the JSON interface.
    dynamics_simulator = simulator.create_dynamics_simulator(bodies, propagator_settings)

    assert len(dynamics_simulator.propagation_results.state_history) == 8641
