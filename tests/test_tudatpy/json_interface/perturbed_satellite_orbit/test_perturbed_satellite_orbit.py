from pathlib import Path

from tudatpy.dynamics import environment_setup, simulator
from tudatpy.interface import spice
from tudatpy.json_interface import (
    load_body_list_settings,
    load_translational_propagator_settings,
)

settings_directory = Path(__file__).parent


def test_perturbed_satellite_orbit_from_json():
    spice.load_standard_kernels()

    body_list_settings = load_body_list_settings(settings_directory / "body_list_settings.json")
    bodies = environment_setup.create_system_of_bodies(body_list_settings)

    propagator_settings = load_translational_propagator_settings(
        settings_directory / "propagator_settings.json", bodies
    )

    dynamics_simulator = simulator.create_dynamics_simulator(bodies, propagator_settings)

    assert len(dynamics_simulator.propagation_results.state_history) == 8641
