from tudatpy.numerical_simulation.environment_setup import (
    aerodynamic_coefficients,
    radiation_pressure,
)
from tudatpy.numerical_simulation.propagation_setup import (
    acceleration,
)
from tudatpy.numerical_simulation import propagation_setup
from tudatpy.numerical_simulation import environment_setup
from tudatpy.interface import spice
from tudatpy import constants


def test_basic_environment_setup() -> None:

    spice.load_standard_kernels()

    simulation_start_epoch = 0.0
    simulation_end_epoch = constants.JULIAN_DAY
    bodies_to_create = ["Sun", "Earth", "Moon", "Mars", "Venus"]
    global_frame_origin = "Earth"
    global_frame_orientation = "J2000"

    # Create system of bodies
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, global_frame_origin, global_frame_orientation
    )
    bodies = environment_setup.create_system_of_bodies(body_settings)

    # Add vehicle
    bodies.create_empty_body("Delfi-C3")
    bodies.get("Delfi-C3").mass = 400.0

    # Update vehicle with aerodynamic's interface
    reference_area = 4.0
    drag_coefficient = 1.2
    aero_coefficient_settings = aerodynamic_coefficients.constant(
        reference_area, [drag_coefficient, 0, 0]
    )
    environment_setup.add_aerodynamic_coefficient_interface(
        bodies, "Delfi-C3", aero_coefficient_settings
    )

    # Update vehicle with SRP interface
    radiation_pressure_settings = radiation_pressure.cannonball(
        "Sun", 4.0, 1.2, ["Earth"]
    )
    environment_setup.add_radiation_pressure_interface(
        bodies, "Delfi-C3", radiation_pressure_settings
    )

    # Propagation settings
    bodies_to_propagate = ["Delfi-C3"]
    central_bodies = ["Earth"]
    acceleration_settings = {
        "Delfi-C3": {
            "Sun": [
                # acceleration.cannonball_radiation_pressure(),
                acceleration.point_mass_gravity(),
            ],
            "Earth": [
                acceleration.point_mass_gravity(),
                acceleration.aerodynamic(),
            ],
            "Moon": [acceleration.point_mass_gravity()],
            "Mars": [acceleration.point_mass_gravity()],
            "Venus": [acceleration.point_mass_gravity()],
        }
    }
    acceleration_models = propagation_setup.create_acceleration_models(
        bodies, acceleration_settings, bodies_to_propagate, central_bodies
    )

    spice.clear_kernels()

    return None
