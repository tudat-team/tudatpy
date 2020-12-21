"""
Case Study: Mars 2020
---------------------
Mars 2020 is a Mars rover mission by NASA's Mars Exploration Program that
includes the Perseverance rover and the Ingenuity helicopter drone. It was
launched on 30 July 2020 at 11:50 UTC, and will touch down in Jezero
crater on Mars on 18 February 2021.
"""

###############################################################################
# IMPORT STATEMENTS ###########################################################
###############################################################################
from tudatpy.kernel.interface import spice_interface
from tudatpy.util import result2array
from tudatpy.kernel.simulation import (environment_setup, propagation_setup)
from tudatpy.kernel.astro.two_body_dynamics import (
    LambertTargeterIzzo,
)

import numpy as np
import matplotlib.pyplot as plt
from tudatpy.kernel import constants

if __name__ == "__main__":
    spice_kernels = spice_interface.get_standard_kernels() + [
        "/home/ggarrett/Downloads/de430.bsp",
        "/home/ggarrett/Downloads/de431_part-1.bsp",
        "/home/ggarrett/Downloads/de431_part-2.bsp"]

    for kernel in spice_kernels:
        spice_interface.load_kernel(kernel)

    import matplotlib.pyplot as plt

    # mission system parent body
    parent_body = "Sun"
    parent_body_mu = spice_interface.get_body_gravitational_parameter(parent_body)

    # define reference frame
    frame_origin = parent_body
    frame_orientation = "ECLIPJ2000"

    # Mission planetary departure parameter definition
    departure_date = "July 30, 2020"
    departure_body = "Earth"
    departure_r_limit = 1.048

    # Mission planetary arrival parameter definition
    rendezvous_date = "18 February, 2021"
    rendezvous_body = "Mars"
    rendezvous_r_limit = 1.05

    # Calculate ephemeris time and time of flight parameters
    departure_reference_ephemeris_time = spice_interface.convert_date_string_to_ephemeris_time(departure_date)
    rendezvous_reference_ephemeris_time = spice_interface.convert_date_string_to_ephemeris_time(rendezvous_date)

    # Determine constant departure parameters
    departure_body_mu = spice_interface.get_body_gravitational_parameter(departure_body)
    departure_body_radius = spice_interface.get_average_radius(departure_body)

    # Determine constant arrival geometry
    rendezvous_body_mu = spice_interface.get_body_gravitational_parameter(rendezvous_body)
    rendezvous_body_radius = spice_interface.get_average_radius(rendezvous_body)

    # Define departure body state at departure conic reference
    # NEW DESIGN
    # departure_body_state_at_reference = cartesian_state_from_spice(
    #     body_name=departure_body,
    #     ephemeris_time=departure_reference_ephemeris_time,
    #     gravitational_parameter=parent_body_mu,
    #     frame_origin=frame_origin,
    #     frame_orientation=frame_orientation,
    #     aberration_corrections="none"
    # )
    # CURRENT DESIGN
    #     m.def("get_body_cartesian_state_at_epoch",
    #       &tudat::spice_interface::getBodyCartesianStateAtEpoch,
    #       py::arg("target_body_name"),
    #       py::arg("observer_body_name"),
    #       py::arg("reference_frame_name"),
    #       py::arg("aberration_corrections"),
    #       py::arg("ephemeris_time"),
    #       "Get Cartesian position of a body, as observed from another body.");
    departure_body_cartesian_state_vector = spice_interface.get_body_cartesian_state_at_epoch(
        target_body_name=departure_body,
        ephemeris_time=departure_reference_ephemeris_time,
        observer_body_name=frame_origin,
        reference_frame_name=frame_orientation,
        aberration_corrections="none"
    )

    # NEW DESIGN
    # departure_body_state_at_reference = CartesianState(
    #     position_vector=departure_body_cartesian_state_vector[:3],
    #     velocity_vector=departure_body_cartesian_state_vector[3:],
    #     gravitational_parameter=parent_body_mu,
    # )

    # Determine rendezvous body state at rendezvous conic reference
    # rendezvous_body_state_at_reference = cartesian_state_from_spice(
    #     body_name=rendezvous_body,
    #     ephemeris_time=rendezvous_reference_ephemeris_time,
    #     gravitational_parameter=parent_body_mu,
    #     frame_origin=frame_origin,
    #     frame_orientation=frame_orientation,
    #     aberration_corrections="none")

    rendezvous_body_cartesian_state_vector = spice_interface.get_body_cartesian_state_at_epoch(
        target_body_name=rendezvous_body,
        ephemeris_time=rendezvous_reference_ephemeris_time,
        observer_body_name=frame_origin,
        reference_frame_name=frame_orientation,
        aberration_corrections="none"
    )

    # NEW DESIGN
    # rendezvous_body_state_at_reference = CartesianState(
    #     position_vector=rendezvous_body_cartesian_state_vector[:3],
    #     velocity_vector=rendezvous_body_cartesian_state_vector[3:],
    #     gravitational_parameter=parent_body_mu,
    # )

    # Calculate coast interplanetary leg
    departure_velocity, rendezvous_velocity = LambertTargeterIzzo(
        # departure_position=departure_body_state_at_reference.position_vector, # NEW DESIGN
        departure_position=departure_body_cartesian_state_vector[:3],
        # arrival_position=rendezvous_body_state_at_reference.position_vector, # NEW DESIGN
        arrival_position=rendezvous_body_cartesian_state_vector[:3],
        time_of_flight=rendezvous_reference_ephemeris_time - departure_reference_ephemeris_time,
        gravitational_parameter=parent_body_mu
    ).get_velocity_vectors()

    simulation_start_epoch = departure_reference_ephemeris_time
    simulation_end_epoch = rendezvous_reference_ephemeris_time

    ###########################################################################
    # CREATE ENVIRONMENT AND VEHICLE ##########################################
    ###########################################################################

    # Create default body settings for "Earth"
    bodies_to_create = [departure_body, rendezvous_body, parent_body]

    # Create default body settings for bodies_to_create, with "SSB"/"J2000" as
    # global frame origin and orientation
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, frame_origin, frame_orientation
    )

    # Create system of bodies (in this case only Earth)
    bodies = environment_setup.create_system_of_bodies(body_settings)

    # Add vehicle object to system of bodies
    bodies.create_empty_body("Mars-2020")

    ###########################################################################
    # CREATE ACCELERATIONS ####################################################
    ###########################################################################

    # Define bodies that are propagated.
    bodies_to_propagate = ["Mars-2020"]

    # Define central bodies of propagation.
    central_bodies = [parent_body]

    # Define accelerations acting on Delfi-C3.
    acceleration_settings_mars2020 = {
        # departure_body: [propagation_setup.acceleration.point_mass_gravity()],
        rendezvous_body: [propagation_setup.acceleration.point_mass_gravity()],
        parent_body: [propagation_setup.acceleration.point_mass_gravity()],
    }

    acceleration_settings = {
        "Mars-2020": acceleration_settings_mars2020,
    }

    # Create acceleration models.
    acceleration_models = propagation_setup.create_acceleration_models(
        bodies, acceleration_settings, bodies_to_propagate, central_bodies
    )

    ###########################################################################
    # CREATE PROPAGATION SETTINGS #############################################
    ###########################################################################

    # Set initial conditions for the Mars 2020 spacecraft
    initial_state = np.hstack([
        departure_body_cartesian_state_vector[:3],
        departure_velocity
    ])

    # Define list of dependent variables to save.
    dependent_variables_to_save = [
        propagation_setup.dependent_variable.relative_position("Mars-2020", frame_origin),
        propagation_setup.dependent_variable.relative_position(rendezvous_body, frame_origin),
        propagation_setup.dependent_variable.relative_position(departure_body, frame_origin),
        propagation_setup.dependent_variable.relative_position("Mars-2020", departure_body)
    ]

    # Create propagation settings.
    propagator_settings = propagation_setup.propagator.translational(
        central_bodies,
        acceleration_models,
        bodies_to_propagate,
        initial_state,
        simulation_end_epoch,
        output_variables=dependent_variables_to_save,
    )

    # Create numerical integrator settings.
    fixed_step_size = 300.0
    integrator_settings = propagation_setup.integrator.runge_kutta_4(
        simulation_start_epoch, fixed_step_size
    )

    ###########################################################################
    # PROPAGATE ORBIT #########################################################
    ###########################################################################

    # Create simulation object and propagate dynamics.
    dynamics_simulator = propagation_setup.SingleArcDynamicsSimulator(
        bodies, integrator_settings, propagator_settings, True
    )

    states = result2array(dynamics_simulator.dependent_variable_history)

    scale = constants.AU
    sun_circle = plt.Circle((0, 0), 0.2)
    fig = plt.gcf()
    ax = fig.gca()
    ax.set_aspect('equal', adjustable='box')
    ax.plot(states[:, 1] / scale, states[:, 2] / scale)
    ax.plot(states[:, 4] / scale, states[:, 5] / scale)
    ax.plot(states[:, 7] / scale, states[:, 8] / scale)
    ax.plot([0], [0])
    plt.show()

    spice_kernel_repr = "\n".join(spice_kernels)
    print(
        f"""
Spice Kernels Loaded
====================
{spice_kernel_repr}

""")
