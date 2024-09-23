# Load standard modules
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import pytest

# Load tudatpy modules
from tudatpy.interface import spice
from tudatpy import numerical_simulation
from tudatpy.numerical_simulation import environment_setup, propagation_setup
from tudatpy.astro import element_conversion
from tudatpy import constants
from tudatpy.util import result2array

# Semantic variable history stuff
import pickle
from tudatpy.numerical_simulation.propagation.dependent_variable_dictionary import (
    create_dependent_variable_dictionary,
)


@pytest.mark.skip(reason="Implementation error to be solved in future release")
def test_dependent_variable_dictionary():

    # %% SET UP A PROPAGATION TO TEST THE SEMANTIC VARIABLE HISTORY FUNCTIONALITY
    # ===========================================================================

    # Load spice kernels
    spice.load_standard_kernels()

    # Set simulation start and end epochs
    simulation_start_epoch = 0.0
    simulation_end_epoch = constants.JULIAN_DAY

    # Define string names for bodies to be created from default.
    bodies_to_create = ["Sun", "Earth", "Moon", "Mars", "Venus"]

    # Use "Earth"/"J2000" as global frame origin and orientation.
    global_frame_origin = "Earth"
    global_frame_orientation = "J2000"

    # Create default body settings, usually from `spice`.
    body_settings = environment_setup.get_default_body_settings(
        bodies_to_create, global_frame_origin, global_frame_orientation
    )

    # Create system of selected celestial bodies
    bodies = environment_setup.create_system_of_bodies(body_settings)

    # Create vehicle objects.
    bodies.create_empty_body("Delfi-C3")

    bodies.get("Delfi-C3").mass = 400.0

    # Create aerodynamic coefficient interface settings, and add to vehicle
    reference_area = 4.0
    drag_coefficient = 1.2
    aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(
        reference_area, [drag_coefficient, 0, 0]
    )
    environment_setup.add_aerodynamic_coefficient_interface(
        bodies, "Delfi-C3", aero_coefficient_settings
    )

    # To account for the pressure of the solar radiation on the satellite, let's add another interface. This takes a radiation pressure coefficient of 1.2, and a radiation area of 4m$^2$. This interface also accounts for the variation in pressure cause by the shadow of Earth.

    # Create radiation pressure settings, and add to vehicle
    reference_area_radiation = 4.0
    radiation_pressure_coefficient = 1.2
    occulting_bodies = ["Earth"]
    radiation_pressure_settings = environment_setup.radiation_pressure.cannonball(
        "Sun",
        reference_area_radiation,
        radiation_pressure_coefficient,
        occulting_bodies,
    )
    environment_setup.add_radiation_pressure_interface(
        bodies, "Delfi-C3", radiation_pressure_settings
    )

    # Define bodies that are propagated
    bodies_to_propagate = ["Delfi-C3"]

    # Define central bodies of propagation
    central_bodies = ["Earth"]

    # Define accelerations acting on Delfi-C3 by Sun and Earth.
    accelerations_settings_delfi_c3 = dict(
        Sun=[
            propagation_setup.acceleration.cannonball_radiation_pressure(),
            propagation_setup.acceleration.point_mass_gravity(),
        ],
        Earth=[
            propagation_setup.acceleration.spherical_harmonic_gravity(5, 5),
            propagation_setup.acceleration.aerodynamic(),
        ],
        Moon=[propagation_setup.acceleration.point_mass_gravity()],
        Mars=[propagation_setup.acceleration.point_mass_gravity()],
        Venus=[propagation_setup.acceleration.point_mass_gravity()],
    )

    # Create global accelerations settings dictionary.
    acceleration_settings = {"Delfi-C3": accelerations_settings_delfi_c3}

    # Create acceleration models.
    acceleration_models = propagation_setup.create_acceleration_models(
        bodies, acceleration_settings, bodies_to_propagate, central_bodies
    )

    # Set initial conditions for the satellite that will be
    # propagated in this simulation. The initial conditions are given in
    # Keplerian elements and later on converted to Cartesian elements
    earth_gravitational_parameter = bodies.get("Earth").gravitational_parameter
    initial_state = element_conversion.keplerian_to_cartesian_elementwise(
        gravitational_parameter=earth_gravitational_parameter,
        semi_major_axis=7500.0e3,
        eccentricity=0.1,
        inclination=np.deg2rad(85.3),
        argument_of_periapsis=np.deg2rad(235.7),
        longitude_of_ascending_node=np.deg2rad(23.4),
        true_anomaly=np.deg2rad(139.87),
    )

    # Define list of dependent variables to save
    dependent_variables_to_save = [
        propagation_setup.dependent_variable.total_acceleration("Delfi-C3"),
        propagation_setup.dependent_variable.keplerian_state("Delfi-C3", "Earth"),
        propagation_setup.dependent_variable.latitude("Delfi-C3", "Earth"),
        propagation_setup.dependent_variable.longitude("Delfi-C3", "Earth"),
        propagation_setup.dependent_variable.single_acceleration_norm(
            propagation_setup.acceleration.point_mass_gravity_type, "Delfi-C3", "Sun"
        ),
        propagation_setup.dependent_variable.single_acceleration_norm(
            propagation_setup.acceleration.point_mass_gravity_type, "Delfi-C3", "Moon"
        ),
        propagation_setup.dependent_variable.single_acceleration_norm(
            propagation_setup.acceleration.point_mass_gravity_type, "Delfi-C3", "Mars"
        ),
        propagation_setup.dependent_variable.single_acceleration_norm(
            propagation_setup.acceleration.point_mass_gravity_type, "Delfi-C3", "Venus"
        ),
        propagation_setup.dependent_variable.single_acceleration_norm(
            propagation_setup.acceleration.spherical_harmonic_gravity_type,
            "Delfi-C3",
            "Earth",
        ),
        propagation_setup.dependent_variable.single_acceleration_norm(
            propagation_setup.acceleration.aerodynamic_type, "Delfi-C3", "Earth"
        ),
        propagation_setup.dependent_variable.single_acceleration_norm(
            propagation_setup.acceleration.cannonball_radiation_pressure_type,
            "Delfi-C3",
            "Sun",
        ),
    ]

    # Create termination settings
    termination_condition = propagation_setup.propagator.time_termination(
        simulation_end_epoch
    )

    # Create numerical integrator settings
    fixed_step_size = 10.0
    integrator_settings = propagation_setup.integrator.runge_kutta_4(fixed_step_size)

    # Create propagation settings
    propagator_settings = propagation_setup.propagator.translational(
        central_bodies,
        acceleration_models,
        bodies_to_propagate,
        initial_state,
        simulation_start_epoch,
        integrator_settings,
        termination_condition,
        output_variables=dependent_variables_to_save,
    )

    # Create simulation object and propagate the dynamics
    dynamics_simulator = numerical_simulation.create_dynamics_simulator(
        bodies, propagator_settings
    )

    # %% TEST SEMANTIC VARIABLE HISTORY FUNCTIONALITY
    # ===============================================

    # Create semantic dependent variable history
    dep_vars_dict = create_dependent_variable_dictionary(dynamics_simulator)

    # Retrieve time history
    dependent_variable_history_array = result2array(
        dynamics_simulator.dependent_variable_history
    )

    arrays_are_equal = lambda a, b: np.allclose(a, b, rtol=1e-10, atol=1e-10)

    # TEST 1: Assert that the time histories obtained from newly created dependent variable settings
    #         objects and from those used to set up the propagation are the same.
    assert arrays_are_equal(
        dep_vars_dict.asarray(dependent_variables_to_save[0]),
        dep_vars_dict.asarray(
            propagation_setup.dependent_variable.total_acceleration("Delfi-C3")
        ),
    )

    # TEST 2: Assert that the time history of the dependent variables is the same as the one obtained from the
    #         dynamics simulator.
    assert arrays_are_equal(
        dep_vars_dict.asarray(dependent_variables_to_save[2]),
        dependent_variable_history_array[:, 10],
    )

    # TEST 3: Assert (by virtue of an error not happening) that result2array works as expected on time histories of
    #         scalar dependent variables.
    #         This will fail for vectorial or matrix dependent variables because the shape of the array associated to each
    #         epoch is not flat. It is not recommended to use result2array together with the create_dependent_variable_dictionary due
    #         to this. The test is useful though as the entries of the semantic variable history associated to scalar dependent
    #         variables should be perfectly compatible with result2array.
    delfi_latitude_array = result2array(dep_vars_dict[dependent_variables_to_save[2]])

    # TEST 4: Assert that the shape of the vectorial dependent variable histories is correct
    assert dep_vars_dict.asarray(dependent_variables_to_save[0])[0].shape == (3,)
    assert dep_vars_dict.asarray(dependent_variables_to_save[1])[0].shape == (6,)
