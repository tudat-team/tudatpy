################################################################################
# IMPORT STATEMENTS ############################################################
################################################################################
from tudatpy.kernel import constants
from tudatpy.kernel.interface import spice_interface
from tudatpy.kernel.simulation import environment_setup
from tudatpy.kernel.simulation import propagation_setup

################################################################################
# GENERAL SIMULATION SETUP #####################################################
################################################################################

# Load spice kernels.
spice_interface.load_standard_kernels()

# Set simulation start epoch.
simulation_start_epoch = 1.0E7

# Set numerical integration fixed step size.
fixed_step_size = 3600.0

# Set simulation end epoch.
simulation_end_epoch = 1.0E7 + 5.0 * constants.JULIAN_YEAR

################################################################################
# SETUP ENVIRONMENT ############################################################
################################################################################

# Define bodies in simulation.
bodies_to_create = bodies_to_propagate = ["Moon", "Earth", "Mars", "Venus",
                                          "Mercury", "Sun"]

# Create bodies in simulation.
body_settings = environment_setup.get_default_body_settings(bodies_to_create)
body_system = environment_setup.create_bodies(body_settings)
environment_setup.set_global_frame_body_ephemerides(body_system, "SSB", "ECLIPJ2000")

################################################################################
# SETUP PROPAGATION ############################################################
################################################################################

results = {}
for propagation_variant in ["barycentric", "hierarchical"]:

    ############################################################################
    # SETUP PROPAGATION : CREATE ACCELERATION MODELS ###########################
    ############################################################################

    # Create barycentric body settings
    acceleration_dict = {}
    for body_i in bodies_to_create:
        current_accelerations = {}
        for body_j in bodies_to_create:
            if body_i != body_j:
                current_accelerations[body_j] = [
                    propagation_setup.AccelerationSettings(propagation_setup.AvailableAcceleration.central_gravity)
                ]
        acceleration_dict[body_i] = current_accelerations
    central_bodies = []

    # Barycentric propagation.
    if propagation_variant == "barycentric":
        central_bodies = ["SSB"] * len(bodies_to_create)

    # Hierarchical parent body propagation.
    elif propagation_variant == "hierarchical":
        for body_name in bodies_to_create:
            if body_name == "Moon":
                central_bodies.append("Earth")
            elif body_name == "Sun":
                central_bodies.append("SSB")
            else:
                central_bodies.append("Sun")

    # Convert acceleration mappings into acceleration models.
    acceleration_models = propagation_setup.create_acceleration_models_dict(
        body_system=body_system,
        selected_acceleration_per_body=acceleration_dict,
        bodies_to_propagate=bodies_to_propagate,
        central_bodies=central_bodies
    )

    ############################################################################
    # SETUP PROPAGATION : PROPAGATION SETTINGS #################################
    ############################################################################

    # Get system initial state.
    system_initial_state = propagation_setup.get_initial_state_of_bodies(
        bodies_to_propagate=bodies_to_propagate,
        central_bodies=central_bodies,
        body_system=body_system,
        initial_time=simulation_start_epoch
    )
    # Create propagation settings.
    propagator_settings = propagation_setup.TranslationalStatePropagatorSettings(
        central_bodies,
        acceleration_models,
        bodies_to_propagate,
        system_initial_state,
        simulation_end_epoch
    )
    # Create numerical integrator settings.
    integrator_settings = propagation_setup.IntegratorSettings(
        propagation_setup.AvailableIntegrators.rk4,
        simulation_start_epoch,
        fixed_step_size
    )

    ############################################################################
    # PROPAGATE ################################################################
    ############################################################################

    # Instantiate the dynamics simulator.
    dynamics_simulator = propagation_setup.SingleArcDynamicsSimulator(
        body_system, integrator_settings, propagator_settings, True)

    # Propagate and store results to outer loop results dictionary.
    results[propagation_variant] = dynamics_simulator.get_equations_of_motion_numerical_solution()

################################################################################
# VISUALISATION / OUTPUT / PRELIMINARY ANALYSIS ################################
################################################################################
