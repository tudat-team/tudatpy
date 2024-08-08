###############################################################################
# IMPORT STATEMENTS ###########################################################
###############################################################################
import numpy as np
import pandas as pd
from tudatpy import elements
from tudatpy.kernel import constants
from tudatpy.kernel import numerical_integrators
from tudatpy.kernel import spice_interface
from tudatpy.kernel import propagators
from tudatpy.kernel import simulation_setup


###############################################################################
# FUNCTION DEFINITION OF TUTORIAL #############################################
###############################################################################
def single_satellite_propagator(
        start_epoch=0.0,
        fixed_step_size=10.0,
        end_epoch=constants.JULIAN_DAY,
        parent_body="Earth",
        frame_orientation="ECLIPJ2000",
        frame_origin="SSB",
        satellite_name="Delfi-C3",
        sat_sma=7500.0E3,
        sat_ecc=0.1,
        sat_inc=np.deg2rad(85.3),
        sat_raan=np.deg2rad(23.4),
        sat_argp=np.deg2rad(235.7),
        sat_nu=np.deg2rad(139.87),
        return_output=False,
        print_output=True,
        output_type='markdown'  # [array, markdown, latex, string]
):
    """
    Tutorial 1 of the tudatpy library in function form.

    Warning
    -------
    This function is not suited to mass generation of output data, as it
    is creating the environment and simulation from scratch for each call.

    Parameters
    ----------
    start_epoch
    fixed_step_size
    end_epoch
    parent_body
    frame_orientation
    frame_origin
    satellite_name
    sat_sma
    sat_ecc
    sat_inc
    sat_raan
    sat_argp
    sat_nu
    return_output
    print_output
    output_type

    Returns
    -------

    """
    # Load spice kernels.
    spice_interface.load_standard_kernels()

    # Set simulation start epoch.
    simulation_start_epoch = start_epoch

    # Set numerical integration fixed step size.
    fixed_step_size = fixed_step_size

    # Set simulation end epoch.
    simulation_end_epoch = end_epoch

    ###########################################################################
    # CREATE ENVIRONMENT ######################################################
    ###########################################################################

    # Create body objects.
    bodies_to_create = [parent_body]

    body_settings = simulation_setup.get_default_body_settings(bodies_to_create)

    body_settings[parent_body].ephemeris_settings = simulation_setup.ConstantEphemerisSettings(
        np.zeros(6))

    body_settings[parent_body].rotation_model_settings.reset_original_frame(frame_orientation)

    # Create Earth Object.
    bodies = simulation_setup.create_bodies(body_settings)

    ###########################################################################
    # CREATE VEHICLE ##########################################################
    ###########################################################################

    # Create vehicle objects.
    bodies[satellite_name] = simulation_setup.Body()

    ###########################################################################
    # FINALIZE BODIES #########################################################
    ###########################################################################

    simulation_setup.set_global_frame_body_ephemerides(bodies, frame_origin,
                                                       frame_orientation)

    ###########################################################################
    # CREATE ACCELERATIONS ####################################################
    ###########################################################################

    # Define bodies that are propagated.
    bodies_to_propagate = [satellite_name]

    # Define central bodies.
    central_bodies = [parent_body]

    # Define accelerations acting on Delfi-C3.
    accelerations_of_satellite = {
        parent_body: [simulation_setup.Acceleration.point_mass_gravity()]
    }

    # Create global accelerations dictionary.
    accelerations = {satellite_name: accelerations_of_satellite}

    # Create acceleration models.
    acceleration_models = simulation_setup.create_acceleration_models_dict(
        bodies, accelerations, bodies_to_propagate, central_bodies)

    ###########################################################################
    # CREATE PROPAGATION SETTINGS #############################################
    ###########################################################################

    # Set initial conditions for the Asterix satellite that will be
    # propagated in this simulation. The initial conditions are given in
    # Keplerian elements and later on converted to Cartesian elements.

    # Set Keplerian elements for satellite
    parent_gravitational_parameter = bodies[
        parent_body].gravity_field_model.get_gravitational_parameter()

    # REVISED CONTEMPORARY DESIGN.
    system_initial_state = elements.keplerian2cartesian(
        mu=parent_gravitational_parameter,
        a=sat_sma,
        ecc=sat_ecc,
        inc=sat_inc,
        raan=sat_raan,
        argp=sat_argp,
        nu=sat_nu)

    # Create propagation settings.
    propagator_settings = propagators.TranslationalStatePropagatorSettings(
        central_bodies,
        acceleration_models,
        bodies_to_propagate,
        system_initial_state,
        simulation_end_epoch
    )
    # Create numerical integrator settings.
    integrator_settings = numerical_integrators.IntegratorSettings(
        numerical_integrators.AvailableIntegrators.rungeKutta4,
        simulation_start_epoch,
        fixed_step_size
    )

    ###########################################################################
    # PROPAGATE ORBIT #########################################################
    ###########################################################################

    # Create simulation object and propagate dynamics.
    dynamics_simulator = propagators.SingleArcDynamicsSimulator(
        bodies, integrator_settings, propagator_settings, True)
    result = dynamics_simulator.get_equations_of_motion_numerical_solution()

    ###########################################################################
    # OUTPUT HANDLING #########################################################
    ###########################################################################

    if output_type != 'array':
        df = pd.DataFrame(data=np.vstack((
            result[simulation_start_epoch] / 1E3,
            result[simulation_end_epoch] / 1E3,
        )), columns=['R_x [km]', 'R_y [km]', 'R_z [km]', 'V_x [km/s]', 'V_y [km/s]', 'V_z [km/s]']
        )
        df.index = ['t_0', 't_f']

    # Output type
    if output_type == 'markdown':
        output = df.to_markdown()

    elif output_type == 'latex':
        output = df.to_latex(escape=False)

    elif output_type == 'df':
        output = df

    else:
        output = np.vstack((
            result[simulation_start_epoch] / 1E3,
            result[simulation_end_epoch] / 1E3,
        ))

    # Print output
    if print_output:
        print(output)

    # Return output
    if return_output:
        return output
    else:
        return None
