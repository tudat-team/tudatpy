from tudatpy import numerical_simulation
from tudatpy.numerical_simulation import environment_setup, propagation_setup
from tudatpy.astro import time_conversion
from tudatpy.interface import spice
from datetime import timedelta
import numpy as np

class CreateEphemerisSettings:
    def __init__(self, SpaceTrackQuery, GetAccelerationSettingsPerRegime):
        self.SpaceTrackQuery = SpaceTrackQuery
        self.GetAccelerationSettingsPerRegime = GetAccelerationSettingsPerRegime

    def set_object_ephemeris_settings(
            self,
            object_name,
            tle_line_1,
            tle_line_2,
            bodies_to_create = ["Sun", "Earth", "Moon"],
            frame_origin = "Earth",
            frame_orientation = "J2000",
            dynamical_model = 'SGP4'
    ):

        spice.load_standard_kernels()
        body_settings = environment_setup.get_default_body_settings(
            bodies_to_create,
            frame_origin,
            frame_orientation)

        body_settings.add_empty_settings(f'{object_name}')

        if dynamical_model == 'SGP4':
            body_settings.get(f'{object_name}').ephemeris_settings = environment_setup.ephemeris.sgp4(tle_line_1, tle_line_2,frame_origin,frame_orientation)
            return environment_setup.ephemeris.sgp4(tle_line_1, tle_line_2,frame_origin,frame_orientation)

        # depending for instance on the dynamical model, we can give different ephemeris types (custom, keplerian, etc...)
        elif dynamical_model == 'LEO-REGIME':
            reference_epoch_utc = self.SpaceTrackQuery.TleUtils.get_tle_reference_epoch(self.SpaceTrackQuery.TleUtils, tle_line_1)  #exact module location TBD`
            end_epoch_utc = reference_epoch_utc + timedelta(hours=5)

            simulation_start_epoch = time_conversion.datetime_to_tudat(reference_epoch_utc).epoch().to_float()
            simulation_end_epoch = time_conversion.datetime_to_tudat(end_epoch_utc).epoch().to_float()
            propagation_times = np.arange(simulation_start_epoch, simulation_end_epoch, 60)
            initial_state = self.SpaceTrackQuery.TleUtils.tle_to_TleEphemeris_object(self.SpaceTrackQuery.TleUtils, tle_line_1, tle_line_2).cartesian_state(simulation_start_epoch)

            bodies = environment_setup.create_system_of_bodies(body_settings)

            bodies_to_propagate = [object_name]
            central_bodies = ["Earth"]

            acceleration_settings = self.GetAccelerationSettingsPerRegime.get_LEO_acceleration_settings()

            ### MAKE THIS INTO A FUNCTION BASED ON REGIME ###
            acceleration_settings = {object_name: acceleration_settings}
            acceleration_models = propagation_setup.create_acceleration_models(
                bodies, acceleration_settings, bodies_to_propagate, central_bodies
            )
            ### MAKE THIS INTO A FUNCTION BASED ON REGIME ###

            termination_settings = propagation_setup.propagator.time_termination(simulation_end_epoch)
            # Create numerical integrator settings
            fixed_step_size = 10.0
            integrator_settings = propagation_setup.integrator.runge_kutta_fixed_step(
                time_step = numerical_simulation.Time(fixed_step_size),
                coefficient_set = propagation_setup.integrator.rk_4)

            # Create propagation settings
            propagator_settings = propagation_setup.propagator.translational(
                central_bodies,
                acceleration_models,
                bodies_to_propagate,
                initial_state,
                numerical_simulation.Time(simulation_start_epoch),
                integrator_settings,
                termination_settings
            )

            # Create simulation object and propagate the dynamics
            dynamics_simulator = numerical_simulation.create_dynamics_simulator(
                bodies, propagator_settings
            )

            # Extract the resulting state history and convert it to an ndarray
            state_history = dynamics_simulator.propagation_results.state_history_float

            print(state_history)

            tabulated_ephemeris = environment_setup.ephemeris.tabulated(
                body_state_history=state_history,
                frame_origin='Earth',
                frame_orientation=frame_orientation
            )

            body_settings.get(object_name).ephemeris_settings = tabulated_ephemeris

            return tabulated_ephemeris