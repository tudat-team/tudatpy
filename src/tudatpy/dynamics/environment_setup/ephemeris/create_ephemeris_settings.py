from tudatpy import dynamics
from tudatpy.dynamics import environment_setup, propagation_setup
from tudatpy.astro import time_representation
from tudatpy.interface import spice
from datetime import timedelta
import numpy as np
from tudatpy.data.spacetrack import SpaceTrackQuery
from tudatpy.data.discos import DiscosQuery
import warnings



class CreateEphemerisSettings:
    def __init__(self, SpaceTrackQuery, GetAccelerationSettingsPerRegime):
        self.spactrack_request = SpaceTrackQuery
        self.GetAccelerationSettingsPerRegime = GetAccelerationSettingsPerRegime
        self.tle_query = self.spactrack_request.DownloadTle(self.spactrack_request)
        self.omm_utils = self.spactrack_request.OMMUtils(self.tle_query)

    def set_object_ephemeris_settings(
            self,
            object_name,
            tle_line_1,
            tle_line_2,
            simulation_start_epoch = None, # float
            simulation_end_epoch = None, # float
            timestep_global = 5, # seconds
            bodies_to_create = ["Sun", "Earth", "Moon"],
            frame_origin = "Earth",
            frame_orientation = "J2000",
            orbital_regime = None,
            sgp4_flag = None,
            aerodynamic = True,
            srp = True,
            mass = 260,
            reference_area = 20,
            drag_coefficient = 1.2,
            radiation_coefficient = 1.2
    ):

        supported_regimes = self.omm_utils.supported_orbital_regimes
        # Checking start and end of simulation
        reference_epoch_tle = time_representation.datetime_to_tudat(self.omm_utils.get_tle_reference_epoch(tle_line_1)).epoch()
        if not simulation_start_epoch and not simulation_end_epoch:
            simulation_start_epoch = reference_epoch_tle
            simulation_start_epoch_datetime = time_representation.DateTime.to_python_datetime(time_representation.DateTime.from_epoch(simulation_start_epoch))
            simulation_end_epoch = time_representation.datetime_to_tudat(simulation_start_epoch_datetime + timedelta(hours=5)).epoch()
            warnings.warn(
                "No simulation start nor end epoch provided.\n"
                "Starting simulation at TLE reference epoch.\n"
                "Ending simulation at TLE reference epoch + 5 hours.",
                UserWarning
            )
        elif not simulation_start_epoch and simulation_end_epoch:
            simulation_end_epoch = time_representation.datetime_to_tudat(simulation_end_epoch).epoch()
            simulation_start_epoch = reference_epoch_tle
            warnings.warn('No simulation start epoch provided. Starting simulation at TLE reference epoch.', UserWarning)
        elif simulation_start_epoch and not simulation_end_epoch:
            simulation_start_epoch_datetime = time_representation.DateTime.to_python_datetime(time_representation.DateTime.from_epoch(simulation_start_epoch))
            simulation_end_epoch = time_representation.datetime_to_tudat(simulation_start_epoch_datetime + timedelta(hours=5)).epoch()
            warnings.warn('No simulation end epoch provided. Ending simulation at TLE reference epoch + 5 hours.')

        if orbital_regime and orbital_regime not in supported_regimes:
            warnings.warn(f'User-defined orbital regime {orbital_regime} is not supported.\nSupported orbital regimes: {supported_regimes}')
            exit()

        spice.load_standard_kernels()
        body_settings = environment_setup.get_default_body_settings(
            bodies_to_create,
            frame_origin,
            frame_orientation)

        body_settings.get("Earth").rotation_model_settings = environment_setup.rotation_model.gcrs_to_itrs(
            environment_setup.rotation_model.iau_2006,
            frame_orientation )
        body_settings.get("Earth").gravity_field_settings.associated_reference_frame = "ITRS"

        body_settings.add_empty_settings(f'{object_name}')
        body_settings.get(object_name).constant_mass = mass

        if aerodynamic:
            #body_settings.get( "Earth" ).atmosphere_settings = environment_setup.atmosphere.nrlmsise00()
            reference_area = 20  # Average projection area of a 3U CubeSat
            drag_coefficient = 1.2
            aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(
                reference_area, [drag_coefficient, 0.0, 0.0]
            )
            # Add the aerodynamic interface to the environment
            body_settings.get(object_name).aerodynamic_coefficient_settings = aero_coefficient_settings

        if srp:
            # Create radiation pressure settings
            reference_area_radiation = 20  # Average projection area of a 3U CubeSat
            radiation_pressure_coefficient = 1.2
            occulting_bodies = dict()
            occulting_bodies["Sun"] = ["Earth"]
            radiation_pressure_settings = environment_setup.radiation_pressure.cannonball_radiation_target(
                reference_area_radiation, radiation_pressure_coefficient, occulting_bodies)
            # Add the radiation pressure interface to the environment
            body_settings.get(object_name).radiation_pressure_target_settings = radiation_pressure_settings

        if not orbital_regime:
            if sgp4_flag:
                sgp4_ephemeris = environment_setup.ephemeris.sgp4(tle_line_1, tle_line_2,frame_origin,frame_orientation)
                tabulated_ephemeris =  environment_setup.ephemeris.tabulated_from_existing(
                    sgp4_ephemeris,
                    simulation_start_epoch,
                    simulation_end_epoch,
                    timestep_global)
            else:
                raise ValueError(f'No dynamical model defined. orbital_regime is {orbital_regime}, and sgp4_flag is {sgp4_flag}')

        else:
            initial_state = self.omm_utils.tle_to_TleEphemeris_object(tle_line_1, tle_line_2).cartesian_state(simulation_start_epoch)
            temp_bodies = environment_setup.create_system_of_bodies(body_settings)

            bodies_to_propagate = [object_name]
            central_bodies = ["Earth"]

            # get acceleration by orbital regime
            acceleration_settings = self.GetAccelerationSettingsPerRegime.get_acceleration_settings(orbital_regime)

            acceleration_settings = {object_name: acceleration_settings}
            acceleration_models = propagation_setup.create_acceleration_models(
                temp_bodies, acceleration_settings, bodies_to_propagate, central_bodies
            )

            termination_settings = propagation_setup.propagator.time_termination(simulation_end_epoch)
            # Create numerical integrator settings
            integrator_settings = propagation_setup.integrator. \
                runge_kutta_fixed_step_size(initial_time_step= time_representation.Time(timestep_global),
                                            coefficient_set=propagation_setup.integrator.CoefficientSets.rkdp_87)

            # Create propagation settings
            propagator_settings = propagation_setup.propagator.translational(
                central_bodies,
                acceleration_models,
                bodies_to_propagate,
                initial_state,
                simulation_start_epoch,
                integrator_settings,
                termination_settings
            )

            # Create simulation object and propagate the dynamics
            dynamics_simulator = dynamics.simulator.create_dynamics_simulator(
                temp_bodies, propagator_settings
            )

            # Extract the resulting state history and convert it to an ndarray
            state_history = dynamics_simulator.propagation_results.state_history_float

            tabulated_ephemeris = environment_setup.ephemeris.tabulated(
                body_state_history=state_history,
                frame_origin='Earth',
                frame_orientation=frame_orientation
            )

        body_settings.get(object_name).ephemeris_settings = tabulated_ephemeris
        bodies = environment_setup.create_system_of_bodies(body_settings)

        if sgp4_flag:
            return tabulated_ephemeris, bodies
        else:
            return tabulated_ephemeris, bodies, integrator_settings, acceleration_settings, propagator_settings, body_settings