###############################################################################
# IMPORT STATEMENTS ###########################################################
###############################################################################
import numpy as np
from tudatpy import constants
from tudatpy import elements
from tudatpy import io
from tudatpy import ephemerides
from tudatpy import interpolators
from tudatpy import numerical_integrators
from tudatpy import spice_interface
from tudatpy import basic_astrodynamics
# from tudatpy import orbital_element_conversions # LEGACY MODULE
from tudatpy import propagators
from tudatpy import aerodynamics
from tudatpy import simulation_setup
# from tudatpy.prototype import Environment

from tudatpy import spice_interface
from tudatpy import simulation_setup
from simulation_setup import Body as _Body
import warnings
import numpy as np
import time


class Body(_Body):

    def __init__(self, environment=None, state=np.zeros(6)):
        """

        Parameters
        ----------
        state
        """
        super().__init__(state)
        self.environment = environment

    @property
    def state(self):
        """
        Current state of the body.

        Returns
        -------
        np.ndarray[dim=6]
        """
        return self.get_state()

    @state.setter
    def state(self, x):
        self.set_state(x)

    @property
    def inertia_tensor(self):
        return self.get_body_inertia_tensor()

    @inertia_tensor.setter
    def inertia_tensor(self, x):
        self.set_body_inertia_tensor(x)

    @property
    def atmosphere_model(self):
        return self.get_atmosphere_model()

    @atmosphere_model.setter
    def atmosphere_model(self, x):
        self.set_atmosphere_model(x)

    def _set_radiation_pressure_interface(self,
                                         reference_area,
                                         coefficient,
                                         occulting,
                                         source,
                                         model="cannon_ball"):
        """

        Parameters
        ----------
        target_body_name : str
            Body to be interfaced to the environment.
        reference_area : float
            Radiation pressure reference area (m).
        coefficient : float or ndarray[ndim=3]
            Radiation pressure coefficient (m).
        occulting : list[str]
            Occulting bodies between source and target_body.
        source : str
            Body emitting the radiation pressure.
        model : {"cannon_ball"}
            Radiation pressure model to be used.

        Returns
        -------
        None

        """
        if model == "cannon_ball":
            t0_0 = time.time()

            rad_press_settings = simulation_setup.CannonBallRadiationPressureInterfaceSettings(
                source, reference_area, coefficient, occulting)
            # Create and set radiation pressure settings
            super(Body, self).set_radiation_pressure_interface(
                source, simulation_setup.create_radiation_pressure_interface(
                    rad_press_settings,
                    self.environment.__inv_dict__()[self],
                    self.environment.__dict__()))
            t0_1 = time.time()
            print('set_radiation_pressure_interface: ', t0_1 - t0_0)
        else:
            raise NotImplementedError(f"Model {model} not implemented.")

    def _set_aerodynamic_coefficient_interface(self,
                                              reference_area,
                                              coefficients,
                                              coefficients_in_aerodynamic_frame=True,
                                              coefficients_in_negative_axis_direction=True,
                                              model="constant_coefficients"
                                              ):
        """

        Parameters
        ----------
        target_body_name
        reference_area
        coefficients
        coefficients_in_aerodynamic_frame
        coefficients_in_negative_axis_direction
        model

        Returns
        -------

        """
        if model is "constant_coefficients":
            t0_0 = time.time()

            aero_c_settings = simulation_setup.ConstantAerodynamicCoefficientSettings(
                reference_area,
                coefficients,
                coefficients_in_aerodynamic_frame,
                coefficients_in_negative_axis_direction
            )
            # Create and set aerodynamic coefficients object
            super(Body, self).set_aerodynamic_coefficient_interface(
                simulation_setup.create_aerodynamic_coefficient_interface(
                    aero_c_settings,
                    self.environment.__inv_dict__()[self])
            )
            t0_1 = time.time()
            print('set_aerodynamic_coefficient_interface: ', t0_1 - t0_0)
        else:
            raise NotImplementedError(f"Model {model} is not implemented.")


class Environment:

    def __init__(self,
                 tudatpy_bodies=None,
                 custom_bodies=None,
                 frame_origin=None,
                 frame_orientation=None,
                 start_epoch=None,
                 end_epoch=None,
                 epoch_margin=300.0):
        """

        Parameters
        ----------
        tudatpy_bodies : list[str]
            List of bodies to be added to the environment that are part
            of the available bodies cataloged in the tudatBundle.
        custom_bodies : list[str]
            List of user-defined bodies to be added to the environment
            which are not part of the available bodies in the tudatBundle.
        frame_origin : str
            Frame origin for simulation output results.
        frame_orientation : str
            Frame orientation for simulation output results.

        Examples
        --------


        """
        # Sets as [] if tudatpy_bodies is None.
        self._tudatpy_bodies = tudatpy_bodies if tudatpy_bodies else []

        # Sets as [] if custom_bodies is None.
        self._custom_bodies = custom_bodies if custom_bodies else []

        # Warns ands sets default origin.
        self._frame_origin = frame_origin or self._default_origin()

        # Warns and sets default orientation.
        self._frame_orientation = frame_orientation or self._default_orient()

        self._start_epoch = start_epoch
        self._end_epoch = end_epoch
        self._epoch_margin = abs(epoch_margin)

        self._body_settings = None
        self._environment_finalized = False

        if self._start_epoch and self._end_epoch:

            # Get default body settings for cataloged bodies in tudatpy.
            self._body_settings = simulation_setup.get_default_body_settings(
                self._tudatpy_bodies,
                self._start_epoch - self._epoch_margin,
                self._end_epoch + self._epoch_margin)

        else:

            # Get default body settings for cataloged bodies in tudatpy.
            self._body_settings = simulation_setup.get_default_body_settings(
                self._tudatpy_bodies)

        # Reset all bodies ephemeris and rotation model orientation frame.
        for body in self._tudatpy_bodies:
            self._body_settings[body].ephemeris_settings.reset_frame_orientation(
                self._frame_orientation)
            self._body_settings[body].rotation_model_settings.reset_original_frame(
                self._frame_orientation)

        # Create cataloged bodies.
        self._bodies = simulation_setup.create_bodies(self._body_settings)

        # Create all custom bodies.
        for body in self._custom_bodies:
            self._bodies[body] = simulation_setup.Body()

    def __dict__(self):
        return self._bodies

    def __inv_dict__(self):
        return {v: k for k, v in self._bodies.items()}

    @staticmethod
    def _default_orient(orientation="ECLIPJ2000"):
        warnings.warn("frame_orientation was not explicitly set, so the "
                      f"default of {orientation} is being used.")
        return orientation

    @staticmethod
    def _default_origin(origin="SSB"):
        warnings.warn("frame_origin was not explicitly set, so the default of "
                      f"{origin} is being used.")
        return origin

    def _get_body(self, body):
        # Check that environment has been finalized.
        # try:
        #     assert self._environment_finalized
        # except AssertionError:
        #     raise EnvironmentError(
        #         "Environment bodies must be finalized before adding "
        #         "interfaces or modifying existing bodies."
        #     )

        # Check that body exists in custom or given tudatpy bodies.
        try:
            assert body in self._tudatpy_bodies + self._custom_bodies
        except AssertionError:
            raise EnvironmentError(
                f"The body named: {body}, does not exist in the environment.")

        return self._bodies[body]

    def _set_body(self, key, value):

        # Check that body being defined is not already existing.
        try:
            assert key not in self._custom_bodies + self._tudatpy_bodies
        except AssertionError:
            raise KeyError("Key already exists in the environment bodies.")

        # Check that value being set is a Body type.
        try:
            assert isinstance(value, simulation_setup.Body) is True
        except AssertionError:
            raise AssertionError(
                "Only bodies can be set as direct attributes of the "
                "environment class.")

        self._bodies[key] = value

    # def __setattr__(self, key, value):
    #     self._set_body(key, value)

    def __str__(self):
        return (
            f"""
Environment details
===================
frame_orientation: {self._frame_orientation}
frame_origin: {self._frame_origin}
tudatpy_bodies: {", ".join(self._tudatpy_bodies)}
custom_bodies: {", ".join(self._custom_bodies)}   
        """
        )

    def __setitem__(self, key, value):
        self._set_body(key, value)
        self._custom_bodies.append(key)

    def __getitem__(self, body):
        return self._get_body(body)

    def __getattr__(self, body):
        return self._get_body(body)

    def finalize_environment(self):
        self._environment_finalized = True
        simulation_setup.set_global_frame_body_ephemerides(
            self._bodies, self._frame_origin, self._frame_orientation)

    #


def main():
    # Load spice kernels.
    spice_interface.load_standard_spice_kernels()

    # Set simulation start epoch.
    simulation_start_epoch = 0.0

    # Set numerical integration fixed step size.
    fixed_step_size = 1.0

    # Set simulation end epoch.
    simulation_end_epoch = constants.JULIAN_DAY

    ###########################################################################
    # CREATE ENVIRONMENT ######################################################
    ###########################################################################

    bodies_to_create = ["Sun", "Earth", "Moon", "Mars", "Venus"]

    env = Environment(
        tudatpy_bodies=bodies_to_create,
        custom_bodies=None,
        frame_origin="SSB",
        frame_orientation="ECLIPJ2000",
        start_epoch=simulation_start_epoch,
        end_epoch=simulation_end_epoch,
        epoch_margin=300.0)

    ###########################################################################
    # CREATE VEHICLE ##########################################################
    ###########################################################################

    env["Delfi-C3"] = Body(env)

    env["Delfi-C3"].set_constant_body_mass(400.0)

    ###########################################################################
    # CREATE VEHICLE - ENVIRONMENT INTERFACE ##################################
    ###########################################################################

    # Create aerodynamic coefficient interface settings
    reference_area = 4.0
    aerodynamic_coefficient = 1.2
    aero_c_settings = simulation_setup.ConstantAerodynamicCoefficientSettings(
        reference_area,
        aerodynamic_coefficient * np.ones(3),
        are_coefficients_in_aerodynamic_frame=True,
        are_coefficients_in_negative_axis_direction=True
    )
    # Create and set aerodynamic coefficients object
    env._bodies["Delfi-C3"].set_aerodynamic_coefficient_interface(
        simulation_setup.create_aerodynamic_coefficient_interface(
            aero_c_settings,
            "Delfi-C3")
    )
    # TODO: Simplify (post 1.0.0 work)

    # Create radiation pressure settings
    reference_area_radiation = 4.0
    radiation_pressure_coefficient = 1.2
    occulting_bodies = ["Earth"]
    rad_press_settings = simulation_setup.CannonBallRadiationPressureInterfaceSettings(
        "Sun", reference_area_radiation, radiation_pressure_coefficient, occulting_bodies)

    # Create and set radiation pressure settings
    env._bodies["Delfi-C3"].set_radiation_pressure_interface(
        "Sun", simulation_setup.create_radiation_pressure_interface(
            rad_press_settings, "Delfi-C3", env._bodies))

    # # Create aerodynamic interface settings for Delfi-C3.
    # env["Delfi-C3"]._set_aerodynamic_coefficient_interface(
    #     reference_area=4.0,
    #     coefficients=1.2 * np.ones(3),
    #     coefficients_in_aerodynamic_frame=True,
    #     coefficients_in_negative_axis_direction=True,
    #     model="constant_coefficients")
    #
    # # Create radiation pressure interface settings for Delfi-C3.
    # env["Delfi-C3"]._set_radiation_pressure_interface(
    #     reference_area=4.0,
    #     coefficient=1.2,
    #     occulting=["Earth"],
    #     source="Sun",
    #     model="cannon_ball")

    ###########################################################################
    # FINALIZE BODIES #########################################################
    ###########################################################################

    # Finalize body creation.
    env.finalize_environment()

    ###########################################################################
    # CREATE ACCELERATIONS ####################################################
    ###########################################################################

    # Define bodies that are propagated.
    bodies_to_propagate = ["Delfi-C3"]

    # Define central bodies.
    central_bodies = ["Earth"]

    # Define unique (Sun, Earth) accelerations acting on Delfi-C3.
    accelerations_of_delfi_c3 = dict(
        Sun=
        [
            simulation_setup.Acceleration.canon_ball_radiation_pressure()
            # AccelerationSettings(AvailableAcceleration.cannon_ball_radiation_pressure) # LEGACY DESIGN.
        ],
        Earth=
        [
            simulation_setup.Acceleration.spherical_harmonic_gravity(5, 5),
            # SphericalHarmonicAccelerationSettings(5, 5), # LEGACY DESIGN.

            simulation_setup.Acceleration.aerodynamic()
            # AccelerationSettings(AvailableAcceleration.aerodynamic) # LEGACY DESIGN.
        ])

    # Define other point mass accelerations acting on Delfi-C3.
    for other in set(bodies_to_create).difference({"Sun", "Earth"}):
        accelerations_of_delfi_c3[other] = [
            simulation_setup.Acceleration.point_mass_gravity()]

    # Create global accelerations dictionary.
    acceleration_dict = {"Delfi-C3": accelerations_of_delfi_c3}

    # Create acceleration models.

    acceleration_models = simulation_setup.create_acceleration_models_dict(
        env._bodies,
        acceleration_dict,
        bodies_to_propagate,
        central_bodies)

    ###########################################################################
    # CREATE PROPAGATION SETTINGS #############################################
    ###########################################################################

    # Set initial conditions for the Asterix satellite that will be
    # propagated in this simulation. The initial conditions are given in
    # Keplerian elements and later on converted to Cartesian elements.

    # Set Keplerian elements for Delfi-C3
    earth_gravitational_parameter = env[
        "Earth"].gravity_field_model.get_gravitational_parameter()

    # REVISED CONTEMPORARY DESIGN.
    system_initial_state = elements.keplerian2cartesian(
        mu=earth_gravitational_parameter,
        a=7500.0E3,
        ecc=0.1,
        inc=np.deg2rad(85.3),
        raan=np.deg2rad(23.4),
        argp=np.deg2rad(235.7),
        nu=np.deg2rad(139.87))

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

    t0 = time.time()
    # Create simulation object and propagate dynamics.
    print(env._bodies)
    dynamics_simulator = propagators.SingleArcDynamicsSimulator(
        env._bodies, integrator_settings, propagator_settings, True)
    print("dynamics_simulator: ", time.time() - t0)
    result = dynamics_simulator.get_equations_of_motion_numerical_solution()

    ###########################################################################
    # PRINT INITIAL AND FINAL STATES ##########################################
    ###########################################################################

    print(
        f"""
Single Earth-Orbiting Satellite Example.
The initial position vector of Delfi-C3 is [km]: \n{
        result[simulation_start_epoch][:3] / 1E3}
The initial velocity vector of Delfi-C3 is [km/s]: \n{
        result[simulation_start_epoch][3:] / 1E3}
After {simulation_end_epoch} seconds the position vector of Delfi-C3 is [km]: \n{
        result[simulation_end_epoch][:3] / 1E3}
And the velocity vector of Delfi-C3 is [km/s]: \n{
        result[simulation_start_epoch][3:] / 1E3}
        """
    )

    ###########################################################################
    # SAVE RESULTS ############################################################
    ###########################################################################

    io.save2txt(
        solution=result,
        filename="singlePerturbedSatellitePropagationHistory.dat",
        directory="./tutorial_2_prototype",
    )

    # Final statement (not required, though good practice in a __main__).
    return 0


if __name__ == "__main__":
    main()
