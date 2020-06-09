# from .core import _spice_interface as spice_interface
# from .core import _simulation_setup as simulation_setup

from tudatpy.kernel import spice_interface
from tudatpy.kernel import simulation_setup
from simulation_setup import Body as _Body
import warnings
import numpy as np


class Body(_Body):

    def __init__(self, state=np.zeros(6)):
        """

        Parameters
        ----------
        state
        """
        super().__init__(state)





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

    def set_radiation_pressure_interface(self,
                                         target_body_name,
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
            rad_press_settings = simulation_setup.CannonBallRadiationPressureInterfaceSettings(
                source, reference_area, coefficient, occulting)

            # Create and set radiation pressure settings
            self.set_radiation_pressure_interface(
                source, simulation_setup.create_radiation_pressure_interface(
                    rad_press_settings, target_body_name, self.bodies))
        else:
            raise NotImplementedError(f"Model {model} not implemented.")

    def set_aerodynamic_coefficient_interface(self,
                                              target_body_name,
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
            aero_c_settings = simulation_setup.ConstantAerodynamicCoefficientSettings(
                reference_area,
                coefficients,
                coefficients_in_aerodynamic_frame,
                coefficients_in_negative_axis_direction
            )
            # Create and set aerodynamic coefficients object
            self.bodies[target_body_name].set_aerodynamic_coefficient_interface(
                simulation_setup.create_aerodynamic_coefficient_interface(
                    aero_c_settings,
                    target_body_name)
            )
        else:
            raise NotImplementedError(f"Model {model} is not implemented.")


class Environment:

    def __init__(self,
                 tudatpy_bodies,
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
        self._tudatpy_bodies = [] or tudatpy_bodies

        # Sets as [] if custom_bodies is None.
        self._custom_bodies = [] or custom_bodies

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
        try:
            assert self._environment_finalized
        except AssertionError:
            raise EnvironmentError(
                "Environment bodies must be finalized before adding "
                "interfaces or modifying existing bodies."
            )

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
            print(key)
            print(self._custom_bodies + self._tudatpy_bodies)
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


if __name__ == "__main__":
    spice_interface.load_standard_spice_kernels()

    env = Environment(tudatpy_bodies=["Earth"],
                      custom_bodies=["Delfi"],
                      frame_origin="SSB",
                      frame_orientation="ECLIPJ2000",
                      start_epoch=None,
                      end_epoch=None,
                      epoch_margin=300.0)

    # After environment has been finalized, then member functions can be called.
    env.finalize_environment()

    # Method 1 of retrieving bodies.
    state = env["Delfi"].get_state()

    # Alternative way of retrieving bodies.
    state = env.Delfi.get_state()

    # Method 2 of setting custom bodies.
    env["DWelfi_n3xt"] = simulation_setup.Body()

    env.Delfi_one = simulation_setup.Body()

    # print(env)


