from .core import _spice_interface as spice_interface
from .core import _simulation_setup as simulation_setup
import warnings


class Environment:

    def __init__(self,
                 available_bodies,
                 user_defined_bodies=None,
                 frame_origin="SSB",
                 frame_orientation="ECLIPJ2000",
                 custom_body_settings=None,
                 alternative_kernels=None,
                 start_epoch=None,
                 end_epoch=None,
                 fixed_dt=None,
                 epoch_margin=300.0):
        """

        Parameters
        ----------
        available_bodies : list[str]
            List of bodies to be added to the environment that are part
            of the available bodies cataloged in the tudatBundle.
        user_defined_bodies : list[str]
            List of user-defined bodies to be added to the environment
            which are not part of the available bodies in the tudatBundle.
        frame_origin : str
            Frame origin for simulation output results.
        frame_orientation : str
            Frame orientation for simulation output results.
        custom_body_settings : dict, optional
            Custom body settings which update in place of the default.
        alternative_kernels : str, optional

        Examples
        --------
        available_bodies= ["Sun", "Earth", "Moon", "Mars", "Venus"]
        user_defined_bodies = ["Delfi-C3"]
        custom_body_settings = {"Delfi-C3":dict(constant_mass=400.0)}
        tutorial_2_environment = Environment(
           available_bodies=available_bodies,
           user_defined_bodies=user_defined_bodies,
           custom_body_settings=custom_body_settings)

        """

        self._available_bodies = available_bodies
        self._user_defined_bodies = user_defined_bodies
        self._frame_origin = frame_origin
        self._frame_orientation = frame_orientation
        self._custom_body_settings = custom_body_settings
        self._alternative_kernels = alternative_kernels
        self._start_epoch = start_epoch
        self._end_epoch = end_epoch
        self._epoch_margin = abs(epoch_margin)
        self._fixed_dt = fixed_dt
        self._reset_kernels = True
        self._reset_bodies = True
        self._reset_vehicles = True
        self._body_settings = None
        self._environment_finalized = False

    def _check_initialize_environment(self):

        if self._reset_kernels:

            # Load standard spice kernels or alternative?
            if self._alternative_kernels:
                # Alternative if defined.
                spice_interface.load_standard_spice_kernels(
                    self._alternative_kernels)
            else:
                # Standard if alternative is not defined (None).
                spice_interface.load_standard_spice_kernels()

            # If kernels are reset, bodies should be too.
            self._reset_bodies = True

        if self._reset_bodies:

            if (self._start_epoch is not None) and (self._end_epoch is not None) and (self._fixed_dt is not None):
                self._body_settings = simulation_setup.get_default_body_settings(
                    self._available_bodies,
                    self._start_epoch - self._epoch_margin,
                    self._end_epoch + self._epoch_margin,
                    self._fixed_dt)
            else:
                if [self._start_epoch, self._end_epoch, self._fixed_dt].count(1) >= 1:
                    warnings.warn(
                        "Ambiguous definition of simulation epoch parameters."
                        "Default body settings are being retrieved without"
                        "definition of start_epoch, end_epoch and fixed_dt, "
                        f"  start_epoch = {self._start_epoch},"
                        f"  end_epoch = {self._end_epoch},"
                        f"  fixed_dt = {self._fixed_dt}")

                # Get default body settings.
                self._body_settings = simulation_setup.get_default_body_settings(
                    self._available_bodies)

            for body in self._available_bodies:
                self._body_settings[body].ephemeris_settings.reset_frame_orientation(self._frame_orientation)
                self._body_settings[body].rotation_model_settings.reset_original_frame(self._frame_orientation)

            # If custom body settings are to be used over the default.
            if self._custom_body_settings:

                # Cycle through bodies in the custom settings.
                for body_str in self._custom_body_settings.keys():

                    # Check if body is in the available bodies.
                    if body_str in self._available_bodies:

                        custom_settings = self._custom_body_settings[body_str]

                        for settings_str in custom_settings.keys():

                            try:
                                setattr(
                                    self._body_settings[body_str],
                                    settings_str,
                                    custom_settings[settings_str])
                            except AttributeError:
                                raise AttributeError(
                                    f"The settings attribute of the BodySettings "
                                    f"class does not exist. Available found at:"
                                    f"https://ggarrett13.github.io/tudatpy/_build/_src_modules/simulation_setup.html#tudatpy.core._simulation_setup.BodySettings"
                                )

                    # Not sure if this is correct for custom bodies,
                    # tutorial examples suggest it may be.
                    elif body_str in self._user_defined_bodies:
                        pass

                    # Body settings provided doesn't correspond to any
                    # defined in union[available, user_defined].
                    else:
                        raise EnvironmentError(
                            f"The body {body_str} contained in the custom body "
                            f"settings does not exist in the available bodies. "
                            f"{body_str} not in {self._available_bodies}")

            # If bodies are reset, vehicles must be re-added.
            self._reset_vehicles = True
            self._bodies = simulation_setup.create_bodies(self._body_settings)

        if self._reset_vehicles:
            for vehicle in self._user_defined_bodies:
                self._bodies[vehicle] = simulation_setup.Body()

        self._reset_kernels = False
        self._reset_bodies = False
        self._reset_vehicles = False

    @property
    def bodies(self):
        self._check_initialize_environment()
        # Finalize body creation.
        if not self._environment_finalized:
            simulation_setup.set_global_frame_body_ephemerides(
                self._bodies, self._frame_origin, self._frame_orientation)

        return self._bodies

    @property
    def vehicles(self):
        return self._user_defined_bodies

    @vehicles.setter
    def vehicles(self, x):
        self._reset_vehicles = True
        self._environment_finalized = False
        self._user_defined_bodies = x

    def finalize_environment(self,
                             frame_origin,
                             frame_orientation):
        self._environment_finalized = True
        simulation_setup.set_global_frame_body_ephemerides(
            self._bodies, self._frame_origin, self._frame_orientation)

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
            self.bodies[target_body_name].set_radiation_pressure_interface(
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
        coefficint_vector
        coefficients_in_aerodynmic_frame
        coefficients_in_negative_axis_direction
        model

        Returns
        -------
        None

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

