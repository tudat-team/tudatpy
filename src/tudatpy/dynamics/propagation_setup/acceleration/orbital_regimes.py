from tudatpy.dynamics import propagation_setup, environment_setup


class GetAccelerationSettingsPerRegime:
    """
    A utility class for generating and modifying TudatPy acceleration settings
    for spacecraft in different orbital regimes (LEO, MEO, GEO, or OTHER).

    This class provides:
    - Default acceleration settings for standard regimes.
    - The ability to add custom accelerations and required body settings.
    - The ability to remove specific accelerations.

    Attributes
    ----------
    required_body_settings : dict
        Mapping of acceleration keywords to a list of required body settings
        (e.g., 'gravity_field_settings', 'ephemeris_settings') needed to
        support that acceleration.
    required_accelerations : dict
        Mapping of acceleration keywords to TudatPy acceleration settings objects.
    """

    def __init__(self):
        self.required_body_settings = {
            'central_gravity': ['ephemeris_settings', 'gravity_field_settings']
        }

        self.required_accelerations = {
            'central_gravity': propagation_setup.acceleration.point_mass_gravity()
        }

    def get_acceleration_settings(self,
                                  bodies,
                                  body_settings,
                                  object_name,
                                  orbital_regime,
                                  aerodynamics=True,
                                  radiation_pressure=True,
                                  add_forces=None,
                                  drop_forces=None):
        """
        Retrieve the acceleration settings for a given orbital regime, with optional
        additions or removals of forces.

        Parameters
        ----------
        bodies : tudatpy.kernel.simulation.environment.SystemOfBodies
            Current environment bodies object.
        body_settings : tudatpy.kernel.simulation.environment_setup.BodyListSettings
            Settings for all bodies in the simulation.
        object_name : str
            Name of the propagated object (spacecraft).
        orbital_regime : str
            Regime identifier. Must be one of:
            'LEO_REGIME', 'MEO_REGIME', 'GEO_REGIME', or 'OTHER'.
        aerodynamics : bool, optional
            If True, include aerodynamic acceleration from Earth (default=True).
        radiation_pressure : bool, optional
            If True, include solar radiation pressure from Sun (default=True).
        add_forces : dict, optional
            Dictionary mapping body names to lists of acceleration keywords to add.
            Example: {"Jupiter": ["central_gravity"]}.
        drop_forces : dict, optional
            Dictionary mapping body names to lists of keywords to remove.
            Example: {"Earth": ["aerodynamic"]}.

        Returns
        -------
        tuple
            acceleration_settings : dict
                Mapping of body names to lists of TudatPy acceleration settings.
            bodies : tudatpy.kernel.simulation.environment.SystemOfBodies
                Updated environment bodies object.
        """

        # Default per regime
        if orbital_regime == 'LEO_REGIME':
            acceleration_settings = self.get_LEO_acceleration_settings(aerodynamics, radiation_pressure)
        elif orbital_regime == 'MEO_REGIME':
            acceleration_settings = self.get_MEO_acceleration_settings(aerodynamics, radiation_pressure)
        elif orbital_regime == 'GEO_REGIME':
            acceleration_settings = self.get_GEO_acceleration_settings(aerodynamics, radiation_pressure)
        elif orbital_regime == 'OTHER':
            acceleration_settings = self.get_GEO_acceleration_settings(aerodynamics, radiation_pressure)
        else:
            raise ValueError(f"Unknown orbital regime: {orbital_regime}")

        # Apply additions
        if add_forces:
            for body, accelerations_to_add in add_forces.items():
                for acceleration_to_add in accelerations_to_add:
                    if not bodies.does_body_exist(body):
                        body_settings.add_empty_settings(body)
                    if acceleration_to_add in self.required_body_settings.keys():
                        body_settings = self.add_required_body_settings(
                            body,
                            self.required_body_settings[acceleration_to_add],
                            body_settings
                        )
                self.append_required_accelerations(body, acceleration_settings, accelerations_to_add)

        # Create updated system of bodies
        bodies = environment_setup.create_system_of_bodies(body_settings)

        # Apply removals
        if drop_forces:
            for body, drop_list in drop_forces.items():
                if body in acceleration_settings:
                    acceleration_settings[body] = [
                        a for a in acceleration_settings[body]
                        if not any(keyword in str(a) for keyword in drop_list)
                    ]

        return acceleration_settings, bodies

    def add_required_body_settings(self, body, body_settings_strings_to_add, body_settings):
        """
        Add required settings (e.g., ephemeris, gravity field) for a given body.

        Parameters
        ----------
        body : str
            Name of the body to modify.
        body_settings_strings_to_add : list
            List of settings keywords to add (from `self.required_body_settings`).
        body_settings : tudatpy.kernel.simulation.environment_setup.BodyListSettings
            Settings object for all bodies.

        Returns
        -------
        tudatpy.kernel.simulation.environment_setup.BodyListSettings
            Updated body settings.
        """
        for body_settings_string in body_settings_strings_to_add:
            if body_settings_string == 'gravity_field_settings':
                body_settings.get(body).gravity_field_settings = environment_setup.gravity_field.central_spice(body)
            elif body_settings_string == 'ephemeris_settings':
                body_settings.get(body).ephemeris_settings = environment_setup.ephemeris.direct_spice(
                    frame_origin='Earth',
                    frame_orientation='J2000'
                )
        return body_settings

    def append_required_accelerations(self, body, acceleration_settings, accelerations_to_add):
        """
        Append specified accelerations to the acceleration settings for a body.

        Parameters
        ----------
        body : str
            Name of the body exerting the acceleration.
        acceleration_settings : dict
            Mapping of body names to lists of acceleration settings.
        accelerations_to_add : list
            List of acceleration keywords to add.
        """
        for acceleration in accelerations_to_add:
            new_accel = self.required_accelerations[acceleration]
            if body in acceleration_settings:
                acceleration_settings[body].append(new_accel)
            else:
                acceleration_settings[body] = [new_accel]

    def get_LEO_acceleration_settings(self, aerodynamics, radiation_pressure):
        """
        Default accelerations for LEO regime.

        Parameters
        ----------
        aerodynamics : bool
            Whether to include aerodynamic acceleration from Earth.
        radiation_pressure : bool
            Whether to include solar radiation pressure from Sun.

        Returns
        -------
        dict
            Mapping of body names to lists of acceleration settings.
        """
        LEO_acceleration_settings_dict = dict(
            Earth=[propagation_setup.acceleration.spherical_harmonic_gravity(5, 5)],
            Sun=[propagation_setup.acceleration.point_mass_gravity()],
            Moon=[propagation_setup.acceleration.point_mass_gravity()],
        )
        if aerodynamics:
            LEO_acceleration_settings_dict['Earth'].append(propagation_setup.acceleration.aerodynamic())
        if radiation_pressure:
            LEO_acceleration_settings_dict['Sun'].append(propagation_setup.acceleration.radiation_pressure())
        return LEO_acceleration_settings_dict

    def get_MEO_acceleration_settings(self, aerodynamics, radiation_pressure):
        """
        Default accelerations for MEO regime.

        Parameters
        ----------
        aerodynamics : bool
            Whether to include aerodynamic acceleration from Earth.
        radiation_pressure : bool
            Whether to include solar radiation pressure from Sun.

        Returns
        -------
        dict
            Mapping of body names to lists of acceleration settings.
        """
        MEO_acceleration_settings_dict = dict(
            Earth=[propagation_setup.acceleration.spherical_harmonic_gravity(5, 5)],
            Sun=[propagation_setup.acceleration.point_mass_gravity()],
            Moon=[propagation_setup.acceleration.spherical_harmonic_gravity(5, 5)],
        )
        if aerodynamics:
            MEO_acceleration_settings_dict['Earth'].append(propagation_setup.acceleration.aerodynamic())
        if radiation_pressure:
            MEO_acceleration_settings_dict['Sun'].append(propagation_setup.acceleration.radiation_pressure())
        return MEO_acceleration_settings_dict

    def get_GEO_acceleration_settings(self, aerodynamics, radiation_pressure):
        """
        Default accelerations for GEO regime.

        Parameters
        ----------
        aerodynamics : bool
            Whether to include aerodynamic acceleration from Earth.
        radiation_pressure : bool
            Whether to include solar radiation pressure from Sun.

        Returns
        -------
        dict
            Mapping of body names to lists of acceleration settings.
        """
        GEO_acceleration_settings_dict = dict(
            Earth=[propagation_setup.acceleration.point_mass_gravity()],
            Sun=[propagation_setup.acceleration.point_mass_gravity()],
            Moon=[propagation_setup.acceleration.spherical_harmonic_gravity(5, 5)],
            Jupiter=[propagation_setup.acceleration.point_mass_gravity()],
            Mars=[propagation_setup.acceleration.point_mass_gravity()],
            Venus=[propagation_setup.acceleration.point_mass_gravity()],
            Saturn=[propagation_setup.acceleration.point_mass_gravity()]
        )
        if aerodynamics:
            GEO_acceleration_settings_dict['Earth'].append(propagation_setup.acceleration.aerodynamic())
        if radiation_pressure:
            GEO_acceleration_settings_dict['Sun'].append(propagation_setup.acceleration.radiation_pressure())
        return GEO_acceleration_settings_dict
