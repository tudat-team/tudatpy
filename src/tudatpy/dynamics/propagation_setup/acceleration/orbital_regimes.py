from mypyc.crash import catch_errors
from tudatpy.dynamics import propagation_setup, environment_setup
from xxlimited import Error


class GetAccelerationSettingsPerRegime:
    def __init__(self):
        pass

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
        add_forces: dict with format {"Body": [acceleration_settings,...]}
        drop_forces: dict with format {"Body": ["force_type_keyword",...]}
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

        ## Apply additions
        if add_forces:
            for body, new_accels in add_forces.items():
                if not bodies.does_body_exist(body):
                    body_settings.add_empty_settings(body)
                    body_settings.get(body).gravity_field_settings = environment_setup.gravity_field.central_spice(body)
                    body_settings.get(body).ephemeris_settings = environment_setup.ephemeris.direct_spice(frame_origin = 'Earth', frame_orientation = 'J2000')
                    bodies = environment_setup.create_system_of_bodies(body_settings)

                if body not in acceleration_settings:
                    acceleration_settings[body] = [] # this is not enough, since we also need to define settings in the first place
                acceleration_settings[body].extend(new_accels)

        # Apply removals
        if drop_forces:
            for body, drop_list in drop_forces.items():
                if body in acceleration_settings:
                    acceleration_settings[body] = [
                        a for a in acceleration_settings[body]
                        if not any(keyword in str(a) for keyword in drop_list)
                    ]

        return acceleration_settings, bodies

    def get_LEO_acceleration_settings(self, aerodynamics, radiation_pressure):


        LEO_acceleration_settings_dict = dict(
            Earth=[propagation_setup.acceleration.spherical_harmonic_gravity(5, 5)],
            Sun = [propagation_setup.acceleration.point_mass_gravity()],
            Moon = [propagation_setup.acceleration.point_mass_gravity()],
        )

        if aerodynamics:
            LEO_acceleration_settings_dict['Earth'].append(propagation_setup.acceleration.aerodynamic())

        if radiation_pressure:
            LEO_acceleration_settings_dict['Sun'].append(propagation_setup.acceleration.radiation_pressure())

        return LEO_acceleration_settings_dict

    def get_MEO_acceleration_settings(self, aerodynamics, radiation_pressure):

        MEO_acceleration_settings_dict = dict(
            Earth=[
                propagation_setup.acceleration.spherical_harmonic_gravity(5, 5)],
            Sun = [propagation_setup.acceleration.point_mass_gravity()],
            Moon = [propagation_setup.acceleration.spherical_harmonic_gravity(5, 5)],
        )

        if aerodynamics:
            MEO_acceleration_settings_dict['Earth'].append(propagation_setup.acceleration.aerodynamic())

        if radiation_pressure:
            MEO_acceleration_settings_dict['Sun'].append(propagation_setup.acceleration.radiation_pressure())

        return MEO_acceleration_settings_dict

    def get_GEO_acceleration_settings(self, aerodynamics, radiation_pressure):

        GEO_acceleration_settings_dict = dict(
            Earth=[propagation_setup.acceleration.point_mass_gravity()],
            Sun = [propagation_setup.acceleration.point_mass_gravity()],
            Moon = [propagation_setup.acceleration.spherical_harmonic_gravity(5, 5)],
            Jupiter = [propagation_setup.acceleration.point_mass_gravity()],
            Mars = [propagation_setup.acceleration.point_mass_gravity()],
            Venus = [propagation_setup.acceleration.point_mass_gravity()],
            Saturn = [propagation_setup.acceleration.point_mass_gravity()]
        )

        if aerodynamics:
            GEO_acceleration_settings_dict['Earth'].append(propagation_setup.acceleration.aerodynamic())

        if radiation_pressure:
            GEO_acceleration_settings_dict['Sun'].append(propagation_setup.acceleration.radiation_pressure())

        return GEO_acceleration_settings_dict