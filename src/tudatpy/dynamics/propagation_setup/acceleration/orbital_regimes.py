from tudatpy.dynamics import propagation_setup

class GetAccelerationSettingsPerRegime:
    def __init__(self, CreateEphemeris_Settings):

        self.CreateEphemeris_Settings = CreateEphemeris_Settings
        pass

    def get_LEO_acceleration_settings(self):

        LEO_acceleration_settings_dict = dict(
            Earth=[
                propagation_setup.acceleration.spherical_harmonic_gravity(5, 5)]
            ,
            Sun = [propagation_setup.acceleration.point_mass_gravity(),
                   propagation_setup.acceleration.radiation_pressure()
                   ],
            Moon = [propagation_setup.acceleration.point_mass_gravity()],
        )
        return LEO_acceleration_settings_dict

    def get_MEO_acceleration_settings(self):

        MEO_acceleration_settings_dict = dict(
            Earth=[
                propagation_setup.acceleration.spherical_harmonic_gravity(5, 5),
                propagation_setup.acceleration.aerodynamic()
            ],
            Sun = [propagation_setup.acceleration.point_mass_gravity(),
                   propagation_setup.acceleration.radiation_pressure()
                   ],
            Moon = [propagation_setup.acceleration.spherical_harmonic_gravity(5, 5)],
        )
        return MEO_acceleration_settings_dict

    def get_GEO_acceleration_settings(self):
        GEO_acceleration_settings_dict = dict(
            Earth=[propagation_setup.acceleration.point_mass_gravity()],
            Sun = [propagation_setup.acceleration.point_mass_gravity(), propagation_setup.acceleration.radiation_pressure()],
            Moon = [propagation_setup.acceleration.spherical_harmonic_gravity(5, 5)],
            Jupiter = [propagation_setup.acceleration.point_mass_gravity()],
            Mars = [propagation_setup.acceleration.point_mass_gravity()],
            Venus = [propagation_setup.acceleration.point_mass_gravity()],
            Saturn = [propagation_setup.acceleration.point_mass_gravity()]
        )
        return GEO_acceleration_settings_dict