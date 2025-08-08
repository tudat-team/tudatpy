from mypyc.crash import catch_errors
from tudatpy.dynamics import propagation_setup
from xxlimited import Error


class GetAccelerationSettingsPerRegime:
    def __init__(self):
        pass

    def get_acceleration_settings(self,
                                  bodies,
                                  object_name,
                                  orbital_regime,
                                  aerodynamics = True,
                                  radiation_pressure = True
                                  ):

        if orbital_regime == 'LEO_REGIME':
            acceleration_settings = self.get_LEO_acceleration_settings(aerodynamics, radiation_pressure)
        elif orbital_regime == 'MEO_REGIME':
            acceleration_settings = self.get_MEO_acceleration_settings(aerodynamics, radiation_pressure)
        elif orbital_regime == 'GEO_REGIME':
            acceleration_settings = self.get_GEO_acceleration_settings(aerodynamics, radiation_pressure)
        elif orbital_regime == 'OTHER':
            acceleration_settings = self.get_GEO_acceleration_settings(aerodynamics, radiation_pressure)

        return acceleration_settings

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