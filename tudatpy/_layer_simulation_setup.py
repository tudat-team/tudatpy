from .core._basic_astrodynamics import AvailableAcceleration
from .core._simulation_setup import AccelerationSettings
from .core._simulation_setup import SphericalHarmonicAccelerationSettings


def modify_simulation_setup(simulation_setup):
    class Acceleration:
        """

        """

        @staticmethod
        def spherical_harmonic_gravity(maximum_degree=8, maximum_order=8):
            """

            Parameters
            ----------
            maximum_degree
            maximum_order

            Returns
            -------

            """
            return SphericalHarmonicAccelerationSettings(maximum_degree, maximum_order)

        @staticmethod
        def aerodynamic():
            return AccelerationSettings(AvailableAcceleration.aerodynamic)

        @staticmethod
        def canon_ball_radiation_pressure():
            return AccelerationSettings(AvailableAcceleration.cannon_ball_radiation_pressure)

        @staticmethod
        def point_mass_gravity():
            return AccelerationSettings(AvailableAcceleration.point_mass_gravity)

        # undefined_acceleration
        # point_mass_gravity
        # central_gravity
        # aerodynamic
        # cannon_ball_radiation_pressure
        # spherical_harmonic_gravity
        # mutual_spherical_harmonic_gravity
        # third_body_point_mass_gravity
        # third_body_central_gravity
        # third_body_spherical_harmonic_gravity
        # third_body_mutual_spherical_harmonic_gravity
        # thrust_acceleration
        # relativistic_correction_acceleration
        # empirical_acceleration
        # direct_tidal_dissipation_in_central_body_acceleration
        # direct_tidal_dissipation_in_orbiting_body_acceleration

    setattr(simulation_setup, "Acceleration", Acceleration())
