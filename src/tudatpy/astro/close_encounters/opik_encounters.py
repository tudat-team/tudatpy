import numpy as np
import math
from tudatpy import constants


class OpikEncounter:
    """
    Idealized Opik encounter. Given (a, e, i), the planet's mass/radius/orbit,
    and this encounter's b and psi, everything else -- theta, phi, U, gamma,
    theta', phi', a', e', i' -- follows analytically. No ephemerides needed.

    b and psi describe the specific encounter (Carusi et al. 1990, Section 3):
    they are given directly, or sampled for Monte Carlo use (psi is uniform
    over [0, 2*pi) in Opik's own statistical treatment -- see sample_psi()).
    """

    GRAVITATIONAL_CONSTANT = constants.GRAVITATIONAL_CONSTANT  # [m^3 kg^-1 s^-2]
    SOLAR_MASS = 1.988475e30  # [KG]

    def __init__(
        self,
        small_body_orbital_elements,  # Opik-normalized: [a (dimensionless), e, i (rad)]
        central_body_mass,  # [KG] planet physical mass
        central_body_radius,  # [m] planet's physical radius
        central_body_orbital_radius,  # [m] planet's heliocentric semi-major axis
        b=None,  # [m] this encounter's impact parameter
        psi=None,  # [rad] this encounter's orientation angle
    ):
        self.a, self.e, self.i = small_body_orbital_elements[0:3]
        self.central_body_mass = central_body_mass
        self.central_body_radius = central_body_radius
        self.central_body_orbital_radius = central_body_orbital_radius
        self.solar_mass = self.SOLAR_MASS

        self._b = b
        self._psi = psi

    @staticmethod
    def sample_psi(rng=None):
        """Draws psi ~ Uniform(0, 2*pi), Opik's own statistical assumption
        for the encounter orientation angle."""
        rng = rng or np.random.default_rng()
        return rng.uniform(0.0, 2 * np.pi)

    @property
    def sun_gravitational_parameter(self):
        """mu of the Sun [m^3/s^2]."""
        return self.GRAVITATIONAL_CONSTANT * self.solar_mass

    @property
    def velocity_scale(self):
        """Planet's circular velocity [m/s] -- the unit U is normalized by."""
        return math.sqrt(self.sun_gravitational_parameter / self.central_body_orbital_radius)

    @property
    def a_physical(self):
        """Small body's semi-major axis, meters."""
        return self.a * self.central_body_orbital_radius

    @property
    def tisserand_parameter(self):
        return 1 / self.a + 2 * math.sqrt(self.a * (1 - self.e**2)) * np.cos(self.i)

    @property
    def planetocentric_velocity_magnitude_normalized(self):
        return math.sqrt(3 - self.tisserand_parameter)

    @property
    def planetocentric_velocity_magnitude_physical(self):
        """|U|, m/s."""
        return self.planetocentric_velocity_magnitude_normalized * self.velocity_scale

    @property
    def planetocentric_velocity_vector_normalized(self):
        Ux = math.sqrt(2 - 1 / self.a - self.a * (1 - self.e**2))
        Uy = math.sqrt(self.a * (1 - self.e**2)) * np.cos(self.i) - 1
        Uz = math.sqrt(self.a * (1 - self.e**2)) * np.sin(self.i)
        return [Ux, Uy, Uz]

    @property
    def planetocentric_velocity_vector_physical(self):
        """Ux, Uy, Uz, m/s."""
        return [
            component * self.velocity_scale
            for component in self.planetocentric_velocity_vector_normalized
        ]

    @property
    def incoming_asymptote_angles_vector(self):
        Ux, Uy, Uz = self.planetocentric_velocity_vector_normalized
        U = self.planetocentric_velocity_magnitude_normalized
        theta = np.arccos(Uy / U)
        phi = np.arctan2(Ux, Uz)
        return [theta, phi]

    @property
    def collision_limiting_impact_parameter_normalized(self):
        """b0: the largest b that would still graze the planet's physical
        surface (gravitational focusing). This is NOT the specific
        encounter's b."""
        Q = self.central_body_radius / self.central_body_orbital_radius  # both in meters
        U = self.planetocentric_velocity_magnitude_normalized
        S = math.sqrt(2 * self.central_body_mass / (self.solar_mass * Q))  # mass -> ratio
        return Q * math.sqrt(1 + S**2 / U**2)

    @property
    def collision_limiting_impact_parameter_physical(self):
        """b0, meters."""
        return (
            self.collision_limiting_impact_parameter_normalized * self.central_body_orbital_radius
        )

    @property
    def b(self):
        """This encounter's impact parameter, Opik-normalized."""
        if self._b is None:
            raise ValueError(
                "b not set: pass one to the constructor, or use EphemerisOpikEncounter."
            )
        return self._b / self.central_body_orbital_radius

    @property
    def b_physical(self):
        """This encounter's impact parameter, meters."""
        return self.b * self.central_body_orbital_radius

    @property
    def psi(self):
        """This encounter's orientation angle. Give one directly, or draw
        one with sample_psi() for Monte Carlo use."""
        if self._psi is None:
            raise ValueError(
                "psi not set: pass one to the constructor, or use EphemerisOpikEncounter."
            )
        return self._psi

    @property
    def deflection_angle(self):
        """
        Deflection angle of the encounter ($\gamma$, as named in Carusi et al., 1990)
        """
        U = self.planetocentric_velocity_magnitude_normalized
        if self.tisserand_parameter > 3:
            raise ValueError(
                f"Tisserand Parameter value is {self.tisserand_parameter}. "
                f"Opik Theory is not valid in this regime."
            )
        mass_ratio = self.central_body_mass / self.solar_mass  # eq. 10's m is a mass ratio, not kg
        return 2 * np.arctan(mass_ratio / (self.b * U**2))

    @property
    def outgoing_asymptote_angles_vector(self):
        """
        Post-Encounter asymptote $\theta$ and $phi$ as a function of the incoming ones, $\gamma$ and $\psi$
        """
        theta, phi = self.incoming_asymptote_angles_vector
        gamma = self.deflection_angle
        psi = self.psi

        cos_theta_out = np.cos(theta) * np.cos(gamma) + np.sin(theta) * np.sin(gamma) * np.cos(psi)
        theta_out = np.arccos(cos_theta_out)

        sin_chi_num = np.sin(gamma) * np.sin(psi)
        cos_chi_den = np.sin(theta) * np.cos(gamma) - np.cos(theta) * np.sin(gamma) * np.cos(psi)
        chi = np.arctan2(sin_chi_num, cos_chi_den)

        phi_out = phi - chi
        return [theta_out, phi_out]

    @property
    def outgoing_orbital_elements(self):
        """
        Post-encounter (a', e', i'), Opik-normalized, same convention as
        self.a/e/i (Carusi et al. 1990). Inverts planetocentric_velocity_vector_normalized:

            Ux = U sin(theta) sin(phi)
            Uy = U cos(theta)
            Uz = U sin(theta) cos(phi)

            a  = 1 / (1 - U^2 - 2*Uy)                      [matches eq. (13)'s denominator]
            e  = sqrt(1 - [(1+Uy)^2 + Uz^2] / a)
            i  = atan2(Uz, 1+Uy)
        """
        theta_out, phi_out = self.outgoing_asymptote_angles_vector
        U = self.planetocentric_velocity_magnitude_normalized

        # Ux_out not needed -- a', e', i' only depend on Uy_out and Uz_out
        Uy_out = U * np.cos(theta_out)
        Uz_out = U * np.sin(theta_out) * np.cos(phi_out)

        a_out = 1 / (1 - U**2 - 2 * Uy_out)
        e_out = math.sqrt(max(0.0, 1 - ((1 + Uy_out) ** 2 + Uz_out**2) / a_out))
        i_out = np.arctan2(Uz_out, 1 + Uy_out)

        return a_out, e_out, i_out

    @property
    def outgoing_orbital_elements_physical(self):
        """(a', e', i') with a' in meters -- e' and i' are already
        dimensionless/angular, so only a' needs converting."""
        a_out, e_out, i_out = self.outgoing_orbital_elements
        return a_out * self.central_body_orbital_radius, e_out, i_out


class EphemerisOpikEncounter(OpikEncounter):
    """
    Real Opik encounter. theta, phi, U, b, and psi are all pulled from a pair
    of tudatpy Ephemeris objects at encounter_epoch, instead of the idealized
    (a, e, i)-only formulas -- deflection_angle, outgoing_asymptote_angles_vector,
    and outgoing_orbital_elements are inherited unchanged from OpikEncounter
    and run on these real values automatically.

    small_body_ephemeris and planet_ephemeris can be ANY tudatpy Ephemeris --
    Kepler, direct_spice, tabulated, a numerically propagated one, whatever --
    as long as it supports .cartesian_state(epoch). Use from_keplerian_elements()
    for the common case of building plain two-body Kepler ephemerides.
    """

    def __init__(
        self,
        small_body_ephemeris,  # tudatpy Ephemeris
        planet_ephemeris,  # tudatpy Ephemeris
        central_body_mass,  # [KG] planet physical mass
        central_body_radius,  # [m] planet's physical radius
        encounter_epoch,
    ):
        from tudatpy.astro import element_conversion

        self.small_body_ephemeris = small_body_ephemeris
        self.planet_ephemeris = planet_ephemeris
        self.encounter_epoch = encounter_epoch

        mu_sun = self.GRAVITATIONAL_CONSTANT * self.SOLAR_MASS

        planet_state = np.asarray(planet_ephemeris.cartesian_state(encounter_epoch)).flatten()
        planet_keplerian = element_conversion.cartesian_to_keplerian(planet_state, mu_sun)
        central_body_orbital_radius = planet_keplerian[0]  # planet's a, meters

        body_state = np.asarray(small_body_ephemeris.cartesian_state(encounter_epoch)).flatten()
        body_keplerian = element_conversion.cartesian_to_keplerian(body_state, mu_sun)
        a = body_keplerian[0] / central_body_orbital_radius  # Opik-normalized
        e, i = body_keplerian[1], body_keplerian[2]

        super().__init__(
            [a, e, i], central_body_mass, central_body_radius, central_body_orbital_radius
        )

    @classmethod
    def from_keplerian_elements(
        cls,
        small_body_orbital_elements,  # tudat: [a(m), e, i(rad), argp(rad), raan(rad), true_anomaly(rad)]
        planet_orbital_elements,  # tudat: [a(m), e, i(rad), argp(rad), raan(rad), true_anomaly(rad)]
        central_body_mass,  # [KG] planet physical mass
        central_body_radius,  # [m] planet's physical radius
        small_body_epoch=0.0,
        planet_epoch=0.0,
        encounter_epoch=None,
        frame_origin="Sun",
        frame_orientation="ECLIPJ2000",
    ):
        """Convenience constructor for the common case: builds plain two-body
        Kepler ephemerides from raw orbital elements, then delegates to the
        main constructor."""
        from tudatpy.dynamics import environment_setup

        mu_sun = cls.GRAVITATIONAL_CONSTANT * cls.SOLAR_MASS

        small_body_settings = environment_setup.ephemeris.keplerian(
            small_body_orbital_elements, small_body_epoch, mu_sun, frame_origin, frame_orientation
        )
        small_body_ephemeris = environment_setup.create_body_ephemeris(
            small_body_settings, "SmallBody"
        )

        planet_settings = environment_setup.ephemeris.keplerian(
            planet_orbital_elements, planet_epoch, mu_sun, frame_origin, frame_orientation
        )
        planet_ephemeris = environment_setup.create_body_ephemeris(planet_settings, "Planet")

        return cls(
            small_body_ephemeris,
            planet_ephemeris,
            central_body_mass,
            central_body_radius,
            encounter_epoch,
        )

    @property
    def encounter_state_vectors(self):
        """r_body, v_body, r_planet, v_planet at self.encounter_epoch (the
        MOID time), evaluated from the two ephemerides above."""
        body_state = np.asarray(
            self.small_body_ephemeris.cartesian_state(self.encounter_epoch)
        ).flatten()
        planet_state = np.asarray(
            self.planet_ephemeris.cartesian_state(self.encounter_epoch)
        ).flatten()

        r_body, v_body = body_state[:3], body_state[3:6]
        r_planet, v_planet = planet_state[:3], planet_state[3:6]
        return r_body, v_body, r_planet, v_planet

    @property
    def _local_frame_state(self):
        """r_rel, v_rel at encounter_epoch, projected onto the planet's local
        frame (X=radial, Y=along-track, Z=orbit-normal)."""
        r_body, v_body, r_planet, v_planet = map(np.array, self.encounter_state_vectors)

        X = r_planet / np.linalg.norm(r_planet)  # X axis
        Z = np.cross(r_planet, v_planet)  # Z axis
        Z /= np.linalg.norm(Z)
        Y = np.cross(Z, X)  # Y axis

        to_frame = lambda v: np.array([v @ X, v @ Y, v @ Z])  # projects onto the planet's frame
        r_rel = to_frame(r_body - r_planet)  # relative position
        v_rel = to_frame(v_body - v_planet)  # relative velocity [Ux, Uy, Uz]
        return r_rel, v_rel

    @property
    def planetocentric_velocity_vector_normalized(self):
        """v_rel in the local frame, normalized by the planet's own circular
        velocity -- everything downstream (deflection_angle, a'/e'/i', ...)
        assumes U is dimensionless, same convention as the idealized model."""
        _, v_rel = self._local_frame_state
        v_planet_circular = math.sqrt(
            self.sun_gravitational_parameter / self.central_body_orbital_radius
        )
        return list(v_rel / v_planet_circular)

    @property
    def planetocentric_velocity_magnitude_normalized(self):
        return float(np.linalg.norm(self.planetocentric_velocity_vector_normalized))

    @property
    def b(self):
        """This encounter's real impact parameter: |h| / |v_rel|."""
        r_rel, v_rel = self._local_frame_state
        h = np.cross(r_rel, v_rel)  # angular momentum vector
        b_physical = np.linalg.norm(h) / np.linalg.norm(v_rel)
        return b_physical / self.central_body_orbital_radius

    @property
    def psi(self):
        """Psi in the tangent plane at the incoming asymptote direction
        (tangent plane of U), from the real position/velocity vectors.
        """
        theta, phi = self.incoming_asymptote_angles_vector
        r_rel, v_rel = self._local_frame_state

        U_hat = v_rel / np.linalg.norm(v_rel)  # relative velocity unit vector
        h_hat = np.cross(r_rel, v_rel)
        h_hat /= np.linalg.norm(h_hat)  # angular momentum unit vector
        b_hat = np.cross(h_hat, U_hat)  # impact parameter unit vector, tangent to U space

        # decompose b in theta and phi components
        du_dtheta = np.array(
            [np.cos(theta) * np.sin(phi), -np.sin(theta), np.cos(theta) * np.cos(phi)]
        )
        du_dphi = np.array([np.cos(phi), 0, -np.sin(phi)])

        # compute arctan2(dU/dphi, -dU/dtheta) = psi
        return np.arctan2(b_hat @ du_dphi, -(b_hat @ du_dtheta))
