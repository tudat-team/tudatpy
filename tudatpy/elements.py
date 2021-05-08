from .kernel.astro import conversion as _oec
import numpy as np


def spherical2cartesian(r,
                        lat,
                        lon,
                        speed,
                        fpa,
                        heading) -> np.ndarray:
    """
    Function to convert spherical state to cartesian.

    Parameters
    ----------
    r : float
        Position vector magnitude (m).
    lat : float
        Latitude (m).
    lon : float
        Longitude (m).
    speed : float
        Magnitude of velocity (m).
    fpa : float
        Flight path angle (rad).
    heading : float
        Heading angle (rad).

    Returns
    -------
    cartesian_state : ndarray
        Cartesian state represented as [Rx, Ry, Rz, Vx, Vz, Vy] with distance in (m) and speed in (m/s).

    """
    spherical_idx = _oec.SphericalOrbitalStateElementIndices
    spherical_state = np.zeros(6)
    spherical_state[int(spherical_idx.radius_index)] = r
    spherical_state[int(spherical_idx.latitude_index)] = lat
    spherical_state[int(spherical_idx.longitude_index)] = lon
    spherical_state[int(spherical_idx.speed_index)] = speed
    spherical_state[int(spherical_idx.flight_path_index)] = fpa
    spherical_state[int(spherical_idx.heading_angle_index)] = heading
    return _oec.convert_spherical_orbital_to_cartesian_state(spherical_state)


def keplerian2cartesian(mu: float = None,
                        sma: float = None,
                        ecc: float = None,
                        inc: float = None,
                        raan: float = None,
                        argp: float = None,
                        theta: float = None,
                        **kwargs) -> np.ndarray:
    """
    Function to convert Keplerian state to cartesian.

    Parameters
    ----------
    mu : float
        Standard gravitational parameter (m^3 / s^2).
        (alias = "gravitational_parameter")
    sma : float
        Semi-major axis (m).
        (alias = "semi_major_axis")
    ecc : float
        Eccentricity (-).
        (alias = "eccentricity")
    inc : float
        Inclination (rad).
        (alias = "inclination")
    raan : float
        Right Ascension of the Ascending Node (rad).
        (alias = "right_ascension_of_the_ascending_node")
    argp : float
        Argument of Perigee (rad).
        (alias = "argument_of_periapsis")
    theta : float
        True anomaly (rad).
        (alias = "true_anomaly")

    Returns
    -------
    cartesian_state : ndarray
        Cartesian state represented as [Rx, Ry, Rz, Vx, Vz, Vy] with distance in (m) and speed in (m/s).

    """
    # TODO: Add epoch overloaded input version.
    _mu = mu if mu else kwargs.get("gravitational_parameter")
    _sma = sma if sma else kwargs.get("semi_major_axis")
    _ecc = ecc if ecc else kwargs.get("eccentricity")
    _inc = inc if inc else kwargs.get("inclination")
    _raan = raan if raan else kwargs.get("right_ascension_of_the_ascending_node")
    _argp = argp if argp else kwargs.get("argument_of_periapsis")
    _theta = theta if theta else kwargs.get("true_anomaly")
    keplerian_idx = _oec.KeplerianElementIndices
    keplerian_state = np.zeros(6)
    keplerian_state[int(keplerian_idx.semi_major_axis_index)] = _sma
    keplerian_state[int(keplerian_idx.eccentricity_index)] = _ecc
    keplerian_state[int(keplerian_idx.inclination_index)] = _inc
    keplerian_state[int(keplerian_idx.longitude_of_ascending_node_index)] = _raan
    keplerian_state[int(keplerian_idx.argument_of_periapsis_index)] = _argp
    keplerian_state[int(keplerian_idx.true_anomaly_index)] = _theta
    return _oec.convert_keplerian_to_cartesian_elements(keplerian_state, _mu)


if __name__ == "__main__":
    res = spherical2cartesian(100, 2, 2, 1000, 1, 1)
    print(res)
