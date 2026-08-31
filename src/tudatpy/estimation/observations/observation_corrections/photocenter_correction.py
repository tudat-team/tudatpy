"""
Functions to calculate photocenter corrections to observations
"""
import numpy as np
from tudatpy.estimation.observations import ObservationCollection
from tudatpy.dynamics.environment_setup import create_ground_station_ephemeris
from tudatpy.dynamics.environment import SystemOfBodies
from numpy.linalg import norm
from ._correction_utils import _offset_vector_to_corrections, _unit, _apply_corrections_to_observation_collection
from tudatpy.constants import SPEED_OF_LIGHT
from collections.abc import Iterable

def _photocenter_offset(diameter : float,
                        body_wrt_observer : np.ndarray,
                        body_wrt_sun : np.ndarray,) -> np.ndarray:
    """Helper function to calculate photocenter offset at one epoch"""
    # Get unit vectors
    body_wrt_observer_unit = _unit(body_wrt_observer)
    body_wrt_sun_unit = _unit(body_wrt_sun)

    # Angle spanned by observer, body and Sun:
    solar_phase_angle = np.arccos(np.dot(body_wrt_observer_unit, body_wrt_sun_unit))

    offset_direction = - _unit((body_wrt_sun_unit - np.dot(body_wrt_sun_unit, body_wrt_observer_unit)
                    * body_wrt_observer_unit))
    # Calculate fraction of radius according to Fuentes-Munoz (2024):
    cot = lambda x: np.cos(x) / np.sin(x)
    num = 2 * (np.sin(solar_phase_angle) + (np.pi - solar_phase_angle) * np.cos(solar_phase_angle))
    denom = 3 * np.pi * (
            cot(solar_phase_angle / 2) - np.sin(solar_phase_angle / 2) * np.log(cot(solar_phase_angle / 4)))
    offset_ratio = num / denom # Fraction of body radius

    offset_magnitude = offset_ratio * (diameter/2) / norm(body_wrt_observer)

    return offset_direction * offset_magnitude


def photocenter_correction_angular_observations_sphere(
        observations: np.ndarray ,
        diameter: float,
        bodies: SystemOfBodies,
        body_name: str,
        observer_body_name: str,
        observer_reference_name: str | None = None
) -> np.ndarray:
    r"""
    Compute photocenter-barycenter corrections assuming a spherical shape model.

    Compute photocenter-barycenter corrections assuming a spherical shape model, according to Fuentes-Munoz
    et al. (2024):

    .. math::

        \frac{d_{L-S}}{R} = \frac{2}{3\pi}\frac{\sin{\alpha} + (\pi-\alpha)\cos{\alpha}}{\cot{\frac{\alpha}{2}}-
        \sin{\frac{\alpha}{2}}\ln{(\cot{\frac{\alpha}{4}})}}

    with :math:`R` the radius of the observed body, and :math:`\alpha` the solar phase angle.

    The SystemOfBodies object must have an apriori ephemeris loaded for the observer, observed body and the sun.
    The corrections from this function should be added to observations (or the apply_* functions should be used).

    Parameters
    ----------
    observations: np.ndarray
        Array of angular observations with columns [time, RA, DEC], with time in s since J2000, and angles in rad
    diameter : float
        Diameter of the body in meter
    bodies : SystemOfBodies
        The SystemOfBodies object that contains ephemerides of observer, body and sun in same reference frame
    body_name : str
        Name of the body being observed, in the SystemOfBodies object
    observer_body_name : str
        Name of the observing (base) body in the SystemOfBodies object
    observer_reference_name : str
        Name of reference point on the observing body. If not given, it is assumed the observer coincides with the origin
            of the observer_body_name body.
    Returns
    -------
    np.ndarray
        Nx2 array of corrections for RA and DEC, to be added to observations
    """
    # Validate inputs
    body_undefined = [not bodies.does_body_exist(name) or bodies.get(name).ephemeris is None
                for name in (body_name, observer_body_name, 'Sun')]
    if any(body_undefined):
        raise ValueError(f'Bodies {body_name}, {observer_body_name}, "Sun" not found in SystemOfBodies object, or their'
                         f'associated ephemerides do not exist')

    if observations.shape[1] != 3:
        raise ValueError(f'Observations must be in shape N x 3 with columns time, RA, DEC')

    if diameter <= 0:
        raise ValueError('Asteroid diameter must be a positive number')

    # Create ephemeris for the reference point if it is given
    if observer_reference_name is not None:
        observer_ephemeris = create_ground_station_ephemeris(
            bodies.get(observer_body_name),
            observer_reference_name,
            bodies
        ) # -> In global frame origin/orientation
    else:
        observer_ephemeris = None

    # Corrections to RA/DEC observations
    corrections = []

    for epoch, ra, dec in observations:

        # Retrieve body, Sun and observer positions in common frame
        body_position = bodies.get(body_name).state_in_base_frame_from_ephemeris(epoch)[:3]
        sun_position = bodies.get('Sun').state_in_base_frame_from_ephemeris(epoch)[:3]

        if observer_ephemeris is not None: # Retrieve from reference point's ephemeris
            observer_position = observer_ephemeris.cartesian_position(epoch)
        else: # Retrieve from observer body ephemeris
            observer_position = bodies.get(observer_body_name).state_in_base_frame_from_ephemeris(epoch)[:3]

        offset_vector = _photocenter_offset(diameter,
                                            (body_position - observer_position),
                                            (body_position - sun_position))

        # Calculate corrections to the astrometry
        corrections.append(
            _offset_vector_to_corrections(offset_vector, ra, dec)
        )

    corrections = np.array(corrections)

    return corrections

def apply_photocenter_correction_to_observation_collection(
        observation_collection: ObservationCollection,
        body_size: float | Iterable[float],
        bodies: SystemOfBodies,
        body_name: str,
        observer_body_name: str,
        observer_reference_name: str | None = None,
        in_place: bool = True,
) -> ObservationCollection | None:
    """
    Computes photocenter corrections and applies them to an observation collection.

    If the body size is given as a float (a diameter in m), the spherical approximation model for the correction,
    :func:`~tudatpy.estimation.observations.observation_corrections.photocenter_correction.photocenter_correction_angular_observations_sphere`, is used.

    If the body size is given as a list of 3 floats (ellipsoid semi-axes) the ellipsoidal approximation for the
    correction is called:
    :func:`~tudatpy.estimation.observations.observation_corrections.photocenter_correction.photocenter_correction_angular_observations_ellipsoid`.

    Parameters
    ----------
    observation_collection : ObservationCollection
        Uncorrected observation collection
    body_size : float | list[float]
        Body size in m. If a float, it is assumed the diameter. If a list of floats, they are assumed to be ellipsoid
        semi-axes.
    bodies : SystemOfBodies
        The SystemOfBodies object that contains ephemerides of observer, body and sun in same reference frame
    body_name : str
        Name of the body being observed, in the SystemOfBodies object
    observer_body_name : str
        Name of the observing (base) body in the SystemOfBodies object
    observer_reference_name : str
        Name of reference point on the observing body. If not given, it is assumed the observer coincides with the origin
            of the observer_body_name body.
    in_place : bool
        Flag if corrections should be applied in-place, or if a new ObservationCollection should be returned. By default,
        True.
    Returns
    -------
    None | ObservationCollection
        Returns None, or the new ObservationCollection depending on the argument 'in_place'
    """
    if isinstance(body_size, Iterable): # Ellipsoidal correction
        return _apply_corrections_to_observation_collection(
            observation_collection=observation_collection,
            bodies=bodies,
            body_name=body_name,
            observer_body_name=observer_body_name,
            observer_reference_name=observer_reference_name,
            in_place=in_place,
            correction_function=photocenter_correction_angular_observations_ellipsoid,
            semi_axes=body_size
        )
    elif isinstance(body_size, (float, int)): # Spherical approximation
        return _apply_corrections_to_observation_collection(
            observation_collection=observation_collection,
            diameter=body_size,
            bodies=bodies,
            body_name=body_name,
            observer_body_name=observer_body_name,
            observer_reference_name=observer_reference_name,
            in_place=in_place,
            correction_function=photocenter_correction_angular_observations_sphere
        )
    else:
        raise TypeError('Body size input type not recognized. Needs to be either float or array of floats')

def _photocenter_correction_ellipsoidal(
        semi_axes: list | np.ndarray,
        e_sun: np.ndarray,
        e_observer : np.ndarray,
) -> np.ndarray:
    """Compute the photocenter-barycenter offset at a single epoch according to Muinonen & Lumme (2015)

    Parameters
    ----------
    semi_axes : list | np.ndarray
        axial parameters of the ellipsoid
    e_sun : np.ndarray
        Unit vector pointing in the direction of the sun, in the principle axis system of the ellipsoid
    e_observer : np.ndarray
        Unit vector pointing in the direction of the observer, in the principle axis system of the ellipsoid

    Returns
    -------
    np.ndarray
        Photocenter-barycenter offset vector in principal axis system of the ellipsoid.
    """
    a, b, c = semi_axes
    c_mat = np.diag((a ** -2, b ** -2, c ** -2))

    # Auxiliary scalar units
    s_sun = np.sqrt(
        e_sun.T @ c_mat @ e_sun
    )
    s_observer = np.sqrt(
        e_observer.T @ c_mat @ e_observer
    )

    cos_alpha_p = e_sun.T @ c_mat @ e_observer / (s_sun * s_observer)
    sin_alpha_p = np.sqrt(1 - cos_alpha_p ** 2)

    s = np.sqrt(s_sun ** 2 + s_observer ** 2 + 2 * s_sun * s_observer * cos_alpha_p)

    cos_lambda_p = (s_sun + s_observer * cos_alpha_p) / s
    sin_lambda_p = s_observer * sin_alpha_p / s

    alpha_p = np.arctan2(sin_alpha_p, cos_alpha_p)
    lambda_p = np.arctan2(sin_lambda_p, cos_lambda_p)

    # Compute dyadic components
    i1 = (
        1/2 * (np.pi - alpha_p) * (np.cos(lambda_p - alpha_p) + np.sin(2*lambda_p) * np.sin(lambda_p - alpha_p)) +
        1/2 * np.cos(lambda_p) * np.sin(alpha_p) - np.log(- np.sin(lambda_p - alpha_p) / np.sin(lambda_p)) *
        np.sin(lambda_p) ** 2 * np.sin(lambda_p - alpha_p)
    )
    i2 = (
        1/2 * (np.pi - alpha_p) * - np.cos(2 * lambda_p) * np.sin(lambda_p - alpha_p) +
        1/2 * np.sin(lambda_p) * np.sin(alpha_p) + np.log(- np.sin(lambda_p - alpha_p) / np.sin(lambda_p)) *
        np.cos(lambda_p) * np.sin(lambda_p) * np.sin(lambda_p - alpha_p)
    )
    i_vec = np.array([i1, i2, 0])

    # Compute rotation matrix K-frame (principal axis) -> K''-frame (computation frame)
    # (note we skip the explicit Euler matrix derivation steps as the angles are not needed)
    n_sun = np.sqrt(c_mat) @ e_sun / (np.sqrt(e_sun.T @ c_mat @ e_sun))
    n_observer = np.sqrt(c_mat) @ e_observer / (np.sqrt(e_observer.T @ c_mat @ e_observer))

    e_x_pp = n_sun
    e_z_pp = np.cross(n_sun, n_observer) / np.linalg.norm(np.cross(n_sun, n_observer))
    e_y_pp = np.cross(e_z_pp, e_x_pp)

    euler_mat = np.vstack((e_x_pp, e_y_pp, e_z_pp))

    cot = lambda x: np.cos(x) / np.sin(x)
    fac = (
        np.cos(lambda_p - alpha_p) + np.cos(lambda_p) + np.sin(lambda_p) * np.sin(lambda_p - alpha_p) *
        np.log(cot(lambda_p / 2) * cot((alpha_p - lambda_p) / 2))
    )
    # Offset in ellipsoid axes frame
    offset = 8 / (3*np.pi * fac) * np.sqrt(np.linalg.inv(c_mat)) @ euler_mat.T @ i_vec

    return offset


def photocenter_correction_angular_observations_ellipsoid(
        observations: np.ndarray ,
        semi_axes: list | np.ndarray,
        bodies: SystemOfBodies,
        body_name: str,
        observer_body_name: str,
        observer_reference_name: str | None = None
) -> np.ndarray:
    r"""
    Compute photocenter-barycenter offset correction assuming an ellipsoidal shape of the observed body, according to
    Muinonen & Lumme (2015).

    Compute photocenter-barycenter offset correction assuming an ellipsoidal shape of the observed body, according to
    Muinonen & Lumme (2015):

    .. math::
        \mathbf{d}(\alpha) = \frac{1}{L(\alpha)}\frac{1}{3}F_0 \tilde{\omega} P(\alpha) abc \frac{S_\odot S_\oplus}{S}
        \sqrt{C^{-1}} R_E^T \begin{bmatrix}
        I_1 \\ I_2 \\ 0 \end{bmatrix}

    An apriori ephemeris should be loaded for the observer, observed body and sun. Furthermore, the observed body must
    have a rotational ephemeris loaded, where the body-fixed X, Y and Z-axes must align with the principal axes of the
    ellipsoid in the same order as the semi_axes parameter.

    Parameters
    ----------
    observations : np.ndarray
        Observations with columns [time, RA, DEC]
    semi_axes : list | np.ndarray
        Ellipsoid semi-axes as [a, b, c]. They must be aligned with the body-fixed (principal axis) frame of the
        ellipsoid, i.e.: a along x-axis, b along y-axis, and c along z-axis.
    bodies : SystemOfBodies
        System of Bodies object
    body_name : str
        Name of body being observed
    observer_body_name : str
        Name of the (base) body acting as observer
    observer_reference_name : str | None
        If applicable, the name of the reference point on the observer body

    Returns
    -------
    np.ndarray
        Corrections to RA, DEC
    """
    # Input validation
    body_undefined = [not bodies.does_body_exist(name) or bodies.get(name).ephemeris is None
                for name in (body_name, observer_body_name, 'Sun')]
    if any(body_undefined):
        raise ValueError(f'Bodies {body_name}, {observer_body_name}, "Sun" not found in SystemOfBodies object, or their'
                         f'associated ephemerides do not exist')

    if bodies.get(body_name).rotation_model is None:
        raise ValueError(f'Body {body_name} needs a rotation model for the photocenter correction computation')

    if observations.shape[1] != 3:
        raise ValueError('Observations must be shaped N x 3')

    if len(semi_axes) != 3:
        raise ValueError('Semi-axes array must be of length 3')

    if any([ax <= 0 for ax in semi_axes]):
        raise ValueError('Semi-axes must be positive')

    # Create ephemeris for the reference point if it is given
    if observer_reference_name is not None:
        observer_ephemeris = create_ground_station_ephemeris(
            bodies.get(observer_body_name),
            observer_reference_name,
            bodies
        )  # -> In global frame origin/orientation
    else:
        observer_ephemeris = None

    corrections = []

    for epoch, ra, dec in observations:

        # Get inertial positions of sun, body and observer
        if observer_ephemeris is not None:
            observer_inertial = observer_ephemeris.cartesian_position(epoch)
        else:
            observer_inertial = bodies.get(observer_body_name).state_in_base_frame_from_ephemeris(epoch)[:3]
        sun_inertial = bodies.get('Sun').state_in_base_frame_from_ephemeris(epoch)[:3]
        body_inertial = bodies.get(body_name).state_in_base_frame_from_ephemeris(epoch)[:3]
        observer_body_distance = np.linalg.norm(observer_inertial - body_inertial)

        # Inertial directions of sun and observer from observed body
        e_sun_inertial = _unit(sun_inertial - body_inertial)
        e_observer_inertial = _unit(observer_inertial - body_inertial)

        # Retrieve rotational state of ellipsoid body at time of emission
        light_time = observer_body_distance / SPEED_OF_LIGHT
        rotation_matrix = bodies.get(body_name).rotation_model.inertial_to_body_fixed_rotation(epoch - light_time)

        e_sun = rotation_matrix @ e_sun_inertial
        e_observer = rotation_matrix @ e_observer_inertial

        # Get offset in axes of ellipsoid
        offset = _photocenter_correction_ellipsoidal(
            semi_axes,
            e_sun,
            e_observer,
        )

        # Rotate offset to inertial frame and retain only its plane-of-sky component
        offset_inertial = rotation_matrix.T @ offset
        offset_sky = offset_inertial - np.dot(offset_inertial, e_observer_inertial) * e_observer_inertial

        # Translate into corrections
        corrections.append(
            _offset_vector_to_corrections(offset_sky / observer_body_distance, ra, dec)
        )

    return np.array(corrections)









