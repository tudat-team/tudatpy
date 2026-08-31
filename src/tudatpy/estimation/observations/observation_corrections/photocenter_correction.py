"""
Functions to calculate photocenter corrections to observations
"""
import numpy as np
from tudatpy.estimation.observations import ObservationCollection
from tudatpy.dynamics.environment_setup import create_ground_station_ephemeris
from tudatpy.dynamics.environment import SystemOfBodies
from ._correction_utils import _offset_vector_to_corrections, _unit, _apply_corrections_to_observation_collection
from tudatpy.constants import SPEED_OF_LIGHT
from collections.abc import Sequence

def apply_photocenter_correction_to_observation_collection(
        observation_collection: ObservationCollection,
        body_dimensions: float | Sequence[float],
        bodies: SystemOfBodies,
        body_name: str,
        observer_body_name: str,
        observer_reference_name: str | None = None,
        in_place: bool = True,
) -> ObservationCollection | None:
    """
    Computes photocenter corrections and applies them to an observation collection.

    Computes photocenter corrections and applies them to an observation collection. Calls the function
    :func:`~tudatpy.estimation.observations.observation_corrections.photocenter_correction.photocenter_correction_angular_observations`,
    and applies the resulting corrections to all angular observations in the :class:`~tudatpy.estimation.observations.ObservationCollection` with the specified
    link-ends.

    Parameters
    ----------
    observation_collection : :class:`~tudatpy.estimation.observations.ObservationCollection`
        ObservationCollection containing angular observations on specified link-ends
    body_dimensions : float | list[float]
        Body size in meters. If a scalar, the body is modeled as a sphere with that radius. If a list of floats, they are
        assumed to be ellipsoid semi-axes.
    bodies : SystemOfBodies
        The SystemOfBodies object
    body_name : str
        Name of the body being observed, in the SystemOfBodies object
    observer_body_name : str
        Name of the observing (base) body in the SystemOfBodies object
    observer_reference_name : str | None
        Name of reference point on the observing body. If not given, it is assumed the observer coincides with the origin
        of the observer_body_name body.
    in_place : bool
        If true, corrections are applied in-place to the ObservationCollection object. If false, a new ObservationCollection
        is returned with the corrections applied. By default, true.

    Returns
    -------
    None | ObservationCollection
        Returns a new observation collection with applied corrections if in_place is False.
    """
    return _apply_corrections_to_observation_collection(
        observation_collection=observation_collection,
        body_name=body_name,
        bodies=bodies,
        observer_body_name=observer_body_name,
        observer_reference_name=observer_reference_name,
        in_place=in_place,
        correction_function=photocenter_correction_angular_observations,
        body_dimensions=body_dimensions,
    )


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


def photocenter_correction_angular_observations(
        observations: np.ndarray ,
        body_dimensions: float | Sequence[float],
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

    When the observed body is a sphere, this collapses into:

    .. math::

        \frac{d_{L-S}}{R} = \frac{2}{3\pi}\frac{\sin{\alpha} + (\pi-\alpha)\cos{\alpha}}{\cot{\frac{\alpha}{2}}-
        \sin{\frac{\alpha}{2}}\ln{(\cot{\frac{\alpha}{4}})}}

    where :math:`R` is the observed body radius, and :math:`\alpha` is the solar phase angle, i.e., the angle spanned
    by the relative position vector to the sun and observer. For the full derivation, the reader is referred to the
    original reference paper.

    An apriori ephemeris should be loaded for the observer, observed body and sun. Furthermore, if the observed body
    is an ellipsoid, it must have a rotational ephemeris loaded, where the body-fixed X, Y and Z-axes must align with
    the principal axes of the ellipsoid in the same order as the body_dimensions parameter. The corrections returned
    by this function should be added to observations (or the apply_* function should be used).

    Parameters
    ----------
    observations : np.ndarray
        Observations with columns [time, RA, DEC] in seconds since J2000 and radians, respectively
    body_dimensions : float | list[float]
        Body size in meters. If a scalar, the body is modeled as a sphere with that radius. If a list of floats, they are
        assumed to be ellipsoid semi-axes.
    bodies : SystemOfBodies
        The SystemOfBodies object
    body_name : str
        Name of the body being observed, in the SystemOfBodies object
    observer_body_name : str
        Name of the observing (base) body in the SystemOfBodies object
    observer_reference_name : str | None
        Name of reference point on the observing body. If not given, it is assumed the observer coincides with the origin
        of the observer_body_name body.

    Returns
    -------
    np.ndarray
        Nx2 array of corrections to RA and DEC in radians. Must be added to observations
    """
    # Check if all bodies are defined in the SystemOfBodies
    body_or_ephemeris_undefined = [not bodies.does_body_exist(name) or bodies.get(name).ephemeris is None
                for name in (body_name, observer_body_name, 'Sun')]
    if any(body_or_ephemeris_undefined):
        raise ValueError(f'Bodies {body_name}, {observer_body_name}, "Sun" not found in SystemOfBodies object, or their'
                         f'associated ephemerides do not exist')

    # Check body shape parameters
    if np.isscalar(body_dimensions): # -> body is a sphere
        if body_dimensions <= 0:
            raise ValueError('Observed body radius must be positive')
        body_dimensions = 3 * [body_dimensions]
        rotation_inertial_to_body_frame = lambda epoch: np.eye(3) # Rotational state irrelevant for a sphere

    else:
        if len(body_dimensions) != 3:
            raise ValueError('body_dimensions must be a float or iterable of length 3')
        if any([axis <= 0 for axis in body_dimensions]):
            raise ValueError('Semi-axes must be positive')
        if bodies.get(body_name).rotation_model is None:
            raise ValueError(f'Body {body_name} needs a rotation model for the photocenter correction')
        rotation_inertial_to_body_frame = bodies.get(body_name).rotation_model.inertial_to_body_fixed_rotation

    if observations.shape[1] != 3:
        raise ValueError('Observations must be shaped N x 3')

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
        rotation_matrix = rotation_inertial_to_body_frame(epoch - light_time)

        e_sun = rotation_matrix @ e_sun_inertial
        e_observer = rotation_matrix @ e_observer_inertial

        # Get offset in axes of ellipsoid
        offset = _photocenter_correction_ellipsoidal(
            body_dimensions,
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









