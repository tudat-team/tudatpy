"""
Functions to calculate photocenter corrections to observations
"""

import numpy as np
from tudatpy.estimation.observations import ObservationCollection
from tudatpy.dynamics.environment_setup import create_ground_station_ephemeris
from tudatpy.dynamics.environment import SystemOfBodies
from ._correction_utils import (
    _offset_vector_to_corrections,
    _unit,
    _apply_corrections_to_observation_collection,
)
from tudatpy.constants import SPEED_OF_LIGHT
from collections.abc import Iterable


def apply_photocenter_correction_to_observation_collection(
    observation_collection: ObservationCollection,
    body_dimensions: float | Iterable[float],
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
        Body size in meters. If a scalar, the body is modeled as a sphere with that radius. If a list of 3 floats, they are
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
    unit_vector_to_sun_body_fixed: np.ndarray,
    unit_vector_to_observer_body_fixed: np.ndarray,
) -> np.ndarray:
    """Compute the photocenter-barycenter offset at a single epoch according to Muinonen & Lumme (2015)

    Parameters
    ----------
    semi_axes : list | np.ndarray
        axial parameters of the ellipsoid
    unit_vector_to_sun_body_fixed : np.ndarray
        Unit vector pointing in the direction of the Sun, in the principal-axis system of the ellipsoid
    unit_vector_to_observer_body_fixed : np.ndarray
        Unit vector pointing in the direction of the observer, in the principal-axis system of the ellipsoid

    Returns
    -------
    np.ndarray
        Photocenter-barycenter offset vector in principal axis system of the ellipsoid.
    """
    # These variables represent the semi-axis symbols a, b, and c in Muinonen & Lumme (2015).
    semi_axis_a, semi_axis_b, semi_axis_c = semi_axes
    inverse_squared_semi_axes_matrix = np.diag(
        (semi_axis_a**-2, semi_axis_b**-2, semi_axis_c**-2)
    )  # Represents C in Muinonen & Lumme (2015).

    # This variable represents S_sun in Muinonen & Lumme (2015).
    scaled_sun_direction_norm = np.sqrt(
        unit_vector_to_sun_body_fixed.T
        @ inverse_squared_semi_axes_matrix
        @ unit_vector_to_sun_body_fixed
    )
    # This variable represents S_observer in Muinonen & Lumme (2015).
    scaled_observer_direction_norm = np.sqrt(
        unit_vector_to_observer_body_fixed.T
        @ inverse_squared_semi_axes_matrix
        @ unit_vector_to_observer_body_fixed
    )

    cosine_scaled_phase_angle = (
        unit_vector_to_sun_body_fixed.T
        @ inverse_squared_semi_axes_matrix
        @ unit_vector_to_observer_body_fixed
        / (scaled_sun_direction_norm * scaled_observer_direction_norm)
    )  # Represents cos(alpha') in Muinonen & Lumme (2015).
    sine_scaled_phase_angle = np.sqrt(
        1 - cosine_scaled_phase_angle**2
    )  # Represents sin(alpha') in Muinonen & Lumme (2015).

    combined_scaled_direction_norm = np.sqrt(
        scaled_sun_direction_norm**2
        + scaled_observer_direction_norm**2
        + 2 * scaled_sun_direction_norm * scaled_observer_direction_norm * cosine_scaled_phase_angle
    )  # Represents S in Muinonen & Lumme (2015).

    cosine_scaled_auxiliary_angle = (
        scaled_sun_direction_norm + scaled_observer_direction_norm * cosine_scaled_phase_angle
    ) / combined_scaled_direction_norm  # Represents cos(lambda') in Muinonen & Lumme (2015).
    sine_scaled_auxiliary_angle = (
        scaled_observer_direction_norm * sine_scaled_phase_angle
    ) / combined_scaled_direction_norm  # Represents sin(lambda') in Muinonen & Lumme (2015).

    scaled_phase_angle = np.arctan2(
        sine_scaled_phase_angle, cosine_scaled_phase_angle
    )  # Represents alpha' in Muinonen & Lumme (2015).
    scaled_auxiliary_angle = np.arctan2(
        sine_scaled_auxiliary_angle, cosine_scaled_auxiliary_angle
    )  # Represents lambda' in Muinonen & Lumme (2015).

    # This variable represents I_1 in Muinonen & Lumme (2015).
    first_dyadic_component = (
        1
        / 2
        * (np.pi - scaled_phase_angle)
        * (
            np.cos(scaled_auxiliary_angle - scaled_phase_angle)
            + np.sin(2 * scaled_auxiliary_angle)
            * np.sin(scaled_auxiliary_angle - scaled_phase_angle)
        )
        + 1 / 2 * np.cos(scaled_auxiliary_angle) * np.sin(scaled_phase_angle)
        - np.log(
            -np.sin(scaled_auxiliary_angle - scaled_phase_angle) / np.sin(scaled_auxiliary_angle)
        )
        * np.sin(scaled_auxiliary_angle) ** 2
        * np.sin(scaled_auxiliary_angle - scaled_phase_angle)
    )
    # This variable represents I_2 in Muinonen & Lumme (2015).
    second_dyadic_component = (
        1
        / 2
        * (np.pi - scaled_phase_angle)
        * -np.cos(2 * scaled_auxiliary_angle)
        * np.sin(scaled_auxiliary_angle - scaled_phase_angle)
        + 1 / 2 * np.sin(scaled_auxiliary_angle) * np.sin(scaled_phase_angle)
        + np.log(
            -np.sin(scaled_auxiliary_angle - scaled_phase_angle) / np.sin(scaled_auxiliary_angle)
        )
        * np.cos(scaled_auxiliary_angle)
        * np.sin(scaled_auxiliary_angle)
        * np.sin(scaled_auxiliary_angle - scaled_phase_angle)
    )
    # This variable represents the vector [I_1, I_2, 0]^T in Muinonen & Lumme (2015).
    dyadic_component_vector = np.array([first_dyadic_component, second_dyadic_component, 0])

    # Compute rotation matrix K-frame (principal axis) -> K''-frame (computation frame)
    # (note we skip the explicit Euler matrix derivation steps as the angles are not needed)
    normalized_scaled_sun_direction = (
        np.sqrt(inverse_squared_semi_axes_matrix)
        @ unit_vector_to_sun_body_fixed
        / scaled_sun_direction_norm
    )  # Represents n_sun in Muinonen & Lumme (2015).
    normalized_scaled_observer_direction = (
        np.sqrt(inverse_squared_semi_axes_matrix)
        @ unit_vector_to_observer_body_fixed
        / scaled_observer_direction_norm
    )  # Represents n_observer in Muinonen & Lumme (2015).

    computation_frame_x_axis = (
        normalized_scaled_sun_direction  # Represents e_x'' in Muinonen & Lumme (2015).
    )
    computation_frame_z_axis = np.cross(
        normalized_scaled_sun_direction, normalized_scaled_observer_direction
    )
    computation_frame_z_axis /= np.linalg.norm(
        computation_frame_z_axis
    )  # Represents e_z'' in Muinonen & Lumme (2015).
    computation_frame_y_axis = np.cross(
        computation_frame_z_axis, computation_frame_x_axis
    )  # Represents e_y'' in Muinonen & Lumme (2015).

    principal_axes_to_computation_frame_rotation = np.vstack(
        (
            computation_frame_x_axis,
            computation_frame_y_axis,
            computation_frame_z_axis,
        )
    )  # Represents R_E in Muinonen & Lumme (2015).

    cotangent = lambda angle: np.cos(angle) / np.sin(angle)
    brightness_phase_factor = (
        np.cos(scaled_auxiliary_angle - scaled_phase_angle)
        + np.cos(scaled_auxiliary_angle)
        + np.sin(scaled_auxiliary_angle)
        * np.sin(scaled_auxiliary_angle - scaled_phase_angle)
        * np.log(
            cotangent(scaled_auxiliary_angle / 2)
            * cotangent((scaled_phase_angle - scaled_auxiliary_angle) / 2)
        )
    )
    # This variable represents d(alpha), expressed in the ellipsoid principal-axis frame.
    photocenter_offset_body_fixed = (
        8
        / (3 * np.pi * brightness_phase_factor)
        * np.sqrt(np.linalg.inv(inverse_squared_semi_axes_matrix))
        @ principal_axes_to_computation_frame_rotation.T
        @ dyadic_component_vector
    )

    return photocenter_offset_body_fixed


def photocenter_correction_angular_observations(
    observations: np.ndarray,
    body_dimensions: float | Iterable[float],
    bodies: SystemOfBodies,
    body_name: str,
    observer_body_name: str,
    observer_reference_name: str | None = None,
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
        Body size in meters. If a scalar, the body is modeled as a sphere with that radius. If a list of 3 floats, they are
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
    body_or_ephemeris_undefined = [
        not bodies.does_body_exist(name) or bodies.get(name).ephemeris is None
        for name in (body_name, observer_body_name, "Sun")
    ]
    if any(body_or_ephemeris_undefined):
        raise ValueError(
            f'Bodies {body_name}, {observer_body_name}, "Sun" not found in SystemOfBodies object, or their'
            f"associated ephemerides do not exist"
        )

    # Check body shape parameters
    if np.isscalar(body_dimensions):  # -> body is a sphere
        if body_dimensions <= 0:
            raise ValueError("Observed body radius must be positive")
        body_dimensions = 3 * [body_dimensions]
        inertial_to_body_fixed_rotation = lambda epoch: np.eye(
            3
        )  # Rotational state irrelevant for a sphere

    else:
        if len(body_dimensions) != 3:
            raise ValueError("body_dimensions must be a float or iterable of length 3")
        if any([axis <= 0 for axis in body_dimensions]):
            raise ValueError("Semi-axes must be positive")
        if bodies.get(body_name).rotation_model is None:
            raise ValueError(
                f"Body {body_name} needs a rotation model for the photocenter correction"
            )
        inertial_to_body_fixed_rotation = bodies.get(
            body_name
        ).rotation_model.inertial_to_body_fixed_rotation

    if observations.shape[1] != 3:
        raise ValueError("Observations must be shaped N x 3")

    # Create ephemeris for the reference point if it is given
    if observer_reference_name is not None:
        observer_ephemeris = create_ground_station_ephemeris(
            bodies.get(observer_body_name), observer_reference_name, bodies
        )  # -> In global frame origin/orientation
    else:
        observer_ephemeris = None

    angular_corrections = []

    for observation_epoch, right_ascension, declination in observations:

        # Get inertial positions of sun, body and observer
        if observer_ephemeris is not None:
            observer_inertial_position = observer_ephemeris.cartesian_position(observation_epoch)
        else:
            observer_inertial_position = bodies.get(
                observer_body_name
            ).state_in_base_frame_from_ephemeris(observation_epoch)[:3]
        sun_inertial_position = bodies.get("Sun").state_in_base_frame_from_ephemeris(
            observation_epoch
        )[:3]
        observed_body_inertial_position = bodies.get(body_name).state_in_base_frame_from_ephemeris(
            observation_epoch
        )[:3]
        observer_to_observed_body_distance = np.linalg.norm(
            observer_inertial_position - observed_body_inertial_position
        )

        # Inertial directions of sun and observer from observed body
        unit_vector_to_sun_inertial = _unit(sun_inertial_position - observed_body_inertial_position)
        unit_vector_to_observer_inertial = _unit(
            observer_inertial_position - observed_body_inertial_position
        )

        # Retrieve rotational state of ellipsoid body at time of emission
        observed_body_to_observer_light_time = observer_to_observed_body_distance / SPEED_OF_LIGHT
        inertial_to_body_fixed_rotation_matrix = inertial_to_body_fixed_rotation(
            observation_epoch - observed_body_to_observer_light_time
        )

        unit_vector_to_sun_body_fixed = (
            inertial_to_body_fixed_rotation_matrix @ unit_vector_to_sun_inertial
        )
        unit_vector_to_observer_body_fixed = (
            inertial_to_body_fixed_rotation_matrix @ unit_vector_to_observer_inertial
        )

        # Get offset in axes of ellipsoid
        photocenter_offset_body_fixed = _photocenter_correction_ellipsoidal(
            body_dimensions,
            unit_vector_to_sun_body_fixed,
            unit_vector_to_observer_body_fixed,
        )

        # Rotate offset to inertial frame and retain only its plane-of-sky component
        photocenter_offset_inertial = (
            inertial_to_body_fixed_rotation_matrix.T @ photocenter_offset_body_fixed
        )
        photocenter_offset_plane_of_sky = (
            photocenter_offset_inertial
            - np.dot(photocenter_offset_inertial, unit_vector_to_observer_inertial)
            * unit_vector_to_observer_inertial
        )

        # Translate into corrections
        angular_corrections.append(
            _offset_vector_to_corrections(
                photocenter_offset_plane_of_sky / observer_to_observed_body_distance,
                right_ascension,
                declination,
            )
        )

    return np.array(angular_corrections)
