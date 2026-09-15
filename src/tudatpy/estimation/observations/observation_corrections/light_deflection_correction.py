"""
Functions to calculate light deflection corrections to observations
"""

import numpy as np
from numpy.linalg import norm
from tudatpy.estimation.observations import ObservationCollection
from tudatpy.dynamics.environment import SystemOfBodies
from tudatpy.dynamics.environment_setup import create_ground_station_ephemeris
from ._correction_utils import (
    _offset_vector_to_corrections,
    _apply_corrections_to_observation_collection,
)
from tudatpy.constants import SPEED_OF_LIGHT
from collections.abc import Iterable


def _light_deflection_single_contribution(
    observer_position: np.ndarray,
    observed_body_position: np.ndarray,
    perturbing_body_position: np.ndarray,
    perturbing_body_gravitational_parameter: float,
) -> np.ndarray:
    """
    Helper function to calculate post-Newtonian light-bending contribution from one body (delta k_pn in Klioner (2003))

    Parameters
    ----------
    observer_position : np.ndarray
        Position vector of observer at time of observation
    observed_body_position : np.ndarray
        Position vector of body being observed at time of emission
    perturbing_body_position : np.ndarray
        Position vector of perturber at reference time
    perturbing_body_gravitational_parameter : float
        Gravitational parameter of perturber

    Returns
    -------
    np.ndarray
        Individual contribution to total light deflection
    """
    # Relative-position vectors from Klioner (2003), Eq. 70: R, r_eA, and r_oA.
    position_observer_with_respect_to_observed_body = (
        observer_position - observed_body_position
    )  # Represents R in Klioner (2003), Eq. 70.
    position_observed_body_with_respect_to_perturbing_body = (
        observed_body_position - perturbing_body_position
    )  # Represents r_eA in Klioner (2003), Eq. 70.
    position_observer_with_respect_to_perturbing_body = (
        observer_position - perturbing_body_position
    )  # Represents r_oA in Klioner (2003), Eq. 70.

    # Calculate contribution of deflection from perturber body
    relativistic_strength_factor = (2 * perturbing_body_gravitational_parameter) / SPEED_OF_LIGHT**2
    cross_product_numerator = np.cross(
        position_observer_with_respect_to_observed_body,
        np.cross(
            position_observed_body_with_respect_to_perturbing_body,
            position_observer_with_respect_to_perturbing_body,
        ),
    )
    normalization_denominator = (
        norm(position_observer_with_respect_to_observed_body)
        * norm(position_observer_with_respect_to_perturbing_body)
        * (
            norm(position_observed_body_with_respect_to_perturbing_body)
            * norm(position_observer_with_respect_to_perturbing_body)
            + np.dot(
                position_observer_with_respect_to_perturbing_body,
                position_observed_body_with_respect_to_perturbing_body,
            )
        )
    )
    light_deflection_offset = (
        relativistic_strength_factor * cross_product_numerator / normalization_denominator
    )

    return light_deflection_offset


def light_deflection_correction_angular_observations(
    observations: np.ndarray,
    bodies: SystemOfBodies,
    body_name: str,
    observer_body_name: str,
    observer_reference_name: str | None = None,
    perturbing_bodies_list: Iterable[str] = ("Sun",),
) -> np.ndarray:
    r"""
    Compute corrections to angular observations for the relativistic deflection of light around massive bodies.

    Compute corrections to angular observations for the relativistic deflection of light around massive bodies, according to
    Klioner (2003) equation 70 (post-Newtonian component only):

    .. math::

        \delta \mathbf{k}_{pN} = - \sum_A \left[\frac{(1+\gamma)GM_A}{c^2} \frac{\mathbf{R} \times
        (\mathbf{r}_{eA} \times \mathbf{r}_oA)}{|\mathbf{R}||\mathbf{r}_{oA}|(|\mathbf{r}_{eA}||\mathbf{r}_{oA}|+
        \mathbf{r}_{oA}\cdot \mathbf{r}_eA)}\right]

    where :math:`\mathbf{R}` is the position of the observer w.r.t the observed body, :math:`\mathbf{r}_{eA}` is the
    position of the observed body w.r.t the perturber, and :math:`\mathbf{r}_{oA}` is the position of
    the observer w.r.t. the perturber.

    A reference ephemeris of the observed body and all the perturbing bodies should be
    loaded in the SystemOfBodies object. The corrections from this function should be added to observations (or the
    apply_* functions should be used).

    Parameters
    ----------
    observations : np.ndarray
        Observations with columns [time, RA, DEC] in seconds since J2000 and radians, respectively
    bodies : SystemOfBodies
        The SystemOfBodies object
    body_name : str
        Name of the body being observed, in the SystemOfBodies object
    observer_body_name : str
        Name of the observing (base) body in the SystemOfBodies object
    observer_reference_name : str | None
        Name of reference point on the observing body. If not given, it is assumed the observer coincides with the origin
        of the observer_body_name body.
    perturbing_bodies_list : list[str]
        Names of the bodies that light-deflection contribution should be computed for, default = 'Sun'

    Returns
    -------
    np.ndarray
        Nx2 array of corrections to RA and DEC in radians. Must be added to observations
    """
    # Input validation
    if observations.shape[1] != 3:
        raise ValueError(f"Observations must be in shape N x 3 with columns time, RA, DEC")
    body_or_ephemeris_undefined = [
        not bodies.does_body_exist(name) or bodies.get(name).ephemeris is None
        for name in list(perturbing_bodies_list) + [body_name, observer_body_name]
    ]
    if any(body_or_ephemeris_undefined):
        raise ValueError(
            "Some or all included bodies in the relativistic light deflection computation are missing"
            "from SystemOfBodies or their associated ephemerides are not specified."
        )

    # Get observer reference point ephemeris if it exists
    if observer_reference_name is not None:
        observer_ephemeris = create_ground_station_ephemeris(
            bodies.get(observer_body_name), observer_reference_name, bodies
        )  # -> in global frame origin/orientation
    else:
        observer_ephemeris = None

    angular_corrections = []

    # Loop over observations
    for observation_epoch, right_ascension, declination in observations:
        # Position of observer at current epoch:
        if observer_ephemeris is not None:
            observer_position = observer_ephemeris.cartesian_position(observation_epoch)
        else:
            observer_position = bodies.get(observer_body_name).state_in_base_frame_from_ephemeris(
                observation_epoch
            )[:3]

        observed_body_position = bodies.get(body_name).state_in_base_frame_from_ephemeris(
            observation_epoch
        )[:3]

        # One iteration of light-time calculation for light from observed body to observer
        observed_body_to_observer_light_time = (
            np.linalg.norm(observer_position - observed_body_position) / SPEED_OF_LIGHT
        )

        # Observed body position at time of emission of light:
        observed_body_position = bodies.get(body_name).state_in_base_frame_from_ephemeris(
            observation_epoch - observed_body_to_observer_light_time
        )[:3]

        total_light_deflection_offset = np.zeros(3)

        # Loop over light-deflecting bodies
        for perturbing_body_name in perturbing_bodies_list:

            perturbing_body_gravitational_parameter = bodies.get(
                perturbing_body_name
            ).gravitational_parameter
            perturbing_body_position = bodies.get(
                perturbing_body_name
            ).state_in_base_frame_from_ephemeris(observation_epoch)[:3]

            # Compute one iteration of light-time from observer to perturbing body
            perturbing_body_to_observer_light_time = (
                np.linalg.norm(observer_position - perturbing_body_position) / SPEED_OF_LIGHT
            )

            # Perturber position at reference time
            perturbing_body_position = bodies.get(
                perturbing_body_name
            ).state_in_base_frame_from_ephemeris(
                observation_epoch - perturbing_body_to_observer_light_time
            )[
                :3
            ]

            perturbing_body_light_deflection_offset = _light_deflection_single_contribution(
                observer_position=observer_position,
                observed_body_position=observed_body_position,
                perturbing_body_position=perturbing_body_position,
                perturbing_body_gravitational_parameter=perturbing_body_gravitational_parameter,
            )
            # Omit the minus sign in the original equation because of the correction convention.
            total_light_deflection_offset += perturbing_body_light_deflection_offset

        angular_corrections.append(
            _offset_vector_to_corrections(
                total_light_deflection_offset,
                right_ascension,
                declination,
            )
        )

    return np.array(angular_corrections)


def apply_light_deflection_correction_to_observation_collection(
    observation_collection: ObservationCollection,
    bodies: SystemOfBodies,
    body_name: str,
    observer_body_name: str,
    observer_reference_name: str | None = None,
    perturbing_bodies_list: Iterable[str] = ("Sun",),
    in_place: bool = True,
) -> ObservationCollection | None:
    """
    Computes relativistic light-deflection corrections and applies them to an observation collection.

    Computes relativistic light-deflection corrections and applies them to an observation collection. Calls the function
    :func:`~tudatpy.estimation.observations.observation_corrections.light_deflection_correction.light_deflection_correction_angular_observations`,
    and applies the resulting corrections to all angular observations in the :class:`~tudatpy.estimation.observations.ObservationCollection` with the specified
    link-ends.

    Parameters
    ----------
    observation_collection : :class:`~tudatpy.estimation.observations.ObservationCollection`
        ObservationCollection containing angular observations on specified link-ends
    bodies : SystemOfBodies
        The SystemOfBodies object
    body_name : str
        Name of the body being observed, in the SystemOfBodies object
    observer_body_name : str
        Name of the observing (base) body in the SystemOfBodies object
    observer_reference_name : str | None
        Name of reference point on the observing body. If not given, it is assumed the observer coincides with the origin
        of the observer_body_name body.
    perturbing_bodies_list : list[str]
        Names of the bodies that light-deflection contribution should be computed for, default = 'Sun'
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
        correction_function=light_deflection_correction_angular_observations,
        in_place=in_place,
        perturbing_bodies_list=perturbing_bodies_list,
    )
