"""
Functions to calculate light deflection corrections to observations
"""
import numpy as np
from numpy.linalg import norm
from tudatpy.estimation.observations import ObservationCollection
from tudatpy.dynamics.environment import SystemOfBodies
from tudatpy.dynamics.environment_setup import create_ground_station_ephemeris
from ._correction_utils import _offset_vector_to_corrections, _apply_corrections_to_observation_collection
from tudatpy.constants import SPEED_OF_LIGHT
from collections.abc import Iterable

def _light_deflection_single_contribution(
        observer_position: np.ndarray,
        body_position: np.ndarray,
        perturber_position: np.ndarray,
        mu_perturber: float
) -> np.ndarray:
    """
    Helper function to calculate post-Newtonian light-bending contribution from one body (delta k_pn in Klioner (2003))

    Parameters
    ----------
    observer_position : np.ndarray
        Position vector of observer at time of observation
    body_position : np.ndarray
        Position vector of body being observed at time of emission
    perturber_position : np.ndarray
        Position vector of perturber at reference time
    mu_perturber : float
        Gravitational parameter of perturber

    Returns
    -------
    np.ndarray
        Individual contribution to total light deflection
    """
    # Define vectors from Klioner 2003
    position_observer_wrt_body = observer_position - body_position  # R
    position_body_wrt_perturber = body_position - perturber_position # r_ea
    position_observer_wrt_perturber = observer_position - perturber_position # r_oA

    # Calculate contribution of deflection from perturber body
    strength_factor = (2 * mu_perturber) / (SPEED_OF_LIGHT ** 2)
    num = np.cross(position_observer_wrt_body, np.cross(position_body_wrt_perturber, position_observer_wrt_perturber))
    denom = norm(position_observer_wrt_body) * norm(position_observer_wrt_perturber) * (
            norm(position_body_wrt_perturber) * norm(position_observer_wrt_perturber) +
            np.dot(position_observer_wrt_perturber, position_body_wrt_perturber))
    offset_vec = strength_factor * num / denom

    return offset_vec


def light_deflection_correction_angular_observations(
        observations: np.ndarray,
        bodies: SystemOfBodies,
        body_name: str,
        observer_body_name: str,
        observer_reference_name: str | None = None,
        perturbing_bodies_list: Iterable[str] = ('Sun',),
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
        Array of angular observations with columns [time, RA, DEC], with time in s since J2000, and angles in rad
    bodies : SystemOfBodies
        SystemOfBodies object
    body_name : str
        Name of the observed body object
    observer_body_name : str
        Name of the observer body
    observer_reference_name : str
        Name of the reference point on the observer body. If not given, it is assumed that the observer location
            coincides with the 'observer_body_name' body center.
    perturbing_bodies_list : Iterable[str]
        Names of the bodies that light-deflection contribution should be computed for, default = 'Sun'

    Returns
    -------
    np.ndarray
        Nx2 array of corrections in RA and DEC. Must be added to observations
    """
    # Input validation
    if observations.shape[1] != 3:
        raise ValueError(f'Observations must be in shape N x 3 with columns time, RA, DEC')
    body_undefined = [not bodies.does_body_exist(name) or bodies.get(name).ephemeris is None
                for name in list(perturbing_bodies_list) + [body_name, observer_body_name]]
    if any(body_undefined):
        raise ValueError('Some or all included bodies in the relativistic light deflection computation are missing'
                         'from SystemOfBodies or their associated ephemerides are not specified.')

    # Get observer reference point ephemeris if it exists
    if observer_reference_name is not None:
        observer_ephemeris = create_ground_station_ephemeris(
            bodies.get(observer_body_name),
            observer_reference_name,
            bodies
        ) # -> in global frame origin/orientation
    else:
        observer_ephemeris = None

    corrections = []

    # Loop over observations
    for epoch, ra, dec in observations:
        # Position of observer at current epoch:
        if observer_ephemeris is not None:
            observer_position = observer_ephemeris.cartesian_position(epoch)
        else:
            observer_position = bodies.get(observer_body_name).state_in_base_frame_from_ephemeris(epoch)[:3]

        body_position = bodies.get(body_name).state_in_base_frame_from_ephemeris(epoch)[:3]

        # One iteration of light-time calculation for light from observed body to observer
        light_time_body_to_observer = np.linalg.norm(
            observer_position - body_position
        ) / SPEED_OF_LIGHT

        # Observed body position at time of emission of light:
        body_position = bodies.get(body_name).state_in_base_frame_from_ephemeris(epoch - light_time_body_to_observer)[:3]

        total_offset = np.zeros(3) # Total deflection accumulated from all bodies

        # Loop over light-deflecting bodies
        for perturber_name in perturbing_bodies_list:

            mu_perturber = bodies.get(perturber_name).gravitational_parameter
            perturber_position = bodies.get(perturber_name).state_in_base_frame_from_ephemeris(epoch)[:3]

            # Compute one iteration of light-time from observer to perturbing body
            light_time_perturber_to_observer = np.linalg.norm(
                observer_position - perturber_position
            ) / SPEED_OF_LIGHT

            # Perturber position at reference time
            perturber_position = bodies.get(perturber_name).state_in_base_frame_from_ephemeris(
                epoch - light_time_perturber_to_observer)[:3]

            offset_from_perturber = _light_deflection_single_contribution( # Single contribution to light deflection
                observer_position= observer_position,
                body_position= body_position,
                perturber_position= perturber_position,
                mu_perturber = mu_perturber,
            )
            total_offset += offset_from_perturber # Note: omit minus sign of original eq because of convention


        corrections.append(
            _offset_vector_to_corrections(total_offset, ra, dec)
        )

    return np.array(corrections)


def apply_light_deflection_correction_to_observation_collection(
        observation_collection: ObservationCollection,
        bodies: SystemOfBodies,
        body_name: str,
        observer_body_name: str,
        observer_reference_name: str | None = None,
        perturbing_bodies_list: Iterable[str] = ('Sun',),
        in_place: bool = True
) -> ObservationCollection | None:
    """
    Computes corrections using
    :func:`~tudatpy.estimation.observations.observation_corrections.light_deflection_correction.light_deflection_correction_angular_observations`
    and applies them to an :class:`~tudatpy.estimation.observations.ObservationCollection` object.

    Parameters
    ----------
    observation_collection : :class:`~tudatpy.estimation.observations.ObservationCollection`
        ObservationCollection object containing the angular observations
    bodies : SystemOfBodies
        SystemOfBodies object
    body_name : str
        Name of the observed body object
    observer_body_name : str
        Name of the observer body
    observer_reference_name : str
        Name of the reference point on the observer body. If not given, it is assumed that the observer location
            coincides with the 'observer_body_name' body center.
    perturbing_bodies_list : Iterable[str]
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
        perturbing_bodies_list=perturbing_bodies_list
    )








