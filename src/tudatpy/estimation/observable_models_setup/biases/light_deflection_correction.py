"""
Functions to calculate light deflection corrections to observations
"""
import numpy as np
from numpy.linalg import norm
from tudatpy.estimation.observations import ObservationCollection, create_new_observation_collection
from tudatpy.estimation.observations.observations_processing import observation_parser
from tudatpy.estimation.observable_models_setup.model_settings import ObservableType
from tudatpy.dynamics.environment import SystemOfBodies
from tudatpy.dynamics.environment_setup import create_ground_station_ephemeris
from tudatpy.estimation.observable_models_setup.biases.observation_correction_utils import _offset_vector_to_corrections
from tudatpy.constants import SPEED_OF_LIGHT

def _calculate_light_deflection(
        observer_pos: np.ndarray,
        body_pos: np.ndarray,
        perturber_pos: np.ndarray,
        mu_perturber: float
) -> np.ndarray:
    """
    Helper function to calculate post-Newtonian light-bending contribution from one body (delta k_pn in Klioner (2003))

    Args:
        observer_pos (np.ndarray): Position vector of the observer
        body_pos (np.ndarray): Position vector of the body being observed
        perturber_pos (np.ndarray): Position vector of the perturber body
        mu_perturber (float): Gravitational parameter of the perturber body

    Returns:
        np.ndarray: Individual contribution to total light deflection
    """
    # Define vectors from Klioner 2003
    r_upper = observer_pos - body_pos  # capital R
    r_ea = body_pos - perturber_pos
    r_oa = observer_pos - perturber_pos

    # Calculate contribution of deflection from perturber body
    factor = (2 * mu_perturber) / (SPEED_OF_LIGHT ** 2)
    num = np.cross(r_upper, np.cross(r_ea, r_oa))
    denom = norm(r_upper) * norm(r_oa) * (norm(r_ea) * norm(r_oa) + np.dot(r_oa, r_ea))
    offset_vec = factor * num / denom

    return offset_vec


def relativistic_light_deflection_from_observations(
        observations: np.ndarray,
        bodies: SystemOfBodies,
        body_name: str,
        observer_body_name: str,
        observer_reference_name: str = None,
        perturbing_bodies_list: list[str] = ('Sun',),
):
    """
    Compute corrections to observations for the relativstic deflection of light around massive bodies.

    Compute corrections to observations for the relativstic deflection of light around massive bodies. Currently,
    the implementation requires that a reference ephemeris for the observed body ('body_name') is provided in the
    bodies object. This can e.g. be retrieved from a high accuracy ephemeris source such as JPL Horizons.

    Parameters
    ----------
    observations : np.ndarray
        Array of angular observations with columns [time, RA, DEC]
    bodies : SystemOfBodies
        SystemOfBodies object
    body_name : str
        Name of the observed body object
    observer_body_name : str
        Name of the observer body
    observer_reference_name : str
        Name of the reference point on the observer body. If not given, it is assumed that the observer location
            coincides with the 'observer_body_name' body center.
    perturbing_bodies_list : list[str]
        Names of the bodies that light-deflection contribution should be computed for, default = 'Sun'

    Returns
    -------

    """
    # Input validation
    if observations.shape[1] != 3:
        raise ValueError(f'Observations must be in shape N x 3 with columns time, RA, DEC')
    body_dne = [not bodies.does_body_exist(name) or bodies.get(name).ephemeris is None
                for name in perturbing_bodies_list + [body_name, observer_body_name]]
    if any(body_dne):
        raise ValueError('Some or all included bodies in the relativistic light deflection computation are missing'
                         'from SystemOfBodies or their associated ephemerides are not specified.')

    # Get observer reference point ephemeris if it exists
    if observer_reference_name is not None:
        observer_ephemeris = create_ground_station_ephemeris(
            bodies.get(observer_body_name),
            observer_reference_name,
            bodies
        ) # -> in global frame origin/orientation

    corrections = []

    # Loop over observations
    for epoch, ra, dec in observations.T:
        # Position of observer at current epoch:
        if observer_reference_name is not None:
            observer_pos = observer_ephemeris.cartesian_position(epoch)
        else:
            observer_pos = bodies.get(observer_body_name).state_in_base_frame_from_ephemeris(epoch)[:3]

        body_pos = bodies.get(body_name).state_in_base_frame_from_ephemeris(epoch)[:3]

        # One iteration of light-time calculation for light from observed body to observer
        light_time_body_to_observer = np.linalg.norm(
            observer_pos - body_pos
        ) / SPEED_OF_LIGHT

        # Observed body position at time of emission of light:
        body_pos = bodies.get(body_name).state_in_base_frame_from_ephemeris(epoch - light_time_body_to_observer)[:3]

        total_offset = np.zeros(3) # Total deflection accumulated from all bodies

        # Loop over light-deflecting bodies
        for perturber_name in perturbing_bodies_list:

            mu_perturber = bodies.get(perturber_name).gravitational_parameter
            perturber_pos = bodies.get(perturber_name).state_in_base_frame_from_ephemeris(epoch)[:3]

            # Compute one iteration of light-time from observer to perturbing body
            light_time_perturber_to_observer = np.linalg.norm(
                observer_pos - perturber_pos
            ) / SPEED_OF_LIGHT

            # Perturber position at reference time
            perturber_pos = bodies.get(perturber_name).state_in_base_frame_from_ephemeris(epoch - light_time_perturber_to_observer)[:3]

            offset_from_perturber = _calculate_light_deflection(
                observer_pos = observer_pos,
                body_pos = body_pos,
                perturber_pos = perturber_pos,
                mu_perturber = mu_perturber,
            )
            total_offset += offset_from_perturber # Note: we omit the minus sign from Klioner eq. 70 since this is just an odd convention in observation directions


        corrections.append(
            _offset_vector_to_corrections(total_offset, ra, dec)
        )

    return np.array(corrections)


def apply_light_deflection_correction_to_observation_collection(
        observation_collection: ObservationCollection,
        bodies: SystemOfBodies,
        body_name: str,
        observer_body_name: str,
        observer_reference_name: str = None,
        perturbing_bodies_list: list[str] = ('Sun',),
        in_place: bool = False
) -> ObservationCollection:
    """
    Compute corrections to observations for the relativstic deflection of light around massive bodies, and apply
    to an observation collection.

    Compute corrections to observations for the relativstic deflection of light around massive bodies, and apply
    to an observation collection. Currently,the implementation requires that a reference ephemeris for the observed body
    ('body_name') is provided in the bodies object. This can e.g. be retrieved from a high accuracy ephemeris source
    such as JPL Horizons.

    Parameters
    ----------
    observation_collection : ObservationCollection
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
    perturbing_bodies_list : list[str]
        Names of the bodies that light-deflection contribution should be computed for, default = 'Sun'
    in_place : bool
        If true, corrections are applied in-place to the Observationcollection object. If false, a new Observationcollection
            is returned with the corrections applied.
    Returns
    -------
    None | ObservationCollection
    """
    # Parser to obtain angular observations for specified observer
    parsers = []
    parsers.append(observation_parser(ObservableType.angular_position_type))
    parsers.append(observation_parser(body_name))
    if observer_reference_name is not None:
        parsers.append(observation_parser(observer_reference_name, is_reference_point=True))
    else:
        parsers.append(observation_parser(observer_body_name))
    parser = observation_parser(parsers, combine_conditions=True)

    observations = np.reshape(
        observation_collection.get_concatenated_observations(parser),
        (-1, 2)
    )
    if len(observations) == 0:
        raise ValueError(f'ObservationCollection does not contain angular observations with specified link-ends.')

    # Compute corrections
    corrections = relativistic_light_deflection_from_observations(
        observations = observations,
        bodies=bodies,
        body_name= body_name,
        observer_body_name= observer_body_name,
        observer_reference_name= observer_reference_name,
        perturbing_bodies_list= perturbing_bodies_list,
    )

    corrected_observations = observations + corrections

    # Wrap RA
    corrected_observations[:,0] = (corrected_observations[:,0] + np.pi) % (2 * np.pi) - np.pi

    if in_place: # Apply to original observation collection
        observation_collection.set_observations(corrected_observations.flatten(), parser)
        return None

    else: # Create new observation collection
        new_observation_collection = create_new_observation_collection(observation_collection)
        new_observation_collection.set_observations(corrected_observations.flatten(), parser)

        return new_observation_collection












