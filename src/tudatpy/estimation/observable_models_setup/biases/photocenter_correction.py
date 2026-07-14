"""
Functions to calculate photocenter corrections to observations
"""
import numpy as np
from tudatpy.estimation.observations import ObservationCollection, create_new_observation_collection
from tudatpy.estimation.observations.observations_processing import observation_parser
from tudatpy.estimation.observable_models_setup.model_settings import ObservableType
from tudatpy.dynamics.environment_setup import create_ground_station_ephemeris
from tudatpy.dynamics.environment import SystemOfBodies
from numpy.linalg import norm
from tudatpy.estimation.observable_models_setup.biases.observation_correction_utils import _offset_vector_to_corrections, _unit


def _photocenter_offset(diameter : float,
                        body_wrt_observer : np.ndarray,
                        body_wrt_sun : np.ndarray,) -> np.ndarray:
    """Helper function to calculate photocenter offset at one epoch"""
    # Get unit vectors
    body_wrt_observer_unit = _unit(body_wrt_observer)
    body_wrt_sun_unit = _unit(body_wrt_sun)

    # Angle spanned by observer, body and Sun:
    solar_phase_angle = np.arccos(np.dot(body_wrt_observer_unit, body_wrt_sun_unit))

    offset_dir = - (body_wrt_sun_unit - np.dot(body_wrt_sun_unit, body_wrt_observer_unit)
                    * body_wrt_observer_unit)
    # Calculate fraction of radius according to Fuentes-Munoz (2024):
    cot = lambda x: np.cos(x) / np.sin(x)
    num = 2 * (np.sin(solar_phase_angle) + (np.pi - solar_phase_angle) * np.cos(solar_phase_angle))
    denom = 3 * np.pi * (
            cot(solar_phase_angle / 2) - np.sin(solar_phase_angle / 2) * np.log(cot(solar_phase_angle / 4)))
    offset_ratio = num / denom # Fraction of body radius

    offset_mag = offset_ratio * (diameter/2) / norm(body_wrt_observer)

    return offset_dir * offset_mag


def photocenter_corrections_from_observations(
        observations: np.ndarray ,
        diameter: float,
        bodies: SystemOfBodies,
        body_name: str,
        observer_body_name: str,
        observer_reference_name: str = None
) -> np.ndarray:
    """
    Calculate corrections to observations to account for the photocenter-barycenter offset.

    Calculate corrections to observations to account for the photocenter-barycenter offset. Currently, this function
    requires that the bodies object contains an ephemeris for the observed body ('body_name') in order to retrieve its
    inertial position as a function of time. The corrections that are output should be added to observations.

    Note that light-time offsets are small and are  neglected accordingly.

    Parameters
    ----------
    observations: np.ndarray
        Observations as [time, RA, DEC] with time in seconds since J2000, RA, DEC as radians
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
        Corrections as [RA, DEC]
    """
    # Validate inputs
    body_dne = [not bodies.does_body_exist(name) or bodies.get(name).ephemeris is None
                for name in (body_name, observer_body_name, 'Sun')]
    if any(body_dne):
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

    # Corrections to RA/DEC observations
    corrections = []

    for epoch, ra, dec in observations:

        # Retrieve body, Sun and observer positions in common frame
        body_pos = bodies.get(body_name).state_in_base_frame_from_ephemeris(epoch)[:3]
        sun_pos = bodies.get('Sun').state_in_base_frame_from_ephemeris(epoch)[:3]

        if observer_reference_name is not None: # Retrieve from reference point's ephemeris
            observer_pos = observer_ephemeris.cartesian_position(epoch)
        else: # Retrieve from observer body ephemeris
            observer_pos = bodies.get(observer_body_name).state_in_base_frame_from_ephemeris(epoch)[:3]

        offset_vector = _photocenter_offset(diameter,
                                            (body_pos - observer_pos),
                                            (body_pos - sun_pos))

        # Calculate corrections to the astrometry
        corrections.append(
            _offset_vector_to_corrections(offset_vector, ra, dec)
        )

    corrections = np.array(corrections)

    return corrections

def apply_photocenter_correction_to_observation_collection(
        observation_collection: ObservationCollection,
        diameter: float,
        bodies: SystemOfBodies,
        body_name: str,
        observer_body_name: str,
        observer_reference_name: str = None,
        in_place: bool = False
) -> ObservationCollection:
    """
    Computes the photocenter-barycenter offset with 'photocenter_corrections_from_observations' and applies the
    corrections to an ObservationCollection.

    Currently, this function requires that the bodies object contains an ephemeris for the observed body ('body_name')
    in order to retrieve its inertial position as a function of time. .

    Parameters
    ----------
    observation_collection : ObservationCollection
        Uncorrected observation collection
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
    None | ObservationCollection
        Returns None, or the new ObservationCollection depending on the argument 'in_place'
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

    # Compute photocenter corrections
    corrections = photocenter_corrections_from_observations(
        observations = observations,
        diameter = diameter,
        bodies = bodies,
        body_name = body_name,
        observer_body_name = observer_body_name,
        observer_reference_name = observer_reference_name,
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



