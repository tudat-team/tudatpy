"""
Utility functions for observation corrections
"""
import numpy as np
from numpy.linalg import norm
from tudatpy.estimation.observations import ObservationCollection, create_new_observation_collection
from tudatpy.estimation.observations.observations_processing import observation_parser
from tudatpy.estimation.observable_models_setup.model_settings import ObservableType
from tudatpy.dynamics.environment import SystemOfBodies
from collections.abc import Callable

def _unit(vector:np.ndarray) -> np.ndarray:
    return vector / norm(vector)


def _offset_vector_to_corrections(offset_vec: np.ndarray,
                                       ra: float,
                                       dec: float):
    """
    Calculate corrections in right ascension and declination from a plane-of-sky offset vector.

    Args:
        offset_vec (np.ndarray): Offset vector such that true dir + offset = observed dir.
        ra: Right ascension (rad)
        dec: Declination (rad)

    Returns:
        float: correction in right ascension
        float: correction in declination
    """
    observed_dir = np.array([np.cos(ra)*np.cos(dec), np.sin(ra)*np.cos(dec), np.sin(dec)])
    # true dir. + offset = observed dir.
    true_dir = observed_dir - offset_vec # Small angle approximation
    true_dir = _unit(true_dir)

    ra_true = np.arctan2(true_dir[1], true_dir[0])
    dec_true = np.arctan2(true_dir[2], np.sqrt(true_dir[0] ** 2 + true_dir[1] ** 2))

    ra_corr = ((ra_true - ra) + np.pi) % (2 * np.pi) - np.pi # Wraps to correct range
    dec_corr = dec_true - dec

    return ra_corr, dec_corr

def _apply_corrections_to_observation_collection(
        observation_collection: ObservationCollection,
        body_name: str,
        bodies: SystemOfBodies,
        observer_body_name: str,
        observer_reference_name: str | None,
        correction_function: Callable,
        in_place: bool = True,
        **kwargs
) -> ObservationCollection | None:
    """
    Helper function that computes the corrections from correction_function, and applies them on an
    ObservationCollection.
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

    observations, times = observation_collection.get_concatenated_observations_and_times(parser)
    if len(observations) == 0:
        raise ValueError(f'ObservationCollection does not contain angular observations with specified link-ends.')

    observations = np.reshape(observations, (-1, 2))
    observations_array = np.column_stack((np.array(times), observations))

    # Compute corrections
    corrections = correction_function(
        observations = observations_array,
        bodies=bodies,
        body_name= body_name,
        observer_body_name= observer_body_name,
        observer_reference_name= observer_reference_name,
        **kwargs
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

