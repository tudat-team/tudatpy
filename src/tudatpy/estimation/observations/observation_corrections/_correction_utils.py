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


def _unit(vector: np.ndarray) -> np.ndarray:
    return vector / norm(vector)


def _offset_vector_to_corrections(
    plane_of_sky_offset: np.ndarray,
    right_ascension: float,
    declination: float,
):
    """
    Convert a plane-of-sky offset vector, i.e. the offset that is perpendicular to the line-of-sight, to
    corrections in right ascension and declination.
    """
    observed_direction = np.array(
        [
            np.cos(right_ascension) * np.cos(declination),
            np.sin(right_ascension) * np.cos(declination),
            np.sin(declination),
        ]
    )
    # True direction + plane-of-sky offset = observed direction.
    true_direction = _unit(observed_direction - plane_of_sky_offset)  # Small-angle approximation

    true_right_ascension = np.arctan2(true_direction[1], true_direction[0])
    true_declination = np.arctan2(
        true_direction[2],
        np.sqrt(true_direction[0] ** 2 + true_direction[1] ** 2),
    )

    right_ascension_correction = ((true_right_ascension - right_ascension) + np.pi) % (
        2 * np.pi
    ) - np.pi
    declination_correction = true_declination - declination

    return right_ascension_correction, declination_correction


def _apply_corrections_to_observation_collection(
    observation_collection: ObservationCollection,
    body_name: str,
    bodies: SystemOfBodies,
    observer_body_name: str,
    observer_reference_name: str | None,
    correction_function: Callable,
    in_place: bool = True,
    **kwargs,
) -> ObservationCollection | None:
    """
    Helper function that computes the corrections from correction_function, and applies them on an
    ObservationCollection.
    """
    # Parser to obtain angular observations for specified observer
    observation_parsers = []
    observation_parsers.append(observation_parser(ObservableType.angular_position_type))
    observation_parsers.append(observation_parser(body_name))
    if observer_reference_name is not None:
        observation_parsers.append(
            observation_parser(observer_reference_name, is_reference_point=True)
        )
    else:
        observation_parsers.append(observation_parser(observer_body_name))
    angular_observation_parser = observation_parser(observation_parsers, combine_conditions=True)

    concatenated_observations, observation_times = (
        observation_collection.get_concatenated_observations_and_times(angular_observation_parser)
    )
    if len(concatenated_observations) == 0:
        raise ValueError(
            f"ObservationCollection does not contain angular observations with specified link-ends."
        )

    angular_observations = np.reshape(concatenated_observations, (-1, 2))
    observations_with_times = np.column_stack((np.array(observation_times), angular_observations))

    # Compute corrections
    angular_corrections = correction_function(
        observations=observations_with_times,
        bodies=bodies,
        body_name=body_name,
        observer_body_name=observer_body_name,
        observer_reference_name=observer_reference_name,
        **kwargs,
    )

    corrected_angular_observations = angular_observations + angular_corrections

    # Wrap RA
    corrected_angular_observations[:, 0] = (corrected_angular_observations[:, 0] + np.pi) % (
        2 * np.pi
    ) - np.pi

    if in_place:  # Apply to original observation collection
        observation_collection.set_observations(
            corrected_angular_observations.flatten(), angular_observation_parser
        )
        return None

    else:  # Create new observation collection
        new_observation_collection = create_new_observation_collection(observation_collection)
        new_observation_collection.set_observations(
            corrected_angular_observations.flatten(), angular_observation_parser
        )

        return new_observation_collection
