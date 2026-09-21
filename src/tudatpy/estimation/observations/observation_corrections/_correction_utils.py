"""
Utility functions for observation corrections
"""

import numpy as np
from numpy.linalg import norm
from tudatpy.estimation.observations import (
    ObservationCollection,
    ObservationDataset,
    ObservationSelectionCondition,
    LinkEndId,
    angular_position,
    create_observation_collection_from_dataset,
    create_observation_dataset_from_collection,
    observation_query,
    receiver,
    transmitter,
)
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


def _epoch_as_float(epoch) -> float:
    if hasattr(epoch, "to_float"):
        return float(epoch.to_float())
    return float(epoch)


def _apply_corrections_to_observation_dataset(
    observation_dataset: ObservationDataset,
    body_name: str,
    bodies: SystemOfBodies,
    observer_body_name: str,
    observer_reference_name: str | None,
    correction_function: Callable,
    in_place: bool = True,
    **kwargs,
) -> ObservationDataset | None:
    """
    Compute corrections first, then apply them atomically to matching dataset rows.
    """
    observer_link_end = LinkEndId(observer_body_name, observer_reference_name or "")
    condition = (
        (observation_query.observable_type == angular_position)
        & (observation_query.transmitter == LinkEndId(body_name, ""))
        & (observation_query.receiver == observer_link_end)
    )
    selected = observation_dataset.get_data(
        condition,
        fields=("times", "observations", "observation_ids", "set_ids"),
        ordering="internal",
    )
    if not selected["observation_ids"]:
        raise ValueError(
            "ObservationDataset does not contain angular observations with the specified link ends."
        )

    for set_id in set(selected["set_ids"]):
        if observation_dataset.get_observation_set_metadata(set_id).reference_link_end != receiver:
            raise ValueError("Angular correction functions require reception-time observations")

    angular_observations = np.asarray(selected["observations"])
    observations_with_times = np.column_stack(
        (np.array([_epoch_as_float(epoch) for epoch in selected["times"]]), angular_observations)
    )

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

    target_dataset = observation_dataset
    if not in_place:
        target_dataset = observation_dataset.create_new_and_keep(
            ObservationSelectionCondition.all()
        )

    replacements_by_id = dict(
        zip(selected["observation_ids"], corrected_angular_observations, strict=True)
    )
    for set_id in set(selected["set_ids"]):
        observation_ids = target_dataset.observation_ids_for_set(set_id)
        observations = target_dataset.observations_for_set(set_id)
        for row_index, observation_id in enumerate(observation_ids):
            if observation_id in replacements_by_id:
                observations[row_index] = replacements_by_id[observation_id]
        target_dataset.set_observations_for_set(set_id, observations)

    return None if in_place else target_dataset


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
    """Compatibility adapter for the dataset correction implementation."""
    observation_dataset = create_observation_dataset_from_collection(observation_collection)
    result = _apply_corrections_to_observation_dataset(
        observation_dataset=observation_dataset,
        body_name=body_name,
        bodies=bodies,
        observer_body_name=observer_body_name,
        observer_reference_name=observer_reference_name,
        correction_function=correction_function,
        in_place=in_place,
        **kwargs,
    )
    if in_place:
        return None

    return create_observation_collection_from_dataset(result)
