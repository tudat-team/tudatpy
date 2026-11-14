"""Sigma rejection uses measurement covariance without modifying residuals or weights."""

import numpy as np

from tudatpy.estimation import observations
from tudatpy.estimation.observable_models_setup import links, model_settings


def angular_dataset(receiver="Earth"):
    dataset = observations.ObservationDataset()
    link = links.LinkDefinition(
        {
            links.transmitter: links.LinkEndId("673", ""),
            links.receiver: links.LinkEndId(receiver, ""),
        }
    )
    dataset.add_observation_set(
        model_settings.angular_position_type, link, [[0.1, 0.2]] * 3, [1, 2, 3], links.receiver
    )
    return dataset


def test_sigma_rejection_selects_whole_event_and_leaves_other_sets_unchanged():
    dataset = angular_dataset()
    other = angular_dataset("Gaia")
    other.set_residuals_for_set(0, [[100, 100]] * 3)
    dataset.add_observation_set_from_dataset(other, 0)
    dataset.set_weight_vector_for_set(0, [4, 4, 4, 4, 4, 4])
    dataset.set_residuals_for_set(0, [[2.5, 0], [0, -2.6], [1, 1]])
    query = observations.observation_query
    selected = query.receiver == links.LinkEndId("Earth", "")
    before = dataset.get_data(fields=("residuals", "weight_matrix"))
    dataset.reject_observations(selected, sigma_limit=5.0, reason="five sigma")
    assert dataset.get_observation_ids(query.rejected) == [1]
    assert dataset.observation_row(1).rejection_reason == "five sigma"
    np.testing.assert_array_equal(dataset.get_residuals(), before["residuals"])
    np.testing.assert_array_equal(
        dataset.get_weight_matrix().toarray(), before["weight_matrix"].toarray()
    )


def test_correlated_weights_use_marginal_covariance_of_complete_set():
    dataset = angular_dataset()
    covariance = np.eye(6)
    covariance[0, 2] = covariance[2, 0] = 0.9
    dataset.set_weight_matrix_for_set(0, np.linalg.inv(covariance))
    dataset.set_residuals_for_set(0, [[3, 0], [0, 0], [0, 6]])
    query = observations.observation_query
    dataset.reject_observations(query.time == 2, "previous rejection")
    dataset.reject_observations(query.active, 5.0)
    # Conditioning on the other CCD would incorrectly reject the 3-sigma row.
    assert dataset.get_observation_ids(query.rejected) == [1, 2]
    assert dataset.observation_row(1).rejection_reason == "previous rejection"


def test_empty_selection_and_existing_selection_only_overload():
    dataset = angular_dataset()
    query = observations.observation_query
    dataset.reject_observations(query.time > 10, sigma_limit=5.0)
    assert dataset.get_observation_ids(query.rejected) == []
    dataset.reject_observations(query.time == 1, "selection only")
    assert dataset.get_observation_ids(query.rejected) == [0]
