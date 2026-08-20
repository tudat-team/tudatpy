import numpy as np
import pytest

from tudatpy.constants import SPEED_OF_LIGHT
from tudatpy.estimation.observable_models_setup.links import (
    body_origin_link_end_id,
    body_reference_point_link_end_id,
    receiver,
    reflector1,
    transmitter,
)
from tudatpy.estimation.observable_models_setup.model_settings import ObservableType
from tudatpy.estimation.observations import ObservationCollection, create_single_observation_set
from tudatpy.estimation.observations.observations_processing import observation_parser
from tudatpy.estimation.observations_setup.ancillary_settings import (
    FrequencyBands,
    ObservationAncillarySimulationSettings,
    ObservationAncillarySimulationVariable,
    dsn_default_turnaround_ratios,
    dsn_n_way_doppler_ancillary_settings,
)


def _n_way_link_ends(station_name):
    return {
        transmitter: body_reference_point_link_end_id("Earth", station_name),
        reflector1: body_origin_link_end_id("Spacecraft"),
        receiver: body_reference_point_link_end_id("Earth", station_name),
    }


def _one_way_link_ends(station_name):
    return {
        transmitter: body_origin_link_end_id("Spacecraft"),
        receiver: body_reference_point_link_end_id("Earth", station_name),
    }


def _scalar_observations(values):
    return [np.array([value]) for value in values]


def test_dsn_n_way_range_residuals_in_state_space():
    range_conversion_factor = 0.15
    ancillary_settings = ObservationAncillarySimulationSettings()
    ancillary_settings.set_float_settings(
        ObservationAncillarySimulationVariable.range_conversion_factor, range_conversion_factor
    )

    residuals = np.array([2.0, -4.0, 0.5])
    observation_set = create_single_observation_set(
        ObservableType.dsn_n_way_range_type,
        _n_way_link_ends("DSS-14"),
        _scalar_observations([10.0, 11.0, 12.0]),
        [1.0, 2.0, 3.0],
        receiver,
        ancillary_settings,
    )
    observation_set.set_residuals(residuals)

    assert observation_set.residual_state_space_conversion_factor == pytest.approx(range_conversion_factor)
    np.testing.assert_allclose(
        observation_set.concatenated_residuals_in_state_space,
        residuals * range_conversion_factor,
    )


def test_dsn_n_way_averaged_doppler_residuals_in_state_space():
    reference_frequency = 7.2e9
    ancillary_settings = dsn_n_way_doppler_ancillary_settings(
        [FrequencyBands.x_band, FrequencyBands.x_band],
        FrequencyBands.x_band,
        reference_frequency,
        60.0,
    )
    expected_factor = SPEED_OF_LIGHT / (
        2.0 * dsn_default_turnaround_ratios(FrequencyBands.x_band, FrequencyBands.x_band) * reference_frequency
    )

    residuals = np.array([0.25, -1.5])
    observation_set = create_single_observation_set(
        ObservableType.dsn_n_way_averaged_doppler_type,
        _n_way_link_ends("DSS-14"),
        _scalar_observations([8.4e9, 8.4e9]),
        [1.0, 2.0],
        receiver,
        ancillary_settings,
    )
    observation_set.set_residuals(residuals)

    assert observation_set.residual_state_space_conversion_factor == pytest.approx(expected_factor)
    np.testing.assert_allclose(
        observation_set.concatenated_residuals_in_state_space,
        residuals * expected_factor,
    )


def test_one_way_doppler_conversion_is_rejected():
    observation_set = create_single_observation_set(
        ObservableType.one_way_instantaneous_doppler_type,
        _one_way_link_ends("DSS-14"),
        _scalar_observations([0.1]),
        [1.0],
        receiver,
    )
    observation_set.set_residuals(np.array([0.2]))

    with pytest.raises(RuntimeError, match="not supported"):
        _ = observation_set.residual_state_space_conversion_factor
    with pytest.raises(RuntimeError, match="not supported"):
        _ = observation_set.concatenated_residuals_in_state_space


def test_collection_parser_order_and_unsupported_rejection():
    range_ancillary_settings = ObservationAncillarySimulationSettings()
    range_ancillary_settings.set_float_settings(ObservationAncillarySimulationVariable.range_conversion_factor, 2.0)
    doppler_ancillary_settings = dsn_n_way_doppler_ancillary_settings(
        [FrequencyBands.x_band, FrequencyBands.x_band],
        FrequencyBands.x_band,
        7.2e9,
        60.0,
    )

    range_set = create_single_observation_set(
        ObservableType.dsn_n_way_range_type,
        _n_way_link_ends("DSS-14"),
        _scalar_observations([1.0, 2.0]),
        [10.0, 11.0],
        receiver,
        range_ancillary_settings,
    )
    range_set.set_residuals(np.array([4.0, 5.0]))

    doppler_set = create_single_observation_set(
        ObservableType.dsn_n_way_averaged_doppler_type,
        _n_way_link_ends("DSS-15"),
        _scalar_observations([8.4e9]),
        [20.0],
        receiver,
        doppler_ancillary_settings,
    )
    doppler_set.set_residuals(np.array([-0.5]))

    unsupported_set = create_single_observation_set(
        ObservableType.one_way_instantaneous_doppler_type,
        _one_way_link_ends("DSS-16"),
        _scalar_observations([0.1]),
        [30.0],
        receiver,
    )
    unsupported_set.set_residuals(np.array([0.3]))

    supported_collection = ObservationCollection([range_set, doppler_set])
    native_residuals = supported_collection.get_concatenated_residuals()
    converted_residuals = supported_collection.get_concatenated_residuals_in_state_space()
    native_blocks = supported_collection.get_residuals()
    converted_blocks = supported_collection.get_residuals_in_state_space()

    assert len(converted_blocks) == len(native_blocks)
    expected = []
    for native_block, observation_set in zip(native_blocks, supported_collection.get_single_observation_sets()):
        expected.append(np.ravel(native_block) * observation_set.residual_state_space_conversion_factor)
    np.testing.assert_allclose(np.ravel(converted_residuals), np.concatenate(expected))
    np.testing.assert_equal(np.ravel(native_residuals).shape, np.ravel(converted_residuals).shape)

    parsed_native = supported_collection.get_concatenated_residuals(
        observation_parser(ObservableType.dsn_n_way_range_type)
    )
    parsed_converted = supported_collection.get_concatenated_residuals_in_state_space(
        observation_parser(ObservableType.dsn_n_way_range_type)
    )
    np.testing.assert_allclose(parsed_converted, parsed_native * 2.0)

    mixed_collection = ObservationCollection([range_set, doppler_set, unsupported_set])
    mixed_collection.get_concatenated_residuals()
    with pytest.raises(RuntimeError, match="not supported"):
        mixed_collection.get_concatenated_residuals_in_state_space()
    mixed_collection.get_concatenated_residuals_in_state_space(
        observation_parser(ObservableType.dsn_n_way_range_type)
    )
