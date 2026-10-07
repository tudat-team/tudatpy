import numpy as np
import pytest

from tudatpy.astro import time_representation as time
from tudatpy.data_input.tracking_data import TrackingData, get_tracking_data_epoch_bounds
from tudatpy.data_input.environment_data import spice


def tracking_data(epochs, scale="TDB"):
    return TrackingData(
        observable_type="AngularPosition",
        link_ends=[(("673", ""), "transmitter"), (("Gaia", ""), "receiver")],
        observations=[[0.1, 0.2] for _ in epochs],
        epochs=epochs,
        reference_link_end="receiver",
        time_scale=scale,
    )


@pytest.mark.parametrize(
    "output_scale", [time.tdb_scale, time.utc_scale, time.tt_scale, time.tai_scale, time.ut1_scale]
)
def test_unsorted_mixed_scales_and_requested_output(output_scale):
    sources = [
        tracking_data([1000.0, -200.0, 0.0], "UTC"),
        tracking_data([500.0, -1000.0, 200.0], "TDB"),
    ]
    converter = time.default_time_scale_converter()
    expected = [
        converter.convert_time(input_scale, output_scale, epoch)
        for epochs, input_scale in [
            ([1000.0, -200.0, 0.0], time.utc_scale),
            ([500.0, -1000.0, 200.0], time.tdb_scale),
        ]
        for epoch in epochs
    ]
    actual = get_tracking_data_epoch_bounds(sources, output_time_scale=output_scale)
    np.testing.assert_allclose(
        [float(value) for value in actual], [min(expected), max(expected)], rtol=0.0, atol=1.0e-9
    )
    np.testing.assert_array_equal(np.asarray(sources[0].epochs, dtype=float), [1000.0, -200.0, 0.0])


def test_default_tdb_preserves_native_time_precision():
    first = time.Time(1.0e9) + time.Time(1.0e-12)
    last = first + time.Time(2.0e-12)
    assert first < last and float(first) == float(last)
    bounds = get_tracking_data_epoch_bounds([tracking_data([last, first])])
    assert bounds == (first, last)
    assert all(isinstance(value, time.Time) for value in bounds)


@pytest.mark.parametrize("sources", [[], [tracking_data([])], [None]])
def test_missing_epochs_raise(sources):
    with pytest.raises(ValueError):
        get_tracking_data_epoch_bounds(sources)


def test_empty_objects_are_skipped():
    assert tuple(
        map(float, get_tracking_data_epoch_bounds([tracking_data([]), tracking_data([2.0, -3.0])]))
    ) == (-3.0, 2.0)


def test_invalid_output_scale_is_rejected():
    with pytest.raises(ValueError, match="support"):
        get_tracking_data_epoch_bounds(
            [tracking_data([0.0])], output_time_scale=time.local_proper_time_scale
        )


@pytest.mark.parametrize(
    "number, identifier",
    [(1, "2000001"), (16, "2000016"), (673, "2000673"), (67125, "2067125"), (999999, "2999999")],
)
def test_original_asteroid_spice_identifier(number, identifier):
    assert spice.asteroid_spice_id(number) == identifier


@pytest.mark.parametrize("number", [-1, 0, 1000000])
def test_invalid_original_asteroid_number_is_rejected(number):
    with pytest.raises(ValueError):
        spice.asteroid_spice_id(number)
