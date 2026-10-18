import numpy as np
import pytest

from tudatpy.data_input.tracking_data.radio_utilities import (
    dual_frequency_plasma_correction,
)


@pytest.mark.parametrize(
    "uplink_frequencies, downlink_frequencies",
    [
        # Frequency diversity on uplink and downlink (e.g. X/X and Ka/Ka), downlink only, and uplink only
        ((7.18e9, 34.3e9), (8.43e9, 32.1e9)),
        ((7.18e9, 7.18e9), (8.43e9, 32.1e9)),
        ((7.18e9, 34.3e9), (8.43e9, 8.43e9)),
    ],
)
def test_dual_frequency_plasma_correction(uplink_frequencies, downlink_frequencies):
    rng = np.random.default_rng(0)
    number_of_observations = 5

    # Inverse squared frequencies, relative to 1 GHz, per (observation, channel, link)
    inverse_squared_frequencies = np.empty((number_of_observations, 2, 2))
    inverse_squared_frequencies[:, :, 0] = (1.0e9 / np.array(uplink_frequencies)) ** 2
    inverse_squared_frequencies[:, :, 1] = (1.0e9 / np.array(downlink_frequencies)) ** 2
    carrier_normalizations = np.tile(downlink_frequencies, (number_of_observations, 1))

    # Non-dispersive signal, common plasma, and different (modelled) uplink and downlink ionospheres
    signal = rng.normal(0.0, 1.0e-8, number_of_observations)
    plasma = rng.normal(0.0, 1.0e-12, number_of_observations)
    uplink_ionosphere = rng.normal(0.0, 1.0e-12, number_of_observations)
    downlink_ionosphere = rng.normal(0.0, 1.0e-12, number_of_observations)

    ionospheric_contributions = (
        inverse_squared_frequencies[:, :, 0] * uplink_ionosphere[:, None]
        + inverse_squared_frequencies[:, :, 1] * downlink_ionosphere[:, None]
    )
    plasma_contributions = np.sum(inverse_squared_frequencies, axis=2) * plasma[:, None]
    normalized_frequencies = signal[:, None] + ionospheric_contributions + plasma_contributions

    correction = dual_frequency_plasma_correction(
        normalized_frequencies,
        inverse_squared_frequencies,
        ionospheric_contributions,
        carrier_normalizations,
    )

    # Doppler observable decreases with received frequency, so the correction is the plasma contribution in Hz
    expected_correction = carrier_normalizations[:, 0] * plasma_contributions[:, 0]
    assert np.allclose(correction, expected_correction, rtol=1.0e-6, atol=0.0)
