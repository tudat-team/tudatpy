"""Utilities for the processing of radiometric tracking data."""

import numpy as np

from tudatpy.kernel.data_input.tracking_data import get_dsn_doppler_frequency_calibration


def dual_frequency_plasma_correction(
    normalized_frequencies,
    inverse_squared_frequencies,
    ionospheric_contributions,
    carrier_normalizations,
):
    """Compute the dual-frequency plasma correction to the Doppler observations of channel 0.

    For each pair of simultaneous Doppler observations at two (sets of) frequencies, the normalized received frequency
    of channel ``i`` is modelled as :math:`z_{i}=s+a_{i}(I_{u,i}+P)+b_{i}(I_{d,i}+P)`, with :math:`s` the
    non-dispersive signal, :math:`a_{i}` and :math:`b_{i}` the uplink and downlink inverse squared frequencies,
    :math:`I_{u,i}` and :math:`I_{d,i}` the (independently modelled) ionospheric contributions, and :math:`P` a plasma
    contribution that is common to uplink and downlink (for instance, the Io plasma torus close to the spacecraft).
    The plasma contribution is estimated from the difference between the two channels (after removing the
    ionospheric contributions), and the correction to the Doppler observable of channel 0 is computed from it.
    The ionospheric contribution is not part of the correction, and should be retained in the observation model.

    Parameters
    ----------
    normalized_frequencies : numpy.ndarray
        Received frequencies, normalized by the carrier normalization, of shape (N, 2) (observation, channel).
    inverse_squared_frequencies : numpy.ndarray
        Uplink and downlink inverse squared frequencies (in any constant unit), of shape (N, 2, 2) (observation,
        channel, link).
    ionospheric_contributions : numpy.ndarray
        Modelled ionospheric contribution to the normalized received frequency, of shape (N, 2).
    carrier_normalizations : numpy.ndarray
        Carrier normalization (turnaround ratio times transmitted frequency) [Hz], of shape (N, 2).

    Returns
    -------
    numpy.ndarray
        Correction [Hz] to be added to the Doppler observations of channel 0, of shape (N,).
    """
    plasma_responses = np.sum(inverse_squared_frequencies, axis=2)
    common_plasma = (
        normalized_frequencies[:, 0]
        - normalized_frequencies[:, 1]
        - ionospheric_contributions[:, 0]
        + ionospheric_contributions[:, 1]
    ) / (plasma_responses[:, 0] - plasma_responses[:, 1])
    return carrier_normalizations[:, 0] * plasma_responses[:, 0] * common_plasma


def compute_dual_frequency_doppler_corrections(
    target,
    calibrator,
    bodies,
    transmitting_frequency_calculators,
    non_dispersive_corrections,
    ionospheric_correction,
):
    """Compute the dual-frequency plasma correction to two-way DSN Doppler tracking data.

    The correction is computed for all observations of ``target`` for which ``calibrator`` contains an observation at
    the same (UTC reception) epoch, using :func:`dual_frequency_plasma_correction`. The two data sets must have the same
    link ends and count interval, but differ in frequency (on the uplink, downlink, or both). The tracking data must be
    receiver-tagged, with a ramped transmitter and constant reception reference frequency (as for DSN ODF data).
    Any corrections already set in the tracking data (e.g. for spacecraft spin) are applied before computing
    the plasma correction, as are link end delays in the ancillary settings. The light time is computed from the
    current ephemerides in ``bodies``, which should be the reference (e.g. reconstructed) ephemerides.

    The correction can be applied to the target data by adding it to the observation corrections (with zero where no
    correction is available), before creating the observation collection with ``apply_corrections=True``.

    Parameters
    ----------
    target : TrackingData
        Doppler data to which the correction is to be applied.
    calibrator : TrackingData
        Doppler data at a different frequency, used to calibrate the target data.
    bodies : SystemOfBodies
        System of bodies, containing the ramp tables of the transmitting stations.
    transmitting_frequency_calculators : list[TransmittingFrequencyCalculator]
        Transmitting frequency calculators for target and calibrator data.
    non_dispersive_corrections : list[LightTimeCorrectionSettings]
        Non-dispersive light-time corrections, used to compute the light time.
    ionospheric_correction : LightTimeCorrectionSettings
        Ionospheric light-time correction, evaluated independently on uplink and downlink.

    Returns
    -------
    numpy.ndarray
        Correction [Hz] to be added to the target observations, NaN where no calibrator observation is available.
    """
    integration_time_key = "Doppler observable integration time"
    if (
        target.link_ends != calibrator.link_ends
        or target.get_ancillary_settings_double()[integration_time_key]
        != calibrator.get_ancillary_settings_double()[integration_time_key]
    ):
        raise ValueError(
            "Error in dual-frequency Doppler correction: target and calibrator data must have the same link ends and count interval"
        )

    calibrator_indices = {epoch: index for index, epoch in enumerate(calibrator.epochs)}
    target_indices = [
        index for index, epoch in enumerate(target.epochs) if epoch in calibrator_indices
    ]
    other_indices = [calibrator_indices[target.epochs[index]] for index in target_indices]

    target_calibration, other_calibration = (
        np.asarray(
            get_dsn_doppler_frequency_calibration(
                data,
                indices,
                bodies,
                frequency_calculator,
                non_dispersive_corrections,
                ionospheric_correction,
            ),
            dtype=np.longdouble,
        )
        for data, frequency_calculator, indices in zip(
            (target, calibrator),
            transmitting_frequency_calculators,
            (target_indices, other_indices),
        )
    )

    # A different emission interval duration (e.g. due to different link end delays) is not a dispersive effect
    normalized_frequencies = np.column_stack(
        (
            target_calibration[:, 0],
            other_calibration[:, 0]
            + (target_calibration[:, 5] - other_calibration[:, 5])
            / target.get_ancillary_settings_double()[integration_time_key],
        )
    )
    corrections = np.full(len(target.epochs), np.nan)
    corrections[target_indices] = dual_frequency_plasma_correction(
        normalized_frequencies,
        np.stack((target_calibration[:, 3:5], other_calibration[:, 3:5]), axis=1),
        np.column_stack((target_calibration[:, 2], other_calibration[:, 2])),
        np.column_stack((target_calibration[:, 1], other_calibration[:, 1])),
    )
    return corrections
