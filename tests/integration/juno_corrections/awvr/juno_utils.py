from collections import Counter

import numpy as np
from tudatpy import constants, data_input, estimation
from tudatpy.astro import time_representation
from tudatpy.astro.time_representation import DateTime
from tudatpy.data_input.environment_data import read_dsn_amc_file, spice
from tudatpy.dynamics import environment
from tudatpy.estimation import observable_models_setup
from tudatpy.estimation.observations.observations_processing import (
    ObservationSetSplitterType,
    observation_parser,
    observation_set_splitter,
)
from tudatpy.estimation.observations_setup import ancillary_settings
from tudatpy.estimation.observations_setup.ancillary_settings import FrequencyBands


def load_juno_default_kernels(unbundled_kernel_path: str):

    default_kernel_names = [
        "naif0012.tls",
        "pck00011.tpc",
        "inpop19a_TDB_m100_p100_spice.tpc",
        "NOE-5-2021.tpc",
    ]

    default_kernel_path = data_input.resource_paths.get_spice_kernel_path()

    # jup380 instead of jup310, de440s since de442s is not available locally
    unbundled_kernel_names = [
        "spk/de440s.bsp",
        "spk/jup380s.bsp",
    ]

    for name in default_kernel_names:
        spice.load_kernel(default_kernel_path + "/" + name)
    for name in unbundled_kernel_names:
        spice.load_kernel(unbundled_kernel_path + "/" + name)


def load_perijove_data(file_path):
    """
    Load perijove data from a text file.

    Parameters
    ----------
    file_path : str
        Path to the text file containing perijove data.

    Returns
    -------
    dict
        Dictionary mapping perijove numbers to their corresponding epochs (TDB).
    """
    raw = np.loadtxt(file_path, delimiter="\t", skiprows=1, dtype=str)

    # The times in the file are UTC (spacecraft event time), convert them to TDB
    time_scale_converter = time_representation.default_time_scale_converter()
    pj_to_epoch = {
        int(row[0]): time_scale_converter.convert_time_object(
            time_representation.utc_scale,
            time_representation.tdb_scale,
            DateTime.from_iso_string(row[1]).to_epoch_time_object(),
        )
        for row in raw
    }

    return pj_to_epoch


def awvr_rejection_reasons(record):
    """Conservative AMC screening; algorithm identifiers are retained metadata."""
    reasons = []
    for flag in (
        "cloud_flag",
        "wvr_flag",
        "mtp_flag",
        "surface_flag",
        "wind_flag",
        "path_delay_flag",
    ):
        if getattr(record, flag) != 0:
            reasons.append(flag)
    for component in ("wet", "dry"):
        delay = getattr(record, f"zenith_{component}_delay")
        if not np.isfinite(delay):
            reasons.append(f"nonfinite_{component}")
        elif delay < 0.0:
            reasons.append(f"negative_{component}")
        elif np.isclose(delay, 9.99, rtol=0.0, atol=1e-12):
            reasons.append(f"missing_{component}")  # AMC missing retrieval: 999 cm
    if not np.isfinite(record.utc_time):
        reasons.append("nonfinite_epoch")
    return reasons


def load_awvr_zenith_delays(amc_directory, unit="01"):
    """Read one redundant AWVR unit at DSS-25; return UTC -> [dry, wet] metres."""
    if unit not in ("01", "02"):
        raise ValueError("AWVR unit must be 01 or 02.")
    files = sorted(amc_directory.glob(f"*_{unit}_*.AMC.txt"))
    if not files:
        raise ValueError(f"No AMC files for unit {unit} in {amc_directory}.")
    delays = {}
    for path in files:
        records = read_dsn_amc_file(path.as_posix()).records
        rejected = Counter()
        algorithms = Counter()
        retained = 0
        for record in records:
            algorithms[record.algorithm_id] += 1
            reasons = awvr_rejection_reasons(record)
            if reasons:
                rejected.update(reasons)
                continue
            if record.utc_time in delays:
                raise ValueError(f"Overlapping AMC epochs in {path}.")
            delays[record.utc_time] = [record.zenith_dry_delay, record.zenith_wet_delay]
            retained += 1
        print(
            f"AWVR {path.name}: {retained}/{len(records)} samples retained; "
            f"algorithms {dict(algorithms)}; rejection reasons (may overlap) {dict(rejected)}"
        )
    if len(delays) < 2:
        raise ValueError("AWVR screening retained fewer than two samples.")
    return dict(sorted(delays.items()))


def awvr_segment_boundaries(zenith_delays, maximum_gap=30.01):
    """UTC endpoints of each uninterrupted group of accepted AMC samples."""
    epochs = sorted(zenith_delays)
    if not epochs:
        return np.array([])
    boundaries = [epochs[0]]
    for first, second in zip(epochs, epochs[1:]):
        if second - first > maximum_gap:
            boundaries.extend((first, second))
    boundaries.append(epochs[-1])
    return np.array(sorted(set(boundaries)))


def crosses_awvr_transition(first_utc, last_utc, boundaries):
    """Whether a station-event interval touches a change of correction source."""
    return np.searchsorted(boundaries, first_utc, side="left") < np.searchsorted(
        boundaries, last_utc, side="right"
    )


def exclude_awvr_source_transition_counts(
    observation_collection, bodies, zenith_delays, time_windows, integration_time=60.0
):
    """Apply the same source-boundary mask in TSAC and AWVR comparisons.

    Reference-SPICE geometric light times locate both Earth events at the count
    endpoints. A one-second guard covers small differences from the fitted link.
    This mask only removes counts crossing a source change; gaps use TSAC.
    """
    boundaries = awvr_segment_boundaries(zenith_delays)
    converter = time_representation.default_time_scale_converter()
    link_types = observable_models_setup.links.LinkEndType
    earth = bodies.get("Earth")
    dss25_position = earth.get_ground_station(
        "DSS-25"
    ).station_state.cartesian_position_at_reference_epoch

    def position(body_name, epoch):
        return spice.get_body_cartesian_state_at_epoch(body_name, "SSB", "J2000", "NONE", epoch)[:3]

    def station_position(station_name, epoch):
        fixed_position = earth.get_ground_station(
            station_name
        ).station_state.cartesian_position_at_reference_epoch
        return (
            position("Earth", epoch)
            + earth.rotation_model.body_fixed_to_inertial_rotation(epoch) @ fixed_position
        )

    def transmission_time(receive_time, transmitter_station, receiver_station):
        receive_position = station_position(receiver_station, receive_time)
        bounce_time = receive_time
        for _ in range(6):
            bounce_time = (
                receive_time
                - np.linalg.norm(position("JUNO", bounce_time) - receive_position)
                / constants.SPEED_OF_LIGHT
            )
        bounce_position = position("JUNO", bounce_time)
        transmit_time = bounce_time
        for _ in range(6):
            transmit_time = (
                bounce_time
                - np.linalg.norm(
                    bounce_position - station_position(transmitter_station, transmit_time)
                )
                / constants.SPEED_OF_LIGHT
            )
        return transmit_time

    retained_sets = []
    removed = 0
    for obs_set in observation_collection.get_single_observation_sets():
        link_ends = obs_set.link_definition.link_ends
        tx_station = link_ends[link_types.transmitter].reference_point
        rx_station = link_ends[link_types.receiver].reference_point
        epochs = obs_set.observation_times
        values = obs_set.list_of_observations
        retained = []
        for index, epoch in enumerate(epochs):
            receive_time = epoch.to_float() if hasattr(epoch, "to_float") else float(epoch)
            reject = False
            if "DSS-25" in (tx_station, rx_station) and any(
                start <= receive_time <= end for start, end in time_windows
            ):
                first_rx = receive_time - integration_time / 2.0 - 1.0
                last_rx = receive_time + integration_time / 2.0 + 1.0
                station_intervals = []
                if rx_station == "DSS-25":
                    station_intervals.append((first_rx, last_rx))
                if tx_station == "DSS-25":
                    station_intervals.append(
                        (
                            transmission_time(first_rx, tx_station, rx_station),
                            transmission_time(last_rx, tx_station, rx_station),
                        )
                    )
                for first_tdb, last_tdb in station_intervals:
                    first_utc = converter.convert_time(
                        time_representation.tdb_scale,
                        time_representation.utc_scale,
                        first_tdb,
                        dss25_position,
                    )
                    last_utc = converter.convert_time(
                        time_representation.tdb_scale,
                        time_representation.utc_scale,
                        last_tdb,
                        dss25_position,
                    )
                    reject |= crosses_awvr_transition(first_utc, last_utc, boundaries)
            if reject:
                removed += 1
            else:
                retained.append(index)
        if retained:
            retained_sets.append(
                estimation.observations.create_single_observation_set(
                    obs_set.observable_type,
                    link_ends,
                    [values[i] for i in retained],
                    [epochs[i] for i in retained],
                    obs_set.reference_link_end,
                    obs_set.ancillary_settings,
                )
            )
    print(f"AWVR source transitions: excluded {removed} counts (same mask for TSAC and AWVR).")
    return estimation.observations.ObservationCollection(retained_sets)


def juno_turnaround_ratio():
    "https://doi.org/10.1007/s11214-017-0428-7"

    return {
        (FrequencyBands.x_band, FrequencyBands.x_band): 880.0 / 749.0,
        (FrequencyBands.ka_band, FrequencyBands.ka_band): 3360.0 / 3599.0,
        (FrequencyBands.x_band, FrequencyBands.ka_band): 3360.0 / 749.0,
    }


def update_dss_25_ka_band_delays(
    observation_collection: estimation.observations.ObservationCollection,
):
    "buccino2021StabilityDSS25Kaband"

    # Define splitting epochs
    splitting_epochs = [DateTime(2020, 11, 24).to_epoch()]

    # Create a splitter to divide sets at specific times
    splitter = observation_set_splitter(
        ObservationSetSplitterType.time_tags_splitter, splitting_epochs
    )

    # Apply the splitter in-place
    observation_collection.split_observation_sets(splitter)
    radiometric_observables = [
        observable_models_setup.model_settings.ObservableType.dsn_n_way_averaged_doppler_type,
        observable_models_setup.model_settings.ObservableType.dsn_n_way_range_type,
    ]
    radiometric_parser = observation_parser(radiometric_observables)
    for single_obs_set in observation_collection.get_single_observation_sets(radiometric_parser):
        if single_obs_set.observable_type not in radiometric_observables:
            raise ValueError(
                f"Unexpected observable type {single_obs_set.observable_type} in observation set."
            )

        obs_set_bands = [
            FrequencyBands(int(x))
            for x in single_obs_set.ancillary_settings.get_float_list_settings(
                ancillary_settings.ObservationAncillarySimulationVariable.frequency_bands
            )
        ]

        if obs_set_bands != [FrequencyBands.ka_band, FrequencyBands.ka_band]:
            continue

        link_definition = single_obs_set.link_definition

        if (
            link_definition.link_end_id(
                observable_models_setup.links.LinkEndType.transmitter
            ).station_name
            != "DSS-25"
            and link_definition.link_end_id(
                observable_models_setup.links.LinkEndType.receiver
            ).station_name
            != "DSS-25"
        ):
            continue

        if (
            single_obs_set.time_bounds[0] < DateTime(2020, 11, 24).to_epoch()
            and single_obs_set.time_bounds[1] > DateTime(2020, 11, 25).to_epoch()
        ):
            raise ValueError("Observation set spans the splitting epoch, which should not happen.")

        original_link_ends_delays = single_obs_set.ancillary_settings.get_float_list_settings(
            ancillary_settings.ObservationAncillarySimulationVariable.link_ends_delays
        )
        updated_link_ends_delays = original_link_ends_delays.copy()

        if single_obs_set.time_bounds[1] < DateTime(2020, 11, 24).to_epoch():
            updated_link_ends_delays[0] = 99.5e-6
            updated_link_ends_delays[2] = 75.5e-6

        elif single_obs_set.time_bounds[0] > DateTime(2020, 11, 24).to_epoch():
            updated_link_ends_delays[0] = 78.5e-6
            updated_link_ends_delays[2] = 78.5e-6

        single_obs_set.ancillary_settings.set_float_list_settings(
            ancillary_settings.ObservationAncillarySimulationVariable.link_ends_delays,
            updated_link_ends_delays,
        )

    return observation_collection


def apply_spin_effect_correction(
    observation_collection: estimation.observations.ObservationCollection,
    bodies: environment.SystemOfBodies,
    apply_positive_correction: bool = True,
):
    """Apply spin effect correction to the observation collection."""

    n_way_doppler_parser = observation_parser(
        [observable_models_setup.model_settings.ObservableType.dsn_n_way_averaged_doppler_type]
    )

    turnaround_ratios = juno_turnaround_ratio()

    for single_obs_set in observation_collection.get_single_observation_sets(n_way_doppler_parser):
        if (
            single_obs_set.observable_type
            != observable_models_setup.model_settings.ObservableType.dsn_n_way_averaged_doppler_type
        ):
            raise ValueError(
                f"Unexpected observable type {single_obs_set.observable_type} in observation set."
            )

        obs_set_bands = [
            FrequencyBands(int(x))
            for x in single_obs_set.ancillary_settings.get_float_list_settings(
                ancillary_settings.ObservationAncillarySimulationVariable.frequency_bands
            )
        ]

        epochs = single_obs_set.observation_times

        # TODO: Use t2 for spin rate
        angular_velocity_magnitudes = [
            np.linalg.norm(
                bodies.get("JUNO").rotation_model.angular_velocity_in_body_fixed_frame(epoch)
            )
            / (2 * np.pi)
            for epoch in epochs
        ]

        # TODO: expose the turnaround ratio from the system model
        # turnaround_ratio = bodies.get("Juno").system_models.get_turnaround_ratio()
        turnaround_ratio = turnaround_ratios[*obs_set_bands]

        sign = 1 if apply_positive_correction else -1

        if obs_set_bands == [FrequencyBands.ka_band, FrequencyBands.ka_band]:
            sign = -sign

        spin_effect_corrections = [
            (1 + sign * turnaround_ratio) * angular_velocity
            for angular_velocity in angular_velocity_magnitudes
        ]

        observations = single_obs_set.list_of_observations

        single_obs_set.set_observations(
            [
                observation + correction
                for observation, correction in zip(observations, spin_effect_corrections)
            ]
        )

    return observation_collection
