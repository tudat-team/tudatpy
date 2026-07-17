import datetime
from pathlib import Path
from bisect import bisect_left, bisect_right
from urllib.request import urlretrieve

import numpy as np
import pytest

from tudatpy.data_input import resource_paths as data_paths
from tudatpy.data_input.environment_data.horizons import HorizonsQuery
from tudatpy.data_input.tracking_data import (
    TrackingData,
    TrackingSupplementaryData,
    TranslationalStateSupplementaryData,
)
from tudatpy.data_input.tracking_data.fdets import FdetDateFormat, read_fdets_data
from tudatpy.data_input.tracking_data.ifms import read_ifms_data
from tudatpy.data_input.tracking_data.jpl_radar import JPLRadarQuery
from tudatpy.data_input.tracking_data.mpc import BatchMPC
from tudatpy.data_input.tracking_data.odf import read_odf_data
from tudatpy.data_input.tracking_data.psf import read_psf_data
from tudatpy.data_input.tracking_data.radar_utilities import (
    DOPPLER_OBSERVABLE,
    RANGE_OBSERVABLE,
    radar_data_to_tracking_data,
)
from tudatpy.data_input.tracking_data.radar_utilities.stations import (
    add_all_radar_ground_stations,
)
from tudatpy.data_input.tracking_data.tnf import read_tnf_data
from tudatpy.astro import time_representation
from tudatpy.estimation.observations import (
    create_observation_collection_from_tracking_data,
    create_compressed_doppler_collection,
    compute_residuals_and_dependent_variables,
    simulate_observations,
    set_tracking_supplementary_data_in_bodies,
)
from tudatpy.kernel import constants
from tudatpy.kernel.dynamics import environment, environment_setup
from tudatpy.kernel.dynamics.environment_setup import (
    ephemeris,
    ground_station,
    rotation_model,
    shape,
    shape_deformation,
)
from tudatpy.kernel.estimation.observable_models_setup import (
    biases,
    light_time_corrections,
    links,
    model_settings,
)
from tudatpy.kernel.estimation.observations_setup import (
    ancillary_settings,
    observations_simulation_settings,
)
from tudatpy.kernel.interface import spice


def _test_data_path() -> Path:
    return Path(data_paths.get_test_data_path())


def _download_file(url: str, directory: Path) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    file_path = directory / url.rsplit("/", 1)[1]
    if not file_path.exists() or file_path.stat().st_size == 0:
        urlretrieve(url, file_path)
    return file_path


def _create_mex_bodies(test_data_path: Path):
    spice.load_standard_kernels()
    spice.load_kernel(
        str(
            test_data_path
            / "estrack_n_way_doppler_observation_model"
            / "MEX_ROB_130101_131231_001_shortened.bsp"
        )
    )

    body_settings = environment_setup.get_default_body_settings(
        ["Earth", "Sun", "Moon", "Mars"], "SSB", "J2000"
    )
    earth_settings = body_settings.get("Earth")
    earth_settings.shape_settings = shape.oblate_spherical_spice()
    earth_settings.rotation_model_settings = rotation_model.gcrs_to_itrs(base_frame="J2000")
    earth_settings.ground_station_settings = ground_station.radio_telescope_stations()
    earth_settings.shape_deformation_settings = [shape_deformation.iers_2010_solid_body_tidal()]

    body_settings.add_empty_settings("MeX")
    body_settings.get("MeX").ephemeris_settings = ephemeris.direct_spice("Mars", "J2000")

    bodies = environment_setup.create_system_of_bodies(body_settings)
    mex_systems = environment.VehicleSystems()
    mex_systems.set_default_transponder_turnaround_ratio_function()
    bodies.get_body("MeX").system_models = mex_systems
    return bodies


def _create_mro_bodies(mro_kernel_path: Path, initial_time: float, final_time: float):
    spice.load_standard_kernels()
    for kernel_file in (
        "mro_sc_psp_120313_120319.bc",
        "mro_hga_psp_120313_120319.bc",
        "mro_sclkscet_00112_65536.tsc",
        "mro_v16.tf",
        "mro_psp22.bsp",
        "mro_struct_v10.bsp",
    ):
        spice.load_kernel(str(mro_kernel_path / kernel_file))

    body_settings = environment_setup.get_default_body_settings(
        ["Earth", "Sun", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Moon"],
        "SSB",
        "J2000",
    )
    earth_settings = body_settings.get("Earth")
    earth_settings.shape_settings = shape.oblate_spherical_spice()
    earth_settings.rotation_model_settings = rotation_model.gcrs_to_itrs(base_frame="J2000")
    earth_settings.ground_station_settings = ground_station.dsn_stations()

    body_settings.add_empty_settings("MRO")
    body_settings.get("MRO").ephemeris_settings = ephemeris.interpolated_spice(
        initial_time,
        final_time,
        10.0,
        "Mars",
        "J2000",
    )
    body_settings.get("MRO").rotation_model_settings = rotation_model.spice(
        "J2000", "MRO_SPACECRAFT", ""
    )

    bodies = environment_setup.create_system_of_bodies(body_settings)
    mro_systems = environment.VehicleSystems()
    mro_systems.set_default_transponder_turnaround_ratio_function()
    bodies.get_body("MRO").system_models = mro_systems
    return bodies


def _download_grail_residual_files(test_data_path: Path):
    grail_data_path = test_data_path / "grail_residuals_download"
    kernel_path = grail_data_path / "kernel_download"
    antenna_path = grail_data_path / "antenna_download"

    for url in (
        "https://naif.jpl.nasa.gov/pub/naif/pds/data/grail-l-spice-6-v1.0/"
        "grlsp_1000/data/sclk/gra_sclkscet_00014.tsc",
        "https://naif.jpl.nasa.gov/pub/naif/pds/data/grail-l-spice-6-v1.0/"
        "grlsp_1000/data/ck/gra_rec_120409_120415.bc",
        "https://naif.jpl.nasa.gov/pub/naif/pds/data/grail-l-spice-6-v1.0/"
        "grlsp_1000/data/spk/grail_120301_120529_sci_v02.bsp",
        "https://naif.jpl.nasa.gov/pub/naif/pds/data/grail-l-spice-6-v1.0/"
        "grlsp_1000/data/fk/grail_v07.tf",
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de440_200625.bpc",
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_de440_250416.tf",
    ):
        _download_file(url, kernel_path)

    antenna_file = _download_file(
        "https://pds-geosciences.wustl.edu/grail/grail-l-lgrs-3-cdr-v1/"
        "grail_0101/level_1b/2012_04_10/vgs1b_2012_04_10_a_04.asc",
        antenna_path,
    )

    return kernel_path, [antenna_file]


def _read_grail_antenna_switch_history(antenna_files, start_time: float, end_time: float):
    switch_times = []
    switch_positions = []
    for antenna_file in antenna_files:
        is_data_block = False
        with Path(antenna_file).open() as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue
                if is_data_block:
                    fields = line.split()
                    if len(fields) != 7:
                        raise RuntimeError(f"Unexpected GRAIL antenna switch line: {line}")
                    antenna_norm = float(fields[2])
                    antenna_direction = np.array(
                        [float(fields[3]), float(fields[4]), float(fields[5])]
                    )
                    switch_times.append(float(fields[0]))
                    switch_positions.append(antenna_norm * antenna_direction)
                elif line == "END OF HEADER":
                    is_data_block = True

    if not switch_positions:
        raise RuntimeError("No GRAIL antenna switch data found.")

    antenna_switch_history = {
        start_time: switch_positions[0],
        end_time: switch_positions[-1],
    }
    antenna_switch_history.update(
        {
            switch_time: switch_position
            for switch_time, switch_position in zip(switch_times, switch_positions)
        }
    )
    return antenna_switch_history


def _create_grail_bodies(test_data_path: Path, initial_time: float, final_time: float):
    kernel_path, _ = _download_grail_residual_files(test_data_path)
    spice.load_standard_kernels()
    for kernel_file in (
        "moon_de440_250416.tf",
        "moon_pa_de440_200625.bpc",
        "grail_v07.tf",
        "gra_rec_120409_120415.bc",
        "gra_sclkscet_00014.tsc",
        "grail_120301_120529_sci_v02.bsp",
    ):
        spice.load_kernel(str(kernel_path / kernel_file))

    body_settings = environment_setup.get_default_body_settings_time_limited(
        ["Earth", "Sun", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Moon"],
        initial_time,
        final_time,
        "SSB",
        "J2000",
    )
    earth_settings = body_settings.get("Earth")
    earth_settings.shape_settings = shape.oblate_spherical_spice()
    earth_settings.rotation_model_settings = rotation_model.gcrs_to_itrs(base_frame="J2000")
    earth_settings.ground_station_settings = ground_station.dsn_stations()
    body_settings.get("Moon").rotation_model_settings = rotation_model.spice(
        "J2000", "MOON_PA_DE440", "MOON_PA_DE440"
    )

    body_settings.add_empty_settings("GRAIL-A")
    body_settings.get("GRAIL-A").ephemeris_settings = ephemeris.interpolated_spice(
        initial_time,
        final_time,
        10.0,
        "Moon",
        "J2000",
    )
    body_settings.get("GRAIL-A").rotation_model_settings = rotation_model.spice(
        "J2000", "GRAIL-A_SPACECRAFT", ""
    )

    bodies = environment_setup.create_system_of_bodies(body_settings)
    grail_systems = environment.VehicleSystems()
    grail_systems.set_default_transponder_turnaround_ratio_function()
    bodies.get_body("GRAIL-A").system_models = grail_systems
    return bodies


def _create_juice_bodies(test_data_path: Path):
    spice.load_standard_kernels()
    spice.load_kernel(str(test_data_path / "juice_orbc_000074_230414_310721_v01.bsp"))

    body_settings = environment_setup.get_default_body_settings(
        ["Earth", "Moon", "Sun", "Jupiter"], "SSB", "J2000"
    )
    earth_settings = body_settings.get("Earth")
    earth_settings.shape_settings = shape.oblate_spherical_spice()
    earth_settings.rotation_model_settings = rotation_model.gcrs_to_itrs(base_frame="J2000")
    earth_settings.shape_deformation_settings = [shape_deformation.iers_2010_solid_body_tidal()]

    station_positions = ground_station.get_radio_telescope_positions()
    earth_settings.ground_station_settings = [
        ground_station.basic_station(
            "NWNORCIA",
            station_positions["NWNORCIA"],
            station_motion_settings=[
                ground_station.linear_station_motion(
                    np.array([-45.00, 10.00, 47.00]) / 1.0e3 / constants.JULIAN_YEAR,
                    0.0,
                )
            ],
        ),
        ground_station.basic_station(
            "YARRAGAD",
            station_positions["YARRAGAD"],
            station_motion_settings=[
                ground_station.linear_station_motion(
                    np.array([-47.45, 9.12, 51.76]) / 1.0e3 / constants.JULIAN_YEAR,
                    0.0,
                )
            ],
        ),
    ]

    body_settings.add_empty_settings("JUICE")
    body_settings.get("JUICE").ephemeris_settings = ephemeris.direct_spice("Earth", "J2000")

    bodies = environment_setup.create_system_of_bodies(body_settings)
    juice_systems = environment.VehicleSystems()
    juice_systems.set_default_transponder_turnaround_ratio_function()
    bodies.get_body("JUICE").system_models = juice_systems
    bodies.get_body("Earth").get_ground_station("NWNORCIA").set_transmitting_frequency_calculator(
        environment.ConstantTransmittingFrequencyCalculator(7180127320.0)
    )
    return bodies


def _create_voyager_triton_psf_bodies(test_data_path: Path, voyager_state, supplementary_data):
    body_settings = environment_setup.BodyListSettings("8", "J2000")
    for body_name in ("8", "VGR2", "TRITON"):
        body_settings.add_empty_settings(body_name)

    body_settings.get("8").ephemeris_settings = ephemeris.direct_spice("SSB", "J2000")
    body_settings.get("VGR2").ephemeris_settings = ephemeris.constant(voyager_state, "8", "J2000")
    body_settings.get("TRITON").ephemeris_settings = ephemeris.direct_spice("8", "J2000", "801")
    body_settings.get("VGR2").rotation_model_settings = rotation_model.constant_rotation_model(
        "J2000", "VGR2_Fixed", np.eye(3)
    )

    bodies = environment_setup.create_system_of_bodies(body_settings)
    set_tracking_supplementary_data_in_bodies(bodies, supplementary_data)
    return bodies


def _keep_first_observations(tracking_data, number_of_observations: int):
    for index in range(len(tracking_data.observations) - 1, number_of_observations - 1, -1):
        tracking_data.remove_single_observation_entry(index)


def _keep_observations_in_time_window(tracking_data, start_time: float, end_time: float):
    filtered_tracking_data = []
    for current_tracking_data in tracking_data:
        epochs = [float(epoch) for epoch in current_tracking_data.epochs]
        first_index_to_keep = bisect_left(epochs, start_time)
        last_index_to_keep = bisect_right(epochs, end_time)
        if first_index_to_keep == last_index_to_keep:
            continue

        for index in range(len(epochs) - 1, last_index_to_keep - 1, -1):
            current_tracking_data.remove_single_observation_entry(index)
        for index in range(first_index_to_keep - 1, -1, -1):
            current_tracking_data.remove_single_observation_entry(index)
        filtered_tracking_data.append(current_tracking_data)

    return filtered_tracking_data


def _set_doppler_integration_time(observation_collection, integration_time: float):
    for observation_sets_per_link in observation_collection.sorted_observation_sets.values():
        for observation_sets in observation_sets_per_link.values():
            for observation_set in observation_sets:
                observation_set.ancillary_settings.set_float_settings(
                    ancillary_settings.doppler_integration_time,
                    integration_time,
                )


def _simulate_dsn_n_way_averaged_doppler(observation_collection, bodies):
    light_time_settings = [
        light_time_corrections.first_order_relativistic_light_time_correction(["Sun"])
    ]
    observation_model_settings = []
    for raw_link_ends_list in observation_collection.link_ends_per_observable_type.values():
        for raw_link_ends in raw_link_ends_list:
            observation_model_settings.append(
                model_settings.dsn_n_way_doppler_averaged(
                    links.link_definition(raw_link_ends),
                    light_time_settings,
                    biases.absolute_bias(np.zeros(1)),
                    light_time_corrections.light_time_convergence_settings(True),
                    False,
                )
            )

    observation_simulators = observations_simulation_settings.create_observation_simulators(
        observation_model_settings, bodies
    )
    observation_simulation_settings = (
        observations_simulation_settings.observation_settings_from_collection(
            observation_collection, bodies
        )
    )
    return simulate_observations(observation_simulation_settings, observation_simulators, bodies)


def _simulate_dsn_n_way_averaged_doppler_with_corrections(
    observation_collection, bodies, light_time_settings
):
    observation_model_settings = []
    for raw_link_ends_list in observation_collection.link_ends_per_observable_type.values():
        for raw_link_ends in raw_link_ends_list:
            observation_model_settings.append(
                model_settings.dsn_n_way_doppler_averaged(
                    links.link_definition(raw_link_ends),
                    light_time_settings,
                    biases.absolute_bias(np.zeros(1)),
                    light_time_corrections.light_time_convergence_settings(True),
                )
            )

    observation_simulators = observations_simulation_settings.create_observation_simulators(
        observation_model_settings, bodies
    )
    observation_simulation_settings = (
        observations_simulation_settings.observation_settings_from_collection(
            observation_collection, bodies
        )
    )
    return simulate_observations(observation_simulation_settings, observation_simulators, bodies)


def _simulate_doppler_measured_frequency(observation_collection, bodies):
    light_time_settings = [
        light_time_corrections.first_order_relativistic_light_time_correction(
            ["Sun", "Moon", "Earth"]
        )
    ]
    observation_model_settings = []
    for raw_link_ends_list in observation_collection.link_ends_per_observable_type.values():
        for raw_link_ends in raw_link_ends_list:
            observation_model_settings.append(
                model_settings.doppler_measured_frequency(
                    links.link_definition(raw_link_ends),
                    light_time_settings,
                )
            )

    observation_simulators = observations_simulation_settings.create_observation_simulators(
        observation_model_settings, bodies
    )
    observation_simulation_settings = (
        observations_simulation_settings.observation_settings_from_collection(
            observation_collection, bodies
        )
    )
    return simulate_observations(observation_simulation_settings, observation_simulators, bodies)


def _simulate_radar_residuals(observation_collection, bodies):
    """Compute radar residuals using the same observable families loaded from the data."""
    light_time_settings = [
        light_time_corrections.first_order_relativistic_light_time_correction(["Sun"])
    ]
    observation_model_settings = []
    for (
        observable_type,
        raw_link_ends_list,
    ) in observation_collection.link_ends_per_observable_type.items():
        for raw_link_ends in raw_link_ends_list:
            link_definition = links.link_definition(raw_link_ends)
            if observable_type == model_settings.n_way_range_type:
                observation_model_settings.append(
                    model_settings.n_way_range(
                        link_definition,
                        light_time_settings,
                        bias_settings=None,
                        time_scale_for_observable=time_representation.utc_scale,
                    )
                )
            elif observable_type == model_settings.doppler_measured_frequency_type:
                observation_model_settings.append(
                    model_settings.doppler_measured_frequency(
                        link_definition,
                        light_time_settings,
                        bias_settings=None,
                    )
                )

    observation_simulators = observations_simulation_settings.create_observation_simulators(
        observation_model_settings, bodies
    )
    compute_residuals_and_dependent_variables(
        observation_collection,
        observation_simulators,
        bodies,
    )
    return np.asarray(observation_collection.get_concatenated_residuals())


def _simulate_angular_position_residuals(observation_collection, bodies):
    """Compute optical angular-position residuals for already-created observations."""
    observation_model_settings = []
    for raw_link_ends_list in observation_collection.link_ends_per_observable_type.values():
        for raw_link_ends in raw_link_ends_list:
            observation_model_settings.append(
                model_settings.angular_position(
                    links.link_definition(raw_link_ends),
                    bias_settings=None,
                )
            )

    observation_simulators = observations_simulation_settings.create_observation_simulators(
        observation_model_settings, bodies
    )
    compute_residuals_and_dependent_variables(
        observation_collection,
        observation_simulators,
        bodies,
    )
    return np.asarray(observation_collection.get_concatenated_residuals())


def _skip_remote_data_error(error: Exception, context: str):
    message = str(error)
    remote_error_markers = [
        "Network is unreachable",
        "Failed to establish a new connection",
        "Max retries exceeded",
        "timed out",
        "Temporary failure",
        "Name or service not known",
        "Connection refused",
        "Service Unavailable",
        "HTTP Error",
    ]
    if any(marker in message for marker in remote_error_markers):
        pytest.skip(f"{context} unavailable: {message}")
    raise error


def _utc_seconds_to_tdb(epochs):
    time_scale_converter = time_representation.default_time_scale_converter()
    return np.array(
        [
            time_scale_converter.convert_time(
                input_scale=time_representation.utc_scale,
                output_scale=time_representation.tdb_scale,
                input_value=float(epoch),
            )
            for epoch in epochs
        ],
        dtype=float,
    )


def _create_itokawa_radar_bodies(radar_table):
    """Create the environment needed to reproduce the selected Itokawa radar arc."""
    spice.load_standard_kernels()
    start_time = float(radar_table["epoch_seconds_UTC"].min()) - 10.0 * constants.JULIAN_DAY
    end_time = float(radar_table["epoch_seconds_UTC"].max()) + 10.0 * constants.JULIAN_DAY

    body_settings = environment_setup.get_default_body_settings(
        ["Sun", "Earth", "Moon"], "SSB", "J2000"
    )
    earth_settings = body_settings.get("Earth")
    earth_settings.shape_settings = shape.oblate_spherical(
        6378137.0,
        1.0 / 298.257223563,
    )
    earth_settings.rotation_model_settings = rotation_model.gcrs_to_itrs(
        rotation_model.iau_2006,
        "J2000",
    )
    earth_settings.gravity_field_settings.associated_reference_frame = "ITRS"

    body_settings.add_empty_settings("101955")
    horizons_states = HorizonsQuery(
        query_id="101955;",
        location="500@10",
        epoch_start=start_time,
        epoch_end=end_time,
        epoch_step="12h",
        extended_query=True,
    ).cartesian(frame_orientation="J2000")
    body_settings.get("101955").ephemeris_settings = ephemeris.tabulated(
        {float(row[0]): np.asarray(row[1:], dtype=float) for row in horizons_states},
        frame_origin="Sun",
        frame_orientation="J2000",
    )

    bodies = environment_setup.create_system_of_bodies(body_settings)
    target_systems = environment.VehicleSystems()
    frequency_bands = [
        ancillary_settings.FrequencyBands.s_band,
        ancillary_settings.FrequencyBands.x_band,
        ancillary_settings.FrequencyBands.ka_band,
        ancillary_settings.FrequencyBands.ku_band,
    ]
    target_systems.set_transponder_turnaround_ratio({(band, band): 1.0 for band in frequency_bands})
    bodies.get_body("101955").system_models = target_systems
    add_all_radar_ground_stations(bodies)
    return bodies


def _itokawa_2005_jpl_radar_arc():
    """Return a compact JPL Itokawa radar arc with both range and Doppler data."""
    return JPLRadarQuery("101955", timeout=60.0).to_radar_data(
        target_body="101955",
        epoch_start=datetime.datetime(2005, 1, 1),
        epoch_end=datetime.datetime(2006, 1, 1),
        target_point="C",
    )


def _radar_sigmas_in_observation_order(observation_collection, radar_table):
    """Match source radar rows to Tudat's concatenated scalar observation order.

    The conversion groups radar rows by link definition and observable type, so
    the concatenated observation vector is not identical to the API row order.
    Matching by scalar value keeps the residual/sigma comparison source-based.
    """
    available_rows = radar_table.reset_index(drop=True).copy()
    sigmas = []
    observable_types = []

    for observed_value in np.asarray(observation_collection.concatenated_observations, dtype=float):
        value_differences = np.abs(available_rows["value"].to_numpy(dtype=float) - observed_value)
        row_index = int(np.argmin(value_differences))
        value_tolerance = max(1.0e-6, abs(observed_value) * 1.0e-12)
        assert value_differences[row_index] < value_tolerance

        matched_row = available_rows.iloc[row_index]
        sigmas.append(float(matched_row["sigma"]))
        observable_types.append(str(matched_row["observable_type"]))
        available_rows = available_rows.drop(available_rows.index[row_index]).reset_index(drop=True)

    return np.asarray(sigmas, dtype=float), np.asarray(observable_types, dtype=str)


def _horizons_cartesian_states(horizons_id, epochs):
    unique_epochs = np.array(sorted(set(float(epoch) for epoch in epochs)))
    return HorizonsQuery(
        query_id=horizons_id,
        location="500@SSB",
        epoch_list=list(unique_epochs),
        extended_query=True,
    ).cartesian(frame_orientation="J2000")


def _pad_state_history_for_lagrange(state_history, minimum_points=6):
    if len(state_history) >= minimum_points:
        return state_history

    epochs = np.array(sorted(state_history), dtype=float)
    states = np.array([state_history[epoch] for epoch in epochs], dtype=float)
    nominal_step = 600.0 if len(epochs) == 1 else np.median(np.diff(epochs))
    if nominal_step <= 0.0:
        nominal_step = 600.0

    center_epoch = np.mean(epochs)
    padded_epochs = (
        center_epoch + (np.arange(minimum_points) - (minimum_points - 1) / 2.0) * nominal_step
    )
    padded_epochs = np.unique(np.concatenate((padded_epochs, epochs)))

    padded_states = {}
    for epoch in padded_epochs:
        if len(epochs) == 1:
            reference_epoch = epochs[0]
            reference_state = states[0]
        elif epoch <= epochs[0]:
            reference_epoch = epochs[0]
            reference_state = states[0]
        elif epoch >= epochs[-1]:
            reference_epoch = epochs[-1]
            reference_state = states[-1]
        else:
            right_index = np.searchsorted(epochs, epoch)
            left_index = right_index - 1
            fraction = (epoch - epochs[left_index]) / (epochs[right_index] - epochs[left_index])
            reference_state = (1.0 - fraction) * states[left_index] + fraction * states[right_index]
            reference_epoch = epoch

        padded_states[float(epoch)] = np.hstack(
            (
                reference_state[:3] + reference_state[3:] * (epoch - reference_epoch),
                reference_state[3:],
            )
        )

    for epoch, state in state_history.items():
        padded_states[float(epoch)] = state

    return dict(sorted(padded_states.items()))


def _local_hst_target_ephemeris_settings(batch, horizons_id):
    table = batch.table.sort_values("epoch_seconds_UTC", kind="stable")
    receive_epochs = _utc_seconds_to_tdb(table["epoch_seconds_UTC"].to_numpy(dtype=float))
    spacecraft_positions = table.loc[
        :,
        ["spacecraft_position_x", "spacecraft_position_y", "spacecraft_position_z"],
    ].to_numpy(dtype=float)
    earth_positions = np.array(
        [
            spice.get_body_cartesian_state_at_epoch(
                "Earth",
                "SSB",
                "J2000",
                "NONE",
                float(epoch),
            )[:3]
            for epoch in receive_epochs
        ]
    )

    target_states_at_receive = _horizons_cartesian_states(horizons_id, receive_epochs)[
        :,
        1:7,
    ]
    receiver_positions = earth_positions + spacecraft_positions
    light_times = (
        np.linalg.norm(target_states_at_receive[:, :3] - receiver_positions, axis=1)
        / constants.SPEED_OF_LIGHT
    )
    transmit_epochs = receive_epochs - light_times

    sample_epochs = []
    for epoch_array in (receive_epochs, transmit_epochs):
        for offset in (-600.0, -300.0, 0.0, 300.0, 600.0):
            sample_epochs.extend(epoch_array + offset)

    target_states = _horizons_cartesian_states(horizons_id, sample_epochs)
    state_history = {float(row[0]): row[1:7] for row in target_states}
    return ephemeris.tabulated(
        _pad_state_history_for_lagrange(state_history),
        frame_origin="SSB",
        frame_orientation="J2000",
    )


def _hubble_space_astrometry_batch(target, observatory):
    batch = BatchMPC()
    batch.get_observations([target], drop_misc_observations=False)
    batch.filter(
        epoch_start=datetime.datetime(1990, 1, 1),
        epoch_end=datetime.datetime(2026, 1, 1),
        observatories=[observatory],
    )
    return batch


def _create_hubble_space_astrometry_bodies(batch, horizons_id, observatory):
    spice.load_standard_kernels()
    target_body = batch.MPC_objects[0]
    body_settings = environment_setup.get_default_body_settings(
        ["Sun", "Earth"],
        "SSB",
        "J2000",
    )
    body_settings.add_empty_settings(target_body)
    body_settings.get(target_body).ephemeris_settings = _local_hst_target_ephemeris_settings(
        batch,
        horizons_id,
    )
    body_settings.add_empty_settings(str(observatory))
    return environment_setup.create_system_of_bodies(body_settings)


def _estimate_voyager_velocity_from_jacobson_references(references, reference_index: int):
    lower_index = reference_index
    upper_index = reference_index
    if reference_index == 0:
        upper_index = 1
    elif reference_index == len(references) - 1:
        lower_index = reference_index - 1
    else:
        lower_index = reference_index - 1
        upper_index = reference_index + 1

    lower_time = spice.convert_julian_date_to_ephemeris_time(references[lower_index][1])
    upper_time = spice.convert_julian_date_to_ephemeris_time(references[upper_index][1])
    return (
        (references[upper_index][2] - references[lower_index][2])
        * 1000.0
        / (upper_time - lower_time)
    )


def _create_jacobson_voyager_state(references, reference_index: int):
    reception_time = spice.convert_julian_date_to_ephemeris_time(references[reference_index][1])
    rotation_from_b1950_to_j2000 = spice.compute_rotation_matrix_between_frames(
        "B1950", "J2000", reception_time
    )
    neptune_state = spice.get_body_cartesian_state_at_epoch(
        "8", "SSB", "J2000", "NONE", reception_time
    )

    voyager_state = np.zeros(6)
    voyager_state[:3] = rotation_from_b1950_to_j2000 @ (1000.0 * references[reference_index][2])
    voyager_state[3:] = np.asarray(neptune_state[3:]) + (
        rotation_from_b1950_to_j2000
        @ _estimate_voyager_velocity_from_jacobson_references(references, reference_index)
    )
    return voyager_state


def _simulate_pixel_coordinates(observation_collection, bodies):
    observation_model_settings = []
    for raw_link_ends_list in observation_collection.link_ends_per_observable_type.values():
        for raw_link_ends in raw_link_ends_list:
            observation_model_settings.append(
                model_settings.pixel_coordinates(
                    links.link_definition(raw_link_ends),
                    [],
                    None,
                    light_time_corrections.light_time_convergence_settings(True),
                    True,
                )
            )

    observation_simulators = observations_simulation_settings.create_observation_simulators(
        observation_model_settings, bodies
    )
    observation_simulation_settings = (
        observations_simulation_settings.observation_settings_from_collection(
            observation_collection, bodies
        )
    )
    return simulate_observations(observation_simulation_settings, observation_simulators, bodies)


def test_ifms_mex_residuals_are_millihertz_level():
    test_data_path = _test_data_path()
    bodies = _create_mex_bodies(test_data_path)

    tracking_data, supplementary_data = read_ifms_data(
        [
            str(
                test_data_path
                / "estrack_n_way_doppler_observation_model"
                / "M32ICL3L02_D2S_133621904_00.TAB.txt"
            )
        ],
        "MeX",
        ["NWNORCIA"],
        "Earth",
        True,
        True,
        [1.0, 0.0],
        0.0,
        0.0,
    )
    set_tracking_supplementary_data_in_bodies(bodies, supplementary_data)
    uncompressed_observations = create_observation_collection_from_tracking_data(
        tracking_data, bodies
    )
    observed_observations = create_compressed_doppler_collection(
        uncompressed_observations,
        60,
        earth_fixed_ground_station_positions=ground_station.get_radio_telescope_positions(),
    )

    computed_observations = _simulate_dsn_n_way_averaged_doppler(observed_observations, bodies)
    residuals = np.asarray(observed_observations.concatenated_observations) - np.asarray(
        computed_observations.concatenated_observations
    )

    assert residuals.size == 321
    assert np.sqrt(np.mean(residuals**2)) < 3.5e-3
    assert abs(np.mean(residuals)) < 1.0e-3


def test_odf_grail_short_arc_residuals_are_millihertz_level():
    test_data_path = _test_data_path()
    grail_data_path = test_data_path / "grail_residuals_download"
    _, antenna_files = _download_grail_residual_files(test_data_path)
    odf_file = _download_file(
        "https://pds-geosciences.wustl.edu/grail/grail-l-rss-2-edr-v1/"
        "grail_0201/odf/gralugf2012_101_0235smmmv1.odf",
        grail_data_path,
    )
    troposphere_file = _download_file(
        "https://pds-geosciences.wustl.edu/grail/grail-l-rss-2-edr-v1/"
        "grail_0201/ancillary/tro/grxlugf2012_092_2012_122.tro",
        grail_data_path,
    )
    ionosphere_file = _download_file(
        "https://pds-geosciences.wustl.edu/grail/grail-l-rss-2-edr-v1/"
        "grail_0201/ancillary/ion/gralugf2012_092_2012_122.ion",
        grail_data_path,
    )

    interval_start = 387314833.0
    interval_end = 387318433.0
    bodies = _create_grail_bodies(test_data_path, interval_start - 3600.0, interval_end + 3600.0)
    tracking_data, supplementary_data = read_odf_data([str(odf_file)], "GRAIL-A", "Earth", False)
    set_tracking_supplementary_data_in_bodies(bodies, supplementary_data)
    doppler_tracking_data = [
        data for data in tracking_data if data.observable_type == "DsnNWayAveragedDoppler"
    ]
    doppler_tracking_data = _keep_observations_in_time_window(
        doppler_tracking_data, interval_start, interval_end
    )

    uncompressed_observations = create_observation_collection_from_tracking_data(
        doppler_tracking_data, bodies
    )
    observed_observations = create_compressed_doppler_collection(uncompressed_observations, 60, 10)
    observed_observations.set_reference_points(
        bodies,
        _read_grail_antenna_switch_history(
            antenna_files, interval_start - 3600.0, interval_end + 3600.0
        ),
        "GRAIL-A",
        links.reflector1,
    )
    light_time_settings = [
        light_time_corrections.first_order_relativistic_light_time_correction(["Sun"]),
        light_time_corrections.dsn_tabulated_tropospheric_light_time_correction(
            [str(troposphere_file)]
        ),
        light_time_corrections.dsn_tabulated_ionospheric_light_time_correction(
            [str(ionosphere_file)], {177: "GRAIL-A"}
        ),
    ]
    computed_observations = _simulate_dsn_n_way_averaged_doppler_with_corrections(
        observed_observations, bodies, light_time_settings
    )
    residuals = np.asarray(observed_observations.concatenated_observations) - np.asarray(
        computed_observations.concatenated_observations
    )

    assert residuals.size == 49
    assert abs(np.mean(residuals)) < 3.0e-3
    assert np.sqrt(np.mean(residuals**2)) < 3.0e-3


def test_tnf_mro_short_arc_residuals_are_low_after_compression():
    test_data_path = _test_data_path()
    mro_data_path = test_data_path / "mro_dsn_observation_model"
    mro_kernel_path = mro_data_path / "kernel_download"
    tnf_file = _download_file(
        "https://pds-geosciences.wustl.edu/mro/mro-m-rss-1-magr-v1/"
        "mrors_0xxx/tnf/mromagr2012_076_0840xmmmv1.tnf",
        mro_data_path / "tnf_download",
    )
    for url in (
        "https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/"
        "mrosp_1000/data/ck/mro_sc_psp_120313_120319.bc",
        "https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/"
        "mrosp_1000/data/ck/mro_hga_psp_120313_120319.bc",
        "https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/"
        "mrosp_1000/data/sclk/mro_sclkscet_00112_65536.tsc",
        "https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/"
        "mrosp_1000/data/fk/mro_v16.tf",
        "https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/"
        "mrosp_1000/data/spk/mro_psp22.bsp",
        "https://naif.jpl.nasa.gov/pub/naif/pds/data/mro-m-spice-6-v1.0/"
        "mrosp_1000/data/spk/mro_struct_v10.bsp",
    ):
        _download_file(url, mro_kernel_path)

    interval_start = 385256700.0
    interval_end = 385260300.0
    bodies = _create_mro_bodies(mro_kernel_path, interval_start - 3600.0, interval_end + 3600.0)
    tracking_data, supplementary_data = read_tnf_data([str(tnf_file)], ["doppler"], "MRO")
    tracking_data = _keep_observations_in_time_window(tracking_data, interval_start, interval_end)
    set_tracking_supplementary_data_in_bodies(bodies, supplementary_data)
    uncompressed_observations = create_observation_collection_from_tracking_data(
        tracking_data, bodies
    )
    observed_observations = create_compressed_doppler_collection(uncompressed_observations, 60, 10)
    observed_observations.set_transponder_delay("MRO", 1.4149e-6)

    mro_center_of_mass_position = np.array([0.0, -1.11, 0.0])
    antenna_position_history = {}
    for observation_times in observed_observations.get_observation_times_objects():
        current_time = observation_times[0].to_float() - 3600.0
        while current_time <= observation_times[-1].to_float() + 3600.0:
            antenna_state = np.zeros(6)
            antenna_state[:3] = (
                spice.get_body_cartesian_position_at_epoch(
                    "-74214",
                    "-74000",
                    "MRO_SPACECRAFT",
                    "none",
                    current_time,
                )
                - mro_center_of_mass_position
            )
            antenna_position_history[current_time] = antenna_state
            current_time += 60.0

    antenna_ephemeris_settings = ephemeris.tabulated(
        antenna_position_history,
        "-74000",
        "MRO_SPACECRAFT",
    )
    antenna_ephemeris = ephemeris.create_ephemeris(antenna_ephemeris_settings, "Antenna")
    observed_observations.set_reference_point(
        bodies,
        antenna_ephemeris,
        "Antenna",
        "MRO",
        links.reflector1,
    )

    light_time_settings = [
        light_time_corrections.approximated_second_order_relativistic_light_time_correction(
            ["Sun"]
        ),
        light_time_corrections.dsn_tabulated_tropospheric_light_time_correction(
            [str(mro_data_path / "mromagr2012_061_2012_092.tro")]
        ),
        light_time_corrections.dsn_tabulated_ionospheric_light_time_correction(
            [str(mro_data_path / "mromagr2012_061_2012_092.ion")], {74: "MRO"}
        ),
    ]
    computed_observations = _simulate_dsn_n_way_averaged_doppler_with_corrections(
        observed_observations, bodies, light_time_settings
    )
    residuals = np.asarray(computed_observations.concatenated_observations) - np.asarray(
        observed_observations.concatenated_observations
    )

    assert residuals.size == 59
    assert abs(np.mean(residuals)) < 1.5e-3
    assert np.sqrt(np.mean(residuals**2)) < 3.0e-3


def test_psf_voyager_triton_pixel_line_residuals_are_subpixel():
    test_data_path = _test_data_path()
    psf_file = test_data_path / "psf" / "psf_vgr2_neptune.txt"
    spice.load_standard_kernels([])
    spice.load_kernel(str(test_data_path / "spice" / "vgr2_nep097.bsp"))

    references = [
        ("2890B+34", 2447480.34066378, np.array([-178334616.8, 346108595.1, 128939829.6])),
        ("4994B+58", 2447550.48720449, np.array([-134196154.2, 260416600.6, 97025687.3])),
        ("5912B+21", 2447581.06662806, np.array([-114954294.2, 223070441.2, 83117675.7])),
        ("7773B+53", 2447643.11764157, np.array([-75908145.9, 147301365.2, 54901517.8])),
        ("8585B+59", 2447670.18760590, np.array([-58873770.6, 114250315.9, 42592968.9])),
        ("9180B+55", 2447690.01869003, np.array([-46394034.5, 90037807.3, 33576066.8])),
    ]

    tracking_data, supplementary_data = read_psf_data(
        [str(psf_file)], "VGR2", {"TRITON": "TRITON"}, False
    )
    all_observations = create_observation_collection_from_tracking_data(
        tracking_data,
        _create_voyager_triton_psf_bodies(test_data_path, np.zeros(6), supplementary_data),
    )
    all_observation_times = np.asarray(all_observations.concatenated_times, dtype=float)

    residuals = []
    for reference_index, reference in enumerate(references):
        reception_time = spice.convert_julian_date_to_ephemeris_time(reference[1])
        component_index = int(np.argmin(np.abs(all_observation_times - reception_time)))
        observation_index = component_index // 2
        assert abs(all_observation_times[component_index] - reception_time) < 1.0e-3

        bodies = _create_voyager_triton_psf_bodies(
            test_data_path,
            _create_jacobson_voyager_state(references, reference_index),
            supplementary_data,
        )
        single_tracking_data = TrackingData(
            tracking_data[0].observable_type,
            tracking_data[0].link_ends,
            [tracking_data[0].observations[observation_index]],
            [tracking_data[0].epochs[observation_index]],
            tracking_data[0].reference_link_end,
            tracking_data[0].time_scale,
        )
        observed_observations = create_observation_collection_from_tracking_data(
            [single_tracking_data], bodies
        )
        computed_observations = _simulate_pixel_coordinates(observed_observations, bodies)
        residuals.append(
            np.asarray(computed_observations.concatenated_observations)
            - np.asarray(observed_observations.concatenated_observations)
        )

    residuals = np.vstack(residuals)
    assert residuals.shape == (6, 2)
    assert np.max(np.linalg.norm(residuals, axis=1)) < 1.0


def test_fdets_juice_short_arc_residual_scatter_is_millihertz_level():
    test_data_path = _test_data_path()
    bodies = _create_juice_bodies(test_data_path)

    tracking_data, supplementary_data = read_fdets_data(
        [str(test_data_path / "Fdets.jui2024.08.20.Yg.r2i.txt")],
        [8422.49e6],
        FdetDateFormat.datetime_string,
        "JUICE",
        ["NWNORCIA"],
        ["YARRAGAD"],
        "Earth",
    )
    assert supplementary_data == []
    _keep_first_observations(tracking_data[0], 120)

    # The FDETS reader currently stores the measured frequency and base frequency.
    # The observation model also needs the link frequency bands.
    tracking_data[0].add_double_vector_ancillary_setting("frequency bands", [1.0, 1.0])
    observed_observations = create_observation_collection_from_tracking_data(tracking_data, bodies)
    computed_observations = _simulate_doppler_measured_frequency(observed_observations, bodies)
    residuals = np.asarray(computed_observations.concatenated_observations) - np.asarray(
        observed_observations.concatenated_observations
    )
    residual_scatter = residuals - np.mean(residuals)

    assert residuals.size == 120
    assert np.sqrt(np.mean(residual_scatter**2)) < 1.0e-2
    assert np.max(np.abs(residual_scatter)) < 2.5e-2


def test_translational_state_supplementary_data_rejects_non_utc_tdb_epochs():
    """Reject supplementary state epochs whose conversion is not explicitly supported."""
    body_settings = environment_setup.BodyListSettings("SSB", "J2000")
    body_settings.add_empty_settings("TestSpacecraft")
    bodies = environment_setup.create_system_of_bodies(body_settings)

    # This supplementary path currently implements the two cases used by the
    # data readers: ephemeris-ready TDB epochs and UTC-tagged source epochs.
    translational_data = TranslationalStateSupplementaryData(
        {0.0: np.zeros(6)},
        "SSB",
        True,
        "TT",
        "J2000",
    )
    supplementary_data = TrackingSupplementaryData("TestSpacecraft", "")
    supplementary_data.translational_state_supplementary_data = translational_data

    with pytest.raises(RuntimeError, match="only TDB and UTC time scales"):
        set_tracking_supplementary_data_in_bodies(bodies, [supplementary_data])


@pytest.mark.remote_data
def test_jpl_itokawa_radar_residuals_are_nonzero_and_bounded():
    """Check JPL Itokawa range and Doppler residuals over a multi-point arc."""
    try:
        radar_table = _itokawa_2005_jpl_radar_arc()
    except Exception as error:
        _skip_remote_data_error(error, "JPL radar data")

    # The selected 2005 arc is intentionally small enough for a remote Horizons
    # ephemeris query, but large enough to cover both radar observable families.
    observation_counts = radar_table["observable_type"].value_counts()
    assert len(radar_table) >= 10
    assert observation_counts[RANGE_OBSERVABLE] >= 10
    assert observation_counts[DOPPLER_OBSERVABLE] >= 3

    tracking_data, supplementary_data = radar_data_to_tracking_data(radar_table)
    try:
        bodies = _create_itokawa_radar_bodies(radar_table)
    except Exception as error:
        _skip_remote_data_error(error, "Horizons Itokawa ephemeris")

    set_tracking_supplementary_data_in_bodies(bodies, supplementary_data)
    observed_observations = create_observation_collection_from_tracking_data(
        tracking_data,
        bodies,
    )

    residuals = _simulate_radar_residuals(observed_observations, bodies)
    sigmas, observable_types = _radar_sigmas_in_observation_order(
        observed_observations,
        radar_table,
    )
    # JPL sigmas are the scale of interest here: the test should fail if the
    # conversion path creates residuals that are large relative to source data.
    normalized_residuals = np.abs(residuals) / sigmas
    range_residuals = residuals[observable_types == RANGE_OBSERVABLE]
    doppler_residuals = residuals[observable_types == DOPPLER_OBSERVABLE]

    assert residuals.size == len(radar_table)
    assert np.all(np.isfinite(residuals))
    assert np.any(np.abs(range_residuals) > 1.0e-6)
    assert np.any(np.abs(doppler_residuals) > 1.0e-9)
    assert np.sqrt(np.mean(normalized_residuals**2)) < 3.0
    assert np.max(normalized_residuals) < 3.5


@pytest.mark.remote_data
def test_hubble_mpc_space_astrometry_residuals_are_sub_arcsecond():
    """Check MPC HST astrometry with spacecraft receiver-state supplementary data."""
    try:
        batch = _hubble_space_astrometry_batch("2003 BF91", "250")
    except Exception as error:
        _skip_remote_data_error(error, "MPC Hubble astrometry")

    assert len(batch.table) > 0
    assert {
        "spacecraft_position_x",
        "spacecraft_position_y",
        "spacecraft_position_z",
    }.issubset(batch.table.columns)
    assert (
        batch.table[["spacecraft_position_x", "spacecraft_position_y", "spacecraft_position_z"]]
        .notna()
        .all(axis=None)
    )

    tracking_data, supplementary_data = batch.to_tracking_dataset(
        add_star_catalog_corrections=False
    )
    try:
        bodies = _create_hubble_space_astrometry_bodies(batch, "2003 BF91;", "250")
    except Exception as error:
        _skip_remote_data_error(error, "Horizons Hubble target ephemeris")

    set_tracking_supplementary_data_in_bodies(bodies, supplementary_data)
    observed_observations = create_observation_collection_from_tracking_data(
        tracking_data,
        bodies,
    )
    residuals = (
        _simulate_angular_position_residuals(observed_observations, bodies)
        .reshape(2, -1, order="F")
        .T
    )
    # RA is periodic; after fixing the HST receiver frame handling, the physical
    # residual can differ from Tudat's raw value by an integer multiple of 2*pi.
    residuals[:, 0] = (residuals[:, 0] + np.pi) % (2.0 * np.pi) - np.pi
    residuals_arcsec = residuals * 180.0 / np.pi * 3600.0

    assert residuals_arcsec.shape[1] == 2
    assert np.all(np.isfinite(residuals_arcsec))
    assert np.max(np.abs(residuals_arcsec)) < 0.1
    assert np.all(np.sqrt(np.mean(residuals_arcsec**2, axis=0)) < 0.05)
