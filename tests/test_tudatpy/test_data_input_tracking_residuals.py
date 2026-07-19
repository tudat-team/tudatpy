from pathlib import Path
from bisect import bisect_left, bisect_right
from urllib.request import urlretrieve

import numpy as np

from tudatpy.data_input import resource_paths as data_paths
from tudatpy.data_input.tracking_data import TrackingData
from tudatpy.data_input.tracking_data.fdets import FdetDateFormat, read_fdets_data
from tudatpy.data_input.tracking_data.ifms import read_ifms_data
from tudatpy.data_input.tracking_data.odf import read_odf_data
from tudatpy.data_input.tracking_data.psf import read_psf_data
from tudatpy.data_input.tracking_data.tnf import read_tnf_data
from tudatpy.estimation.observations import (
    ObservationDataset,
    create_observation_dataset_from_tracking_data,
    create_compressed_doppler_dataset,
    observation_query,
    observation_simulation_settings_from_dataset,
    set_tracking_supplementary_data_in_bodies,
    simulate_observation_dataset,
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
from tudatpy.kernel.estimation.observations_setup import observations_simulation_settings
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


def _dataset_observable_link_definitions(observation_dataset):
    seen_definitions = set()
    for set_id in range(observation_dataset.number_of_observation_sets):
        metadata = observation_dataset.get_observation_set_metadata(set_id)
        definition_key = (metadata.observable_type, metadata.link_definition_id)
        if definition_key in seen_definitions:
            continue
        seen_definitions.add(definition_key)
        yield metadata.observable_type, observation_dataset.link_definition(
            metadata.link_definition_id
        )


def _observation_vector(observation_dataset):
    return np.asarray(observation_dataset.ordered_flattened_observation_data().observation_vector)


def _set_dataset_reference_points_from_switch_history(
    observation_dataset,
    bodies,
    antenna_switch_history,
    spacecraft_name,
    link_end_type,
):
    switch_times = sorted(antenna_switch_history)
    if len(switch_times) < 2:
        raise RuntimeError("Antenna switch history must bound at least one time interval.")

    reference_point_names = {}
    combined_dataset = ObservationDataset()
    for interval_index, (start_time, end_time) in enumerate(
        zip(switch_times[:-1], switch_times[1:])
    ):
        position = np.asarray(antenna_switch_history[start_time])
        position_key = tuple(position)
        if position_key not in reference_point_names:
            reference_point_name = f"Antenna{len(reference_point_names) + 1}"
            bodies.get_body(spacecraft_name).system_models.set_reference_point(
                reference_point_name, position
            )
            reference_point_names[position_key] = reference_point_name

        interval_condition = observation_query.time <= end_time
        if interval_index == 0:
            interval_condition &= observation_query.time >= start_time
        else:
            interval_condition &= observation_query.time > start_time

        interval_dataset = observation_dataset.create_new_and_keep(interval_condition)
        if interval_dataset.number_of_observations == 0:
            continue
        interval_dataset.set_link_end_reference_point(
            spacecraft_name,
            reference_point_names[position_key],
            link_end_type,
        )
        for set_id in range(interval_dataset.number_of_observation_sets):
            combined_dataset.add_observation_set_from_dataset(interval_dataset, set_id)

    return combined_dataset


def _simulate_dsn_n_way_averaged_doppler(observation_dataset, bodies):
    light_time_settings = [
        light_time_corrections.first_order_relativistic_light_time_correction(["Sun"])
    ]
    observation_model_settings = []
    for _, link_definition in _dataset_observable_link_definitions(observation_dataset):
        observation_model_settings.append(
            model_settings.dsn_n_way_doppler_averaged(
                link_definition,
                light_time_settings,
                biases.absolute_bias(np.zeros(1)),
                light_time_corrections.light_time_convergence_settings(True),
                False,
            )
        )

    observation_simulators = observations_simulation_settings.create_observation_simulators(
        observation_model_settings, bodies
    )
    observation_simulation_settings = observation_simulation_settings_from_dataset(
        observation_dataset, bodies
    )
    return simulate_observation_dataset(
        observation_simulation_settings, observation_simulators, bodies
    )


def _simulate_dsn_n_way_averaged_doppler_with_corrections(
    observation_dataset, bodies, light_time_settings
):
    observation_model_settings = []
    for _, link_definition in _dataset_observable_link_definitions(observation_dataset):
        observation_model_settings.append(
            model_settings.dsn_n_way_doppler_averaged(
                link_definition,
                light_time_settings,
                biases.absolute_bias(np.zeros(1)),
                light_time_corrections.light_time_convergence_settings(True),
            )
        )

    observation_simulators = observations_simulation_settings.create_observation_simulators(
        observation_model_settings, bodies
    )
    observation_simulation_settings = observation_simulation_settings_from_dataset(
        observation_dataset, bodies
    )
    return simulate_observation_dataset(
        observation_simulation_settings, observation_simulators, bodies
    )


def _simulate_doppler_measured_frequency(observation_dataset, bodies):
    light_time_settings = [
        light_time_corrections.first_order_relativistic_light_time_correction(
            ["Sun", "Moon", "Earth"]
        )
    ]
    observation_model_settings = []
    for _, link_definition in _dataset_observable_link_definitions(observation_dataset):
        observation_model_settings.append(
            model_settings.doppler_measured_frequency(
                link_definition,
                light_time_settings,
            )
        )

    observation_simulators = observations_simulation_settings.create_observation_simulators(
        observation_model_settings, bodies
    )
    observation_simulation_settings = observation_simulation_settings_from_dataset(
        observation_dataset, bodies
    )
    return simulate_observation_dataset(
        observation_simulation_settings, observation_simulators, bodies
    )


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


def _simulate_pixel_coordinates(observation_dataset, bodies):
    observation_model_settings = []
    for _, link_definition in _dataset_observable_link_definitions(observation_dataset):
        observation_model_settings.append(
            model_settings.pixel_coordinates(
                link_definition,
                [],
                None,
                light_time_corrections.light_time_convergence_settings(True),
                True,
            )
        )

    observation_simulators = observations_simulation_settings.create_observation_simulators(
        observation_model_settings, bodies
    )
    observation_simulation_settings = observation_simulation_settings_from_dataset(
        observation_dataset, bodies
    )
    return simulate_observation_dataset(
        observation_simulation_settings, observation_simulators, bodies
    )


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
    uncompressed_observations = create_observation_dataset_from_tracking_data(tracking_data, bodies)
    observed_observations = create_compressed_doppler_dataset(
        uncompressed_observations,
        60,
        earth_fixed_ground_station_positions=ground_station.get_radio_telescope_positions(),
    )

    computed_observations = _simulate_dsn_n_way_averaged_doppler(observed_observations, bodies)
    residuals = _observation_vector(observed_observations) - _observation_vector(
        computed_observations
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

    uncompressed_observations = create_observation_dataset_from_tracking_data(
        doppler_tracking_data, bodies
    )
    observed_observations = create_compressed_doppler_dataset(uncompressed_observations, 60, 10)
    observed_observations = _set_dataset_reference_points_from_switch_history(
        observed_observations,
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
    residuals = _observation_vector(observed_observations) - _observation_vector(
        computed_observations
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
    uncompressed_observations = create_observation_dataset_from_tracking_data(tracking_data, bodies)
    observed_observations = create_compressed_doppler_dataset(uncompressed_observations, 60, 10)
    bodies.get_body("MRO").system_models.transponder_delay = 1.4149e-6

    mro_center_of_mass_position = np.array([0.0, -1.11, 0.0])
    antenna_position_history = {}
    for set_id in range(observed_observations.number_of_observation_sets):
        observation_times = observed_observations.observation_times_for_set(set_id)
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
    bodies.get_body("MRO").system_models.set_reference_point(
        "Antenna",
        antenna_ephemeris,
    )
    observed_observations.set_link_end_reference_point(
        "MRO",
        "Antenna",
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
    residuals = _observation_vector(computed_observations) - _observation_vector(
        observed_observations
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
    all_observations = create_observation_dataset_from_tracking_data(
        tracking_data,
        _create_voyager_triton_psf_bodies(test_data_path, np.zeros(6), supplementary_data),
    )
    all_observation_times = np.asarray(
        all_observations.ordered_flattened_observation_data().times, dtype=float
    )

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
        observed_observations = create_observation_dataset_from_tracking_data(
            [single_tracking_data], bodies
        )
        computed_observations = _simulate_pixel_coordinates(observed_observations, bodies)
        residuals.append(
            _observation_vector(computed_observations) - _observation_vector(observed_observations)
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
    tracking_data[0].add_string_vector_ancillary_setting("frequency bands", ["X-band", "X-band"])
    observed_observations = create_observation_dataset_from_tracking_data(tracking_data, bodies)
    computed_observations = _simulate_doppler_measured_frequency(observed_observations, bodies)
    residuals = _observation_vector(computed_observations) - _observation_vector(
        observed_observations
    )
    residual_scatter = residuals - np.mean(residuals)

    assert residuals.size == 120
    assert np.sqrt(np.mean(residual_scatter**2)) < 1.0e-2
    assert np.max(np.abs(residual_scatter)) < 2.5e-2
