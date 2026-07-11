from pathlib import Path
from urllib.request import urlretrieve

import numpy as np

from tudatpy.data_access.tracking import (
    create_observation_collection,
    set_tracking_supplementary_data_in_bodies,
)
from tudatpy.data_access.tracking.fdets import FdetDateFormat, read_fdets_files
from tudatpy.data_access.tracking.ifms import read_ifms_files
from tudatpy.data_access.tracking.odf import read_odf_files
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
    observations_simulation_settings,
    observations_wrapper,
)
from tudatpy.kernel.interface import spice


def _test_data_path() -> Path:
    repository_root = Path(__file__).resolve().parents[2]
    for relative_path in (
        "cmake-build-release/tests/data",
        "cmake-build-debug/tests/data",
        "tests/data",
    ):
        data_path = repository_root / relative_path
        if data_path.exists():
            return data_path
    raise RuntimeError("Could not find Tudat test data directory.")


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


def _create_grail_bodies(test_data_path: Path):
    spice.load_standard_kernels()
    spice.load_kernel(str(test_data_path / "grail_shortened.bsp"))

    body_settings = environment_setup.get_default_body_settings(
        ["Earth", "Sun", "Moon"], "SSB", "J2000"
    )
    earth_settings = body_settings.get("Earth")
    earth_settings.shape_settings = shape.oblate_spherical_spice()
    earth_settings.rotation_model_settings = rotation_model.gcrs_to_itrs(base_frame="J2000")
    earth_settings.ground_station_settings = ground_station.dsn_stations()

    body_settings.add_empty_settings("GRAIL-A")
    body_settings.get("GRAIL-A").ephemeris_settings = ephemeris.direct_spice("Moon", "J2000")

    bodies = environment_setup.create_system_of_bodies(body_settings)
    grail_systems = environment.VehicleSystems()
    grail_systems.set_default_transponder_turnaround_ratio_function()
    grail_systems.set_reference_point("Antenna", np.array([-0.082, 0.152, -0.810]))
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


def _keep_first_observations(tracking_data, number_of_observations: int):
    for index in range(len(tracking_data.observations) - 1, number_of_observations - 1, -1):
        tracking_data.remove_single_observation_entry(index)


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
    return observations_wrapper.simulate_observations(
        observation_simulation_settings, observation_simulators, bodies
    )


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
    return observations_wrapper.simulate_observations(
        observation_simulation_settings, observation_simulators, bodies
    )


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
    return observations_wrapper.simulate_observations(
        observation_simulation_settings, observation_simulators, bodies
    )


def test_ifms_mex_residuals_are_millihertz_level():
    test_data_path = _test_data_path()
    bodies = _create_mex_bodies(test_data_path)

    tracking_data, supplementary_data = read_ifms_files(
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
    uncompressed_observations = create_observation_collection(tracking_data, bodies)
    observed_observations = observations_wrapper.create_compressed_doppler_collection(
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

    bodies = _create_grail_bodies(test_data_path)
    tracking_data, supplementary_data = read_odf_files([str(odf_file)], "GRAIL-A", "Earth", False)
    set_tracking_supplementary_data_in_bodies(bodies, supplementary_data)
    doppler_tracking_data = [
        data for data in tracking_data if data.observable_type == "DsnNWayAveragedDoppler"
    ][10000:10120]

    observed_observations = create_observation_collection(doppler_tracking_data, bodies)
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

    assert residuals.size == 120
    assert np.sqrt(np.mean(residuals**2)) < 1.0e-2
    assert np.max(np.abs(residuals)) < 3.0e-2


def test_fdets_juice_short_arc_residual_scatter_is_millihertz_level():
    test_data_path = _test_data_path()
    bodies = _create_juice_bodies(test_data_path)

    tracking_data, supplementary_data = read_fdets_files(
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
    tracking_data[0].add_ancillary_settings("frequency bands", [1.0, 1.0])
    observed_observations = create_observation_collection(tracking_data, bodies)
    computed_observations = _simulate_doppler_measured_frequency(observed_observations, bodies)
    residuals = np.asarray(computed_observations.concatenated_observations) - np.asarray(
        observed_observations.concatenated_observations
    )
    residual_scatter = residuals - np.mean(residuals)

    assert residuals.size == 120
    assert np.sqrt(np.mean(residual_scatter**2)) < 1.0e-2
    assert np.max(np.abs(residual_scatter)) < 2.5e-2
