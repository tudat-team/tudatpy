import warnings
from pathlib import Path

# The ObservationCollection interface used below is deprecated in this build
warnings.filterwarnings("ignore", category=DeprecationWarning)

import matplotlib

matplotlib.use("Agg")
import numpy as np
from matplotlib import pyplot as plt
from tudatpy import constants
from tudatpy.astro import gravitation
from tudatpy.data_input.environment_data import spice
from tudatpy.data_input.tracking_data import (
    RampedFrequencySupplementaryData,
    TrackingSupplementaryData,
)
from tudatpy.data_input.tracking_data.odf import read_odf_data
from tudatpy.dynamics import environment_setup, parameters_setup, propagation_setup
from tudatpy.estimation import (
    estimation_analysis,
    observable_models_setup,
    observations,
    observations_setup,
)
from tudatpy.estimation.observable_models_setup import model_settings
from tudatpy.estimation.observable_models_setup.links import LinkEndType
from tudatpy.estimation.observable_models_setup.model_settings import ObservableType
from tudatpy.estimation.observations.observations_processing import (
    ObservationSetSplitterType,
    observation_parser,
    observation_set_splitter,
)
from tudatpy.estimation.observations_setup.ancillary_settings import (
    FrequencyBands,
    ObservationAncillarySimulationVariable,
)

import juno_utils

sc_name = "JUNO"
sc_id = 61

# One AWVR unit at DSS-25. Set use_awvr=False for TSAC on the same retained counts.
use_awvr = True
awvr_unit = "01"
use_measured_dry_delay = False

# data/output paths relative to this script
mission_file_directory = Path(__file__).parent / "data"
mission_output_directory = Path(__file__).parent / "output"


odf_file_names_per_pj = {
    3: [
        "GRV_JUGR_2016346_1050M25MC001V01.ODF",
        "GRV_JUGR_2016346_1740XMMMC005V01.ODF",
    ],
    6: [
        "GRV_JUGR_2017138_2346M25MC001V01.ODF",
        "GRV_JUGR_2017139_1621XMMMC005V01.ODF",
    ],
}
pjs = list(odf_file_names_per_pj)

# Observations are retained within [perijove - 12 h, perijove + 8 h]
pj_time_window_offsets = (-12 * 3600.0, 8 * 3600.0)
# Observation sets are split at perijove + offset, so that no set is dropped by the time windows
pj_split_offset = 8 * 3600.0

pj_to_epoch = juno_utils.load_perijove_data(mission_file_directory / "perijove-times.txt")
pj_time_windows = {
    pj: tuple(pj_to_epoch[pj] + offset for offset in pj_time_window_offsets) for pj in pjs
}
pj_time_parsers = {
    pj: observation_parser(tuple(t.to_float() for t in window))
    for pj, window in pj_time_windows.items()
}

# Doppler counts are compressed to 60 s integration time (from 1 s science and 5 s navigation data)
doppler_integration_time = 60.0
# Weight of all Doppler observations
doppler_noise = 0.002


def load_kernels():
    # ## Load SPICE Kernels

    kernels_directory = mission_file_directory / "kernels"

    juno_utils.load_juno_default_kernels(kernels_directory.as_posix())

    frames_def_kernel = kernels_directory / "fk" / "juno_v12.tf"
    sclk_kernel = kernels_directory / "sclk" / "JNO_SCLKSCET.00212.tsc"
    orientation_kernel_names = [
        "juno_sc_rec_161210_161211_v01.bc",
        "juno_sc_rec_161211_161217_v01.bc",
        "juno_sc_rec_161218_161224_v01.bc",
        "juno_sc_rec_161225_161231_v01.bc",
        "juno_sc_rec_170101_170107_v01.bc",
        "juno_sc_rec_170108_170114_v01.bc",
        "juno_sc_rec_170115_170121_v01.bc",
        "juno_sc_rec_170122_170128_v01.bc",
        "juno_sc_rec_170129_170204_v01.bc",
        "juno_sc_rec_170514_170520_v01.bc",
        "juno_sc_rec_170518_170518_v01.bc",
        "juno_sc_rec_170518_170519_v02.bc",
    ]
    orientation_kernels = [kernels_directory / "ck" / name for name in orientation_kernel_names]
    trajectory_kernels = [
        # kernels_directory / "spk" / "spk_rec_161115_170106_170113.bsp",
        # kernels_directory / "spk" / "spk_rec_170422_170608_170621.bsp",
        kernels_directory / "spk" / "juno_rec_orbit.bsp",
        kernels_directory / "spk" / "juno_struct_v04.bsp",
    ]

    # Load frame definition file (required for HGA and spacecraft-fixed frames definition)
    spice.load_kernel(frames_def_kernel.as_posix())
    # Load SC clock file (required for conversion between SCLK and UTC, used when retrieving SC orientation)
    spice.load_kernel(sclk_kernel.as_posix())

    # Load orientation kernels
    for orientation_file in orientation_kernels:
        spice.load_kernel(orientation_file.as_posix())
    # Load trajectory kernels
    for trajectory_file in trajectory_kernels:
        spice.load_kernel(trajectory_file.as_posix())


def setup_observations():
    # ## Load Observations
    #

    odf_files = [
        (mission_file_directory / "odf" / name).as_posix()
        for pj in pjs
        for name in odf_file_names_per_pj[pj]
    ]
    tracking_data, tracking_supplementary_data = read_odf_data(
        odf_files,
        spacecraft_name=sc_name,
    )

    # Only the Doppler data is used
    tracking_data = [
        t_data for t_data in tracking_data if t_data.observable_type == "DsnNWayAveragedDoppler"
    ]

    return tracking_data, tracking_supplementary_data


def setup_bodies(tracking_data, tracking_supplementary_data):
    # ## Setup Environment

    # Create default body settings for celestial bodies
    bodies_to_create = [
        "Earth",
        "Sun",
        "Jupiter",
        "Io",
        "Europa",
        "Ganymede",
        "Callisto",
        "Amalthea",
    ]
    global_frame_origin = "SSB"
    # J2000 instead of ECLIPJ2000, since the Jupiter rotation model below is defined w.r.t. J2000
    global_frame_orientation = "J2000"

    start_time = min(min(t_data.epochs) for t_data in tracking_data) - 86400
    end_time = max(max(t_data.epochs) for t_data in tracking_data) + 86400

    body_settings = environment_setup.get_default_body_settings_time_limited(
        bodies_to_create,
        start_time,
        end_time,
        global_frame_origin,
        global_frame_orientation,
        # time_step=
    )

    # Modify Earth default settings
    body_settings.get("Earth").rotation_model_settings = (
        environment_setup.rotation_model.gcrs_to_itrs(
            environment_setup.rotation_model.iau_2006,
            global_frame_orientation,
        )
    )
    # since we defined the Earth rotation model as GCRS to ITRS, we need to set the gravity field reference frame to ITRS as well
    body_settings.get("Earth").gravity_field_settings.associated_reference_frame = "ITRS"

    body_settings.get("Earth").shape_settings = environment_setup.shape.oblate_spherical_spice()
    body_settings.get("Earth").ground_station_settings = (
        environment_setup.ground_station.dsn_stations()
    )

    # DSS-25 transmits X- and Ka-band simultaneously, but a ground station has only a single
    # transmitting frequency. The Ka-band links are therefore defined with the DSN stations
    # of "EarthKa": a copy of the Earth (same ephemeris, rotation and stations), without mass
    body_settings.add_empty_settings("EarthKa")
    body_settings.get("EarthKa").ephemeris_settings = body_settings.get("Earth").ephemeris_settings
    body_settings.get("EarthKa").rotation_model_settings = body_settings.get(
        "Earth"
    ).rotation_model_settings
    # Earth radii from pck00011, as used by oblate_spherical_spice for the Earth
    body_settings.get("EarthKa").shape_settings = environment_setup.shape.oblate_spherical(
        6378136.6, 1.0 - 6356751.9 / 6378136.6
    )
    body_settings.get("EarthKa").ground_station_settings = (
        environment_setup.ground_station.dsn_stations()
    )

    # ### Setup Jupiter Body Settings
    #
    # Gravity field with nominal J2-J12 (Iess et al. 2018) and small degree-2 tesseral terms,
    # pole and rotation rate, and the tide raised by the moons and Sun (k22 only)

    cosine_coefficients = np.zeros((25, 25))
    sine_coefficients = np.zeros((25, 25))
    cosine_coefficients[0, 0] = 1.0
    cosine_coefficients[2:13, 0] = -1e-6 * np.array(
        [
            14696.575,
            -0.048,
            -586.608,
            -0.072,
            34.202,
            0.117,
            -2.421,
            -0.115,
            0.175,
            0.020,
            0.051,
        ]
    )
    cosine_coefficients[2, 1] = -0.006e-6
    sine_coefficients[2, 1] = -0.017e-6
    cosine_coefficients[2, 2] = 0.002e-6
    sine_coefficients[2, 2] = 0.008e-6
    cosine_coefficients, sine_coefficients = gravitation.normalize_spherical_harmonic_coefficients(
        cosine_coefficients, sine_coefficients
    )
    body_settings.get("Jupiter").gravity_field_settings = (
        environment_setup.gravity_field.spherical_harmonic(
            1.26686533e17,
            71492000.0,
            cosine_coefficients,
            sine_coefficients,
            "IAU_JUPITER",
        )
    )
    body_settings.get("Jupiter").rotation_model_settings = (
        environment_setup.rotation_model.iau_rotation_model(
            global_frame_orientation,
            "IAU_JUPITER",
            np.deg2rad(284.95),
            np.deg2rad([268.0570781451589, 64.4958232029158]),
            np.deg2rad(870.536) / 86400.0,
            np.deg2rad([-0.006554328185586419, 0.002476496988122852]) / (36525.0 * 86400.0),
            {},
            {},
            reference_epoch_j2000=0.0,
        )
    )
    body_settings.get("Jupiter").gravity_field_variation_settings = [
        environment_setup.gravity_field_variation.mode_coupled_solid_body_tide(
            ["Io", "Europa", "Ganymede", "Callisto", "Sun"],
            {(2, 2): {(2, 2): 0.589}},
        )
    ]

    # ### Setup JUNO Body Settings

    body_settings.add_empty_settings(sc_name)

    # Retrieve translational ephemeris from SPICE (replaced by the propagated orbit for the estimation)
    body_settings.get(sc_name).ephemeris_settings = environment_setup.ephemeris.interpolated_spice(
        start_time,
        end_time,
        10.0,
        "Jupiter",
        global_frame_orientation,
    )

    # Retrieve rotational ephemeris from SPICE
    body_settings.get(sc_name).rotation_model_settings = environment_setup.rotation_model.spice(
        global_frame_orientation, f"{sc_name}_SPACECRAFT"
    )

    # Point-mass gravity field (as spherical harmonics, needed for the mutual gravity with the moons below)
    body_settings.get(sc_name).gravity_field_settings = (
        environment_setup.gravity_field.spherical_harmonic(
            constants.GRAVITATIONAL_CONSTANT * 2378.0,
            1.0,
            np.ones((1, 1)),
            np.zeros((1, 1)),
            f"{sc_name}_SPACECRAFT",
        )
    )

    # Cannonball solar radiation pressure, occulted by Jupiter
    body_settings.get(sc_name).radiation_pressure_target_settings = (
        environment_setup.radiation_pressure.cannonball_radiation_target(
            60.0, 1.2, {"Sun": ["Jupiter"]}
        )
    )

    # ### Create Bodies
    #
    # We can now create our bodies from the settings and set the default DSN transponder turnaround ratios in the spacecraft object.

    # Create environment
    bodies = environment_setup.create_system_of_bodies(body_settings)

    bodies.get(sc_name).system_models.set_transponder_turnaround_ratio(
        juno_utils.juno_turnaround_ratio()
    )

    # Set the ramp tables: X-band uplink on the Earth stations, Ka-band uplink on the EarthKa stations
    # (the ODF reader adds a non-finite placeholder ramp per station, which is skipped)
    band_supplementary_data = []
    for station_data in tracking_supplementary_data:
        ramps = [
            ramp
            for block in station_data.frequency_supplementary_data
            for ramp in block.frequency_ramps
            if np.isfinite(ramp.start_frequency) and np.isfinite(ramp.frequency_rate)
        ]
        for body_name, band_ramps in (
            ("Earth", [ramp for ramp in ramps if ramp.start_frequency < 20e9]),
            ("EarthKa", [ramp for ramp in ramps if ramp.start_frequency > 20e9]),
        ):
            if len(band_ramps) > 0:
                band_data = TrackingSupplementaryData(body_name, station_data.reference_point_name)
                band_data.set_frequency_supplementary_data(
                    [RampedFrequencySupplementaryData(band_ramps)]
                )
                band_supplementary_data.append(band_data)

    observations.set_tracking_supplementary_data_in_bodies(bodies, band_supplementary_data)

    return bodies


def create_observation_collection(tracking_data, bodies, frequency_band, awvr_zenith_delays=None):
    # ## Create Observation Collection
    #
    # With the bodies set up, the primitive tracking data objects can be turned into an observation collection, which will be used to compute the pre-fit residuals.
    # The ground station locations from the bodies object is used to perform the time tag conversion from UTC to TDB.
    # Before computing the observations, the tracking supplementary data (in this case the ramp tables) needs to be linked with the environment, which is done through the `set_tracking_supplementary_data_in_bodies` function.
    #
    # One collection is created per frequency band (X/X or Ka/Ka).

    band_tracking_data = [
        t_data
        for t_data in tracking_data
        if t_data.get_ancillary_settings_string_vector()["frequency bands"]
        == [frequency_band, frequency_band]
    ]

    # Compress science (1 s) and navigation (5 s) data separately to 60 s counts,
    # without bridging gaps in the data
    compressed_collections = []
    for original_integration_time in (1.0, 5.0):
        integration_time_tracking_data = [
            t_data
            for t_data in band_tracking_data
            if t_data.get_ancillary_settings_double()["Doppler observable integration time"]
            == original_integration_time
        ]
        if len(integration_time_tracking_data) == 0:
            continue

        observation_collection = observations.create_observation_collection_from_tracking_data(
            integration_time_tracking_data, bodies
        )
        compression_ratio = round(doppler_integration_time / original_integration_time)
        compressed_collections.append(
            observations.create_compressed_doppler_collection(
                observation_collection,
                compression_ratio=compression_ratio,
                max_arc_gap=original_integration_time + 0.01,
            )
        )
    compressed_observations = observations.merge_observation_collections(compressed_collections)

    # Keep complete counts on one correction source at both Earth events.
    # This selection is also applied when use_awvr=False for a fair TSAC comparison.
    if awvr_zenith_delays is not None:
        compressed_observations = juno_utils.exclude_awvr_source_transition_counts(
            compressed_observations,
            bodies,
            awvr_zenith_delays,
            [(start.to_float(), end.to_float()) for start, end in pj_time_windows.values()],
            doppler_integration_time,
        )

    # ### Update link-end delays

    compressed_observations = juno_utils.update_dss_25_ka_band_delays(compressed_observations)

    # ### Reference Observations to Spacecraft Antenna
    #
    # The reference point of the observations is set in the `ObservationCollection` to the high gain antenna.
    # We use the integrated spacecraft reference point of the SPICE kernels (NAIF ID `-74`) as frame origin for the antenna ephemeris, as opposed to the reference point of the mechanical assembly (NAIF ID `-74000`), since we have also assigned our spacecraft the ephemeris of `-74`.
    #
    start_time = min(min(t_data.epochs) for t_data in tracking_data) - 86400
    end_time = max(max(t_data.epochs) for t_data in tracking_data) + 86400

    hga_ephemeris_settings = environment_setup.ephemeris.interpolated_spice(
        start_time,
        end_time,
        10.0,
        frame_origin=sc_name,
        frame_orientation=f"{sc_name}_SPACECRAFT",
        body_name_to_use=f"{sc_name}_HGA",
    )

    hga_ephemeris = environment_setup.ephemeris.create_ephemeris(
        hga_ephemeris_settings, f"{sc_name}_HGA"
    )

    # Set the spacecraft's reference point position to that of the antenna (in the spacecraft-fixed frame)
    compressed_observations.set_reference_point(
        bodies,
        hga_ephemeris,
        antenna_name=f"{sc_name}_HGA",
        spacecraft_name=sc_name,
        link_end_type=observable_models_setup.links.LinkEndType.reflector1,
    )

    # Create a parser to select only Doppler observations
    doppler_parser = observation_parser(ObservableType.dsn_n_way_averaged_doppler_type)

    # One time-window parser per perijove, also used to split the observations for plotting

    any_pj_time_parser = observation_parser(
        list(pj_time_parsers.values()), combine_conditions=False
    )

    obs_time_bounds_parser = observation_parser(
        [doppler_parser, any_pj_time_parser], combine_conditions=True
    )

    # Split the observation sets at perijove + offset, so that sets extending beyond
    # the time windows are not dropped as a whole by the time-bounds parser
    pj_split_epochs = sorted((pj_to_epoch[pj] + pj_split_offset).to_float() for pj in pjs)
    compressed_observations.split_observation_sets(
        observation_set_splitter(ObservationSetSplitterType.time_tags_splitter, pj_split_epochs)
    )

    compressed_observations = observations.create_new_observation_collection(
        compressed_observations, obs_time_bounds_parser
    )

    compressed_observations = juno_utils.apply_spin_effect_correction(
        compressed_observations, bodies, apply_positive_correction=True
    )

    # Move the Ka-band links to the EarthKa stations
    if frequency_band == "Ka-band":
        ka_observation_sets = []
        for obs_set in compressed_observations.get_single_observation_sets():
            link_ends = obs_set.link_definition.link_ends
            for link_end_type in (LinkEndType.transmitter, LinkEndType.receiver):
                link_ends[link_end_type] = (
                    observable_models_setup.links.body_reference_point_link_end_id(
                        "EarthKa", link_ends[link_end_type].reference_point
                    )
                )
            ka_observation_sets.append(
                observations.create_single_observation_set(
                    obs_set.observable_type,
                    link_ends,
                    obs_set.list_of_observations,
                    obs_set.observation_times,
                    obs_set.reference_link_end,
                    obs_set.ancillary_settings,
                )
            )
        compressed_observations = observations.ObservationCollection(ka_observation_sets)

    return compressed_observations


def setup_observation_models(compressed_observations, awvr_zenith_delays=None):
    tro_files = [
        mission_file_directory / "tro" / "GRV_JUGR_2016336_2017001_V01.TRO",
        mission_file_directory / "tro" / "GRV_JUGR_2017121_2017152_V01.TRO",
    ]
    ion_files = [
        (mission_file_directory / "ion" / "GRV_JUGR_2016336_2017001_V01.ION").as_posix(),
        (mission_file_directory / "ion" / "GRV_JUGR_2017121_2017152_V01.ION").as_posix(),
    ]
    spacecraft_name_per_id = {sc_id: sc_name}

    # Shapiro delay of Jupiter as well (the link passes close to Jupiter)
    relativistic_lt_correction = observable_models_setup.light_time_corrections.first_order_relativistic_light_time_correction(
        ["Sun", "Earth", "Jupiter"]
    )

    # The media corrections need the Earth body of the link (Earth or EarthKa)
    doppler_lt_corrections = dict()
    for earth_name in ("Earth", "EarthKa"):
        tropospheric_lt_correction = observable_models_setup.light_time_corrections.dsn_tabulated_tropospheric_light_time_correction(
            [f.as_posix() for f in tro_files],
            body_with_atmosphere_name=earth_name,
        )
        if use_awvr:
            if awvr_zenith_delays is None:
                raise ValueError("AWVR correction requires screened AMC samples.")
            # Niell mapping is native to this correction. Measured wet REPLACES
            # TSAC wet; TSAC dry remains active unless explicitly selected above.
            tropospheric_lt_correction = observable_models_setup.light_time_corrections.measured_tropospheric_light_time_correction(
                {"DSS-25": awvr_zenith_delays},
                tropospheric_lt_correction,
                maximum_interpolation_gap=30.01,
                use_measured_dry_delay=use_measured_dry_delay,
                body_with_atmosphere_name=earth_name,
            )
        ionospheric_lt_correction = observable_models_setup.light_time_corrections.dsn_tabulated_ionospheric_light_time_correction(
            ion_files, spacecraft_name_per_id, body_with_atmosphere_name=earth_name
        )

        # solar_corona_lt_correction = observable_models_setup.light_time_corrections.inverse_power_series_solar_corona_light_time_correction(
        #     [1.70e12],
        #     [2.44],
        # )

        doppler_lt_corrections[earth_name] = [
            relativistic_lt_correction,
            tropospheric_lt_correction,
            ionospheric_lt_correction,
            # solar_corona_lt_correction,
        ]

    doppler_link_definitions = compressed_observations.link_definitions_per_observable[
        model_settings.ObservableType.dsn_n_way_averaged_doppler_type
    ]
    observable_model_settings = [
        model_settings.dsn_n_way_doppler_averaged(
            link_definition,
            doppler_lt_corrections[link_definition.link_ends[LinkEndType.receiver].body_name],
        )
        for link_definition in doppler_link_definitions
    ]

    return observable_model_settings


def setup_propagator(bodies):
    # ## Setup Propagation
    #
    # One arc per perijove, propagated forwards and backwards from the perijove

    # Replace the SPICE ephemeris of Juno by a multi-arc tabulated ephemeris, which is set from the propagation
    juno_ephemeris_settings = environment_setup.ephemeris.tabulated({}, "Jupiter", "J2000")
    juno_ephemeris_settings.make_multi_arc_ephemeris = True
    bodies.get(sc_name).ephemeris = environment_setup.ephemeris.create_ephemeris(
        juno_ephemeris_settings, sc_name
    )

    accelerations = {
        "Jupiter": [
            propagation_setup.acceleration.spherical_harmonic_gravity(24, 2),
            propagation_setup.acceleration.relativistic_correction(use_schwarzschild=True),
        ],
        "Sun": [
            propagation_setup.acceleration.point_mass_gravity(),
            propagation_setup.acceleration.radiation_pressure(),
        ],
        "Amalthea": [propagation_setup.acceleration.point_mass_gravity()],
    }
    # Mutual gravity, to include the indirect effect of Jupiter's oblateness (acceleration of Jupiter by the moon)
    for moon in ("Io", "Europa", "Ganymede", "Callisto"):
        accelerations[moon] = [
            propagation_setup.acceleration.mutual_spherical_harmonic_gravity(0, 0, 0, 0, 24, 2)
        ]
    acceleration_models = propagation_setup.create_acceleration_models(
        bodies, {sc_name: accelerations}, [sc_name], ["Jupiter"]
    )

    arc_propagator_settings = []
    for pj in pjs:
        # Arc covers the observation window, plus the light time (< 1 h) and count interval
        arc_start = pj_to_epoch[pj] + pj_time_window_offsets[0] - 2 * 3600.0
        arc_end = pj_to_epoch[pj] + pj_time_window_offsets[1] + 3600.0

        integrator_settings = propagation_setup.integrator.runge_kutta_variable_step(
            5.0,
            propagation_setup.integrator.rkf_78,
            propagation_setup.integrator.step_size_control_elementwise_scalar_tolerance(
                1e-12, 1e-6
            ),
            propagation_setup.integrator.step_size_validation(1e-4, 30.0),
        )
        termination_settings = propagation_setup.propagator.non_sequential_termination(
            propagation_setup.propagator.time_termination(arc_end.to_float()),
            propagation_setup.propagator.time_termination(arc_start.to_float()),
        )
        initial_state = spice.get_body_cartesian_state_at_epoch(
            sc_name, "Jupiter", "J2000", "NONE", pj_to_epoch[pj].to_float()
        )
        arc_propagator_settings.append(
            propagation_setup.propagator.translational(
                ["Jupiter"],
                acceleration_models,
                [sc_name],
                initial_state,
                pj_to_epoch[pj],
                integrator_settings,
                termination_settings,
            )
        )

    return propagation_setup.propagator.multi_arc(arc_propagator_settings)


def setup_parameters(bodies, propagator_settings):
    # ## Setup Estimated Parameters
    #
    # Arc initial states, Jupiter GM, zonals up to degree 24, degree-2 tesserals, pole (and rate), and k22

    parameter_settings = parameters_setup.initial_states(propagator_settings, bodies)
    parameter_settings += [
        parameters_setup.gravitational_parameter("Jupiter"),
        parameters_setup.spherical_harmonics_c_coefficients_block(
            "Jupiter", [(degree, 0) for degree in range(2, 25)] + [(2, 1), (2, 2)]
        ),
        parameters_setup.spherical_harmonics_s_coefficients_block("Jupiter", [(2, 1), (2, 2)]),
        parameters_setup.iau_rotation_model_pole("Jupiter"),
        parameters_setup.iau_rotation_model_pole_rate("Jupiter"),
        parameters_setup.mode_coupled_k_love_numbers(
            "Jupiter",
            {(2, 2): [(2, 2)]},
            ["Io", "Europa", "Ganymede", "Callisto", "Sun"],
        ),
    ]
    parameters = parameters_setup.create_parameter_set(
        parameter_settings, bodies, propagator_settings
    )

    # A priori uncertainties of GM, normalized C20-C240, C21, C22, S21, S22 and k22 (index 12-39 and 44),
    # no a priori constraint on the initial states and pole
    apriori_sigmas = np.zeros(parameters.parameter_set_size)
    apriori_sigmas[12:40] = [
        2.0e9,
        1.7003488324182765e-06,
        1.1183414280158652e-06,
        3.072131943422285e-07,
        5.19961517908383e-07,
        5.360610733932582e-07,
        9.724281614730766e-07,
        1.2165561521800082e-06,
        2.1429809213083833e-06,
        2.9901545509984408e-06,
        4.861187730584302e-06,
        6.912547704001333e-06,
        1.0107220740508835e-05,
        1.3913359148836212e-05,
        1.7907432221431677e-05,
        2.3099748765459737e-05,
        2.5411378209153507e-05,
        2.9997358800358413e-05,
        2.6826896414277985e-05,
        2.840978870558551e-05,
        1.868645607624754e-05,
        1.7462608277408972e-05,
        6.472932895590083e-06,
        5.314218210170195e-06,
        3.5719291877840033e-06,
        3.0218687322870547e-06,
        5.4261426479473155e-06,
        2.182517456661895e-06,
    ]
    apriori_sigmas[44] = 33.928849588331296
    inverse_apriori_covariance = np.diag(
        [1.0 / sigma**2 if sigma > 0 else 0.0 for sigma in apriori_sigmas]
    )

    # The a priori values are the nominal values. The estimation starts from J2, J4 and J6 increased
    # by their a priori sigma (C20, C40, C60 decreased), and with J3, J5 and J7-J12 set to zero
    apriori_parameters = parameters.parameter_vector
    initial_parameters = apriori_parameters.copy()
    initial_parameters[[13, 15, 17]] -= apriori_sigmas[[13, 15, 17]]
    initial_parameters[[14, 16, 18, 19, 20, 21, 22, 23]] = 0.0

    return parameters, initial_parameters, apriori_parameters, inverse_apriori_covariance


def estimate(
    estimator,
    parameters,
    compressed_observations,
    initial_parameters,
    apriori_parameters,
    inverse_apriori_covariance,
):
    parameters.parameter_vector = initial_parameters

    estimation_input = estimation_analysis.EstimationInput(
        compressed_observations,
        inverse_apriori_covariance=inverse_apriori_covariance,
        convergence_checker=estimation_analysis.estimation_convergence_checker(
            maximum_iterations=6
        ),
        apply_final_parameter_correction=False,
        apply_apriori_parameter_deviation=True,
    )
    estimation_output = estimator.perform_estimation(
        estimation_input, apriori_parameter_vector=apriori_parameters
    )

    # Store the residuals of the best iteration in the collection, for plotting
    compressed_observations.set_residuals(estimation_output.final_residuals)

    estimated_parameters = estimation_output.final_parameters
    formal_errors = estimation_output.formal_errors
    for degree in range(2, 13):
        # unnormalized J_n = -sqrt(2n+1) C_n0
        index = 11 + degree
        print(
            f"J{degree} = {-np.sqrt(2 * degree + 1) * estimated_parameters[index] * 1e6:10.3f} "
            f"+/- {np.sqrt(2 * degree + 1) * formal_errors[index] * 1e6:.3f} (x 1e-6)"
        )
    print(f"k22 = {estimated_parameters[44]:.3f} +/- {formal_errors[44]:.3f}")

    return estimation_output


# ## Plot Results
#
# The post-fit residuals are plotted per perijove, in mHz and in mm/s.

# Observations with abs(residual) above this threshold [Hz] are removed as outliers
residual_threshold = 0.5


def get_observation_set_bands(obs_set):
    return tuple(
        FrequencyBands(int(band))
        for band in obs_set.ancillary_settings.get_float_list_settings(
            ObservationAncillarySimulationVariable.frequency_bands
        )
    )


def doppler_hz_to_m_per_s_scale(obs_set):
    """Hz -> two-way range-rate [m/s] for DSN N-way averaged Doppler.

    obs ≈ M_ref f_ref - M f_T (1 - rho_dot_tot / c)  (Moyer 2000, eq. 13-54, DSN sign convention)
    => d(obs)/d(rho_dot_tot) = M f_T / c
    """
    # approximation: replace with the ramp-table transmit frequency if f_ref is not the sky frequency
    f_T = obs_set.ancillary_settings.get_float_settings(
        ObservationAncillarySimulationVariable.doppler_reference_frequency
    )
    turnaround_ratio = juno_utils.juno_turnaround_ratio()[get_observation_set_bands(obs_set)]
    return constants.SPEED_OF_LIGHT / (turnaround_ratio * f_T)


def filter_residuals(observation_collection, threshold):
    # create a new observation collection to keep all original observations
    filtered_observations = observations.create_new_observation_collection(observation_collection)
    filtered_observations.filter_observations(
        observations.observations_processing.observation_filter(
            observations.observations_processing.ObservationFilterType.residual_filtering,
            threshold,
        )
    )
    filtered_observations.remove_empty_observation_sets()
    return filtered_observations


def get_pj_residuals(compressed_observations, pj, band):
    """Time since perijove [h] (spacecraft time) and residuals [mHz] and equivalent one-way [mm/s]."""
    hours, residuals_mhz, residuals_mms = [], [], []
    for obs_set in compressed_observations.get_single_observation_sets(pj_time_parsers[pj]):
        if get_observation_set_bands(obs_set) != (band, band):
            continue
        for epoch, residual in zip(obs_set.observation_times, np.ravel(obs_set.residuals)):
            # receive time minus downlink light time
            light_time = (
                np.linalg.norm(
                    spice.get_body_cartesian_position_at_epoch(
                        sc_name, "Earth", "J2000", "NONE", epoch.to_float()
                    )
                )
                / constants.SPEED_OF_LIGHT
            )
            hours.append(((epoch - pj_to_epoch[pj]).to_float() - light_time) / 3600.0)
            residuals_mhz.append(residual * 1e3)
            residuals_mms.append(residual * doppler_hz_to_m_per_s_scale(obs_set) / 2.0 * 1e3)

    return np.array(hours), np.array(residuals_mhz), np.array(residuals_mms)


def plot_results(x_observations, combined_observations):
    for pj in pjs:
        fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)

        for observation_collection, band, label, color, marker in (
            (x_observations, FrequencyBands.x_band, "X-band only", "tab:gray", "x"),
            (
                combined_observations,
                FrequencyBands.x_band,
                "Combined Ka- and X-band: X-band",
                "tab:orange",
                ".",
            ),
            (
                combined_observations,
                FrequencyBands.ka_band,
                "Combined Ka- and X-band: Ka-band",
                "tab:blue",
                ".",
            ),
        ):
            hours, residuals_mhz, residuals_mms = get_pj_residuals(observation_collection, pj, band)
            for ax, residuals, unit in zip(axes, (residuals_mhz, residuals_mms), ("mHz", "mm/s")):
                rms = np.sqrt(np.mean(residuals**2))
                ax.plot(
                    hours,
                    residuals,
                    linestyle="none",
                    marker=marker,
                    markersize=3,
                    color=color,
                    label=f"{label} (RMS: {rms:.3g} {unit})",
                )

        for ax in axes:
            ax.axvline(0, color="0.4", linewidth=0.8)
            ax.grid(alpha=0.25)
            ax.legend(fontsize=8, markerscale=2)
        axes[0].set_ylabel("Residual [mHz]")
        axes[1].set_ylabel("Residual (equivalent one-way)\n[mm/s]")
        axes[1].set_xlabel("Time since perijove [h]")
        fig.suptitle(f"PJ{pj:02d}")
        fig.tight_layout()

        output_directory = (
            mission_output_directory / "awvr" if use_awvr else mission_output_directory
        )
        file_path = output_directory / f"pj{pj}_postfit_residuals.pdf"
        file_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(file_path)
        plt.close(fig)


def main():
    load_kernels()

    tracking_data, tracking_supplementary_data = setup_observations()
    bodies = setup_bodies(tracking_data, tracking_supplementary_data)

    awvr_zenith_delays = juno_utils.load_awvr_zenith_delays(
        mission_file_directory / "amc", awvr_unit
    )

    x_observations = create_observation_collection(
        tracking_data, bodies, "X-band", awvr_zenith_delays
    )
    ka_observations = create_observation_collection(
        tracking_data, bodies, "Ka-band", awvr_zenith_delays
    )

    # Remove gross outliers, from the pre-fit residuals w.r.t. the SPICE (reconstructed) orbit
    observable_model_settings = setup_observation_models(
        observations.merge_observation_collections([x_observations, ka_observations]),
        awvr_zenith_delays,
    )
    observation_simulators = (
        observations_setup.observations_simulation_settings.create_observation_simulators(
            observable_model_settings, bodies
        )
    )
    observations.compute_residuals_and_dependent_variables(
        x_observations, observation_simulators, bodies
    )
    observations.compute_residuals_and_dependent_variables(
        ka_observations, observation_simulators, bodies
    )
    x_observations = filter_residuals(x_observations, residual_threshold)
    ka_observations = filter_residuals(ka_observations, residual_threshold)

    for pj in pjs:
        for band, observation_collection in (("X", x_observations), ("Ka", ka_observations)):
            for obs_set in observation_collection.get_single_observation_sets(pj_time_parsers[pj]):
                link_ends = obs_set.link_definition.link_ends
                print(
                    f"PJ{pj} {band}-band {link_ends[LinkEndType.transmitter].reference_point}"
                    f" -> {link_ends[LinkEndType.receiver].reference_point}: "
                    f"{len(obs_set.observation_times)} observations"
                )

    # ## Estimation
    #
    # Two fits from the same initial parameters: X-band only, and X- and Ka-band combined

    propagator_settings = setup_propagator(bodies)
    parameters, initial_parameters, apriori_parameters, inverse_apriori_covariance = (
        setup_parameters(bodies, propagator_settings)
    )

    # Copy of the X-band observations, so that the X-band only residuals are not overwritten
    combined_observations = observations.merge_observation_collections(
        [observations.create_new_observation_collection(x_observations), ka_observations]
    )
    for observation_collection in (x_observations, combined_observations):
        observation_collection.set_constant_weight(1.0 / doppler_noise**2)

    estimator = estimation_analysis.Estimator(
        bodies, parameters, observable_model_settings, propagator_settings
    )

    print("\nX-band only estimation")
    estimate(
        estimator,
        parameters,
        x_observations,
        initial_parameters,
        apriori_parameters,
        inverse_apriori_covariance,
    )

    print("\nCombined Ka- and X-band estimation")
    estimate(
        estimator,
        parameters,
        combined_observations,
        initial_parameters,
        apriori_parameters,
        inverse_apriori_covariance,
    )

    plot_results(x_observations, combined_observations)


if __name__ == "__main__":
    main()
