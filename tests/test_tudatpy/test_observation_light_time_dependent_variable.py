import numpy as np
import pytest

from tudatpy import constants
from tudatpy.astro import element_conversion
from tudatpy.interface import spice
from tudatpy.dynamics import environment_setup
from tudatpy.estimation.observable_models_setup import links, model_settings
from tudatpy.estimation.observations_setup import (
    ancillary_settings,
    observations_dependent_variables,
    observations_simulation_settings,
)
from tudatpy.estimation.observations_setup.observations_wrapper import simulate_observations


def _create_bodies():
    spice.load_standard_kernels()
    body_settings = environment_setup.get_default_body_settings(["Earth", "Moon", "Sun"], "Earth")
    body_settings.add_empty_settings("MoonOrbiter")
    kepler_elements = np.zeros(6)
    kepler_elements[0] = 2.0e6
    kepler_elements[1] = 0.1
    kepler_elements[2] = 1.0
    body_settings.get("MoonOrbiter").ephemeris_settings = environment_setup.ephemeris.keplerian(
        kepler_elements,
        0.0,
        spice.get_body_gravitational_parameter("Moon"),
        "Moon",
    )
    body_settings.get("Earth").ground_station_settings = [
        environment_setup.ground_station.basic_station(
            "Station1",
            [0.0, 0.35, 0.0],
            element_conversion.geodetic_position_type,
        ),
        environment_setup.ground_station.basic_station(
            "Station2",
            [0.2, 0.45, 0.0],
            element_conversion.geodetic_position_type,
        ),
    ]
    return environment_setup.create_system_of_bodies(body_settings)


def test_one_way_light_time_dependent_variable():
    bodies = _create_bodies()
    link_ends = {
        links.transmitter: links.body_reference_point_link_end_id("Earth", "Station1"),
        links.receiver: links.body_origin_link_end_id("MoonOrbiter"),
    }
    link_definition = links.LinkDefinition(link_ends)
    observation_time = 5.0e8 + 1000.0
    light_time_settings = (
        observations_dependent_variables.light_time_between_link_ends_dependent_variable(
            links.transmitter, links.receiver
        )
    )

    observation_settings = [model_settings.one_way_range(link_definition)]
    observation_simulators = observations_simulation_settings.create_observation_simulators(
        observation_settings, bodies
    )
    simulation_settings = [
        observations_simulation_settings.tabulated_simulation_settings(
            model_settings.one_way_range_type,
            link_definition,
            [observation_time],
        )
    ]
    observations_dependent_variables.add_dependent_variables_to_observable(
        simulation_settings,
        [light_time_settings],
        bodies,
        model_settings.one_way_range_type,
    )
    collection = simulate_observations(simulation_settings, observation_simulators, bodies)
    history = collection.dependent_variable_history(light_time_settings)
    assert len(history) == 1
    computed_light_time = float(list(history.values())[0][0])
    observed_range = float(collection.concatenated_observations[0])
    np.testing.assert_allclose(
        computed_light_time, observed_range / constants.SPEED_OF_LIGHT, rtol=1.0e-13, atol=0.0
    )


def test_two_way_light_time_excludes_retransmission_delay():
    # Two-way range elapsed time includes the spacecraft delay; this variable must not.
    bodies = _create_bodies()
    link_ends = {
        links.transmitter: links.body_reference_point_link_end_id("Earth", "Station1"),
        links.retransmitter: links.body_origin_link_end_id("MoonOrbiter"),
        links.receiver: links.body_reference_point_link_end_id("Earth", "Station1"),
    }
    link_definition = links.LinkDefinition(link_ends)
    observation_time = 5.0e8 + 1000.0
    retransmission_delay = 1.0e-3
    uplink_settings = (
        observations_dependent_variables.light_time_between_link_ends_dependent_variable(
            links.transmitter, links.retransmitter
        )
    )
    downlink_settings = (
        observations_dependent_variables.light_time_between_link_ends_dependent_variable(
            links.retransmitter, links.receiver
        )
    )
    combined_settings = (
        observations_dependent_variables.light_time_between_link_ends_dependent_variable(
            links.transmitter, links.receiver
        )
    )
    epoch_settings = observations_dependent_variables.link_end_epochs_dependent_variable()

    observation_settings = [model_settings.n_way_range(link_definition)]
    observation_simulators = observations_simulation_settings.create_observation_simulators(
        observation_settings, bodies
    )
    simulation_settings = [
        observations_simulation_settings.tabulated_simulation_settings(
            model_settings.n_way_range_type,
            link_definition,
            [observation_time],
            ancillary_settings=ancillary_settings.two_way_range_ancillary_settings(
                retransmission_delay
            ),
        )
    ]
    observations_dependent_variables.add_dependent_variables_to_observable(
        simulation_settings,
        [uplink_settings, downlink_settings, combined_settings, epoch_settings],
        bodies,
        model_settings.n_way_range_type,
    )
    collection = simulate_observations(simulation_settings, observation_simulators, bodies)

    uplink = float(list(collection.dependent_variable_history(uplink_settings).values())[0][0])
    downlink = float(list(collection.dependent_variable_history(downlink_settings).values())[0][0])
    combined = float(list(collection.dependent_variable_history(combined_settings).values())[0][0])
    epochs = np.array(list(collection.dependent_variable_history(epoch_settings).values())[0])
    observed_range = float(collection.concatenated_observations[0])

    np.testing.assert_allclose(
        combined,
        observed_range / constants.SPEED_OF_LIGHT - retransmission_delay,
        rtol=1.0e-13,
        atol=0.0,
    )
    np.testing.assert_allclose(combined, uplink + downlink, rtol=1.0e-13, atol=0.0)
    np.testing.assert_allclose(epochs[2] - epochs[1], retransmission_delay, rtol=0.0, atol=1.0e-6)
    assert abs((epochs[3] - epochs[0]) - combined) > 0.5 * retransmission_delay


def test_unspecified_light_time_selects_only_the_full_path():
    bodies = _create_bodies()
    link_ends = {
        links.transmitter: links.body_reference_point_link_end_id("Earth", "Station1"),
        links.retransmitter: links.body_origin_link_end_id("MoonOrbiter"),
        links.receiver: links.body_reference_point_link_end_id("Earth", "Station1"),
    }
    link_definition = links.LinkDefinition(link_ends)
    observation_time = 5.0e8 + 1000.0
    default_settings = (
        observations_dependent_variables.light_time_between_link_ends_dependent_variable()
    )
    explicit_settings = (
        observations_dependent_variables.light_time_between_link_ends_dependent_variable(
            links.transmitter, links.receiver
        )
    )
    observation_simulators = observations_simulation_settings.create_observation_simulators(
        [model_settings.n_way_range(link_definition)], bodies
    )
    simulation_settings = [
        observations_simulation_settings.tabulated_simulation_settings(
            model_settings.n_way_range_type,
            link_definition,
            [observation_time],
        )
    ]
    observations_dependent_variables.add_dependent_variables_to_observable(
        simulation_settings,
        [default_settings],
        bodies,
        model_settings.n_way_range_type,
    )
    assert len(simulation_settings[0].dependent_variable_settings_list) == 1

    collection = simulate_observations(simulation_settings, observation_simulators, bodies)
    default_light_time = float(
        list(collection.dependent_variable_history(default_settings).values())[0][0]
    )
    explicit_light_time = float(
        list(collection.dependent_variable_history(explicit_settings).values())[0][0]
    )
    np.testing.assert_allclose(default_light_time, explicit_light_time, rtol=1.0e-13, atol=0.0)


def test_station_specific_light_time_setting_skips_other_links_in_a_batch():
    bodies = _create_bodies()
    station1_link = links.LinkDefinition(
        {
            links.transmitter: links.body_reference_point_link_end_id("Earth", "Station1"),
            links.receiver: links.body_origin_link_end_id("MoonOrbiter"),
        }
    )
    station2_link = links.LinkDefinition(
        {
            links.transmitter: links.body_reference_point_link_end_id("Earth", "Station2"),
            links.receiver: links.body_origin_link_end_id("MoonOrbiter"),
        }
    )
    simulation_settings = [
        observations_simulation_settings.tabulated_simulation_settings(
            model_settings.one_way_range_type, station1_link, [5.0e8]
        ),
        observations_simulation_settings.tabulated_simulation_settings(
            model_settings.one_way_range_type, station2_link, [5.0e8]
        ),
    ]
    station1_light_time = (
        observations_dependent_variables.light_time_between_link_ends_dependent_variable(
            links.transmitter,
            links.receiver,
            links.body_reference_point_link_end_id("Earth", "Station1"),
            links.body_origin_link_end_id("MoonOrbiter"),
        )
    )

    observations_dependent_variables.add_dependent_variables_to_all(
        simulation_settings, [station1_light_time], bodies
    )

    assert len(simulation_settings[0].dependent_variable_settings_list) == 1
    assert len(simulation_settings[1].dependent_variable_settings_list) == 0


def test_reversed_light_time_selection_raises():
    bodies = _create_bodies()
    link_ends = {
        links.transmitter: links.body_reference_point_link_end_id("Earth", "Station1"),
        links.retransmitter: links.body_origin_link_end_id("MoonOrbiter"),
        links.receiver: links.body_reference_point_link_end_id("Earth", "Station1"),
    }
    link_definition = links.LinkDefinition(link_ends)
    simulation_settings = [
        observations_simulation_settings.tabulated_simulation_settings(
            model_settings.n_way_range_type,
            link_definition,
            [1.0e7],
        )
    ]
    with pytest.raises(RuntimeError):
        observations_dependent_variables.add_dependent_variables_to_observable(
            simulation_settings,
            [
                observations_dependent_variables.light_time_between_link_ends_dependent_variable(
                    links.receiver, links.transmitter
                )
            ],
            bodies,
            model_settings.n_way_range_type,
        )
