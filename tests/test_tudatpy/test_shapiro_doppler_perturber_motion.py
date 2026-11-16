"""Instantaneous one-way Doppler with a first-order relativistic (Shapiro) light-time correction.

The Shapiro part of the Doppler must equal minus the time derivative of the Shapiro part of the light time,
including the motion of the perturbing body during the observation. Geometry follows unitTestDopplerModels.cpp:
ground station on Earth -> Mars, Shapiro body Sun (GM x 1e4) or Earth (moving at ~30 km/s w.r.t. the SSB).
Differencing with/without the correction removes the geometric Doppler, leaving a clean Shapiro-only comparison.
"""

import numpy as np
import pytest

from tudatpy.astro.time_representation import Time
from tudatpy.dynamics import environment_setup
from tudatpy.estimation import observable_models_setup, observations_setup
from tudatpy.interface import spice

C = 299792458.0
links = observable_models_setup.links
models = observable_models_setup.model_settings
corrections = observable_models_setup.light_time_corrections
simulation = observations_setup.observations_simulation_settings
wrapper = observations_setup.observations_wrapper


@pytest.mark.parametrize("perturber, gm_scaling", [("Sun", 1.0e4), ("Earth", 1.0)])
def test_shapiro_doppler_includes_perturber_motion(perturber, gm_scaling):
    spice.load_standard_kernels()
    settings = environment_setup.get_default_body_settings(
        ["Earth", "Sun", "Mars"], "SSB", "ECLIPJ2000"
    )
    settings.get(perturber).gravity_field_settings = environment_setup.gravity_field.central(
        spice.get_body_gravitational_parameter(perturber) * gm_scaling
    )
    bodies = environment_setup.create_system_of_bodies(settings)
    environment_setup.add_ground_station(
        bodies.get("Earth"), "Station1", np.array([-1917032.190, -6029782.349, -801376.113])
    )
    link = links.LinkDefinition(
        {
            links.transmitter: links.body_reference_point_link_end_id("Earth", "Station1"),
            links.receiver: links.body_origin_link_end_id("Mars"),
        }
    )

    def simulate(model, observable_type, times):
        sim_settings = simulation.tabulated_simulation_settings(
            observable_type, link, [Time(float(t)) for t in times], links.receiver
        )
        simulators = simulation.create_observation_simulators([model], bodies)
        observations = wrapper.simulate_observations([sim_settings], simulators, bodies)
        return np.asarray(observations.concatenated_observations, float).ravel()

    shapiro = [corrections.first_order_relativistic_light_time_correction([perturber])]
    times = 3.5 * 86400.0 + np.arange(0.0, 86400.0, 7200.0)

    def doppler(corr):  # model output in m/s
        model = models.one_way_doppler_instantaneous(link, corr)
        return simulate(model, models.ObservableType.one_way_instantaneous_doppler_type, times) / C

    def light_time(corr, t):
        return (
            simulate(models.one_way_range(link, corr), models.ObservableType.one_way_range_type, t)
            / C
        )

    shapiro_doppler = doppler(shapiro) - doppler([])
    h = 25.0
    shapiro_light_time_rate = (
        (light_time(shapiro, times + h) - light_time([], times + h))
        - (light_time(shapiro, times - h) - light_time([], times - h))
    ) / (2.0 * h)

    # Without the perturber's motion the error is 6.8e-12 (Sun x 1e4) and 2.6e-13 (Earth); the numerical
    # derivative of light-time differences is good to ~1e-14
    assert np.max(np.abs(shapiro_doppler + shapiro_light_time_rate)) < 3.0e-14
