"""Instantaneous Doppler (doppler_measured_frequency) with the tabulated ionospheric light-time correction.

The two-way model called its one-way legs without the ancillary settings that carry the link frequencies, and the
legs computed their light-time partials with a null pointer. Any ionospheric correction (which scales with
frequency) then failed: on 1.0.0 with an access violation, on develop already at model creation ("transmitted
frequency not defined" for the spacecraft, requested by the downlink leg).

Reference: a TRK-2-23 CHPART card whose S-band range delay grows linearly over the day,
rho(t) = c0 + c1 * x, x = 2 (t - t_start) / (t_end - t_start) - 1, so drho/dt = 2 c1 / 86400 s is constant and the
ionospheric part of the received frequency is exact:
    df = f_rx / c * drho/dt * [(f_ref / f_up)^2 + (f_ref / f_down)^2],   f_ref = 2295 MHz,
(phase advance, so the observed frequency increases while the group delay grows).
"""

import numpy as np
import pytest

from tudatpy.astro.time_representation import DateTime, Time
from tudatpy.dynamics import environment_setup
from tudatpy.dynamics.environment import PiecewiseLinearFrequencyInterpolator
from tudatpy.estimation import observable_models_setup, observations_setup
from tudatpy.interface import spice

C = 299792458.0
F_REF = 2295.0e6
F_UP = 7.153e9
C0, C1 = 100.0, 1000.0  # S-band delay [m]
links = observable_models_setup.links
models = observable_models_setup.model_settings
corrections = observable_models_setup.light_time_corrections
ancillary = observations_setup.ancillary_settings
simulation = observations_setup.observations_simulation_settings
wrapper = observations_setup.observations_wrapper
bands = ancillary.FrequencyBands


@pytest.mark.parametrize(
    "transmitter, receiver, complex_id, downlink_band, turnaround",
    [
        ("DSS-55", "DSS-55", "C60", bands.x_band, 880.0 / 749.0),
        ("DSS-14", "DSS-25", "C10", bands.ka_band, 3344.0 / 749.0),
    ],
)
def test_instantaneous_doppler_with_tabulated_ionosphere(
    tmp_path, transmitter, receiver, complex_id, downlink_band, turnaround
):
    card = tmp_path / f"ion_{complex_id}.csp"
    card.write_text(
        f"ADJUST(DOPRNG)BY NRMPOW({C0}, {C1})\n"
        f"MODEL(CHPART)FROM(16/08/27,00:00)TO(16/08/28,00:00)DSN({complex_id})SCID(61).\n"
    )

    spice.load_standard_kernels()
    epochs = spice.convert_date_string_to_ephemeris_time("2016-08-27 06:00:00 UTC") + np.arange(
        0.0, 12 * 3600.0, 3600.0
    )
    settings = environment_setup.get_default_body_settings(["Sun", "Earth", "Moon"], "SSB", "J2000")
    settings.get("Earth").ground_station_settings = environment_setup.ground_station.dsn_stations()
    settings.add_empty_settings("Probe")
    settings.get("Probe").ephemeris_settings = environment_setup.ephemeris.keplerian(
        [5.2 * 1.495978707e11, 0.05, 0.02, 0.3, 0.4, 1.0],
        epochs[0],
        1.32712440018e20,
        "Sun",
        "J2000",
    )
    bodies = environment_setup.create_system_of_bodies(settings)
    utc0 = DateTime(2016, 8, 26, 0, 0, 0.0).to_epoch()
    bodies.get("Earth").get_ground_station(transmitter).set_transmitting_frequency_calculator(
        PiecewiseLinearFrequencyInterpolator(
            [Time(utc0)], [Time(utc0 + 3 * 86400.0)], [0.0], [F_UP]
        )
    )
    bodies.get("Probe").system_models.set_default_transponder_turnaround_ratio_function()

    link = links.LinkDefinition(
        {
            links.transmitter: links.body_reference_point_link_end_id("Earth", transmitter),
            links.reflector1: links.body_origin_link_end_id("Probe"),
            links.receiver: links.body_reference_point_link_end_id("Earth", receiver),
        }
    )
    ionosphere = [
        corrections.dsn_tabulated_ionospheric_light_time_correction([str(card)], {61: "Probe"})
    ]

    def received_frequency(light_time_corrections):
        model = models.doppler_measured_frequency(link, light_time_corrections)
        sim_settings = simulation.tabulated_simulation_settings(
            models.ObservableType.doppler_measured_frequency_type,
            link,
            [Time(float(t)) for t in epochs],
            ancillary_settings=ancillary.doppler_measured_frequency_ancillary_settings(
                [bands.x_band, downlink_band]
            ),
        )
        simulators = simulation.create_observation_simulators([model], bodies)
        observations = wrapper.simulate_observations([sim_settings], simulators, bodies)
        return np.asarray(observations.concatenated_observations, float).ravel()

    ionospheric_doppler = received_frequency(ionosphere) - received_frequency([])

    f_down = F_UP * turnaround
    rate = 2.0 * C1 / 86400.0
    expected = f_down / C * rate * ((F_REF / F_UP) ** 2 + (F_REF / f_down) ** 2)
    # expected: 115.2 mHz (X/X), 266.6 mHz (X/Ka); limit 0.5 mHz
    assert np.max(np.abs(ionospheric_doppler - expected)) < 0.5e-3
