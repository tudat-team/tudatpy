"""Frozen Matera2022 CRD input and native SLRF2020 history for the shared LLR model."""

import numpy as np

from tudatpy import constants
from tudatpy.astro import time_representation as time
from tudatpy.data_input.environment_data import station_history as sinex
from tudatpy.data_input.tracking_data import slr
from tudatpy.dynamics import environment


def read_points(path):
    """Keep both calibrated transmit-tagged normal points and their measured weather."""

    points = []
    required = {
        "cdp_pad_id": 7941,
        "system_number": 77,
        "occupancy": 1,
        "station_epoch_time_scale": 4,
        "range_type": 2,
        "tropospheric_correction_applied": 0,
        "center_of_mass_correction_applied": 0,
        "receive_amplitude_correction_applied": 1,
        "station_system_delay_applied": 1,
        "spacecraft_system_delay_applied": 0,
        "data_quality_alert": 0,
    }

    for session in slr.read_slr_data([str(path)]):
        config = session.configuration
        if config.target_name.lower() != "apollo15" or any(
            getattr(config, key) != value for key, value in required.items()
        ):
            raise ValueError("Unexpected Matera target, station or CRD correction flags")
        if (config.start_year, config.start_month, config.start_day) != (2022, 8, 21):
            raise ValueError("Unexpected date in frozen Matera case")
        if (
            len(session.data.normal_point_data),
            len(session.data.meteorological_data),
            len(session.data.calibration_data),
        ) != (1, 1, 1):
            raise ValueError("Expected one point, measured weather and calibration per session")

        record = session.data.normal_point_data[0]
        meteo = session.data.meteorological_data[0]
        calibration = session.data.calibration_data[0]
        wavelength = config.system_configurations[
            record.system_configuration_id
        ].transmit_wavelength_nm
        if (
            record.epoch_event != 2
            or meteo.origin_of_values != 0
            or wavelength != 532.0
            or record.second_of_day != meteo.second_of_day
        ):
            raise ValueError("Expected transmit UTC tag, measured weather and 532 nm")
        if calibration.system_delay_ps != 108680.1:
            raise ValueError("Unexpected already-applied calibration metadata")

        midnight = time.DateTime(2022, 8, 21, 0, 0, 0.0).to_epoch_time_object()
        hour = int(record.second_of_day // 3600)
        minute = int(record.second_of_day % 3600 // 60)
        seconds = record.second_of_day % 60
        points.append(
            dict(
                index=len(points),
                station="MATERA",
                utc=midnight + record.second_of_day,
                utc_label=f"2022-08-21T{hour:02}:{minute:02}:{seconds:015.12f}",
                observed_rtt_s=record.two_way_time_of_flight,
                reflector=3,
                photons=record.number_of_returns,
                sigma_rtt_s=None,
                quality=str(config.data_quality_alert),
                pressure_pa=meteo.pressure * 100,
                temperature_k=meteo.temperature,
                humidity=meteo.humidity * 0.01,
                wavelength_m=wavelength * 1e-9,
                duration_s=record.normal_point_window_length,
                bin_rms_ps=record.bin_rms,
                epoch_event=record.epoch_event,
                station_system_delay_applied=config.station_system_delay_applied,
                system_delay_ps=calibration.system_delay_ps,
            )
        )

    if len(points) != 2:
        raise ValueError("The frozen Matera case must retain both normal points")
    return points


def station_history(data_dir, dataset):
    """Bind the parsed solution to native nominal-civil-UTC station motion."""

    products = {
        name: sinex.read_sinex_station_history(str(data_dir / record["file"]))
        for name, record in dataset["station_inputs"].items()
    }
    selected = [
        key
        for key in sinex.get_station_coordinate_solution_keys(products["solution"])
        if key.site.site_code == "7941" and key.site.point_code == "A"
    ]
    if len(selected) != 1 or selected[0].solution_id != "1":
        raise ValueError("Expected one Matera A/1 station solution")
    key = selected[0]

    def matches(candidate):
        """Identify the frozen Matera A/1 solution without relying on object identity."""

        return (candidate.site.site_code, candidate.site.point_code, candidate.solution_id) == (
            "7941",
            "A",
            "1",
        )

    source = sinex.get_station_coordinate_solution(products["solution"], key, True)
    if source.site.domes_id != "12734S008" or any(
        epoch.original != source.position_reference_epoch.original
        for epoch in source.velocity_reference_epochs
    ):
        raise ValueError("Unexpected Matera monument or station velocity reference")

    # The official companion files explicitly establish that no Matera PSD/discontinuity term is supplied.
    if any(p.key.site.site_code == "7941" for p in products["post_seismic"].parameters.values()):
        raise ValueError("Matera post-seismic terms require a declared nonzero history model")
    if any(d.key.site.site_code == "7941" for d in products["discontinuities"].discontinuities):
        raise ValueError("Matera discontinuities require explicit solution interval handling")

    spans = [span for span in products["solution"].observation_spans if matches(span.key)]
    if len(spans) != 1:
        raise ValueError("Expected one source observation span for the frozen station solution")

    def tai(epoch):
        """Convert a finite SINEX civil UTC epoch to TAI, preserving open ends."""

        civil = epoch.civil_utc
        return (
            None
            if civil is None
            else environment.sinex_civil_utc_to_tai(
                civil.year, civil.day_of_year, civil.seconds_of_day
            )
        )

    validity = environment.StationHistoryInterval()
    validity.start_tai, validity.end_tai = tai(spans[0].start), tai(spans[0].end)
    eccentricities = [
        e
        for e in products["eccentricity"].eccentricities
        if matches(e.key) and e.basis == sinex.SinexEccentricityBasis.xyz
    ]
    first_tai = environment.sinex_civil_utc_to_tai(2022, 233, 0)
    last_tai = first_tai + 86400.0
    if len(eccentricities) != 1 or any(eccentricities[0].components_metres):
        raise ValueError("Expected the published zero XYZ optical-axis eccentricity")
    offset = eccentricities[0]
    if (tai(offset.start) is not None and tai(offset.start) > first_tai) or (
        tai(offset.end) is not None and tai(offset.end) < last_tai
    ):
        raise ValueError("Zero eccentricity does not cover the complete frozen day")

    civil = source.position_reference_epoch.civil_utc
    linear = environment.StationHistoryLinearSolution()
    linear.source, linear.validity = source, validity
    linear.reference_time = environment.sinex_civil_utc_to_tai(
        civil.year, civil.day_of_year, civil.seconds_of_day
    )
    linear.position_metres = np.asarray(source.position_metres)
    linear.velocity_metres_per_julian_year = (
        np.asarray(source.velocity_metres_per_second) * constants.JULIAN_YEAR
    )
    metadata = environment.StationHistoryMetadata()
    metadata.site, metadata.domes_id = key.site, source.site.domes_id
    metadata.frame, metadata.release = "ITRF2020", "SLRF2020 250205"
    metadata.origin, metadata.tide_system = "whole-Earth centre of mass", "tide-free"
    metadata.reference_point_kind = (
        environment.StationReferencePointKind.calibrated_instrument_point
    )
    metadata.reference_point_name = "MATERA optical axis intersection"
    history = environment.StationCoordinateHistory([linear], [], [], metadata)
    ledger = dict(
        dataset["station"],
        reference_utc=source.position_reference_epoch.original,
        position_m=list(source.position_metres),
        velocity_m_per_s=list(source.velocity_metres_per_second),
        zero_eccentricity_start=offset.start.original,
        zero_eccentricity_end=offset.end.original,
        post_seismic="No Matera terms in pinned official ITRF2020 SLR PSD product",
    )
    return history, ledger
