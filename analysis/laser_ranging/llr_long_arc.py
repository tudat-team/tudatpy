"""Reduce the frozen March-April 2017 arc with one explicitly requested source exclusion."""

import argparse
from collections import defaultdict
import csv
from datetime import date
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from tudatpy import constants, kernel
from tudatpy.astro import time_representation as time
from tudatpy.data_input.tracking_data import slr

import llr_residuals as model

# This interval and both source releases were selected before any long-arc residual.
MANIFEST = Path(__file__).with_name("llr_long_arc_dataset.json")


def is_requested_exclusion(station, source_file, source_line, utc_launch, manifest):
    """Match only the exact source identity excluded at Dominic's explicit request."""

    excluded = manifest["excluded_normal_point"]
    return (station, source_file, str(source_line), utc_launch) == (
        excluded["station"],
        excluded["source_file"],
        str(excluded["source_line"]),
        excluded["utc_launch"],
    )


def read_groups(root, manifest):
    """Keep all frozen rows, original identities and measured weather in bounded groups."""

    apollo = manifest["apollo"]

    points = model.read_points(
        root / apollo["file"],
        apollo["sha256"],
        (date(2017, 3, 1), date(2017, 5, 1)),
        apollo["source_count"],
    )

    nights = defaultdict(list)

    for point, line in zip(points, apollo["source_line_numbers"]):
        point.update(
            station="APOLLO",
            source_file=apollo["file"],
            source_line=line,
            source_id="70610",
            cdp="7045",
            domes="49447S001",
            source_session=point["utc_label"][:10],
            epoch_event=2,
            system_configuration="Group-e",
            bin_rms_ps=None,
        )
        nights[point["utc_label"][:10]].append(point)

    groups = [
        (values, None, apollo["station_profile"], {"station": "APOLLO", "group": night})
        for night, values in sorted(nights.items())
    ]

    grasse = manifest["grasse"]

    inventory_path = root / grasse["inventory"]["file"]

    if model.sha256(inventory_path) != grasse["inventory"]["sha256"]:
        raise ValueError("Grasse inventory differs from its frozen audit")
    with inventory_path.open(newline="") as stream:
        inventory = list(csv.DictReader(stream))

    source_rows = {
        (
            r["source_file"],
            r["utc_launch"][:10],
            float(r["seconds_of_day_printed"]),
            float(r["two_way_time_of_flight_s"]),
        ): r
        for r in inventory
    }

    seen = set()

    required = dict(
        cdp_pad_id=7845,
        system_number=78,
        occupancy=1,
        station_epoch_time_scale=4,
        range_type=2,
        data_quality_alert=0,
        tropospheric_correction_applied=0,
        center_of_mass_correction_applied=0,
        receive_amplitude_correction_applied=0,
        station_system_delay_applied=1,
        spacecraft_system_delay_applied=0,
    )

    for product in grasse["files"]:
        path = root / product["file"]

        if model.sha256(path) != product["sha256"]:
            raise ValueError(f"Grasse provider input changed: {path.name}")
        for session in slr.read_slr_data([str(path)]):
            config = session.configuration

            if any(getattr(config, key) != value for key, value in required.items()):
                raise ValueError("Grasse header differs from the frozen metadata audit")
            night = date(config.start_year, config.start_month, config.start_day).isoformat()

            midnight = time.DateTime(
                config.start_year, config.start_month, config.start_day, 0, 0, 0.0
            ).to_epoch_time_object()

            weather = {}

            for record in session.data.meteorological_data:
                if record.origin_of_values != 0:
                    raise ValueError("Expected measured Grasse meteorological records")
                epoch = float(midnight + record.second_of_day)

                value = [record.temperature, record.pressure * 100.0, record.humidity * 0.01]

                if epoch in weather and weather[epoch] != value:
                    raise ValueError("Conflicting same-epoch measured meteorology")
                weather[epoch] = value

            if not weather:
                raise ValueError("Measured Grasse weather is required")
            values = []

            for record in session.data.normal_point_data:
                key = (path.name, night, record.second_of_day, record.two_way_time_of_flight)

                source = source_rows[key]

                if key in seen or record.epoch_event != 2:
                    raise ValueError("Unexpected duplicate or non-transmit Grasse record")
                seen.add(key)
                utc = midnight + record.second_of_day

                nearest = min(abs(float(utc) - epoch) for epoch in weather)

                if float(utc) < min(weather) - 1000 or float(utc) + 5 > max(weather) + 1000:
                    raise ValueError(
                        "Grasse event exceeds declared measured-weather endpoint holding"
                    )
                wavelength = config.system_configurations[
                    record.system_configuration_id
                ].transmit_wavelength_nm

                if wavelength != float(source["wavelength_nm"]):
                    raise ValueError("Grasse wavelength differs from its printed configuration")
                values.append(
                    dict(
                        index=len(seen) - 1,
                        station="GRASSE",
                        source_file=path.name,
                        source_line=int(source["source_line"]),
                        source_id="GRSM",
                        cdp="7845",
                        domes="10002S002",
                        source_session=source["session_id"],
                        utc=utc,
                        utc_label=source["utc_launch"],
                        observed_rtt_s=record.two_way_time_of_flight,
                        reflector=int(source["reflector_code"]),
                        photons=record.number_of_returns,
                        sigma_rtt_s=None,
                        quality=str(config.data_quality_alert),
                        wavelength_m=wavelength * 1e-9,
                        duration_s=record.normal_point_window_length,
                        bin_rms_ps=record.bin_rms,
                        epoch_event=record.epoch_event,
                        system_configuration=record.system_configuration_id,
                        nearest_weather_s=nearest,
                        weather_mode=(
                            "held-first"
                            if float(utc) < min(weather)
                            else "held-last" if float(utc) > max(weather) else "interpolated"
                        ),
                    )
                )
            if not values:
                raise ValueError("The frozen Grasse session unexpectedly has no normal points")
            groups.append(
                (
                    sorted(values, key=lambda p: float(p["utc"])),
                    weather,
                    grasse["station_profile"],
                    dict(
                        station="GRASSE",
                        group=f"{path.name}:{values[0]['source_session']}",
                        flags={key: getattr(config, key) for key in required},
                        weather_count=len(weather),
                        weather_start_utc_s=min(weather),
                        weather_end_utc_s=max(weather),
                        weather_endpoint_hold_limit_s=1000,
                        calibration_records=len(session.data.calibration_data),
                        calibration_reapplied_m=0.0,
                    ),
                )
            )
    if len(seen) != grasse["source_count"] or len(inventory) != len(seen):
        raise ValueError("The long arc must retain all audited Grasse normal points")
    return groups


def geometry(point, bodies, diagnostics):
    """Describe solved line-of-sight and lunar illumination without fitting a correction."""

    transmit, bounce, _, receive = diagnostics["link_end_epochs_tdb_s"]

    moon = bodies.get("Moon")

    earth = bodies.get("Earth")

    rotation = moon.rotation_model.body_fixed_to_inertial_rotation(bounce)

    moon_centre = np.asarray(moon.state_in_base_frame_from_ephemeris(time.Time(bounce))).reshape(6)[
        :3
    ]

    reflector = model.REFLECTORS[point["reflector"]][0]

    local = np.asarray(
        moon.get_ground_station(reflector).station_state.get_cartesian_position(bounce, "Moon")
    ).reshape(3)

    target = moon_centre + rotation @ local

    directions = []

    for epoch in (transmit, receive):
        station = earth.get_ground_station(point["station"]).station_state

        body_fixed = np.asarray(station.get_cartesian_position(epoch, "Earth")).reshape(3)

        observer = (
            np.asarray(earth.state_in_base_frame_from_ephemeris(time.Time(epoch))).reshape(6)[:3]
            + earth.rotation_model.body_fixed_to_inertial_rotation(epoch) @ body_fixed
        )

        line = observer - target

        directions.append(rotation.T @ (line / np.linalg.norm(line)))
    earth_direction = (
        np.asarray(earth.state_in_base_frame_from_ephemeris(time.Time(bounce))).reshape(6)[:3]
        - moon_centre
    )

    sun_direction = (
        np.asarray(bodies.get("Sun").state_in_base_frame_from_ephemeris(time.Time(bounce))).reshape(
            6
        )[:3]
        - moon_centre
    )

    earth_direction /= np.linalg.norm(earth_direction)
    sun_direction /= np.linalg.norm(sun_direction)
    result = {
        "lunar_phase_angle_deg": float(
            np.rad2deg(np.arccos(np.clip(np.dot(earth_direction, sun_direction), -1, 1)))
        ),
        "transmit_elevation_deg": float(np.rad2deg(diagnostics["transmit_elevation_rad"][0])),
        "receive_elevation_deg": float(np.rad2deg(diagnostics["receive_elevation_rad"][0])),
    }

    for label, vector in [
        ("earth_direction_pa", rotation.T @ earth_direction),
        ("sun_direction_pa", rotation.T @ sun_direction),
        ("transmit_los_pa", directions[0]),
        ("receive_los_pa", directions[1]),
    ]:
        result.update({f"{label}_{axis}": float(value) for axis, value in zip("xyz", vector)})
    return result


def statistics(rows):
    """Report descriptive raw statistics without changing any residual."""

    values = np.array([row["one_way_residual_m"] for row in rows])

    return dict(
        count=len(rows),
        mean_m=float(np.mean(values)),
        rms_about_zero_m=float(np.sqrt(np.mean(values**2))),
        scatter_about_mean_m=float(np.std(values)),
    )


def write_csv(path, rows):
    """Write the declared row fields without residual editing."""

    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))

        writer.writeheader()
        writer.writerows(rows)


def event_vectors(point, bodies, diagnostic):
    """Reconstruct station vectors with BCRS spatial terms at the diagnostic event doubles."""

    rows = []
    reflector = model.REFLECTORS[point["reflector"]][0]

    for body_name, name, epoch in zip(
        ("Earth", "Moon", "Moon", "Earth"),
        (point["station"], reflector, reflector, point["station"]),
        diagnostic["link_end_epochs_tdb_s"],
    ):
        body = bodies.get(body_name)
        local = body.get_ground_station(name).station_state.get_cartesian_position(epoch, "SSB")
        rotation = body.rotation_model.body_fixed_to_inertial_rotation(epoch)
        rows.append(
            dict(
                body=body_name,
                station=name,
                epoch_tdb_s=epoch,
                origin_for_station_motion="SSB",
                body_fixed_position_m=np.asarray(local).ravel().tolist(),
                rotation_to_icrf=np.asarray(rotation).tolist(),
                body_states={
                    name: np.asarray(
                        bodies.get(name).state_in_base_frame_from_ephemeris(time.Time(epoch))
                    )
                    .ravel()
                    .tolist()
                    for name in ("Earth", "Moon", "Sun")
                },
            )
        )
        if body_name == "Moon":
            rows[-1]["intrinsic_reference_radius_position_m"] = (
                np.asarray(
                    body.get_ground_station(
                        name + "_intrinsic"
                    ).station_state.get_cartesian_position(epoch, "Moon")
                )
                .ravel()
                .tolist()
            )

    return dict(
        utc_launch=point["utc_label"],
        station=point["station"],
        reflector=reflector,
        launch_utc_split=point["utc"].__getstate__(),
        diagnostic=diagnostic,
        vectors=rows,
    )


def main():
    """Run only the frozen primary station pair, preserving all inputs and raw results."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--kernel-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--station", choices=["APOLLO", "GRASSE"])
    parser.add_argument("--march9-control", action="store_true")
    parser.add_argument("--ephemeris", choices=("de440", "inpop17a"), default="de440")
    parser.add_argument("--glass-convention", choices=("jpl-nh", "none"))
    parser.add_argument("--control-twelve", action="store_true")
    args = parser.parse_args()

    if args.ephemeris == "inpop17a" and args.glass_convention is None:
        parser.error("INPOP prism convention is unresolved: explicitly declare jpl-nh or none")
    if args.control_twelve and (args.march9_control or args.station is not None):
        parser.error("The predeclared twelve-point control requires both stations")

    started = perf_counter()

    manifest = json.loads(MANIFEST.read_text())

    ephemeris_path = Path(model.__file__).with_name("llr_dataset.json")

    ephemeris = json.loads(ephemeris_path.read_text())
    frame_alias = ephemeris["analysis_frame_alias"]

    if args.ephemeris == "inpop17a":
        ephemeris_path = Path(model.__file__).with_name("inpop17a_dataset.json")
        ephemeris = json.loads(ephemeris_path.read_text())

    kernel_hashes, _ = model.load_ephemeris_kernels(args.kernel_dir, ephemeris, frame_alias)
    pool_count = model.spice.get_total_count_of_kernels_loaded()

    if pool_count != len(kernel_hashes):
        raise RuntimeError("The cleared kernel pool differs from the declared load inventory")

    groups = read_groups(args.data_root, manifest)
    excluded_count = sum(
        is_requested_exclusion(
            point["station"],
            point["source_file"],
            point["source_line"],
            point["utc_label"],
            manifest,
        )
        for points, _, _, _ in groups
        for point in points
    )
    if excluded_count != 1:
        raise ValueError("The frozen source must contain the requested exclusion exactly once")

    if args.station:
        groups = [group for group in groups if group[3]["station"] == args.station]
    expected_count = sum(
        manifest[station.lower()]["retained_count"]
        for station in ("APOLLO", "GRASSE")
        if args.station is None or station == args.station
    )

    if args.march9_control:
        if args.station == "GRASSE":
            parser.error("The March 9 control contains APOLLO observations")
        groups = [
            group
            for group in groups
            if group[3]["station"] == "APOLLO" and group[3]["group"] == "2017-03-09"
        ]
        expected_count = 11

    if args.control_twelve:
        selected = []

        for points, weather, profile, ledger in groups:
            days = (
                ("2017-03-22", "2017-04-07")
                if ledger["station"] == "APOLLO"
                else ("2017-03-22", "2017-04-06")
            )
            if points[0]["utc_label"][:10] not in days:
                continue
            values = [
                next(point for point in points if point["reflector"] == code)
                for code in (0, 2, 3)
                if any(point["reflector"] == code for point in points)
            ]
            if values:
                selected.append(
                    (
                        sorted(values, key=lambda point: float(point["utc"])),
                        weather,
                        profile,
                        ledger,
                    )
                )
        groups = selected
        expected_count = 12

    args.output.mkdir(parents=True, exist_ok=False)
    rows, group_ledger, events = [], [], []

    for points, weather, profile, ledger in groups:
        setup_start = perf_counter()

        bodies, position, _ = model.create_environment(
            points,
            station_profile=profile,
            meteorological_history=weather,
            lunar_frame=ephemeris["reflector_coordinates"]["frame"],
            reflector_coordinates=ephemeris["reflector_coordinates"].get("coordinates_m"),
            gravitational_parameters=ephemeris.get("gravitational_parameters_m3_s2"),
            lunar_model=ephemeris.get("lunar_model", "de440"),
        )

        ledger["setup_seconds"] = perf_counter() - setup_start
        converter = bodies.get("Earth").rotation_model.angles_calculator.time_scale_converter

        solver_seconds = 0.0
        first_row = len(rows)

        for point in points:
            if is_requested_exclusion(
                point["station"],
                point["source_file"],
                point["source_line"],
                point["utc_label"],
                manifest,
            ):
                continue
            epoch = converter.convert_time_object(
                time.utc_scale, time.tdb_scale, point["utc"], position
            )

            correction = model.complete_light_time_corrections(
                point, explicit_glass=args.glass_convention != "none"
            )

            solver_start = perf_counter()

            computed, diagnostic = model.computed_range(point, epoch, bodies, correction)
            solver_seconds += perf_counter() - solver_start
            residual = 0.5 * (constants.SPEED_OF_LIGHT * point["observed_rtt_s"] - computed)

            if args.control_twelve:
                events.append(event_vectors(point, bodies, diagnostic))

            rows.append(
                dict(
                    utc_launch=point["utc_label"],
                    station=point["station"],
                    reflector=model.REFLECTORS[point["reflector"]][0],
                    one_way_residual_m=residual,
                    one_way_sigma_m=(
                        0.5 * constants.SPEED_OF_LIGHT * point["sigma_rtt_s"]
                        if point["sigma_rtt_s"] is not None
                        else None
                    ),
                    source_file=point["source_file"],
                    source_line=point["source_line"],
                    source_session=point["source_session"],
                    source_station_id=point["source_id"],
                    station_cdp=point["cdp"],
                    station_domes=point["domes"],
                    station_profile=profile,
                    quality=point["quality"],
                    epoch_event=point["epoch_event"],
                    system_configuration=point["system_configuration"],
                    wavelength_nm=point["wavelength_m"] * 1e9,
                    observed_rtt_s=point["observed_rtt_s"],
                    computed_full_two_way_range_m=computed,
                    photons=point["photons"],
                    bin_rms_ps_metadata_only=point["bin_rms_ps"],
                    nearest_measured_weather_s=point.get("nearest_weather_s", 0.0),
                    weather_mode=point.get("weather_mode", "measured-at-NP"),
                    station_system_delay_applied=(1 if point["station"] == "GRASSE" else None),
                    calibration_reapplied_m=0.0,
                    **geometry(point, bodies, diagnostic),
                )
            )
        ledger["solver_seconds"] = solver_seconds
        ledger["source_count"] = len(points)
        ledger["count"] = len(rows) - first_row
        group_ledger.append(ledger)
        print(json.dumps(ledger), flush=True)
    rows.sort(key=lambda row: (row["utc_launch"], row["station"], row["source_line"]))
    if len(rows) != expected_count:
        raise RuntimeError("The declared retained arc count differs from the explicit selection")
    write_csv(args.output / "residuals.csv", rows)
    nights = defaultdict(list)

    for row in rows:
        nights[(row["utc_launch"][:10], row["station"], row["reflector"])].append(row)
    nightly = [
        dict(
            utc_night=night,
            station=station,
            reflector=reflector,
            first_utc=values[0]["utc_launch"],
            last_utc=values[-1]["utc_launch"],
            **statistics(values),
        )
        for (night, station, reflector), values in sorted(nights.items())
    ]

    write_csv(args.output / "nightly_raw_statistics.csv", nightly)
    if args.control_twelve:
        (args.output / "event_diagnostics.json").write_text(json.dumps(events, indent=2) + "\n")
    summary = dict(
        count=len(rows),
        source_count=sum(g["source_count"] for g in group_ledger),
        excluded_count=sum(g["source_count"] - g["count"] for g in group_ledger),
        excluded_normal_point=manifest["excluded_normal_point"],
        station_selection=args.station or "APOLLO and GRASSE",
        station_statistics={
            station: statistics([r for r in rows if r["station"] == station])
            for station in sorted({r["station"] for r in rows})
        },
        station_reflector_statistics={
            f"{station}/{reflector}": statistics(values)
            for (station, reflector), values in (
                (key, [r for r in rows if (r["station"], r["reflector"]) == key])
                for key in sorted({(r["station"], r["reflector"]) for r in rows})
            )
        },
        manifest_sha256=model.sha256(MANIFEST),
        ephemeris_manifest_sha256=model.sha256(ephemeris_path),
        kernel_binary_sha256=model.sha256(Path(kernel.__file__)),
        kernel_sha256=kernel_hashes,
        kernel_pool=dict(load_order=list(kernel_hashes), actual_native_count=pool_count),
        source_sha256={
            path.name: model.sha256(path)
            for path in [Path(__file__), Path(model.__file__), model.STATION_PROFILE_PATH]
        },
        earth_orientation_resource_sha256={
            path.name: model.sha256(path)
            for path in sorted(
                Path(model.resource_paths.get_earth_orientation_path()).glob("*.txt")
            )
        },
        station_profiles={name: model.STATION_PROFILES[name] for name in {g[2] for g in groups}},
        groups=group_ledger,
        setup_seconds=sum(g["setup_seconds"] for g in group_ledger),
        solver_seconds=sum(g["solver_seconds"] for g in group_ledger),
        total_seconds=perf_counter() - started,
        physics=f"{args.ephemeris}/ICRF, 30 s Time cubic tables, no loading, 50 iterations and fractional tolerance 1e-13",
        explicit_glass_convention=args.glass_convention or "jpl-nh",
        reduction_limit=ephemeris.get("optical_limit", "Published JPL front-face n*h convention"),
        control_selection=(
            "First chronological A11/A14/A15 on March22 and APOLLOApr7/GrasseApr6"
            if args.control_twelve
            else None
        ),
        fitting="No fitted bias, demeaning, clipping, wavelength correction or empirical libration term",
        geometry_precision="Geometry uses native diagnostic event doubles; solver retains Time internally",
        uncertainty="APOLLO supplied NP RTT uncertainty only; Grasse bin RMS remains metadata",
    )

    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(
        json.dumps(
            {
                key: summary[key]
                for key in [
                    "count",
                    "station_statistics",
                    "setup_seconds",
                    "solver_seconds",
                    "total_seconds",
                ]
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
