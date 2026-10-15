"""Reduce the pinned April 2024 LAGEOS-1 normal points against a fixed SP3 orbit.

No propagation, estimation, residual editing or fitted bias is performed. See the
adjacent README for input provenance, coordinate conventions and limitations.
"""

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys
import tarfile

import numpy as np
import tudatpy
from tudatpy import kernel

from tudatpy import constants
from tudatpy.astro import time_representation as time
from tudatpy.data_input import environment_data, resource_paths
from tudatpy.data_input.environment_data import earth_orientation, ilrs, station_history
from tudatpy.data_input.tracking_data import slr
from tudatpy.dynamics import environment, environment_setup
from tudatpy.dynamics.environment_setup import ephemeris, ground_station
from tudatpy.dynamics.environment_setup import rotation_model, shape_deformation
from tudatpy.estimation import observable_models_setup, observations, observations_setup
from tudatpy.interface import spice
from tudatpy.math import interpolators

# Fixed selections precede any residual calculation; values are not fitted here.
STATIONS = {"7839": "11001S002", "7840": "13212S001"}
SPATIAL_TT_TO_TDB_SCALE = 1.48082686741e-8
C = constants.SPEED_OF_LIGHT


def file_record(path):
    """Identify an actual local input by its bytes, rather than its filename alone."""
    path = Path(path)
    digest = hashlib.sha256()

    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)

    return {
        "path": str(path.resolve()),
        "bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
    }


def unpack_input(path, directory):
    """Decompress one provider gzip without changing the original download."""
    path = Path(path)

    if path.suffix != ".gz":
        return path

    target = directory / path.stem

    with gzip.open(path, "rb") as source, target.open("wb") as destination:
        shutil.copyfileobj(source, destination)

    return target


def calendar_epoch(year, month, day):
    """Construct midnight in a declared time coordinate without a float epoch."""
    return time.date_time_components_to_epoch_time_object(year, month, day, 0, 0, 0.0)


def sinex_bound(token, is_end):
    """Interpret the DHF's zero sentinel as an open interval endpoint."""
    if token == "00:000:00000":
        return math.inf if is_end else -math.inf

    return ilrs.convert_sinex_datetime_to_seconds_since_epoch(token)


def read_corrections(data_root):
    """Read the pinned LAGEOS-1 CoM release and the dated ILRS data-handling file."""
    com_path = data_root / "metadata/com6.220915.tar.gz"
    dhf_path = data_root / "metadata/ILRS_Data_Handling_File_2024.02.13.snx"
    com_rows = []
    dhf_rows = []

    with tarfile.open(com_path, "r:gz") as archive:
        text = archive.extractfile("com_lg1.dat").read().decode("ascii")

    for number, line in enumerate(text.splitlines(), 1):
        fields = line.split()

        if not fields or fields[0] not in STATIONS:
            continue

        start = calendar_epoch(int(fields[3]), int(fields[2]), int(fields[1]))
        # com_6multi.for compares integer MJDs inclusively, including the full printed end date.
        end = calendar_epoch(int(fields[6]), int(fields[5]), int(fields[4])) + 86400.0
        com_rows.append((fields[0], start, end, int(fields[7]), float(fields[8]) * 1e-3, number))

    block = ""

    for number, line in enumerate(dhf_path.read_text().splitlines(), 1):
        if line.startswith("+"):
            block = line[1:].strip()
        elif line.startswith("-"):
            block = ""
        elif block in ("MODEL/RANGE_BIAS", "SOLUTION/DATA_HANDLING"):
            fields = line.split()

            if fields and fields[0] in STATIONS and fields[1] in ("51", "--"):
                dhf_rows.append(
                    (
                        fields[0],
                        sinex_bound(fields[4], False),
                        sinex_bound(fields[5], True),
                        fields[6],
                        fields[7:],
                        block,
                        number,
                    )
                )

    return com_rows, dhf_rows


def corrections_at(station, utc, wavelength_nm, com_rows, dhf_rows):
    """Select published corrections in their stated intervals, rejecting other DHF instructions."""
    candidates = [
        row
        for row in com_rows
        if row[0] == station and row[1] <= utc < row[2] and abs(row[3] - wavelength_nm) <= 50
    ]

    if len(candidates) != 1:
        raise ValueError(f"Expected one CoM row for {station}, found {len(candidates)}")

    active = [row for row in dhf_rows if row[0] == station and row[1] <= float(utc) < row[2]]
    biases = []

    for row in active:
        if row[3] != "R" or row[5] != "MODEL/RANGE_BIAS" or "mm" not in row[4]:
            raise ValueError(f"Unsupported active DHF instruction at line {row[6]}: {row}")

        biases.append((float(row[4][0]) * 1e-3, row[6]))

    if len(biases) != 1:
        raise ValueError(f"Expected one published range bias for {station}")

    # The CoM length is added; the DHF bias is subtracted from one-way observed range.
    return candidates[0][4], biases[0][0], candidates[0][5], biases[0][1]


class FixedOrbit:
    """Interpolate one UTC SLRF2020 orbit and supply its BCRS state to Tudat."""

    def __init__(self, data, earth, sun, converter, points):
        """Keep local TT seconds for interpolation and reject unsupported product labels."""
        if data.time_system != "UTC" or data.coordinate_system != "SLR20":
            raise ValueError("This analysis requires the selected UTC SLRF2020 provider product")

        if not any("SLRF2020" in comment for comment in data.comments):
            raise ValueError("The SP3 must explicitly identify SLRF2020")

        if len(data.satellite_ids) != 1:
            raise ValueError("Select a single-satellite LAGEOS-1 SP3")

        records = data.records[data.satellite_ids[0]]

        if any(
            not r.position_available or not r.velocity_available or r.orbit_predicted or r.maneuver
            for r in records
        ):
            raise ValueError(
                "Missing, predicted or maneuver SP3 records require a separate analysis"
            )

        self.earth, self.sun, self.converter = earth, sun, converter
        self.rotation = earth.rotation_model
        self.origin = converter.convert_time_object(time.utc_scale, time.tt_scale, records[0].epoch)
        self.epochs = [
            converter.convert_time_object(time.utc_scale, time.tt_scale, r.epoch) for r in records
        ]
        self.seconds = np.array([float(t - self.origin) for t in self.epochs])
        self.states = np.array([np.r_[r.position, r.velocity].reshape(6) for r in records])

        if not np.allclose(np.diff(self.seconds), data.epoch_interval, rtol=0, atol=1e-9):
            raise ValueError("SP3 gaps or a changing cadence are unsupported")

        # Both orders use the same interior domain; the polynomial is never extrapolated.
        self.margin = max(points, 10) // 2
        self.lower, self.upper = self.seconds[self.margin], self.seconds[-self.margin - 1]
        self.interpolators = {
            n: self.make_interpolator(self.seconds, self.states, n) for n in {points, 8, 10}
        }
        self.points = points
        self.maximum_callback_half_ulp_seconds = 0.0

    @staticmethod
    def make_interpolator(epochs, states, points):
        """Use Tudat's existing Lagrange interpolator on small epoch differences."""
        settings = interpolators.lagrange_interpolation(
            points,
            boundary_interpolation=interpolators.throw_exception_at_boundary,
            lagrange_boundary_handling=interpolators.lagrange_no_boundary_interpolation,
        )

        return interpolators.create_one_dimensional_vector_interpolator(
            {float(t): state for t, state in zip(epochs, states)}, settings
        )

    def terrestrial_state(self, tt, points=None):
        """Return a state only inside the common, fully supported interpolation interval."""
        seconds = float(tt - self.origin)

        if not self.lower <= seconds <= self.upper:
            raise ValueError("Epoch outside the SP3 interpolation interior")

        return np.asarray(self.interpolators[points or self.points].interpolate(seconds)).reshape(6)

    def relative_position(self, tdb):
        """Apply IERS 11.19 spatial scaling and the topocentric simultaneity term."""
        earth_state = np.asarray(self.earth.state_in_base_frame_from_ephemeris(tdb)).reshape(6)
        sun_state = np.asarray(self.sun.state_in_base_frame_from_ephemeris(tdb)).reshape(6)
        tt_geocentre = self.converter.convert_time_object(
            time.tdb_scale, time.tt_scale, time.Time(tdb)
        )
        tt_site = tt_geocentre

        # Iterate the position-dependent offset before evaluating and rotating the fixed orbit.
        for _ in range(4):
            rotation_time = self.converter.convert_time_object(
                time.tt_scale, time.tdb_scale, tt_site
            )
            rotation = self.rotation.body_fixed_to_inertial_rotation(float(rotation_time))
            position = rotation @ self.terrestrial_state(tt_site)[:3]
            tt_site = tt_geocentre - float(np.dot(earth_state[3:], position) / C**2)

        rotation_time = self.converter.convert_time_object(time.tt_scale, time.tdb_scale, tt_site)
        position = (
            self.rotation.body_fixed_to_inertial_rotation(float(rotation_time))
            @ self.terrestrial_state(tt_site)[:3]
        )
        potential = self.sun.gravity_field_model.gravitational_parameter / np.linalg.norm(
            sun_state[:3] - earth_state[:3]
        )

        # Match the native station model's TT-compatible to TDB-compatible spatial convention.
        position -= (potential / C**2 + SPATIAL_TT_TO_TDB_SCALE) * position + (
            0.5 * earth_state[3:] * np.dot(earth_state[3:], position) / C**2
        )

        return position

    def state(self, tdb):
        """Supply a fixed-orbit BCRS state at the custom ephemeris's double TDB boundary."""
        self.maximum_callback_half_ulp_seconds = max(
            self.maximum_callback_half_ulp_seconds, abs(np.spacing(tdb)) / 2
        )
        earth_state = np.asarray(self.earth.state_in_base_frame_from_ephemeris(tdb)).reshape(6)
        position = self.relative_position(tdb)

        # The observable uses positions; this symmetric derivative supplies a consistent velocity.
        step = 0.5
        velocity = (self.relative_position(tdb + step) - self.relative_position(tdb - step)) / (
            2 * step
        )

        return earth_state + np.r_[position, velocity]

    def diagnostics(self):
        """Report withheld source-node errors without using the observed residuals."""
        errors = {8: [], 10: []}

        # Omit one central sample in each local window; do not compare a fitted node to itself.
        for index in range(16, len(self.seconds) - 16, 41):
            indices = [j for j in range(index - 8, index + 9) if j != index]
            origin = self.seconds[index]

            for points in errors:
                model = self.make_interpolator(
                    self.seconds[indices] - origin, self.states[indices], points
                )
                estimate = np.asarray(model.interpolate(0.0)).reshape(6)
                errors[points].append(float(np.linalg.norm(estimate[:3] - self.states[index, :3])))

        return {
            str(n): {
                "count": len(values),
                "rms_m": float(np.sqrt(np.mean(np.square(values)))),
                "max_m": max(values),
            }
            for n, values in errors.items()
        }


def station_simultaneity(earth, nominal_position):
    """Return the small rigid-station retiming correction absent from the native spatial model."""
    rotation = earth.rotation_model

    def position(tdb):
        """Retard the station orientation by v_E dot r divided by c squared."""
        matrix = rotation.body_fixed_to_inertial_rotation(tdb)
        derivative = rotation.time_derivative_body_fixed_to_inertial_rotation(tdb)
        velocity = np.asarray(earth.state_in_base_frame_from_ephemeris(tdb)).reshape(6)[3:]
        offset = float(np.dot(velocity, matrix @ nominal_position) / C**2)

        return -offset * (matrix.T @ derivative @ nominal_position)

    return position


def load_stations(station_file, eccentricity_file, converter, require_velocity):
    """Select the two unique monument solutions without projecting unrelated SINEX pieces."""
    history = station_history.read_sinex_station_history(str(station_file))
    keys = station_history.get_station_coordinate_solution_keys(history)
    eccentricities = ilrs.read_sinex_station_eccentricities(str(eccentricity_file))
    settings, ledger = [], []
    week_start = float(calendar_epoch(2024, 3, 31))
    week_end = float(calendar_epoch(2024, 4, 7))

    for code, domes in STATIONS.items():
        selected = [key for key in keys if key.site.site_code == code]

        if len(selected) != 1:
            raise ValueError(f"Expected one selected monument solution for {code}")

        solution = station_history.get_station_coordinate_solution(
            history, selected[0], require_velocity
        )

        if solution.site.domes_id != domes:
            raise ValueError("Selected station code and DOMES identity disagree")

        # The archived combined solution explicitly supplies zero XYZ telescope eccentricities for this week.
        offsets = [
            row
            for row in eccentricities.get(domes, [])
            if row.start_epoch <= week_start and (row.has_open_end or row.end_epoch >= week_end)
        ]

        if len(offsets) != 1 or np.any(np.asarray(offsets[0].eccentricity) != 0):
            raise ValueError(
                "This analysis requires the documented zero XYZ eccentricity over the whole week"
            )

        civil = solution.position_reference_epoch.civil_utc
        tai = environment.sinex_civil_utc_to_tai(
            civil.year, civil.day_of_year, civil.seconds_of_day
        )
        reference = converter.convert_time_object(time.tai_scale, time.tdb_scale, tai)
        velocity = solution.velocity_metres_per_second
        motions = []

        if all(value is not None for value in velocity):
            if any(
                epoch.original != solution.position_reference_epoch.original
                for epoch in solution.velocity_reference_epochs
            ):
                raise ValueError(
                    "Distinct station position/velocity origins require another declared model"
                )

            motions.append(
                ground_station.linear_station_motion(np.array(velocity), float(reference))
            )
        elif any(value is not None for value in velocity):
            raise ValueError("A partial station velocity cannot define this experiment")

        settings.append(
            ground_station.basic_station(
                domes, np.array(solution.position_metres), station_motion_settings=motions
            )
        )
        ledger.append(
            {
                "station": code,
                "domes": domes,
                "point": selected[0].site.point_code,
                "solution": selected[0].solution_id,
                "reference_utc": solution.position_reference_epoch.original,
                "position_m": list(solution.position_metres),
                "velocity_m_per_s": list(velocity),
                "eccentricity_m": [0.0, 0.0, 0.0],
                "velocity_interpretation": (
                    "SI linear velocity"
                    if motions
                    else "weekly position snapshot, no velocity supplied"
                ),
            }
        )

    return settings, ledger


def set_weather(station, session):
    """Attach measured CRD weather in UTC, rejecting duplicate conflicting samples."""
    config = session.configuration
    start = calendar_epoch(config.start_year, config.start_month, config.start_day)
    end = calendar_epoch(config.end_year, config.end_month, config.end_day)
    midpoint = (
        float(end - start)
        + config.end_hour * 3600
        + config.end_minute * 60
        + config.end_second
        + config.start_hour * 3600
        + config.start_minute * 60
        + config.start_second
    ) / 2
    values = {}

    # Weather records may be unsorted. The nearest day to this short session determines midnight rollover.
    for record in session.data.meteorological_data:
        if record.origin_of_values != 0:
            raise ValueError("This analysis requires measured weather records")

        local = record.second_of_day + 86400 * math.floor(
            (midpoint - record.second_of_day) / 86400 + 0.5
        )
        utc = float(start + local)
        value = np.array([record.temperature, record.pressure * 100, record.humidity * 0.01])

        if (
            not np.all(np.isfinite(value))
            or value[0] <= 0
            or value[1] <= 0
            or not 0 <= value[2] <= 1
        ):
            raise ValueError("Missing or invalid measured CRD weather")

        if utc in values and not np.array_equal(values[utc], value):
            raise ValueError("Conflicting weather samples at the same UTC epoch")

        values[utc] = value

    if len(values) < 2:
        raise ValueError("At least two weather samples are required")

    model = interpolators.create_one_dimensional_vector_interpolator(
        dict(sorted(values.items())),
        interpolators.linear_interpolation(
            boundary_interpolation=interpolators.throw_exception_at_boundary
        ),
    )
    station.set_station_meteo_data(
        environment.ContinuousInterpolatedMeteoData(
            model,
            {
                environment.temperature_meteo_data: 0,
                environment.pressure_meteo_data: 1,
                environment.relative_humidity_meteo_data: 2,
            },
        )
    )

    return min(values), max(values)


def statistics(rows):
    """Keep the signed mean, uncentred RMS and centred scatter distinct."""
    values = np.array([row["observed_minus_computed_m"] for row in rows])

    return (
        {
            "count": len(values),
            "mean_m": float(np.mean(values)),
            "rms_m": float(np.sqrt(np.mean(values**2))),
            "std_m": float(np.std(values)),
            "median_m": float(np.median(values)),
        }
        if len(values)
        else {"count": 0}
    )


def reduce_session(session, session_index, bodies, orbit, correction_rows):
    """Evaluate every supported normal point, preserving duplicates and recording exclusions."""
    config = session.configuration
    station_code = str(config.cdp_pad_id)
    station = bodies.get("Earth").get_ground_station(STATIONS[station_code])
    converter = orbit.converter
    rows, excluded = [], []

    # These selected passes contain raw two-way ranges and already calibrated station delays.
    required = {
        "range_type": 2,
        "tropospheric_correction_applied": 0,
        "center_of_mass_correction_applied": 0,
        "receive_amplitude_correction_applied": 0,
        "station_system_delay_applied": 1,
        "spacecraft_system_delay_applied": 0,
        "data_quality_alert": 0,
    }

    if config.ilrs_satellite_id != "7603901" or config.station_epoch_time_scale not in (3, 4, 7):
        raise ValueError("Unexpected target or CRD station time reference")

    for name, expected in required.items():
        if getattr(config, name) != expected:
            raise ValueError(f"Unsupported CRD {name}: {getattr(config, name)}")

    weather_start, weather_end = set_weather(station, session)
    extracted = slr.extract_normal_point_measurements(session)
    nominal = np.asarray(station.station_state.cartesian_position_at_reference_epoch).reshape(3)
    links = observable_models_setup.links
    model_settings = observable_models_setup.model_settings
    corrections = observable_models_setup.light_time_corrections
    simulation = observations_setup.observations_simulation_settings
    link = links.LinkDefinition(
        {
            links.transmitter: links.body_reference_point_link_end_id(
                "Earth", STATIONS[station_code]
            ),
            links.reflector1: links.body_origin_link_end_id("LAGEOS1"),
            links.receiver: links.body_reference_point_link_end_id("Earth", STATIONS[station_code]),
        }
    )

    for index, ((utc, measured_range), record) in enumerate(
        zip(extracted, session.data.normal_point_data)
    ):
        if record.epoch_event != 2:
            raise ValueError("This selected analysis requires transmit-epoch normal points")

        wavelength_nm = config.system_configurations[
            record.system_configuration_id
        ].transmit_wavelength_nm
        com, bias, com_line, dhf_line = corrections_at(
            station_code, utc, wavelength_nm, *correction_rows
        )
        transmit = converter.convert_time_object(time.utc_scale, time.tdb_scale, utc, nominal)
        receive = converter.convert_time_object(
            time.utc_scale, time.tdb_scale, utc + record.two_way_time_of_flight, nominal
        )
        reason = None

        # Selection uses coverage and elevation only, never the size or sign of a residual.
        if (
            not weather_start <= float(utc)
            or float(utc + record.two_way_time_of_flight) > weather_end
        ):
            reason = "outside measured weather interval"

        try:
            elevations = []

            for endpoint in (transmit, receive):
                orbit.terrestrial_state(
                    converter.convert_time_object(time.tdb_scale, time.tt_scale, endpoint) - 2.0
                )
                orbit.terrestrial_state(
                    converter.convert_time_object(time.tdb_scale, time.tt_scale, endpoint) + 2.0
                )
                tdb = float(endpoint)
                local = np.asarray(
                    station.station_state.get_cartesian_position(endpoint, "SSB")
                ).reshape(3)
                station_bcrs = (
                    np.asarray(bodies.get("Earth").state_in_base_frame_from_ephemeris(tdb)).reshape(
                        6
                    )[:3]
                    + bodies.get("Earth").rotation_model.body_fixed_to_inertial_rotation(tdb)
                    @ local
                )
                bounce = float(transmit + record.two_way_time_of_flight / 2)
                direction = orbit.state(bounce)[:3] - station_bcrs
                elevations.append(
                    station.pointing_angles_calculator.calculate_elevation_angle(direction, tdb)
                )

            if min(elevations) < np.deg2rad(20):
                reason = "below provider 20 degree elevation cutoff"
        except ValueError as error:
            if str(error) != "Epoch outside the SP3 interpolation interior":
                raise

            reason = str(error)

        if reason:
            excluded.append(
                {
                    "session": session_index,
                    "normal_point": index,
                    "station": station_code,
                    "reason": reason,
                }
            )
            continue

        # The third positional input is the native fractional light-time tolerance, despite its Python argument name.
        convergence = corrections.light_time_convergence_settings(
            True, 50, 1e-12, corrections.throw_exception
        )
        settings = model_settings.two_way_range(
            link,
            [
                corrections.first_order_relativistic_light_time_correction(["Sun", "Earth"]),
                corrections.mendes_pavlis_tropospheric_light_time_correction(wavelength_nm * 1e-9),
            ],
            light_time_convergence_settings=convergence,
            time_scale_for_observable=time.tai_scale,
        )
        simulators = simulation.create_observation_simulators([settings], bodies)
        schedule = simulation.tabulated_simulation_settings(
            model_settings.n_way_range_type, link, [transmit], links.transmitter
        )
        result = observations.simulate_observations([schedule], simulators, bodies)
        computed = np.asarray(result.concatenated_observations).reshape(-1)

        if len(computed) != 1 or not np.isfinite(computed[0]):
            raise RuntimeError(
                "A selected normal point did not produce exactly one finite observable"
            )

        hours, seconds = utc.__getstate__()
        tt = converter.convert_time_object(time.utc_scale, time.tt_scale, utc)
        difference = orbit.terrestrial_state(tt, 8)[:3] - orbit.terrestrial_state(tt, 10)[:3]
        rows.append(
            {
                "station": station_code,
                "domes": STATIONS[station_code],
                "session": session_index,
                "normal_point": index,
                "source": config.source,
                "utc_hours_since_j2000": hours,
                "utc_seconds_into_hour": seconds,
                "epoch_event": record.epoch_event,
                "system_configuration": record.system_configuration_id,
                "wavelength_nm": wavelength_nm,
                "tof_seconds": record.two_way_time_of_flight,
                "observed_reflector_range_m": measured_range,
                "com_added_m": com,
                "dhf_subtracted_m": bias,
                "com_source_line": com_line,
                "dhf_source_line": dhf_line,
                "calibration_reapplied_m": 0.0,
                "computed_full_two_way_range_m": float(computed[0]),
                "observed_minus_computed_m": measured_range + com - bias - 0.5 * computed[0],
                "minimum_elevation_deg": float(np.rad2deg(min(elevations))),
                "np_bin_rms_ps": record.bin_rms,
                "sp3_order8_minus10_at_transmit_norm_m": float(np.linalg.norm(difference)),
            }
        )

    return rows, excluded


def main():
    """Assemble the fixed experiment and write raw residuals plus a reproducible correction ledger."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--orbit", choices=("combined", "asi"), default="combined")
    parser.add_argument(
        "--station-solution", choices=("slrf2020", "asi-weekly"), default="slrf2020"
    )
    parser.add_argument("--got410c-blq", type=Path)
    parser.add_argument("--points", choices=(8, 10), type=int, default=10)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    input_directory = args.output / "uncompressed-inputs"
    input_directory.mkdir()
    root = args.data_root
    orbit_name = (
        "ilrsa.orb.lageos1.240406.v81.sp3.gz"
        if args.orbit == "combined"
        else "asi.orb.lageos1.240406.v80.sp3.gz"
    )
    station_path = (
        root / "metadata/SLRF2020_POS+VEL_2023.10.02.snx"
        if args.station_solution == "slrf2020"
        else root / "2024-04/asi.pos+eop.240406.v80.snx.gz"
    )
    eop_path = Path(resource_paths.get_earth_orientation_path()) / "eopc04_20_IAU2000.62-now.txt"
    crd_files = sorted((root / "2024-04").glob("lageos1_20240[34]*.np2"))

    if len(crd_files) != 7:
        raise ValueError("The pinned week requires all seven daily CRDv2 files")

    inputs = [root / "2024-04" / orbit_name, station_path, eop_path, *crd_files]
    inputs.extend(
        root / path
        for path in (
            "2024-04/ilrsa.pos+eop.240406.v80.snx.gz",
            "metadata/com6.220915.tar.gz",
            "metadata/ILRS_Data_Handling_File_2024.02.13.snx",
        )
    )

    if args.got410c_blq:
        inputs.append(args.got410c_blq)

    # Pin every provider input against the download manifest before reading models.
    manifest = json.loads(Path(__file__).with_name("dataset.json").read_text())
    known = {Path(item["path"]).name: item["sha256"] for item in manifest}
    provenance = [file_record(path) for path in inputs]

    for record in provenance:
        expected = known.get(Path(record["path"]).name)

        if expected and record["sha256"] != expected:
            raise ValueError(f"Downloaded input changed: {record['path']}")

    # Load only the shared Sun/Earth/Moon ephemeris, constants, body shapes and leap-second kernel.
    kernel_directory = Path(resource_paths.get_spice_kernel_path())
    spice_inputs = [
        kernel_directory / name
        for name in (
            "pck00010.tpc",
            "inpop19a_TDB_m100_p100_spice.tpc",
            "inpop19a_TDB_m100_p100_spice.bsp",
            "naif0012.tls",
        )
    ]
    spice.clear_kernels()

    for path in spice_inputs:
        spice.load_kernel(str(path))

    provenance.extend(file_record(path) for path in spice_inputs)
    settings = environment_setup.get_default_body_settings(["Earth", "Sun", "Moon"], "SSB", "J2000")
    settings.get("Earth").rotation_model_settings = rotation_model.gcrs_to_itrs(
        base_frame="J2000", eop_data=earth_orientation.read_earth_orientation_data(str(eop_path))
    )
    deformations = [shape_deformation.iers_2010_solid_body_tidal(), shape_deformation.pole_tidal()]

    if args.got410c_blq:
        deformations.append(shape_deformation.ocean_tidal([str(args.got410c_blq)]))

    settings.get("Earth").shape_deformation_settings = deformations
    bodies = environment_setup.create_system_of_bodies(settings)
    earth = bodies.get("Earth")
    converter = earth.rotation_model.angles_calculator.time_scale_converter
    orbit = FixedOrbit(
        environment_data.read_sp3_file(str(unpack_input(inputs[0], input_directory))),
        earth,
        bodies.get("Sun"),
        converter,
        args.points,
    )
    bodies.create_empty_body("LAGEOS1")
    bodies.get("LAGEOS1").ephemeris = environment_setup.create_body_ephemeris(
        ephemeris.custom_ephemeris(orbit.state, "SSB", "J2000"), "LAGEOS1"
    )

    # Refresh the existing frame callbacks after assigning the newly added body's ephemeris.
    bodies.add_body(bodies.get("LAGEOS1"), "LAGEOS1")
    station_file = unpack_input(station_path, input_directory)

    # Both selected stations have zero XYZ eccentricities in the pinned weekly combined solution.
    eccentricity_file = unpack_input(
        root / "2024-04/ilrsa.pos+eop.240406.v80.snx.gz", input_directory
    )
    station_settings, station_ledger = load_stations(
        station_file, eccentricity_file, converter, args.station_solution == "slrf2020"
    )

    for station in station_settings:
        nominal = np.asarray(station.station_position).reshape(3)
        station.add_station_motion_settings(ground_station.body_deformation_station_motion())
        station.add_station_motion_settings(
            ground_station.custom_station_motion(station_simultaneity(earth, nominal))
        )
        station.add_station_motion_settings(ground_station.relativistic_station_motion(["Sun"]))
        environment_setup.add_ground_station(earth, station, bodies)

    com_rows, dhf_rows = read_corrections(root)
    rows, excluded, session_metadata = [], [], []
    sessions = slr.read_slr_data([str(path) for path in crd_files])

    for index, session in enumerate(sessions):
        if str(session.configuration.cdp_pad_id) not in STATIONS:
            continue

        config = session.configuration
        session_metadata.append(
            {
                "session": index,
                "source": config.source,
                "station": config.cdp_pad_id,
                "system_number": config.system_number,
                "occupancy": config.occupancy,
                "station_epoch_time_scale": config.station_epoch_time_scale,
                "data_release": config.data_release,
                "flags": {
                    name: getattr(config, name)
                    for name in (
                        "tropospheric_correction_applied",
                        "center_of_mass_correction_applied",
                        "receive_amplitude_correction_applied",
                        "station_system_delay_applied",
                        "spacecraft_system_delay_applied",
                        "range_type",
                        "data_quality_alert",
                    )
                },
                "calibration_records": len(session.data.calibration_data),
                "range_supplement_records": len(session.data.range_supplement_data),
            }
        )
        values, omissions = reduce_session(session, index, bodies, orbit, (com_rows, dhf_rows))
        rows.extend(values)
        excluded.extend(omissions)

    if not rows:
        raise RuntimeError(
            "No selected normal points survived the declared physical coverage checks"
        )

    with (args.output / "residuals.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "orbit": args.orbit,
        "station_solution": args.station_solution,
        "inputs": provenance,
        "software": {
            "interpreter": sys.executable,
            "tudatpy": tudatpy.__file__,
            "kernel": file_record(kernel.__file__),
            "script": file_record(__file__),
        },
        "station_coordinates": station_ledger,
        "session_metadata": session_metadata,
        "all": statistics(rows),
        "stations": {
            code: statistics([r for r in rows if r["station"] == code]) for code in STATIONS
        },
        "sessions": {
            str(i): statistics([r for r in rows if r["session"] == i])
            for i in sorted({r["session"] for r in rows})
        },
        "excluded": excluded,
        "sp3_holdout": orbit.diagnostics(),
        "maximum_custom_callback_half_ulp_seconds": orbit.maximum_callback_half_ulp_seconds,
        "relative_callback_position_quantization_bound_m": orbit.maximum_callback_half_ulp_seconds
        * float(
            np.max(
                np.linalg.norm(orbit.states[:, 3:], axis=1)
                + 7.3e-5 * np.linalg.norm(orbit.states[:, :3], axis=1)
            )
        ),
        "sp3_position_quantization_m": 0.001,
        "ocean_loading": (
            "user-declared GOT4.10c BLQ"
            if args.got410c_blq
            else "OMITTED: exact GOT4.10c unavailable"
        ),
        "models": {
            "rotation": "IERS 2010, IAU 2006, selected C04 file",
            "observable_time_scale": "TAI",
            "optical": "Mendes-Pavlis + FCULa; CRD measured weather; C0 wavelength; CO2 375 ppm",
            "station_tides": "IERS2010 solid and pole; native shared mean-pole convention",
            "range_convention": "one-way O-C = c*TOF/2 + CoM - DHF_bias - computed_two_way/2",
            "station_retiming": "first-order rigid nominal orientation; deformation/tectonic retiming omitted",
            "spatial": "IERS11.19 TT/TDB scale and solar potential; native station model shares this choice",
            "relativity": "Sun and Earth light time; no fitted biases or estimation",
        },
    }
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(json.dumps(summary["all"], indent=2))


if __name__ == "__main__":
    main()
