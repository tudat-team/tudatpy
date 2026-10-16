"""Audit three absolute launch epochs independently of native TDB conversion."""

from pathlib import Path
from decimal import Decimal
from datetime import date
import hashlib
import json
import math
import subprocess
import sys

import erfa
import numpy as np
import tudatpy.kernel
from tudatpy.astro import time_representation as time

WORK = Path("/home/dominic/Tudat/laser-ranging-work")
SOURCE = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-signature")
sys.path.insert(0, str(SOURCE / "analysis/laser_ranging"))
import llr_residuals as model
import llr_long_arc as arc


def sha(path):
    """Identify the local source or evidence without altering its bytes."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def reference_calendar(label):
    """Parse the saved ISO tag independently, retaining its printed decimal fraction."""
    day, clock = label.rstrip("Z").split("T")
    year, month, day = map(int, day.split("-"))
    hour, minute, seconds = clock.split(":")
    second = Decimal(seconds)
    integer = (date(year, month, day) - date(2000, 1, 1)).days * 86400 - 43200
    integer += int(hour) * 3600 + int(minute) * 60 + int(second)
    epoch = time.Time(integer) + float(second - int(second))
    utc1, utc2 = erfa.dtf2d("UTC", year, month, day, int(hour), int(minute), float(second))
    tai1, tai2 = erfa.utctai(utc1, utc2)
    tt1, tt2 = erfa.taitt(tai1, tai2)
    tt_minus_utc = float(((tt1 - utc1) + (tt2 - utc2)) * 86400)
    return epoch, (tt1, tt2), tt_minus_utc


manifest = json.loads(arc.MANIFEST.read_text())
ephemeris = json.loads((Path(model.__file__).with_name("llr_dataset.json")).read_text())
model.load_ephemeris_kernels(
    WORK / "data/llr/kernels", ephemeris, ephemeris["analysis_frame_alias"]
)
groups = arc.read_groups(WORK / "research/longarc-20170304", manifest)
apollo = [
    (group, point) for group in groups if group[3]["station"] == "APOLLO" for point in group[0]
]
grasse = next(
    (group, point)
    for group in groups
    for point in group[0]
    if point["utc_label"].startswith("2017-03-09T18:34:15.297622")
)
selected = [apollo[0], grasse, apollo[-1]]
rows = []
for group, point in selected:
    points, weather, profile, ledger = group
    bodies, nominal, _ = model.create_environment(
        points, station_profile=profile, meteorological_history=weather
    )
    converter = bodies.get("Earth").rotation_model.angles_calculator.time_scale_converter
    utc = point["utc"]
    launch = converter.convert_time_object(time.utc_scale, time.tdb_scale, utc, nominal)
    tt = converter.convert_time_object(time.utc_scale, time.tt_scale, utc)
    geocentric = converter.convert_time_object(time.utc_scale, time.tdb_scale, utc)
    ut1 = converter.convert_time_object(time.utc_scale, time.ut1_scale, utc)
    independent_utc, (tt1, tt2), tt_offset = reference_calendar(point["utc_label"])
    # UT1 is deliberately the selected native EOP result; only the TDB/calendar route is independent.
    hours, seconds = ut1.__getstate__()
    ut_fraction = (((hours + 12) % 24) * 3600 + seconds) / 86400
    longitude = math.atan2(nominal[1], nominal[0])
    u, v = math.hypot(nominal[0], nominal[1]) / 1000, nominal[2] / 1000
    ref_geo = float(erfa.dtdb(tt1, tt2, 0.0, 0.0, 0.0, 0.0))
    ref_topo = float(erfa.dtdb(tt1, tt2, ut_fraction, longitude, u, v))
    reference = independent_utc + tt_offset + ref_topo
    ref_geo_epoch = independent_utc + tt_offset + ref_geo
    saved = {}
    original_settings = model.sim.tabulated_simulation_settings
    original_simulate = model.observations.simulate_observations

    def capture_settings(*args, **kwargs):
        """Inspect the existing native Time-valued tabulated epoch, without replacing it."""
        settings = original_settings(*args, **kwargs)
        saved["simulation_epoch"] = settings.simulation_times[0]
        return settings

    def capture_result(*args, **kwargs):
        """Inspect the native result's retained Time-valued reference epoch."""
        result = original_simulate(*args, **kwargs)
        saved["result_epoch"] = result.concatenated_times_objects[0]
        return result

    model.sim.tabulated_simulation_settings = capture_settings
    model.observations.simulate_observations = capture_result
    try:
        computed, diagnostics = model.computed_range(
            point, launch, bodies, model.complete_light_time_corrections(point)
        )
    finally:
        model.sim.tabulated_simulation_settings = original_settings
        model.observations.simulate_observations = original_simulate
    earth = np.asarray(bodies.get("Earth").state_in_base_frame_from_ephemeris(launch)).reshape(6)
    moon = np.asarray(bodies.get("Moon").state_in_base_frame_from_ephemeris(launch)).reshape(6)
    relative = moon - earth
    range_rate = float(relative[:3] @ relative[3:] / np.linalg.norm(relative[:3]))
    row = dict(
        station=point["station"],
        utc=point["utc_label"],
        source_file=point["source_file"],
        source_line=point["source_line"],
        station_profile=profile,
        nominal_xyz_m=nominal.tolist(),
        utc_time_split=utc.__getstate__(),
        launch_tdb_time_split=launch.__getstate__(),
        selected_ut1_time_split=ut1.__getstate__(),
        utc_parser_minus_independent_label_s=float(utc - independent_utc),
        native_tt_minus_utc_s=float(tt - utc),
        erfa_tt_minus_utc_s=tt_offset,
        native_topocentric_tdb_minus_tt_s=float(launch - tt),
        erfa_topocentric_tdb_minus_tt_s=ref_topo,
        native_geocentric_tdb_minus_tt_s=float(geocentric - tt),
        erfa_geocentric_tdb_minus_tt_s=ref_geo,
        absolute_native_minus_independent_tdb_s=float(launch - reference),
        absolute_geocentric_minus_independent_tdb_s=float(geocentric - ref_geo_epoch),
        simulation_epoch_minus_launch_s=float(saved["simulation_epoch"] - launch),
        retained_result_epoch_minus_launch_s=float(saved["result_epoch"] - launch),
        dependent_transmit_double_equals_launch_double=diagnostics["link_end_epochs_tdb_s"][0]
        == float(launch),
        dependent_transmit_double_roundoff_s=float(
            time.Time(diagnostics["link_end_epochs_tdb_s"][0]) - launch
        ),
        native_computed_atomic_range_m=computed,
        geocentric_range_rate_m_per_s=range_rate,
        hypothetical_range_change_from_one_extra_geocentric_tdb_minus_tt_m=range_rate * ref_geo,
    )
    assert (
        row["simulation_epoch_minus_launch_s"] == row["retained_result_epoch_minus_launch_s"] == 0
    )
    assert row["dependent_transmit_double_equals_launch_double"]
    assert abs(row["absolute_native_minus_independent_tdb_s"]) < 1e-9
    rows.append(row)

header = WORK / "research/llr-physical-terms/header.440"
triples = [
    list(map(int, line.split()))
    for line in header.read_text().split("GROUP   1050")[1].split("GROUP")[0].splitlines()
    if line.strip()
]
assert len(triples) == 3 and all(len(line) == 15 for line in triples)
assert triples[1][-1] == triples[2][-1] == 0
output = dict(
    scope="Read-only absolute epoch audit; three original observations, no changed clock/tolerance or full arc rerun.",
    source_revision=subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=SOURCE, text=True
    ).strip(),
    kernel_sha256=sha(tudatpy.kernel.__file__),
    erfa_version=erfa.__version__,
    model_source_sha256=sha(model.__file__),
    driver_source_sha256=sha(arc.__file__),
    audit_source_sha256=sha(__file__),
    independent_reference="ERFA2.0.1.5 SOFA-derived dtf2d/utctai/taitt/dtdb; independent Python route and calendar parser. Native selected UT1 supplied only for topocentric term.",
    primary_reference="https://www.iausofa.org/s/sofa_ts_c.pdf sections4.3.4/4.4; localSOFA dtdb.c notes1,6,7.",
    sofa_source="/home/dominic/Tudat/sofa-cmake/src/dtdb.c",
    sofa_source_sha256=sha("/home/dominic/Tudat/sofa-cmake/src/dtdb.c"),
    de440_header_sha256=sha(header),
    de440_group1050=triples,
    de440_tt_minus_tdb_available=False,
    de440_time_ephemeris_note="Column15 TT−TDB has zero coefficients/subintervals in original cached DE440 ASCII; no independent DE440 time-ephemeris comparison is possible with this file. Fairhead/SOFA is not claimed exactDE440.",
    implementation_trace=[
        "llr_long_arc.py converts UTC→TDB once with convert_time_object and selected station/converter.",
        "computed_range passes this Time into tabulated_simulation_settings with reference transmitter.",
        "simulateObservationWithCheck passes observationTime directly into computeObservationsWithLinkEndData.",
        "NWayRange forwards input Time directly to full-link light-time calculator; TAI selection adjusts only final RX−TX interval.",
    ],
    numerical_limit="Independent ERFA uses splitJD, native splitTime comparisons avoid subtracting large epoch floats. Diagnostic link-end epochs are double, explicitly measured separately; native solver input and retained collection epoch preserve Time exactly. SOFA dtdb note7 specifies better than±3ns over1950–2050 relative to DE405 numerical time ephemerides, not a verified bound versusDE440. At100m/s even10ns is1µm; independent analytic agreement rules out missing~1.6ms but is not an independent DE440 time integration.",
    rows=rows,
    max_abs_absolute_epoch_difference_s=max(
        abs(row["absolute_native_minus_independent_tdb_s"]) for row in rows
    ),
    conclusion="No missed/doubled absolute TT↔TDB conversion found. Positive~1.6ms periodic correction is present exactlyonce before the native solver; atomic output clock does not remap the launch.",
)
out = WORK / "logs/llr-absolute-launch-clock-audit.json"
out.write_text(json.dumps(output, indent=2) + "\n")
print(json.dumps(output, indent=2))
