"""Screen original annual lunar NPs with directDE440 geometry, no atmospheric/range solves."""

import csv
from datetime import date
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, "/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-llr/analysis/laser_ranging")
import llr_residuals as model

ROOT = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-llr/analysis/laser_ranging")
DATA = Path("/home/dominic/Tudat/laser-ranging-work/research/llr-year-2017")
OUT = Path("/home/dominic/Tudat/laser-ranging-work/results/llr-year-2017-geometry-preflight")
manifest = json.loads(
    (
        Path(
            "/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-data/analysis/laser_ranging/llr_year_2017_dataset.json"
        )
    ).read_text()
)
ephemeris = json.loads((ROOT / "llr_dataset.json").read_text())
kernels, _ = model.load_ephemeris_kernels(
    Path("/home/dominic/Tudat/laser-ranging-work/data/llr/kernels"),
    ephemeris,
    ephemeris["analysis_frame_alias"],
)
settings = model.setup.BodyListSettings("SSB", "ICRF")
settings.add_empty_settings("Earth")
earth = settings.get("Earth")
earth.ephemeris_settings = model.setup.ephemeris.direct_spice("SSB", "ICRF")
earth.rotation_model_settings = model.setup.rotation_model.gcrs_to_itrs(base_frame="ICRF")
earth.shape_settings = model.setup.shape.oblate_spherical(6378137.0, 1 / 298.257223563)
bodies = model.setup.create_system_of_bodies(settings)
rotation = bodies.get("Earth").rotation_model
converter = rotation.angles_calculator.time_scale_converter
profiles = {"APOLLO": "inpop17a-gir", "GRASSE": "inpop17a-grasse-gir"}
positions = {
    name: np.asarray(model.STATION_PROFILES[profile]["position_m"])
    for name, profile in profiles.items()
}
normals = {}
for name, position in positions.items():
    geodetic = model.element_conversion.convert_position_elements(
        position,
        model.element_conversion.cartesian_position_type,
        model.element_conversion.geodetic_position_type,
        bodies.get("Earth").shape_model,
        1e-6,
    )
    h, lat, lon = geodetic
    normals[name] = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])

apollo = manifest["apollo"]
points = model.read_points(
    DATA / apollo["file"], apollo["sha256"], (date(2017, 1, 1), date(2018, 1, 1)), 311
)
for point, line in zip(points, apollo["source_line_numbers"]):
    point.update(station="APOLLO", source_file=apollo["file"], source_line=line)
for row in csv.DictReader((DATA / manifest["grasse"]["inventory"]["file"]).open(newline="")):
    day, clock = row["utc_launch"].split("T")
    yy, mm, dd = map(int, day.split("-"))
    hh, minute, sec = clock.split(":")
    points.append(
        dict(
            station="GRASSE",
            source_file=row["source_file"],
            source_line=int(row["source_line"]),
            utc_label=row["utc_launch"],
            utc=model.time.DateTime(
                yy, mm, dd, int(hh), int(minute), float(sec)
            ).to_epoch_time_object(),
            observed_rtt_s=float(row["two_way_time_of_flight_s"]),
            reflector=int(row["reflector_code"]),
        )
    )
assert len(points) == 1582
rows = []
for point in points:
    station = point["station"]
    position = positions[station]
    transmit = converter.convert_time_object(
        model.time.utc_scale, model.time.tdb_scale, point["utc"], position
    )
    bounce = float(transmit) + 0.5 * point["observed_rtt_s"]
    target = model.barycentric_state("Moon", bounce)[:3]
    elevations = []
    for epoch in (float(transmit), float(transmit) + point["observed_rtt_s"]):
        matrix = rotation.body_fixed_to_inertial_rotation(epoch)
        observer = model.barycentric_state("Earth", epoch)[:3] + matrix @ position
        direction = target - observer
        elevations.append(
            float(
                np.rad2deg(
                    np.arcsin(
                        np.dot(matrix @ normals[station], direction) / np.linalg.norm(direction)
                    )
                )
            )
        )
    rows.append(
        dict(
            station=station,
            source_file=point["source_file"],
            source_line=point["source_line"],
            utc_launch=point["utc_label"],
            reflector=model.REFLECTORS[point["reflector"]][0],
            observed_rtt_s=point["observed_rtt_s"],
            transmit_tdb_s=float(transmit),
            approximate_bounce_tdb_s=bounce,
            approximate_transmit_center_elevation_deg=elevations[0],
            approximate_receive_center_elevation_deg=elevations[1],
            minimum_center_elevation_deg=min(elevations),
            candidate_at_or_below4deg=min(elevations) <= 4,
        )
    )
OUT.mkdir(exist_ok=False)
with (OUT / "all_source_geometry.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
candidates = [r for r in rows if r["candidate_at_or_below4deg"]]
summary = dict(
    source_count=len(rows),
    station_source_counts={s: sum(r["station"] == s for r in rows) for s in profiles},
    candidates=candidates,
    minimum_rows=sorted(rows, key=lambda r: r["minimum_center_elevation_deg"])[:8],
    assumptions="DirectDE440 coherentkernels, existingnative GCRS-to-ITRS ICRF rotation and WGS84 geodeticUp; publishednominalJ2000station XYZ, no deformation/drift/spatialmap. Bounce=transmit+observedRTT/2; receive=transmit+observedRTT. Mooncenteronly;4deg threshold leaves1deg margin above FCUL3deg domain. No observations solved or atmosphericcorrection evaluated.",
    kernel_binary_sha256=model.sha256(Path(model.kernel.__file__)),
    kernel_sha256=kernels,
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
)
(OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
