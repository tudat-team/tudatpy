"""Compare original ASCII Eq8 PA-to-ICRF rotations with two native/SPICE matrices."""

from pathlib import Path
import hashlib
import json
import sys
import mpmath as mp
import numpy as np
import spiceypy as spice
from jplephem.pck import PCK

root = Path("/home/dominic/Tudat/laser-ranging-work")
source = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-data/analysis/laser_ranging")
sys.path.insert(0, str(source))
import llr_residuals as model

mp.mp.dps = 80
manifest = json.loads((source / "llr_dataset.json").read_text())
model.load_ephemeris_kernels(root / "data/llr/kernels", manifest, manifest["analysis_frame_alias"])
spice.kclear()
for item in manifest["kernels"]:
    spice.furnsh(str(root / "data/llr/kernels" / item["file"]))
spice.furnsh(str(source / "icrf_spice_alias.tf"))
points = model.read_points(
    root / "data/llr" / manifest["normal_points"]["file"], manifest["normal_points"]["sha256"]
)
bodies, _, _ = model.create_environment(points, station_profile="inpop17a-gir")
selected = []
for case, directory in [
    ("apollo", "apollo-20170309-survey-stable-cubic30"),
    ("matera", "matera-20220821-stable-cubic30"),
]:
    event = json.loads((root / "results" / directory / "event_diagnostics.json").read_text())[0]
    hours, seconds = event["launch_tdb_time_split"]
    epoch = float(mp.mpf(hours) * 3600 + mp.mpf(seconds))
    selected.append((case, epoch))
header = (root / "data/llr/header.440").read_text()
indices, counts, intervals = [
    list(map(int, row.split()))
    for row in header.split("GROUP   1050")[1].split("GROUP")[0].splitlines()
    if row.strip()
]
records = []
ascii_path = root / "research/de440-coefficients/ascp01950.440"
with ascii_path.open() as stream:
    while line := stream.readline():
        if not line.strip():
            continue
        number, length = map(int, line.split())
        fields = "".join(stream.readline() for _ in range((length + 2) // 3)).split()[:length]
        first, last = [mp.mpf(v.replace("D", "e")) for v in fields[:2]]
        if any(first <= mp.mpf("2451545") + mp.mpf(epoch) / 86400 < last for _, epoch in selected):
            records.append([mp.mpf(v.replace("D", "e")) for v in fields])


def clenshaw(coefficients, argument):
    """Evaluate ascending Chebyshev coefficients with 80-digit arithmetic."""
    prior = earlier = mp.mpf(0)
    for coefficient in coefficients[:0:-1]:
        value = coefficient + 2 * argument * prior - earlier
        earlier, prior = prior, value
    return coefficients[0] + argument * prior - earlier


def ascii_angles(epoch):
    """Select quantity13 using the original DE440 ASCII header pointers."""
    day = mp.mpf("2451545") + mp.mpf(epoch) / 86400
    record = next(row for row in records if row[0] <= day < row[1])
    duration = (record[1] - record[0]) / intervals[12]
    subinterval = int(mp.floor((day - record[0]) / duration))
    argument = 2 * (day - record[0] - subinterval * duration) / duration - 1
    start = indices[12] - 1 + 3 * counts[12] * subinterval
    return [
        clenshaw(record[start + axis * counts[12] : start + (axis + 1) * counts[12]], argument)
        for axis in range(3)
    ]


def rotation(angle, axis):
    """Use the passive right-handed matrices defined by Park equations9/11."""
    cosine, sine = mp.cos(angle), mp.sin(angle)
    if axis == 3:
        return mp.matrix([[cosine, sine, 0], [-sine, cosine, 0], [0, 0, 1]])
    return mp.matrix([[1, 0, 0], [0, cosine, sine], [0, -sine, cosine]])


def difference_metrics(first, second):
    """Report matrix error and equivalent small-angle norm, without fitting axes."""
    delta = first - second
    frobenius = mp.sqrt(sum(value * value for value in delta))
    return dict(
        maximum_element=float(max(abs(value) for value in delta)),
        frobenius=float(frobenius),
        small_angle_norm_rad=float(frobenius / mp.sqrt(2)),
        surface_bound_m=float(frobenius * 1738000),
    )


pck = PCK.open(str(root / "data/llr/kernels/moon_pa_de440_200625.bpc"))
report = dict(
    kernel_sha256=model.sha256(
        Path("/home/dominic/Tudat/laser-ranging-build/src/tudatpy/kernel.so")
    ),
    ascii_sha256=model.sha256(ascii_path),
    header_sha256=model.sha256(root / "data/llr/header.440"),
    pck_segment_frames=[
        dict(frame_class_id=s.body, base_frame_id=s.frame, type=s.data_type) for s in pck.segments
    ],
    native_type=type(bodies.get("Moon").rotation_model).__name__,
    cases=[],
)
for case, epoch in selected:
    phi, theta, psi = ascii_angles(epoch)
    reference = rotation(-phi, 3) * rotation(-theta, 1) * rotation(-psi, 3)
    native = mp.matrix(
        bodies.get("Moon").rotation_model.body_fixed_to_inertial_rotation(epoch).tolist()
    )
    direct = mp.matrix(spice.pxform("MOON_PA_DE440", "J2000", epoch).tolist())
    alias = mp.matrix(spice.pxform("MOON_PA_DE440", "ICRF", epoch).tolist())
    report["cases"].append(
        dict(
            case=case,
            epoch_tdb_s=epoch,
            ascii_angles_rad=list(map(float, (phi, theta, psi))),
            ascii_eq8_vs_native=difference_metrics(reference, native),
            ascii_eq8_vs_spice=difference_metrics(reference, direct),
            native_vs_spice=difference_metrics(native, direct),
            alias_vs_j2000=difference_metrics(alias, direct),
            native_orthogonality_max=float(
                np.max(
                    np.abs(
                        np.asarray(native.tolist(), dtype=float).T
                        @ np.asarray(native.tolist(), dtype=float)
                        - np.eye(3)
                    )
                )
            ),
        )
    )
report["probe_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(root / "research/llr-physical-terms/full_orientation_audit.json").write_text(
    json.dumps(report, indent=2) + "\n"
)
print(json.dumps(report, indent=2))
