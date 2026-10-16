"""Fixed-event Shapiro source/roundoff audit; no residual solves or model changes."""

from pathlib import Path
import importlib.util
import json
import math
import subprocess
import sys
import mpmath as mp
import numpy as np
import tudatpy.kernel
from tudatpy.astro import time_representation as time

WORK = Path("/home/dominic/Tudat/laser-ranging-work")
SOURCE = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-ranging")
SCRIPT = SOURCE / "analysis/laser_ranging/llr_residuals.py"
sys.path.insert(0, str(SCRIPT.parent))
spec = importlib.util.spec_from_file_location("llr_shapiro_audit", SCRIPT)
llr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(llr)
manifest = json.loads(SCRIPT.with_name("llr_dataset.json").read_text())
llr.spice.clear_kernels()
for item in manifest["kernels"]:
    path = WORK / "data/llr/kernels" / item["file"]
    assert llr.sha256(path) == item["sha256"]
    llr.spice.load_kernel(str(path))
points = llr.read_points(
    WORK / "data/llr" / manifest["normal_points"]["file"], manifest["normal_points"]["sha256"]
)
bodies, nominal, loading = llr.create_environment(points, 30, "inpop17a-gir")
assert loading is None
EVENTS = WORK / "results/apollo-20170309-inpop17a-gir-cubic30/event_diagnostics.json"
events = json.loads(EVENTS.read_text())
C = llr.constants.SPEED_OF_LIGHT
mp.mp.dps = 80


def endpoint(body, station, epoch):
    """Reconstruct the same deformed surface reference point at a saved event."""
    obj = bodies.get(body)
    fixed = np.asarray(
        obj.get_ground_station(station).station_state.get_cartesian_position(time.Time(epoch), body)
    ).reshape(3)
    return (
        np.asarray(obj.state_in_base_frame_from_ephemeris(time.Time(epoch))).reshape(6)[:3]
        + obj.rotation_model.body_fixed_to_inertial_rotation(epoch) @ fixed
    )


def delay(x1, x2, center, gm, extended=False):
    """Evaluate the documented logarithm with double or 80-digit exact-input arithmetic."""
    if extended:
        x1, x2, center = [[mp.mpf(float(v)) for v in vector] for vector in (x1, x2, center)]
        norm = lambda vector: mp.sqrt(sum(v * v for v in vector))
        r1 = norm([a - b for a, b in zip(x1, center)])
        r2 = norm([a - b for a, b in zip(x2, center)])
        rho = norm([a - b for a, b in zip(x2, x1)])
        return float(2 * mp.mpf(gm) / mp.mpf(C) ** 3 * mp.log((r1 + r2 + rho) / (r1 + r2 - rho)))
    r1 = np.linalg.norm(x1 - center)
    r2 = np.linalg.norm(x2 - center)
    rho = np.linalg.norm(x2 - x1)
    return 2 * gm / C**3 * math.log((r1 + r2 + rho) / (r1 + r2 - rho))


rows = []
for point, event in zip(points, events):
    epochs = event["complete"]["link_end_epochs_tdb_s"]
    reflector = llr.REFLECTORS[point["reflector"]][0]
    xyz = [
        endpoint("Earth", "APOLLO", epochs[0]),
        endpoint("Moon", reflector, epochs[1]),
        endpoint("Moon", reflector, epochs[2]),
        endpoint("Earth", "APOLLO", epochs[3]),
    ]
    legs = []
    for leg, (i, j, tx, rx) in enumerate([(0, 1, "Earth", "Moon"), (2, 3, "Moon", "Earth")]):
        t1, t2 = epochs[i], epochs[j]
        components = {}
        for body in ["Sun", "Earth", "Moon"]:
            fraction = 0 if body == tx else 1 if body == rx else 0.5
            epoch = t1 + fraction * (t2 - t1)
            state = np.asarray(
                bodies.get(body).state_in_base_frame_from_ephemeris(time.Time(epoch))
            ).reshape(6)
            gm = bodies.get(body).gravitational_parameter
            value = delay(xyz[i], xyz[j], state[:3], gm)
            precise = delay(xyz[i], xyz[j], state[:3], gm, True)
            alternatives = {}
            for label, f in [("tx", 0), ("mid", 0.5), ("rx", 1)]:
                center = np.asarray(
                    bodies.get(body).state_in_base_frame_from_ephemeris(
                        time.Time(t1 + f * (t2 - t1))
                    )
                ).reshape(6)[:3]
                alternatives[label] = delay(xyz[i], xyz[j], center, gm) - value
            components[body] = {
                "native_body_epoch_tdb_s": epoch,
                "epoch_fraction": fraction,
                "delay_range_m": value * C,
                "double_minus_80digit_m": (value - precise) * C,
                "alternative_frozen_body_epoch_minus_native_m": {
                    key: v * C for key, v in alternatives.items()
                },
                "velocity_over_c_times_delay_m": np.linalg.norm(state[3:]) / C * abs(value * C),
                "body_velocity_m_s": state[3:].tolist(),
            }
        total = sum(item["delay_range_m"] for item in components.values())
        saved = (
            event["complete"]["uplink_corrections_s" if leg == 0 else "downlink_corrections_s"][0]
            * C
        )
        # Distant planetary point-mass terms: meaningful missing-mm bound, not a fitted correction.
        planets = {}
        for body in ["1", "2", "4", "5", "6", "7", "8", "9"]:
            center = llr.barycentric_state(body, (t1 + t2) / 2)[:3]
            gm = llr.spice.get_body_gravitational_parameter(body)
            planets[body] = delay(xyz[i], xyz[j], center, gm) * C
        legs.append(
            {
                "leg": "uplink" if leg == 0 else "downlink",
                "components": components,
                "reconstructed_minus_saved_native_m": total - saved,
                "total_delay_range_m": total,
                "omitted_planet_delay_range_m": planets,
            }
        )
    rows.append({"index": point["index"], "reflector": reflector, "legs": legs})
all_components = [
    comp for row in rows for leg in row["legs"] for comp in leg["components"].values()
]
result = {
    "source_revision": subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=SOURCE, text=True
    ).strip(),
    "kernel_sha256": llr.sha256(Path(tudatpy.kernel.__file__)),
    "audit_script_sha256": llr.sha256(Path(__file__)),
    "environment_script_sha256": llr.sha256(SCRIPT),
    "saved_events_sha256": llr.sha256(EVENTS),
    "primary_sources": [
        "https://ssd.jpl.nasa.gov/ftp/eph/planets/ioms/LLR_Model_2020_DR.pdf section3.17 Eq19",
        "https://iers-conventions.obspm.fr/content/chapter11/tn36_c11.pdf section11.2 Eq11.17",
    ],
    "native_route": "Each leg gets actual deformed station and reflector positions at its separate events. Sun frozen at leg midpoint; Earth frozen at Earth endpoint; Moon frozen at reflection event. GM from loaded DE440. No moving-body velocity/retarded correction. Bending disabled.",
    "interpretation": "Selected body epochs are the physically closest endpoints for Earth/Moon. Midpoint is suitable for distant Sun on this short segment. Endpoint-epoch alternatives and v/c scales are diagnostics, not a claimed exact moving-monopole solution or proposed correction. IERS closest-approach VLBI prescription11.3 is for an incoming infinite ray and cannot be substituted blindly for finite LLR. JPL2020 does not specify body epochs in Eq19. Omitted planets and leading moving-body scales are below centimetres here.",
    "max_abs_double_minus_80digit_m": max(abs(x["double_minus_80digit_m"]) for x in all_components),
    "max_abs_saved_reconstruction_m": max(
        abs(leg["reconstructed_minus_saved_native_m"]) for row in rows for leg in row["legs"]
    ),
    "max_alternative_epoch_change_m": max(
        abs(v)
        for x in all_components
        for v in x["alternative_frozen_body_epoch_minus_native_m"].values()
    ),
    "max_sum_velocity_over_c_scale_m": max(
        sum(x["velocity_over_c_times_delay_m"] for x in leg["components"].values())
        for row in rows
        for leg in row["legs"]
    ),
    "max_sum_omitted_planets_one_way_m": max(
        sum(leg["omitted_planet_delay_range_m"].values()) for row in rows for leg in row["legs"]
    ),
    "rows": rows,
}
(WORK / "logs/llr-shapiro-fixed-event-audit.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({key: value for key, value in result.items() if key != "rows"}, indent=2))
for name in ["Sun", "Earth", "Moon"]:
    values = [leg["components"][name] for row in rows for leg in row["legs"]]
    print(
        name,
        "delay_range_m",
        [min(x["delay_range_m"] for x in values), max(x["delay_range_m"] for x in values)],
        "epoch_envelope_m",
        max(
            abs(v)
            for x in values
            for v in x["alternative_frozen_body_epoch_minus_native_m"].values()
        ),
    )
