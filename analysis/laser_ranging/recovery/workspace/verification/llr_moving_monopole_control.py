"""Finite-endpoint Hees2014 first-PM constant-velocity control at frozen LLR events."""

from pathlib import Path
import json
import sys
import subprocess
import mpmath as mp
import numpy as np
from tudatpy import kernel
from tudatpy.astro import time_representation as time

TREE = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-data")
WORK = Path("/home/dominic/Tudat/laser-ranging-work")
sys.path.insert(0, str(TREE / "analysis/laser_ranging"))
import llr_residuals as model
import llr_long_arc as arc

mp.mp.dps = 80
C = mp.mpf(model.constants.SPEED_OF_LIGHT)


def dot(a, b):
    """Retain arbitrary precision during short finite-path vector operations."""

    return sum(x * y for x, y in zip(a, b))


def norm(a):
    """Evaluate an ordinary Euclidean vector norm."""

    return mp.sqrt(dot(a, a))


def vec(values):
    """Represent the exact supplied binary doubles at the selected decimal precision."""

    return [mp.mpf(float(v)) for v in values]


def static(xa, xb, centre, gm):
    """Evaluate the native static monopole logarithm as a length, gammaPPN=1."""

    a = norm([x - y for x, y in zip(xa, centre)])
    b = norm([x - y for x, y in zip(xb, centre)])
    r = norm([x - y for x, y in zip(xb, xa)])
    return 2 * gm / C**2 * mp.log((a + b + r) / (a + b - r))


def moving(xa, xb, centre, velocity, tb_minus_t0, gm, quadrature=False):
    """Use Hees Eqs26/30/45 and, independently, the null-segment potential integral23."""

    displacement = [b - a for a, b in zip(xa, xb)]
    r = norm(displacement)
    n = [x / r for x in displacement]
    beta = [v / C for v in velocity]
    lorentz = 1 / mp.sqrt(1 - dot(beta, beta))
    factor = lorentz * (1 - dot(beta, n))
    k = lorentz**2 / (1 + lorentz)
    vectors = []
    for x in [xa, xb]:
        relative = [a - b for a, b in zip(x, centre)]
        vectors.append(
            [
                relative[i]
                + k * beta[i] * dot(beta, relative)
                - lorentz * velocity[i] * tb_minus_t0
                for i in range(3)
            ]
        )
    transformed_a = [vectors[0][i] + lorentz * beta[i] * r for i in range(3)]
    transformed_b = vectors[1]
    ar = norm(transformed_a)
    br = norm(transformed_b)
    transformed_r = r * factor
    result = 2 * gm / C**2 * factor * mp.log((ar + br + transformed_r) / (ar + br - transformed_r))
    integral_error = None
    if quadrature:
        g = [b - a for a, b in zip(transformed_a, transformed_b)]
        integral = mp.quad(
            lambda s: 1 / norm([transformed_b[i] - s * g[i] for i in range(3)]),
            [0, mp.mpf(".01"), mp.mpf(".1"), mp.mpf(".5"), mp.mpf(".9"), mp.mpf(".99"), 1],
        )
        direct = 2 * r * gm / C**2 * factor**2 * integral
        integral_error = direct - result
    return result, factor, integral_error


def endpoint(bodies, body, station, epoch):
    """Reconstruct the complete native SSB reference point, including terminal relativistic mapping."""

    obj = bodies.get(body)
    fixed = np.asarray(
        obj.get_ground_station(station).station_state.get_cartesian_position(
            time.Time(epoch), "SSB"
        )
    ).reshape(3)
    centre = np.asarray(obj.state_in_base_frame_from_ephemeris(time.Time(epoch))).reshape(6)[:3]
    return centre + obj.rotation_model.body_fixed_to_inertial_rotation(time.Time(epoch)) @ fixed


manifest = json.loads(Path(model.__file__).with_name("llr_dataset.json").read_text())
model.load_ephemeris_kernels(WORK / "data/llr/kernels", manifest, manifest["analysis_frame_alias"])
points = model.read_points(
    WORK / "data/llr" / manifest["normal_points"]["file"], manifest["normal_points"]["sha256"]
)
apollo_events_path = WORK / "results/llr-six-de440-inpop17a-attempt1/event_diagnostics.json"
apollo_events = json.loads(apollo_events_path.read_text())
apollo_bodies, _, _ = model.create_environment(points, station_profile="inpop17a-gir")
cases = []
for index in [0, 10]:
    cases.append(
        dict(
            station="APOLLO",
            point=points[index],
            bodies=apollo_bodies,
            epochs=apollo_events[index]["complete"]["link_end_epochs_tdb_s"],
            saved=apollo_events[index]["complete"],
        )
    )
groups = arc.read_groups(WORK / "research/longarc-20170304", json.loads(arc.MANIFEST.read_text()))
points, weather, profile, _ = next(
    group for group in groups if group[3]["group"] == "apollo15_201703.npt:29"
)
grasse_bodies, _, _ = model.create_environment(
    points, station_profile=profile, meteorological_history=weather
)
grasse_events_path = (
    WORK / "results/llr-longarc-20170301-20170430-attempt2/extreme_precision_control.json"
)
grasse_events = json.loads(grasse_events_path.read_text())
point = next(p for p in points if p["source_line"] == 171)
event = next(
    r for r in grasse_events["rows"] if r["source_line"] == 171 and r["spice_knot_step_s"] == 30
)
cases.append(
    dict(
        station="GRASSE",
        point=point,
        bodies=grasse_bodies,
        epochs=event["native_event_epochs_tdb_s"],
        saved=None,
    )
)
rows = []
legs = []
for case in cases:
    bodies = case["bodies"]
    reflector = model.REFLECTORS[case["point"]["reflector"]][0]
    epochs = case["epochs"]
    states = [
        endpoint(bodies, "Earth", case["station"], epochs[0]),
        endpoint(bodies, "Moon", reflector, epochs[1]),
        endpoint(bodies, "Moon", reflector, epochs[2]),
        endpoint(bodies, "Earth", case["station"], epochs[3]),
    ]
    for leg, (i, j, tx, rx) in enumerate([(0, 1, "Earth", "Moon"), (2, 3, "Moon", "Earth")]):
        xa, xb = vec(states[i]), vec(states[j])
        ta, tb = vec([epochs[i], epochs[j]])
        subset = []
        for body in ["Sun", "Earth", "Moon"]:
            fraction = 0 if body == tx else 1 if body == rx else 0.5
            t0 = epochs[i] + fraction * (epochs[j] - epochs[i])
            state = np.asarray(
                bodies.get(body).state_in_base_frame_from_ephemeris(time.Time(t0))
            ).reshape(6)
            centre, velocity = vec(state[:3]), vec(state[3:])
            gm = mp.mpf(bodies.get(body).gravitational_parameter)
            frozen = static(xa, xb, centre, gm)
            exact, factor, integral_error = moving(
                xa, xb, centre, velocity, tb - mp.mpf(t0), gm, True
            )
            zero, _, _ = moving(xa, xb, centre, [mp.mpf(0)] * 3, tb - mp.mpf(t0), gm)
            with mp.workdps(50):
                fifty, _, _ = moving(xa, xb, centre, velocity, tb - mp.mpf(t0), gm)
            row = dict(
                station=case["station"],
                source_line=case["point"].get("source_line", case["point"]["index"]),
                reflector=reflector,
                leg="uplink" if leg == 0 else "downlink",
                body=body,
                endpoint_epochs_tdb_s=[epochs[i], epochs[j]],
                endpoint_bcrs_m=[states[i].tolist(), states[j].tolist()],
                body_reference_epoch_tdb_s=t0,
                body_reference_bcrs_state=state.tolist(),
                gm_si=float(gm),
                static_delay_length_m=float(frozen),
                moving_delay_length_m=float(exact),
                moving_minus_static_m=float(exact - frozen),
                prefactor_only_minus_static_m=float((factor - 1) * frozen),
                transformed_endpoint_remainder_m=float(exact - factor * frozen),
                heuristic_speed_over_c_times_static_m=float(norm(velocity) / C * abs(frozen)),
                quadrature_minus_closed_form_m=float(integral_error),
                zero_velocity_minus_static_m=float(zero - frozen),
                precision50_minus80_m=float(fifty - exact),
                recorded_dt_minus_null_geometric_time_s=float(
                    tb - ta - norm([b - a for a, b in zip(xa, xb)]) / C
                ),
            )
            rows.append(row)
            subset.append(row)
        total = sum(r["static_delay_length_m"] for r in subset)
        saved = (
            None
            if case["saved"] is None
            else case["saved"]["uplink_corrections_s" if leg == 0 else "downlink_corrections_s"][0]
            * float(C)
        )
        legs.append(
            dict(
                station=case["station"],
                reflector=reflector,
                source_line=case["point"].get("source_line", case["point"]["index"]),
                leg="uplink" if leg == 0 else "downlink",
                static_delay_length_m=total,
                moving_minus_static_length_m=sum(r["moving_minus_static_m"] for r in subset),
                reconstructed_minus_saved_native_length_m=None if saved is None else total - saved,
            )
        )
result = dict(
    scope="Three frozen NPs, six existing solved legs, no new observation solve/model correction.",
    source_revision=subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=TREE, text=True
    ).strip(),
    kernel_sha256=model.sha256(Path(kernel.__file__)),
    probe_sha256=model.sha256(Path(__file__)),
    environment_script_sha256=model.sha256(Path(model.__file__)),
    event_inputs={
        str(p.relative_to(WORK)): model.sha256(p) for p in [apollo_events_path, grasse_events_path]
    },
    primary_paper={
        "url": "https://arxiv.org/pdf/1406.6600",
        "sha256": model.sha256(WORK / "research/llr-physical-plumbing/Hees-2014-moving-TTF.pdf"),
        "equations": [22, 23, 26, 30, 34, 44, 45],
    },
    arithmetic_decimal_digits=80,
    speed_of_light_m_s=float(C),
    ppn_gamma=1,
    body_motion="Constant velocity tangent to native body ephemeris at original native static-evaluation epoch (Earth endpoint, Moon reflection, Sun leg midpoint).",
    assumptions="First post-Minkowskian monopole in GR, all velocity orders, spin/multipoles/acceleration/higherG omitted. Null reference segment z0=ctB-lambdaR and z=xB-lambda(xB-xA); actual saved eventdt not imposed inside the Lorentz transform. Delta is a length, delayseconds=Delta/c. This is a fixed-endpoint diagnostic, not an iterated two-way observable change.",
    rows=rows,
    legs=legs,
)
(WORK / "logs/llr-moving-monopole-control.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(legs, indent=2))
print("max_closed_integral_error_m", max(abs(r["quadrature_minus_closed_form_m"]) for r in rows))
