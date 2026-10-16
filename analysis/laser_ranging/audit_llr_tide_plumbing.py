"""Check actual native lunar displacement parameters without solving new observations."""

import argparse
import json
from pathlib import Path
import subprocess
from unittest.mock import patch
from decimal import Decimal, localcontext

import numpy as np
from numpy.polynomial.legendre import Legendre
from tudatpy import constants, kernel
from tudatpy.astro import time_representation as time
from tudatpy.interface import spice

import llr_residuals as model

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument("--data-dir", type=Path, required=True)
parser.add_argument("--kernel-dir", type=Path, required=True)
parser.add_argument("--event-diagnostics", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()

TREE = Path(__file__).resolve().parents[2]


def original_constants(path):
    """Derive SI gravitational parameters and reference radius from the original JPL header."""

    text = path.read_text()

    names = text.split("GROUP   1040")[1].split("GROUP")[0].split()[1:]

    values = text.split("GROUP   1041")[1].split("GROUP")[0].split()[1:]

    source = dict(zip(names, [Decimal(v.replace("D", "e")) for v in values]))

    with localcontext() as context:
        context.prec = 60

        conversion = (source["AU"] * 1000) ** 3 / Decimal(86400) ** 2

        earth = source["GMB"] * conversion / (1 + 1 / source["EMRAT"])

        gm = dict(
            Sun=float(source["GMS"] * conversion),
            Earth=float(earth),
            Moon=float(earth / source["EMRAT"]),
        )

    return (
        gm,
        float(source["AM"] * 1000),
        {key: str(source[key]) for key in ["GMS", "GMB", "EMRAT", "AU", "AM", "RE"]},
    )


def independent_terms(nominal, epoch, gm, radius):
    """Use JPL Eq10 and the degree-three potential gradient with primary constants."""

    unit = nominal / np.linalg.norm(nominal)

    terms = {}

    for name in ["Earth", "Sun"]:
        vector = np.asarray(
            spice.get_body_cartesian_state_at_epoch(
                name, "Moon", "MOON_PA_DE440", "NONE", float(epoch)
            )
        ).reshape(6)[:3]

        distance = np.linalg.norm(vector)

        direction = vector / distance

        cosine = unit @ direction

        for degree, h, l in [(2, 0.0423, 0.0107), (3, 0.0234, 0.0030)]:
            polynomial = Legendre.basis(degree)

            amplitude = gm[name] / gm["Moon"] * radius * (radius / distance) ** (degree + 1)

            terms[f"{name}/{degree}"] = amplitude * (
                h * polynomial(cosine) * unit
                + l * polynomial.deriv()(cosine) * (direction - cosine * unit)
            )

    return terms


manifest = json.loads(Path(model.__file__).with_name("llr_dataset.json").read_text())

model.load_ephemeris_kernels(args.kernel_dir, manifest, manifest["analysis_frame_alias"])
points = model.read_points(
    args.data_dir / manifest["normal_points"]["file"], manifest["normal_points"]["sha256"]
)

event_path = args.event_diagnostics

events = json.loads(event_path.read_text())

epochs = [time.Time(events[index]["complete"]["link_end_epochs_tdb_s"][1]) for index in [0, 5, 10]]

gm, radius, header_constants = original_constants(args.data_dir / "header.440")

factory = model.setup.shape_deformation.basic_solid_body_tidal

shape = model.setup.shape.spherical

variants = [
    ("actual", None, 1, False, 1),
    ("degree2", 2, 1, False, 1),
    ("degree3", 3, 1, False, 1),
    ("degree2_radius2", 2, 2, False, 1),
    ("degree3_radius2", 3, 2, False, 1),
    ("shape2_explicit", 2, 1, False, 2),
    ("shape2_fallback", 2, 1, True, 2),
]

output = {}

for label, degree, radius_scale, fallback, shape_scale in variants:
    calls = []

    def settings(bodies, love, reference):
        """Record supplied settings and vary one declared native-factory parameter only."""

        selected = love if degree is None else {degree: love[degree]}

        actual_radius = float("nan") if fallback else reference * radius_scale

        calls.append(
            dict(
                bodies=bodies,
                love_numbers={str(k): list(v) for k, v in selected.items()},
                reference_radius_m=None if fallback else actual_radius,
            )
        )
        return factory(bodies, selected, actual_radius)

    with patch.object(
        model.setup.shape_deformation, "basic_solid_body_tidal", settings
    ), patch.object(model.setup.shape, "spherical", lambda r: shape(shape_scale * r)):
        bodies, _, _ = model.create_environment(points, station_profile="inpop17a-gir")

    moon = bodies.get("Moon")

    rows = []

    for code, (name, xyz) in sorted(model.REFLECTORS.items()):
        xyz = np.asarray(xyz)

        model.setup.add_ground_station(
            moon,
            model.setup.ground_station.basic_station(
                "probe_" + name,
                xyz,
                station_motion_settings=[
                    model.setup.ground_station.body_deformation_station_motion()
                ],
            ),
            bodies,
        )
        station = moon.get_ground_station("probe_" + name).station_state

        for index, epoch in enumerate(epochs):
            native = np.asarray(station.get_cartesian_position(epoch, "Moon")).reshape(3) - xyz

            components = independent_terms(xyz, epoch, gm, radius)

            expected = sum(
                v
                for key, v in components.items()
                if degree is None or key.endswith("/" + str(degree))
            )

            row = dict(
                reflector=name,
                epoch_index=index,
                native_displacement_m=native.tolist(),
                independent_displacement_m=expected.tolist(),
                components_m={key: v.tolist() for key, v in components.items()},
                native_minus_independent_m=(native - expected).tolist(),
            )

            if label == "actual" and code in {p["reflector"] for p in points}:
                intrinsic = np.asarray(
                    moon.get_ground_station(
                        name + "_intrinsic"
                    ).station_state.get_cartesian_position(epoch, "Moon")
                ).reshape(3)

                rotation = moon.rotation_model.body_fixed_to_inertial_rotation(epoch)

                state = model.barycentric_state("Moon", epoch)

                solar_potential = gm["Sun"] / np.linalg.norm(
                    model.barycentric_state("Sun", epoch)[:3] - state[:3]
                )

                inertial = rotation @ intrinsic

                spatial = (
                    -(1.4825e-8 + solar_potential / constants.SPEED_OF_LIGHT**2) * inertial
                    - 0.5 * (state[3:] @ inertial) * state[3:] / constants.SPEED_OF_LIGHT**2
                )

                final = np.asarray(
                    moon.get_ground_station(name).station_state.get_cartesian_position(
                        epoch, "Moon"
                    )
                ).reshape(3)

                row.update(
                    probe_minus_actual_intrinsic_m=(xyz + native - intrinsic).tolist(),
                    final_minus_one_tide_and_spatial_m=(
                        final - intrinsic - rotation.T @ spatial
                    ).tolist(),
                )
            rows.append(row)
    output[label] = dict(
        factory_calls=calls,
        actual_moon_shape_radius_m=moon.shape_model.average_radius,
        gm_si={n: bodies.get(n).gravitational_parameter for n in gm},
        rows=rows,
    )

# Compare native controls directly: these isolate parameter plumbing from the independent formula.
checks = {}

for label, base, factor in [
    ("degree2_radius2", "degree2", 16),
    ("degree3_radius2", "degree3", 32),
    ("shape2_explicit", "degree2", 1),
    ("shape2_fallback", "degree2", 16),
]:
    errors = [
        np.linalg.norm(
            np.asarray(a["native_displacement_m"]) - factor * np.asarray(b["native_displacement_m"])
        )
        for a, b in zip(output[label]["rows"], output[base]["rows"])
    ]

    checks[label] = dict(expected_displacement_factor=factor, max_scaling_error_m=max(errors))

checks["actual_degree_sum_max_error_m"] = max(
    np.linalg.norm(
        np.asarray(a["native_displacement_m"])
        - np.asarray(b["native_displacement_m"])
        - np.asarray(c["native_displacement_m"])
    )
    for a, b, c in zip(
        output["actual"]["rows"], output["degree2"]["rows"], output["degree3"]["rows"]
    )
)

checks["actual_independent_max_error_m"] = max(
    np.linalg.norm(row["native_minus_independent_m"]) for row in output["actual"]["rows"]
)

checks["actual_probe_intrinsic_max_error_m"] = max(
    np.linalg.norm(row["probe_minus_actual_intrinsic_m"])
    for row in output["actual"]["rows"]
    if "probe_minus_actual_intrinsic_m" in row
)

checks["actual_once_only_max_error_m"] = max(
    np.linalg.norm(row["final_minus_one_tide_and_spatial_m"])
    for row in output["actual"]["rows"]
    if "final_minus_one_tide_and_spatial_m" in row
)

result = dict(
    analysis_source_revision=subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=TREE, text=True
    ).strip(),
    native_kernel_sha256=model.sha256(Path(kernel.__file__)),
    probe_sha256=model.sha256(Path(__file__)),
    header_sha256=model.sha256(args.data_dir / "header.440"),
    production_script_sha256=model.sha256(Path(model.__file__)),
    event_source_sha256=model.sha256(event_path),
    selected_event_indices=[0, 5, 10],
    original_header_constants=header_constants,
    original_gm_si=gm,
    original_AM_radius_m=radius,
    checks=checks,
    variants=output,
)

args.output.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({k: v for k, v in result.items() if k != "variants"}, indent=2))
