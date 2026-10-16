"""Audit existing lunar elastic tides at all eleven saved native reflection events.

Independent Legendre-polynomial displacement: JPL LLR Model 2020 Eq.10,
section 3.12; degree three is the corresponding degree-n tidal potential.
This script changes no residual, station coordinate, tide convention or fit.
"""

from pathlib import Path
import hashlib
import importlib.util
import json
import subprocess

import numpy as np
from numpy.polynomial.legendre import Legendre
import tudatpy.kernel
from tudatpy import constants
from tudatpy.astro import time_representation as time
from tudatpy.dynamics import environment
from tudatpy.interface import spice

WORK = Path("/home/dominic/Tudat/laser-ranging-work")
SOURCE = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-ranging")
SCRIPT = SOURCE / "analysis/laser_ranging/llr_residuals.py"
spec = importlib.util.spec_from_file_location("llr_audit_environment", SCRIPT)
llr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(llr)
manifest = json.loads(SCRIPT.with_name("llr_dataset.json").read_text())
points = llr.read_points(
    WORK / "data/llr" / manifest["normal_points"]["file"], manifest["normal_points"]["sha256"]
)
spice.clear_kernels()
loaded = []
for record in manifest["kernels"]:
    path = WORK / "data/llr/kernels" / record["file"]
    assert llr.sha256(path) == record["sha256"]
    spice.load_kernel(str(path))
    loaded.append({"file": str(path), "sha256": record["sha256"]})

factory_calls = []
factory = llr.setup.shape_deformation.basic_solid_body_tidal


def capture_factory(raising_bodies, love_numbers, radius):
    """Record the actual five supplied constants while calling the existing factory."""
    factory_calls.append(
        {"bodies": raising_bodies, "love_numbers": love_numbers, "radius_m": radius}
    )
    return factory(raising_bodies, love_numbers, radius)


llr.setup.shape_deformation.basic_solid_body_tidal = capture_factory
try:
    bodies, earth_station = llr.create_environment(points, 30)
finally:
    llr.setup.shape_deformation.basic_solid_body_tidal = factory
assert factory_calls == [
    {
        "bodies": ["Earth", "Sun"],
        "love_numbers": {2: (0.0423, 0.0107), 3: (0.0234, 0.0030)},
        "radius_m": 1738000.0,
    }
]
gm = {name: bodies.get(name).gravitational_parameter for name in ["Moon", "Earth", "Sun"]}
assert all(gm[name] == spice.get_body_gravitational_parameter(name) for name in gm)
moon = bodies.get("Moon")
events_path = WORK / "results/apollo-20170309-inpop17a-gir-cubic30/event_diagnostics.json"
events = json.loads(events_path.read_text())
assert len(events) == len(points) == 11


def independent_displacement(nominal, epoch):
    """Sum radial/tangential gradients of P_n using independently evaluated polynomials."""
    unit = nominal / np.linalg.norm(nominal)
    rotation = np.asarray(
        spice.compute_rotation_matrix_between_frames("J2000", "MOON_PA_DE440", float(epoch))
    )
    moon_state = np.asarray(moon.state_in_base_frame_from_ephemeris(epoch)).reshape(6)
    contributions = {}
    for raising in ["Earth", "Sun"]:
        vector = rotation @ (
            np.asarray(bodies.get(raising).state_in_base_frame_from_ephemeris(epoch)).reshape(6)[:3]
            - moon_state[:3]
        )
        distance = np.linalg.norm(vector)
        direction = vector / distance
        cosine = float(unit @ direction)
        for degree, h, l in [(2, 0.0423, 0.0107), (3, 0.0234, 0.0030)]:
            polynomial = Legendre.basis(degree)
            scale = gm[raising] / gm["Moon"] * 1738000.0 * (1738000.0 / distance) ** (degree + 1)
            contributions[f"{raising}_degree_{degree}"] = scale * (
                h * polynomial(cosine) * unit
                + l * polynomial.deriv()(cosine) * (direction - cosine * unit)
            )
    return sum(contributions.values()), contributions, rotation


rows = []
for point, event in zip(points, events):
    assert point["index"] == event["index"]
    name, xyz = llr.REFLECTORS[point["reflector"]]
    xyz = np.asarray(xyz)
    t = event["complete"]["link_end_epochs_tdb_s"]
    assert t[1] == t[2]
    epoch = time.Time(t[1])
    expected, components, rotation = independent_displacement(xyz, epoch)
    intrinsic = moon.get_ground_station(name + "_intrinsic").station_state
    np.testing.assert_array_equal(intrinsic.cartesian_position_at_reference_epoch, xyz)
    native = np.asarray(intrinsic.get_cartesian_position(epoch, "Moon")).reshape(3) - xyz
    rotation_difference = (
        np.asarray(moon.rotation_model.inertial_to_body_fixed_rotation(epoch)) - rotation
    )
    np.testing.assert_allclose(native, expected, rtol=0, atol=3e-10)
    np.testing.assert_allclose(rotation_difference, 0, rtol=0, atol=1e-14)

    # Check the final reflector includes this tide once, followed by the stated spatial terms.
    to_inertial = rotation.T
    deformed_inertial = to_inertial @ (xyz + expected)
    moon_direct = llr.barycentric_state("Moon", epoch)
    solar_distance = np.linalg.norm(llr.barycentric_state("Sun", epoch)[:3] - moon_direct[:3])
    spatial = (
        -(1.4825e-8 + gm["Sun"] / solar_distance / constants.SPEED_OF_LIGHT**2) * deformed_inertial
    )
    spatial -= (
        0.5 * (moon_direct[3:] @ deformed_inertial) * moon_direct[3:] / constants.SPEED_OF_LIGHT**2
    )
    full_expected = xyz + expected + rotation @ spatial
    full_native = np.asarray(
        moon.get_ground_station(name).station_state.get_cartesian_position(epoch, "Moon")
    ).reshape(3)
    np.testing.assert_allclose(full_native, full_expected, rtol=0, atol=5e-10)
    # Half an event-output ULP cannot materially change this slowly varying displacement.
    half_ulp = float(np.spacing(t[1]) / 2)
    minus = independent_displacement(xyz, epoch - half_ulp)[0]
    plus = independent_displacement(xyz, epoch + half_ulp)[0]
    rows.append(
        {
            "index": point["index"],
            "reflector": name,
            "reflection_tdb_seconds": t[1],
            "nominal_pa_xyz_m": xyz.tolist(),
            "native_tide_pa_m": native.tolist(),
            "independent_tide_pa_m": expected.tolist(),
            "difference_pa_m": (native - expected).tolist(),
            "components_pa_m": {key: value.tolist() for key, value in components.items()},
            "native_total_position_error_m": (full_native - full_expected).tolist(),
            "rotation_max_abs_difference": float(np.max(abs(rotation_difference))),
            "event_output_half_ulp_s": half_ulp,
            "tide_change_across_full_event_ulp_m": float(np.linalg.norm(plus - minus)),
        }
    )

# Trace the actual native deformation callback's Time routing without changing its states.
trace = []
for name in ["Moon", "Earth", "Sun"]:
    body = bodies.get(name)
    old = body.ephemeris

    def forwarding(epoch, original=old, body_name=name):
        """Observe incoming native Time while forwarding unchanged ephemeris states."""
        trace.append((body_name, epoch))
        return original.cartesian_state(epoch)

    body.ephemeris = environment.CustomEphemerisTimeObject(forwarding, "SSB", "J2000")
    bodies.add_body(body, name)
probe = time.Time(rows[0]["reflection_tdb_seconds"]) + 20e-9
moon.get_ground_station(rows[0]["reflector"] + "_intrinsic").station_state.get_cartesian_position(
    probe, "Moon"
)
assert trace and {name for name, _ in trace} == {"Moon", "Earth", "Sun"}
assert all(isinstance(epoch, time.Time) and epoch == probe for _, epoch in trace)

result = {
    "source_revision": subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=SOURCE, text=True
    ).strip(),
    "native_kernel_sha256": llr.sha256(Path(tudatpy.kernel.__file__)),
    "analysis_script_sha256": llr.sha256(SCRIPT),
    "audit_script_sha256": llr.sha256(Path(__file__)),
    "event_input": {
        "file": str(events_path),
        "sha256": llr.sha256(events_path),
        "precision": "Saved native events are doubles. Evaluated at Time(double); full-ULP tidal sensitivity is recorded.",
    },
    "primary_source": {
        "url": "https://ssd.jpl.nasa.gov/ftp/eph/planets/ioms/LLR_Model_2020_DR.pdf",
        "cached_pdf_sha256": llr.sha256(WORK / "research/JPL-LLR-Model-2020.pdf"),
        "sections": "3.9 Eq10 constant part unchanged; 3.11 T2 in TDB; 3.12 lunar elastic tides; 3.13 thermal; Appendix A2 constant tide",
    },
    "loaded_kernels": loaded,
    "factory_calls": factory_calls,
    "gm_si": gm,
    "formula": "sum_n,b (GM_b/GM_M) R (R/d_b)^(n+1) [h_n P_n(c) u + l_n P_n_prime(c) (q_b-c u)]",
    "native_routing_trace": [
        {"body": name, "time_split": epoch.__getstate__(), "equals_precise_probe": epoch == probe}
        for name, epoch in trace
    ],
    "rows": rows,
    "max_native_tide_vector_error_m": max(
        float(np.linalg.norm(row["difference_pa_m"])) for row in rows
    ),
    "max_final_reflector_position_error_m": max(
        float(np.linalg.norm(row["native_total_position_error_m"])) for row in rows
    ),
    "max_event_output_full_ulp_tide_change_m": max(
        row["tide_change_across_full_event_ulp_m"] for row in rows
    ),
    "constant_tide": "Full equilibrium degree2/3 displacement retained, no mean subtraction; final reflector matches one tide plus spatial terms.",
    "omitted_terms_source_sizes": {
        "dissipation": "about2mm or less, JPL3.12",
        "Apollo_thermal": "1–1.5mm, JPL3.13",
        "Lunokhod_thermal": "about5mm vertical with horizontal structure, JPL3.13",
        "lunar_pole_tide": "<0.2mm and omitted by JPL3.12",
    },
    "concrete_bug_found": False,
}
result["component_summaries"] = {}
for name in sorted({row["reflector"] for row in rows}):
    selected = [row for row in rows if row["reflector"] == name]
    summary = {}
    for term in selected[0]["components_pa_m"]:
        norms = [1000 * np.linalg.norm(row["components_pa_m"][term]) for row in selected]
        summary[term] = {"norm_mm_min": min(norms), "norm_mm_max": max(norms)}
    result["component_summaries"][name] = {"count": len(selected), "components": summary}
out = WORK / "logs/lunar-elastic-tide-audit.json"
out.write_text(json.dumps(result, indent=2) + "\n")
print(
    json.dumps(
        {
            key: value
            for key, value in result.items()
            if key not in ["rows", "loaded_kernels", "native_routing_trace"]
        },
        indent=2,
    )
)
