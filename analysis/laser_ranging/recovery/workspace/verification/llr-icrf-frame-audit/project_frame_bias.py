"""Project only the Earth rotation frame-bias change at frozen LLR event epochs."""

from pathlib import Path
import hashlib
import importlib.util
import json
import sys
import numpy as np
from tudatpy.astro import time_representation as time
from tudatpy.interface import spice

ROOT = Path("/home/dominic/Tudat/laser-ranging-work")
SOURCE = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-ranging/analysis/laser_ranging")
DIAGNOSTIC = Path(__file__).parent
sys.path.insert(0, str(SOURCE))


def module(name, path):
    """Load the two exact source versions without modifying either one."""
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


old = module("frozen_llr", SOURCE / "llr_residuals.py")
new = module("icrf_llr", DIAGNOSTIC / "llr_residuals.py")
manifest = json.loads((SOURCE / "llr_dataset.json").read_text())
spice.clear_kernels()
for record in manifest["kernels"]:
    path = ROOT / "data/llr/kernels" / record["file"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == record["sha256"]
    spice.load_kernel(str(path))
assert spice.get_total_count_of_kernels_loaded() == 6
spice.load_kernel(str(DIAGNOSTIC / "icrf_spice_alias.tf"))
assert spice.get_total_count_of_kernels_loaded() == 7
results = {}
for case, profile, directory in [
    ("apollo_survey", "survey2005", "apollo-20170309-survey-stable-cubic30"),
    ("apollo_inpop17a", "inpop17a-gir", "apollo-20170309-inpop17a-stable-cubic30"),
    ("matera", "survey2005", "matera-20220821-stable-cubic30"),
]:
    if case == "matera":
        import matera_case

        data = json.loads((SOURCE / "matera_dataset.json").read_text())
        points = matera_case.read_points(ROOT / "data/llr" / data["normal_points"]["file"])
        history, _ = matera_case.station_history(ROOT / "data/llr", data)
    else:
        points = old.read_points(
            ROOT / "data/llr" / manifest["normal_points"]["file"],
            manifest["normal_points"]["sha256"],
        )
        history = None
    original, xyz, _ = old.create_environment(
        points, station_profile=profile, station_history=history
    )
    revised, _, _ = new.create_environment(points, station_profile=profile, station_history=history)
    events = json.loads((ROOT / "results" / directory / "event_diagnostics.json").read_text())
    result = []
    for point, event in zip(points, events):
        transmit, bounce, _, receive = event["complete"]["link_end_epochs_tdb_s"]
        correction = []
        station_changes = []
        bias_matrices = []
        for epoch in [transmit, receive]:
            a = original.get("Earth").rotation_model.body_fixed_to_inertial_rotation(epoch)
            b = revised.get("Earth").rotation_model.body_fixed_to_inertial_rotation(epoch)
            intrinsic = (
                original.get("Earth")
                .get_ground_station(point.get("station", "APOLLO") + "_intrinsic")
                .station_state
            )
            position = np.asarray(intrinsic.get_cartesian_position(epoch, "Earth")).reshape(3)
            station_change = (b - a) @ position
            station = (
                np.asarray(
                    original.get("Earth").state_in_base_frame_from_ephemeris(time.Time(epoch))
                ).reshape(6)[:3]
                + a @ position
            )
            moon = original.get("Moon")
            moon_state = np.asarray(
                moon.state_in_base_frame_from_ephemeris(time.Time(bounce))
            ).reshape(6)[:3]
            reflector = (
                moon.rotation_model.body_fixed_to_inertial_rotation(bounce)
                @ old.REFLECTORS[point["reflector"]][1]
            )
            direction = moon_state + reflector - station
            direction /= np.linalg.norm(direction)
            correction.append(-0.5 * float(np.dot(direction, station_change)))
            station_changes.append(station_change.tolist())
            bias_matrices.append((a @ b.T).tolist())
            for body in ["Earth", "Moon", "Sun"]:
                assert np.array_equal(
                    original.get(body).state_in_base_frame_from_ephemeris(time.Time(epoch)),
                    revised.get(body).state_in_base_frame_from_ephemeris(time.Time(epoch)),
                )
            assert np.array_equal(
                original.get("Moon").rotation_model.body_fixed_to_inertial_rotation(epoch),
                revised.get("Moon").rotation_model.body_fixed_to_inertial_rotation(epoch),
            )
        result.append(
            {
                "index": point["index"],
                "reflector": old.REFLECTORS[point["reflector"]][0],
                "predicted_computed_half_range_change_m": sum(correction),
                "station_change_transmit_receive_m": station_changes,
                "old_from_new_matrix_transmit_receive": bias_matrices,
            }
        )
    results[case] = result
    print(
        case,
        "projected change m",
        [r["predicted_computed_half_range_change_m"] for r in result],
        flush=True,
    )
output = {
    "kernel_binary_sha256": hashlib.sha256(
        Path(__import__("tudatpy").kernel.__file__).read_bytes()
    ).hexdigest(),
    "fixed_event_approximation": "Project rotation-only change of full intrinsic station position at original event doubles; subsequent full re-solves retain all physical corrections.",
    "translation_and_moon_pa_rotation_bitwise_identical_at_all_probed_epochs": True,
    "kernel_load_count": 7,
    "cases": results,
}
(DIAGNOSTIC / "projected_frame_bias.json").write_text(json.dumps(output, indent=2) + "\n")
