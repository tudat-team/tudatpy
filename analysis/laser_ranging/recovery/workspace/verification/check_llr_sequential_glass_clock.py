"""Separate endpoint clock arithmetic from archived sequential glass increments."""

from pathlib import Path
import csv
import importlib.util
import json
import sys

import numpy as np
import tudatpy.kernel
from tudatpy.astro import time_representation as time

WORK = Path("/home/dominic/Tudat/laser-ranging-work")
SCRIPT = Path(
    "/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-ranging/analysis/laser_ranging/llr_residuals.py"
)
sys.path.insert(0, str(SCRIPT.parent))
spec = importlib.util.spec_from_file_location("llr_glass_clock", SCRIPT)
llr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(llr)
manifest = json.loads(SCRIPT.with_name("llr_dataset.json").read_text())
llr.spice.clear_kernels()
for record in manifest["kernels"]:
    llr.spice.load_kernel(str(WORK / "data/llr/kernels" / record["file"]))
points = llr.read_points(
    WORK / "data/llr" / manifest["normal_points"]["file"], manifest["normal_points"]["sha256"]
)
bodies, position, _ = llr.create_environment(points, 30, "survey2005")
converter = bodies.get("Earth").rotation_model.angles_calculator.time_scale_converter
source = WORK / "results/apollo-20170309-cubic30"
summary = json.loads((source / "summary.json").read_text())
np.testing.assert_array_equal(position, summary["station_xyz_m"])
events = json.loads((source / "event_diagnostics.json").read_text())
rows = list(csv.DictReader((source / "normal_points.csv").open(newline="")))
C = llr.constants.SPEED_OF_LIGHT


def clock_interval(event):
    """Match native nWayRange double endpoint inputs with precise offset subtraction."""
    tx, rx = event["link_end_epochs_tdb_s"][0], event["link_end_epochs_tdb_s"][-1]
    return float(
        converter.get_time_object_difference(
            time.tdb_scale, time.tai_scale, time.Time(rx), position
        )
        - converter.get_time_object_difference(
            time.tdb_scale, time.tai_scale, time.Time(tx), position
        )
    )


checks = []
for event, row in zip(events, rows):
    if row["reflector"] != "A14":
        continue
    a, b = event["atmospheric"], event["complete"]
    fixed = C / 2 * (b["uplink_corrections_s"][-1] + b["downlink_corrections_s"][-1])
    sequential = float(row["glass_increment_m"])
    clock_difference = (clock_interval(b) - clock_interval(a)) * C / 2
    checks.append(
        {
            "index": event["index"],
            "sequential_glass_increment_m": sequential,
            "fixed_event_glass_half_range_m": fixed,
            "sequential_minus_fixed_m": sequential - fixed,
            "native_clock_increment_change_m": clock_difference,
            "remainder_after_clock_increment_change_m": sequential - fixed - clock_difference,
            "atmospheric_event_doubles": a["link_end_epochs_tdb_s"],
            "complete_event_doubles": b["link_end_epochs_tdb_s"],
        }
    )
result = {
    "kernel_sha256": llr.sha256(Path(tudatpy.kernel.__file__)),
    "station_xyz_m": position.tolist(),
    "event_source": str(source / "event_diagnostics.json"),
    "rows": checks,
    "scope": "Archived A14 event/correction evaluation only; no range re-solves, tolerance changes or residual adjustment.",
    "explanation": "Native nWayRange uses saved double endpoint times for clock offsets. The finalA14 receive double shifts oneULP between sequential solves. Native1ps inverse stopping changes its clock offset discontinuously; measured clockincrement change is compared to the sequentialglass discrepancy. Remaining geometry/arithmetic is bounded by the displayed remainder.",
}
(WORK / "logs/llr-glass-clock-check.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
