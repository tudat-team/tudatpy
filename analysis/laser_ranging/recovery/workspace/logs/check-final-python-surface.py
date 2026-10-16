"""Resolve edited documentation names against the built kernel, then exercise binding routes."""

from pathlib import Path
import datetime
import hashlib
import importlib
import json
import re
import subprocess
import sys

import numpy as np
import tudatpy
import tudatpy.kernel
from tudatpy.astro import time_representation as tr
from tudatpy.data_input import environment_data as ed
from tudatpy.data_input.tracking_data import slr
from tudatpy.dynamics import environment_setup as es
from tudatpy.estimation.observable_models_setup import light_time_corrections as lc

source = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-ranging")
logs = Path("/home/dominic/Tudat/laser-ranging-work/logs")
revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip()
changed = subprocess.check_output(
    ["git", "diff", "--name-only", "ebe3115ec..HEAD", "--", "*.rst"], cwd=source, text=True
).splitlines()
resolved = []
for relative in changed:
    module = None
    autosummary = False
    for number, line in enumerate((source / relative).read_text().splitlines(), 1):
        current = re.match(r"^\.\. currentmodule::\s*(\S+)", line)
        if current:
            module = current.group(1)
        if line.startswith(".. autosummary::"):
            autosummary = True
            continue
        if autosummary and line and not line[0].isspace():
            autosummary = False
        if not autosummary or not line.strip() or line.lstrip().startswith(":"):
            continue
        name = line.strip().lstrip("~")
        name = re.sub(r".*<([^>]+)>$", r"\1", name)
        full_name = name if name.startswith("tudatpy.") else module + "." + name
        parts = full_name.split(".")
        for split in range(len(parts), 0, -1):
            try:
                value = importlib.import_module(".".join(parts[:split]))
            except ModuleNotFoundError:
                continue
            for part in parts[split:]:
                value = getattr(value, part)
            break
        else:
            raise ImportError(full_name)
        resolved.append({"rst": relative, "line": number, "name": full_name})

scalar = lc.LaserTargetSignature(0.25)
response = scalar.calculate_response([1e7, 0, 0], [1e7, 0, 0], 0.0)
assert isinstance(response, lc.LaserTargetResponse)
assert abs(response.mean_range_correction + 0.25) < 1e-15
array = lc.LaserTargetSignature([[0, 0, 0], [0.1, 0, 0]], [1, 1], 0.01)
array_response = array.calculate_response([1e7, 0, 0], [1e7, 0, 0], 0.0)
assert np.isfinite(array_response.range_spread)
assert isinstance(
    lc.laser_target_light_time_correction("LAGEOS1", scalar),
    lc.LaserTargetLightTimeCorrectionSettings,
)
assert isinstance(
    lc.mendes_pavlis_tropospheric_light_time_correction(532e-9), lc.LightTimeCorrectionSettings
)
motion = es.ground_station.body_deformation_station_motion()
assert isinstance(motion, es.ground_station.BodyDeformationStationMotionSettings)

sp3_path = Path(
    "/home/dominic/Tudat/laser-ranging-work/results/slr-baseline-combined-slrf2020-20261016-attempt1/uncompressed-inputs/ilrsa.orb.lageos1.240406.v81.sp3"
)
sp3 = ed.read_sp3_file(str(sp3_path))
assert isinstance(sp3, ed.Sp3Data)
records = next(iter(sp3.records.values()))
assert isinstance(records[0], ed.Sp3Record) and isinstance(records[0].epoch, tr.Time)
crd_path = Path("/home/dominic/Tudat/laser-ranging-work/data/slr/2024-04/lageos1_20240406.np2")
passes = slr.read_slr_data([str(crd_path)])
assert passes and isinstance(passes[0], slr.CrdPass)
points = slr.extract_normal_point_measurements(passes[0])
assert points and isinstance(points[0][0], tr.Time)

settings = es.BodyListSettings("SSB", "GCRS")
settings.add_empty_settings("Earth")
settings.get("Earth").rotation_model_settings = es.rotation_model.gcrs_to_itrs()
bodies = es.create_system_of_bodies(settings)
converter = bodies.get("Earth").rotation_model.angles_calculator.time_scale_converter
epoch = tr.Time(765201600)
tdb = converter.convert_time_object(tr.utc_scale, tr.tdb_scale, epoch)
roundtrip = converter.convert_time_object(tr.tdb_scale, tr.utc_scale, tdb)
assert isinstance(tdb, tr.Time) and abs(float(roundtrip - epoch)) < 1e-8

kernel = Path(tudatpy.kernel.__file__)
result = {
    "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "source_revision": revision,
    "documentation_baseline": "ebe3115ec",
    "interpreter": sys.executable,
    "package_origin": tudatpy.__file__,
    "kernel_origin": str(kernel),
    "kernel_sha256": hashlib.sha256(kernel.read_bytes()).hexdigest(),
    "kernel_bytes": kernel.stat().st_size,
    "autosummary_names": resolved,
    "autosummary_count": len(resolved),
    "runtime_smoke": {
        "scalar_target": True,
        "array_target": True,
        "target_factory": True,
        "optical_factory": True,
        "body_deformation_station_motion": True,
        "sp3_reader": {
            "path": str(sp3_path),
            "sha256": hashlib.sha256(sp3_path.read_bytes()).hexdigest(),
            "satellites": len(sp3.records),
        },
        "crd_reader": {
            "path": str(crd_path),
            "sha256": hashlib.sha256(crd_path.read_bytes()).hexdigest(),
            "passes": len(passes),
            "first_pass_points": len(points),
        },
        "selected_earth_rotation_converter": {
            "utc_tdb_utc_error_seconds": float(roundtrip - epoch)
        },
    },
}
(logs / "final-python-surface.json").write_text(json.dumps(result, indent=2) + "\n")
print(
    json.dumps(
        {
            "autosummary_count": len(resolved),
            "runtime_smoke": result["runtime_smoke"],
            "kernel_sha256": result["kernel_sha256"],
        },
        indent=2,
    )
)
