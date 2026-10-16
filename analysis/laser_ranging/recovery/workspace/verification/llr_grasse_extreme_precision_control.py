"""Control three unchanged Grasse NPs in their original full-session environment."""

import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from tudatpy import constants, kernel
from tudatpy.astro import time_representation as time

sys.path.insert(0, "/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-data/analysis/laser_ranging")
import llr_long_arc as arc
import llr_residuals as model

work = Path("/home/dominic/Tudat/laser-ranging-work")

manifest = json.loads(arc.MANIFEST.read_text())

ephemeris = json.loads((Path(model.__file__).parent / "llr_dataset.json").read_text())

model.load_ephemeris_kernels(
    work / "data/llr/kernels", ephemeris, ephemeris["analysis_frame_alias"]
)
points, weather, profile, ledger = next(
    group
    for group in arc.read_groups(work / "research/longarc-20170304", manifest)
    if group[3]["group"] == "apollo15_201703.npt:29"
)

selected = [point for point in points if point["source_line"] in (168, 171, 174)]

original = {
    int(row["source_line"]): row
    for row in csv.DictReader(
        (work / "results/llr-longarc-20170301-20170430-attempt2/residuals.csv").open()
    )
    if row["source_file"] == "apollo15_201703.npt"
}

rows = []

for step in (30, 15):
    bodies, position, _ = model.create_environment(
        points, spice_knot_step=step, station_profile=profile, meteorological_history=weather
    )

    converter = bodies.get("Earth").rotation_model.angles_calculator.time_scale_converter

    for point in selected:
        epoch = converter.convert_time_object(
            time.utc_scale, time.tdb_scale, point["utc"], position
        )

        computed, diagnostic = model.computed_range(
            point, epoch, bodies, model.complete_light_time_corrections(point)
        )

        residual = 0.5 * (constants.SPEED_OF_LIGHT * point["observed_rtt_s"] - computed)

        rows.append(
            dict(
                source_line=point["source_line"],
                utc_launch=point["utc_label"],
                spice_knot_step_s=step,
                observed_rtt_s=point["observed_rtt_s"],
                computed_two_way_range_m=computed,
                one_way_residual_m=residual,
                half_range_delta_from_frozen30_m=0.5
                * (
                    computed
                    - float(original[point["source_line"]]["computed_full_two_way_range_m"])
                ),
                native_event_epochs_tdb_s=diagnostic["link_end_epochs_tdb_s"],
            )
        )
        if step == 30:
            assert computed == float(
                original[point["source_line"]]["computed_full_two_way_range_m"]
            )
output = work / "results/llr-longarc-20170301-20170430-attempt2/extreme_precision_control.json"

result = dict(
    analysis_source_revision="6e7d4b807",
    native_source_revision="a366f436c6f6ff03aa57f8b1610075784196bdf8",
    kernel_binary_sha256=model.sha256(Path(kernel.__file__)),
    control_script_sha256=model.sha256(Path(__file__)),
    source_sha256={
        path.name: model.sha256(path)
        for path in [Path(model.__file__), Path(arc.__file__), model.STATION_PROFILE_PATH]
    },
    manifest_sha256=model.sha256(arc.MANIFEST),
    control="Full original fourteen-NP session and fourteen measured weather records determine both environments; solve only the selected three NPs, unchanged50 iterations/fractional1e-13 tolerance, no loading/fits/cuts.",
    profile=profile,
    source_session=ledger["group"],
    source_file_sha256=model.sha256(work / "research/longarc-20170304/crd/apollo15_201703.npt"),
    full_session_np_count=len(points),
    measured_weather_count=len(weather),
    rows=rows,
)

output.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(rows, indent=2))
