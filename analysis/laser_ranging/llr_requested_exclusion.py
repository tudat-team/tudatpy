"""Create current305-point packages from saved ranges, excluding only the requested source row."""

import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from llr_long_arc import MANIFEST, is_requested_exclusion, statistics, write_csv
import plot_llr_residuals as plotting


def raw_statistics(rows):
    """Reuse raw range statistics and retain absolute order statistics without centering."""

    values = np.array([float(row["one_way_residual_m"]) for row in rows])
    return dict(
        **statistics([dict(one_way_residual_m=value) for value in values]),
        median_absolute_residual_m=float(np.median(abs(values))),
        percentile95_absolute_residual_m=float(np.percentile(abs(values), 95)),
    )


def main():
    """Filter the three immutable306-point CSVs byte-for-byte, without computing observations."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    manifest = json.loads(MANIFEST.read_text())
    root = Path(__file__).parent / "results"
    inpop = root / "inpop17a_long_arc_20170301_20170430/full"
    cases = {
        "de440": root / "paired_long_arc_20170301_20170430/residuals.csv",
        "inpop17a_jpl_nh": inpop / "jpl_nh/residuals.csv",
        "inpop17a_zero_explicit_glass": inpop / "zero_explicit_glass/residuals.csv",
    }
    args.output.mkdir(parents=True, exist_ok=False)
    summary = dict(
        excluded_normal_point=manifest["excluded_normal_point"],
        manifest_sha256=hashlib.sha256(MANIFEST.read_bytes()).hexdigest(),
        postprocessing_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        method="One exact source-identity exclusion; no solves, residual threshold, centering or fitting",
        absolute_order_statistics="All retained raw values; NumPy default linear95th percentile",
        cases={},
    )

    for label, source in cases.items():
        raw = source.read_bytes()
        source_hash = hashlib.sha256(raw).hexdigest()
        if label == "de440":
            expected_hash = json.loads((source.parent / "provenance.json").read_text())[
                "output_sha256"
            ]["residuals.csv"]
        else:
            expected_hash = json.loads((inpop.parent / "provenance.json").read_text())["cases"][
                "full/" + label[len("inpop17a_") :]
            ]["residual_csv_sha256"]
        if source_hash != expected_hash:
            raise ValueError("The original computed CSV differs from its pinned provenance")
        lines = raw.splitlines(keepends=True)
        rows = list(csv.DictReader(raw.decode().splitlines()))
        if len(rows) != 306 or len(lines) != len(rows) + 1:
            raise ValueError("Expected the original306 single-line CSV records")
        source_counts = {
            station: sum(row["station"] == station for row in rows)
            for station in ("APOLLO", "GRASSE")
        }
        if source_counts != {"APOLLO": 69, "GRASSE": 237}:
            raise ValueError("Original station source coverage changed")
        selected = [
            is_requested_exclusion(
                row["station"], row["source_file"], row["source_line"], row["utc_launch"], manifest
            )
            for row in rows
        ]
        if sum(selected) != 1 or rows[selected.index(True)]["reflector"] != "A15":
            raise ValueError("The requested source identity must identify exactly one A15 point")
        retained = [row for row, excluded in zip(rows, selected) if not excluded]
        if len(retained) != 305 or sum(row["station"] == "GRASSE" for row in retained) != 236:
            raise ValueError("Expected69 APOLLO and236 retained Grasse points")
        target = args.output / label
        target.mkdir()
        output = lines[0] + b"".join(
            line for line, excluded in zip(lines[1:], selected) if not excluded
        )
        (target / "residuals.csv").write_bytes(output)
        nights = defaultdict(list)

        for row in retained:
            nights[(row["utc_launch"][:10], row["station"], row["reflector"])].append(row)
        write_csv(
            target / "nightly_raw_statistics.csv",
            [
                dict(
                    utc_night=night, station=station, reflector=reflector, **raw_statistics(values)
                )
                for (night, station, reflector), values in sorted(nights.items())
            ],
        )
        figures = target / "figures"
        figures.mkdir()
        points = plotting.read_points(target / "residuals.csv")
        title = label + ";305 retained points; one explicitly requested source exclusion"
        plotting.plot(points, figures / "residuals_full_span", title, False)
        plotting.plot(
            points, figures / "residuals_y_detail", title + "; ±200mm detail", False, (-200, 200)
        )
        summary["cases"][label] = dict(
            source_csv=str(source.relative_to(root)),
            source_csv_sha256=source_hash,
            retained_csv_sha256=hashlib.sha256(output).hexdigest(),
            source_count=306,
            retained_count=len(retained),
            source_station_counts=source_counts,
            all_retained_csv_lines_byte_identical=True,
            station_statistics={
                station: raw_statistics([r for r in retained if r["station"] == station])
                for station in ("APOLLO", "GRASSE")
            },
            station_reflector_statistics={
                station + "/" + reflector: raw_statistics(values)
                for (station, reflector), values in (
                    (key, [r for r in retained if (r["station"], r["reflector"]) == key])
                    for key in sorted({(r["station"], r["reflector"]) for r in retained})
                )
            },
            descriptive_plot_statistics=plotting.describe(points),
        )
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        json.dumps(
            {key: value["station_statistics"] for key, value in summary["cases"].items()}, indent=2
        )
    )


if __name__ == "__main__":
    main()
