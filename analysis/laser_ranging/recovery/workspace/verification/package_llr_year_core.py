"""Package completed annual raw results without recomputing observations."""

from pathlib import Path
import csv
import hashlib
import json
import math
import shutil
import statistics
import subprocess

WORK = Path("/home/dominic/Tudat/laser-ranging-work")
SOURCE = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-signature")
OUTPUT = SOURCE / "analysis/laser_ranging/results/year_2017"
INTEGRATION = "5e4934584"
CASES = ["de440", "inpop17a_jpl_nh", "inpop17a_no_separate_glass"]
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
OUTPUT.mkdir(exist_ok=False)
provenance = dict(
    executed_source_revision="9082ce261",
    integration_source_revision=INTEGRATION,
    native_source_revision="a366f436c6f6ff03aa57f8b1610075784196bdf8",
    kernel_binary_sha256="cab49e738cfcf0ea5f251cb140ab8919d5372dbc3e5ba09314d0d5be4eedfcef",
    numerical_threads=1,
    launch_mode="Three independent concurrent native interpreters; no shared mutable kernel pool",
    source_blobs={},
    cases={},
)
for name in [
    "llr_long_arc.py",
    "llr_residuals.py",
    "llr_year_2017_dataset.json",
    "llr_dataset.json",
    "inpop17a_dataset.json",
    "llr_station_profiles.json",
    "icrf_spice_alias.tf",
    "plot_llr_residuals.py",
]:
    relative = "analysis/laser_ranging/" + name
    blob = subprocess.check_output(["git", "show", f"{INTEGRATION}:{relative}"], cwd=SOURCE)
    assert hashlib.sha256(blob).hexdigest() == sha(SOURCE / relative)
    provenance["source_blobs"][name] = dict(
        sha256=sha(SOURCE / relative),
        git_blob=subprocess.check_output(
            ["git", "rev-parse", f"{INTEGRATION}:{relative}"], cwd=SOURCE, text=True
        ).strip(),
    )
summary = {}
table = []
for case in CASES:
    src = WORK / "results/llr-year-2017-attempt2/full" / case
    data = json.loads((src / "summary.json").read_text())
    assert data["count"] == 1579 and data["source_count"] == 1582 and data["excluded_count"] == 3
    assert (
        data["station_statistics"]["APOLLO"]["count"] == 311
        and data["station_statistics"]["GRASSE"]["count"] == 1268
    )
    rows = list(csv.DictReader((src / "residuals.csv").open(newline="")))
    assert len(rows) == 1579
    target = OUTPUT / case
    target.mkdir()
    shutil.copy2(src / "residuals.csv", target / "residuals.csv")
    summary[case] = {key: value for key, value in data.items() if key != "groups"}
    record = json.loads(
        (WORK / "logs" / f"llr-year-2017-attempt2-full-{case}-provenance.json").read_text()
    )
    record["integration_source_revision"] = INTEGRATION
    record["root_equivalence"] = (
        "All listed model, driver, manifest, station and plotting source blobs verified against root integration revision after execution."
    )
    record["external_summary"] = str(src / "summary.json")
    record["external_summary_sha256"] = sha(src / "summary.json")
    record["partial_csv_sha256"] = sha(src / "residuals.partial.csv")
    boundary = WORK / "results/llr-year-2017-attempt1/boundary" / case / "residuals.csv"
    controls = list(csv.DictReader(boundary.open(newline="")))
    lookup = {(r["station"], r["source_file"], r["source_line"], r["utc_launch"]): r for r in rows}
    for row in controls:
        assert (
            lookup[(row["station"], row["source_file"], row["source_line"], row["utc_launch"])]
            == row
        )
    record["first_last_control"] = dict(
        csv_sha256=sha(boundary),
        count=4,
        all_fields_exact_against_full=True,
        columns=len(rows[0]),
        source_note="Before the two source-invalid February records were identified; these unchanged four boundary observations are reproduced by the final eligible arc.",
    )
    provenance["cases"][case] = record
    for station in ["APOLLO", "GRASSE"]:
        values = [1000 * float(r["one_way_residual_m"]) for r in rows if r["station"] == station]
        n = len(values)
        mean = math.fsum(values) / n
        absolute = sorted(map(abs, values))
        q = 0.95 * (n - 1)
        left = int(q)
        right = math.ceil(q)
        table.append(
            dict(
                case=case,
                station=station,
                count=n,
                mean_mm=mean,
                rms_mm=math.sqrt(math.fsum(v * v for v in values) / n),
                scatter_population_mm=math.sqrt(math.fsum((v - mean) ** 2 for v in values) / n),
                median_absolute_mm=statistics.median(absolute),
                p95_absolute_linear_mm=absolute[left]
                + (q - left) * (absolute[right] - absolute[left]),
            )
        )
(OUTPUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
(OUTPUT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
with (OUTPUT / "raw_statistics.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(table[0]))
    writer.writeheader()
    writer.writerows(table)
print(json.dumps(table, indent=2))
