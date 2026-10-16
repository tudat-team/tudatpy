"""Package unchanged INPOP controls and raw descriptive statistics without any observable solves."""

import csv
import hashlib
import json
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-llr")
WORK = Path("/home/dominic/Tudat/laser-ranging-work")
OUT = ROOT / "analysis/laser_ranging/results/inpop17a_long_arc_20170301_20170430"
REV = "e26fa744306173551e98385684457aa477a1dbb4"


def sha(path):
    """Identify exact retained bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stats(rows):
    """Describe every raw residual without centering, weighting or rejecting points."""
    x = np.array([float(row["one_way_residual_m"]) for row in rows])
    return dict(
        count=len(x),
        mean_m=float(x.mean()),
        rms_about_zero_m=float(np.sqrt(np.mean(x * x))),
        scatter_about_mean_m=float(x.std()),
        median_absolute_residual_m=float(np.median(abs(x))),
        percentile95_absolute_residual_m=float(np.percentile(abs(x), 95)),
        minimum_m=float(x.min()),
        maximum_m=float(x.max()),
    )


OUT.mkdir(exist_ok=False)
cases = {}
provenance = dict(
    analysis_source_revision=REV,
    native_source_revision="a366f436c6f6ff03aa57f8b1610075784196bdf8",
    retrieve_source="git show REV:analysis/laser_ranging/PATH",
    source_blobs={},
    numerical_threads=dict(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1"),
    numpy_version=np.__version__,
    cases={},
)
for name in (
    "llr_residuals.py",
    "llr_long_arc.py",
    "inpop17a_dataset.json",
    "llr_long_arc_dataset.json",
    "llr_station_profiles.json",
    "plot_llr_residuals.py",
):
    path = "analysis/laser_ranging/" + name
    blob = subprocess.check_output(
        ["git", "rev-parse", REV + ":" + path], cwd=ROOT, text=True
    ).strip()
    data = subprocess.check_output(["git", "show", REV + ":" + path], cwd=ROOT)
    provenance["source_blobs"][name] = dict(git_blob=blob, sha256=hashlib.sha256(data).hexdigest())

models = {}
for size in ("twelve", "full"):
    for label, convention in [("jpl_nh", "jpl-nh"), ("zero_explicit_glass", "none")]:
        source = WORK / "results" / f"llr-inpop17a-{size}-{convention}-attempt1"
        key = size + "/" + label
        target = OUT / size / label
        target.mkdir(parents=True)
        shutil.copyfile(source / "residuals.csv", target / "residuals.csv")
        rows = list(csv.DictReader((source / "residuals.csv").open()))
        summary = json.loads((source / "summary.json").read_text())
        bystation = {
            station: stats([r for r in rows if r["station"] == station])
            for station in ("APOLLO", "GRASSE")
        }
        byreflector = {
            station
            + "/"
            + ref: stats([r for r in rows if r["station"] == station and r["reflector"] == ref])
            for station, ref in sorted({(r["station"], r["reflector"]) for r in rows})
        }
        bywave = {
            station
            + "/"
            + wave: stats(
                [r for r in rows if r["station"] == station and r["wavelength_nm"] == wave]
            )
            for station, wave in sorted({(r["station"], r["wavelength_nm"]) for r in rows})
        }
        extreme = max(rows, key=lambda r: abs(float(r["one_way_residual_m"])))
        cases[key] = dict(
            count=len(rows),
            station_statistics=bystation,
            station_reflector_statistics=byreflector,
            station_wavelength_statistics=bywave,
            extreme=extreme,
            timing_seconds={
                k: summary[k] for k in ("setup_seconds", "solver_seconds", "total_seconds")
            },
        )
        provenance["cases"][key] = dict(
            external_directory="results/" + source.name,
            external_summary_sha256=sha(source / "summary.json"),
            residual_csv_sha256=sha(source / "residuals.csv"),
            source_sha256=summary["source_sha256"],
        )
        figures = target / "figures"
        figures.mkdir()
        stems = ["residuals_full_span", "residuals_y_detail"]
        if size == "full":
            stems += ["residuals_apollo_time_detail", "residuals_grasse_time_detail"]
        for stem in stems:
            for suffix in (".png", ".pdf"):
                shutil.copyfile(source / "figures" / (stem + suffix), figures / (stem + suffix))
        plotsummary = json.loads((source / "figures/plot_summary.json").read_text())
        plotsummary["inputs"][0]["path"] = "../residuals.csv"
        (figures / "plot_summary.json").write_text(json.dumps(plotsummary, indent=2) + "\n")
        if size == "full":
            models[label] = rows
        else:
            events = source / "event_diagnostics.json"
            labeldata = dict(
                event_diagnostics_sha256=sha(events),
                vectors="Reconstructed at native diagnostic double epochs; not exact solver Time endpoints",
                earth_station_motion_origin="Earth: before the native relativistic_station_motion BCRS spatial map",
                moon_station_motion_origin="Moon: custom lunar spatial map is already included",
                future_export="Source after the executed revision requests SSB and records origin_for_station_motion; ranges are unchanged",
            )
            (source / "event_vector_conventions.json").write_text(
                json.dumps(labeldata, indent=2) + "\n"
            )
            provenance["cases"][key]["saved_events"] = labeldata

summary0 = json.loads((WORK / "results/llr-inpop17a-full-jpl-nh-attempt1/summary.json").read_text())
for key in (
    "manifest_sha256",
    "ephemeris_manifest_sha256",
    "kernel_binary_sha256",
    "kernel_sha256",
    "kernel_pool",
    "earth_orientation_resource_sha256",
):
    provenance[key] = summary0[key]
provenance["command"] = [
    "python",
    "analysis/laser_ranging/llr_long_arc.py",
    "--data-root",
    "$DATA_ROOT",
    "--kernel-dir",
    "$INPOP_KERNEL_DIR",
    "--ephemeris",
    "inpop17a",
    "--glass-convention",
    "$CONVENTION",
    "--output",
    "$NEW_OUTPUT",
]
provenance["twelve_command_addition"] = ["--control-twelve"]
provenance["current_source_note"] = (
    "Only diagnostic SSB export and release-to-fit caveat changed after this executed revision; the source blob and executed ephemeris manifest remain recoverable in Git."
)
provenance["provider_archive"] = dict(
    url="https://ftp.imcce.fr/pub/ephem/planets/inpop17a/inpop17a_TDB_m100_p100_spice.tar.gz",
    sha256="52afaf265f3d1dd4d03b8e2239574af3d3d8c65c0b57e5a8c56fed0b137b98b2",
)
provenance["original_de440_comparison"] = dict(
    csv="../paired_long_arc_20170301_20170430/residuals.csv",
    sha256=sha(
        ROOT / "analysis/laser_ranging/results/paired_long_arc_20170301_20170430/residuals.csv"
    ),
    limit="Original full DE440 run predates scalar-clock-only lineardrift optimization; one complete11-point default night is bit-exact against the pre-INPOP current core.",
)

baseline = list(
    csv.DictReader(
        (
            ROOT / "analysis/laser_ranging/results/paired_long_arc_20170301_20170430/residuals.csv"
        ).open()
    )
)
models["de440_reference"] = baseline
reference_stats = {
    station: stats([r for r in baseline if r["station"] == station])
    for station in ("APOLLO", "GRASSE")
}
nightly = defaultdict(dict)
for name, rows in models.items():
    groups = defaultdict(list)
    for row in rows:
        groups[(row["utc_launch"][:10], row["station"], row["reflector"])].append(row)
    for key, group in groups.items():
        nightly[key][name] = stats(group)["mean_m"]
with (OUT / "nightly_model_comparison.csv").open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(
        [
            "utc_night",
            "station",
            "reflector",
            "de440_reference_mean_m",
            "inpop17a_jpl_nh_mean_m",
            "inpop17a_zero_explicit_glass_mean_m",
        ]
    )
    for key, values in sorted(nightly.items()):
        writer.writerow(
            [*key, values["de440_reference"], values["jpl_nh"], values["zero_explicit_glass"]]
        )
controls = OUT / "controls"
controls.mkdir()
for source, name in [
    (WORK / "results/llr-inpop17a-default-control.json", "default_equivalence.json"),
    (
        WORK / "results/llr-inpop17a-twelve-jpl-nh-attempt1/tide_radius_control.json",
        "tide_radius_control.json",
    ),
]:
    payload = json.loads(source.read_text())
    (controls / name).write_text(json.dumps(payload, indent=2) + "\n")
    provenance[name + "_external_sha256"] = sha(source)
summary = dict(
    cases=cases,
    de440_reference_station_statistics=reference_stats,
    absolute_order_statistics="Raw |O-C| from every station point, no centering/rejection. NumPy percentile95 uses the default linear order-statistic interpolation.",
    optical_interpretation="Both settings were predeclared before the12-point controls. Neither establishes the published prism reference convention.",
    time_interpretation="Nightly means are unweighted descriptions of unchanged residuals; no fitted/detrended correction or period inferred.",
)
(OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
provenance["package_generator_sha256"] = sha(Path(__file__))
provenance["output_sha256"] = {
    str(path.relative_to(OUT)): sha(path) for path in sorted(OUT.rglob("*")) if path.is_file()
}
(OUT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
print(OUT)
