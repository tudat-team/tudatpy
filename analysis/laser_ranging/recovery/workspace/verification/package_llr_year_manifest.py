from pathlib import Path
import json, csv, collections, hashlib, subprocess

r = Path("/home/dominic/Tudat/laser-ranging-work/research/llr-year-2017")
repo = Path("/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-data/analysis/laser_ranging")
m = json.loads(
    subprocess.check_output(
        ["git", "show", "696856971:analysis/laser_ranging/llr_long_arc_dataset.json"], cwd=repo
    )
)
a = json.loads((r / "apollo_subset_sources.json").read_text())
v = json.loads((r / "native-reader-validation.json").read_text())
sessions = json.loads((r / "crd_sessions.json").read_text())
rows = list(csv.DictReader((r / "grasse_exact_inventory.csv").open()))
g = [q for q in rows if q["station_cdp"] == "7845"]
gs = [q for q in sessions if q["station"][1] == "7845"]
m.update(
    interval_start_utc="2017-01-01",
    interval_stop_exclusive_utc="2018-01-01",
    selection="All calendar-2017 APOLLO Group-e rows and all Grasse 7845 normal points in the 60 official monthly CRDv1 products; only the explicitly identified user-requested Grasse NP is excluded. Other stations in these products are outside the predeclared station scope.",
)
m["apollo"].update(a)
m["apollo"].pop("count", None)
m["apollo"].update(
    source_count=311,
    retained_count=311,
    nights=v["apollo"]["nights"],
    reflector_counts=v["apollo"]["reflector_counts"],
    source_metadata=v["apollo"]["source_metadata"],
    first_utc_launch=v["apollo"]["first_utc"],
    last_utc_launch=v["apollo"]["last_utc"],
)
gm = m["grasse"]
gm.update(
    files=json.loads((r / "SOURCES.json").read_text()),
    source_count=len(g),
    retained_count=len(g) - 1,
    sessions=len(gs),
    nights=len({q["utc_launch"][:10] for q in g}),
    reflector_counts=dict(collections.Counter(q["target"] for q in g)),
    first_utc_launch=min(q["utc_launch"] for q in g),
    last_utc_launch=max(q["utc_launch"] for q in g),
    wavelength_counts=dict(collections.Counter(q["wavelength_nm"] for q in g)),
    maximum_nearest_weather_seconds=max(float(q["nearest_weather_s"]) for q in g),
)
gm.pop("count", None)
gm["inventory"] = {
    "file": "grasse_exact_inventory.csv",
    "sha256": hashlib.sha256((r / "grasse_exact_inventory.csv").read_bytes()).hexdigest(),
}
gm["weather_coverage_audit"] = {
    "status": "Three retained NPs exceed the prior 1000 s session endpoint-hold guard; source audit only, modeling policy pending root review",
    "normal_points_exceeding_prior_hold_guard": [
        {
            k: q[k]
            for k in [
                "source_file",
                "source_line",
                "utc_launch",
                "nearest_weather_s",
                "weather_endpoint_hold_required_s",
            ]
        }
        for q in g
        if float(q["weather_endpoint_hold_required_s"]) > 1000
    ],
    "june22_final_measured_weather": {
        "utc": "2017-06-22T09:29:16",
        "pressure_hPa": 880.13,
        "temperature_K": 295.95,
        "relative_humidity_percent": 65,
        "origin": 0,
    },
    "daily_source_check": "Official June22 daily A15/L1/L2 files repeat the monthly meteorology; no later same-station weather across any selected target. No measurements are invented and no additional NP is excluded.",
}
m.pop("overlap_preflight", None)
m["source_inventory_audit"] = {
    "monthly_products": 60,
    "provider_missing_products": 0,
    "byte_identical_frozen_march_april_products": 10,
    "all_crd_normal_points": sum(v["station_np_counts"].values()),
    "out_of_scope_station_np_counts": {"7941": 33},
    "out_of_scope_selection": "Matera 7941 is outside the predeclared APOLLO/Grasse station scope, irrespective of metadata quality.",
    "native_reader_validation": v,
    "preparation_script": "prepare_llr_year_2017.py",
}
(repo / "llr_year_2017_dataset.json").write_text(json.dumps(m, indent=2) + "\n")
