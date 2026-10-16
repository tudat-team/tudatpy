"""Compare pinned normal point releases by reflector and acquisition window."""

import csv
import hashlib
import json
import re
from decimal import Decimal as D
from pathlib import Path

ROOT = Path("/home/dominic/Tudat/laser-ranging-work")
DATA = ROOT / "data/llr"
OUT = ROOT / "verification/apollo-20170309-release-comparison"
OUT.mkdir(exist_ok=True)
WIDTHS = (2, 4, 2, 2, 2, 2, 9, 14, 1, 5, 3, 6, 3, 1, 6, 4, 2, 5, 1, 4, 1)
TARGETS = {0: "apollo11", 1: "luna17", 2: "apollo14", 3: "apollo15", 4: "luna21"}


def sha(path):
    """Return the SHA256 of the exact downloaded bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def civil_sod(fields):
    """Return the integer UTC seconds of day in an H4 endpoint."""
    return D(fields[3]) * 3600 + D(fields[4]) * 60 + D(fields[5])


current = []
files = []
for target in ("apollo11", "apollo14", "apollo15", "luna17"):
    path = DATA / f"{target}_20170309-current.np2"
    files.append(
        {
            "path": path.name,
            "sha256": sha(path),
            "url": f"https://edc.dgfi.tum.de/pub/slr/data/npt_crd_v2/quarantine/7045/2017/{target}_20170309.np2",
        }
    )
    for block in path.read_text().split("h8"):
        lines = block.strip().splitlines()
        if not any(line.startswith("11 ") for line in lines):
            continue
        by_type = {line.split()[0]: line.split() for line in lines}
        h4 = by_type["h4"]
        np = by_type["11"]
        met = by_type["20"]
        text = "\n".join(lines)
        current.append(
            dict(
                target=target,
                start_sod=civil_sod(h4[2:8]),
                end_sod=civil_sod(h4[8:14]),
                tag_sod=D(np[1]),
                rtt_s=D(np[2]),
                photons=int(np[6]),
                bin_rms_ps=D(np[7]),
                sigma_ps=D(re.search(r"normal point uncertainty of ([\d.]+)", text)[1]),
                pressure_hpa=D(met[2]),
                temperature_k=D(met[3]),
                humidity_percent=D(met[4]),
                duration_s=D(np[5]),
                wavelength_nm=D(by_type["c0"][2]),
                epoch_event=np[4],
                record_flags=h4[14:],
                station_h2=" ".join(by_type["h2"][1:]),
                source=path.name,
                acs_mode=(
                    "photon-by-photon"
                    if "Utilizing photon-by-photon" in text
                    else "historic differential estimate"
                ),
                acs_statement=" ".join(
                    line[3:]
                    for line in lines
                    if (
                        "Applied ACS differential" in line
                        or "ACS corrections applied to range photons" in line
                    )
                ),
                headers={key: " ".join(by_type[key]) for key in ("h1", "h2", "h3", "h4", "h5")},
                calibration_record=" ".join(by_type["40"]),
            )
        )
rows = []
used = set()
legacy = DATA / "apollo_20170309_np2.txt"
for i, line in enumerate(legacy.read_text().splitlines(), 1):
    fields = []
    offset = 0
    for width in WIDTHS:
        fields.append(line[offset : offset + width].strip())
        offset += width
    tag = D(fields[4]) * 3600 + D(fields[5]) * 60 + D(fields[6]) * D("1e-7")
    target = TARGETS[int(fields[8])]
    matches = [
        n
        for n, p in enumerate(current)
        if p["target"] == target and p["start_sod"] <= tag <= p["end_sod"]
    ]
    assert len(matches) == 1, (i, matches)
    j = matches[0]
    assert j not in used
    used.add(j)
    point = current[j]
    rtt = D(fields[7]) * D("1e-13")
    delta_rtt = point["rtt_s"] - rtt
    delta_tag = point["tag_sod"] - tag
    rows.append(
        dict(
            point=i,
            target=target,
            legacy_tag_sod=tag,
            current_tag_sod=point["tag_sod"],
            delta_tag_s=delta_tag,
            legacy_rtt_s=rtt,
            current_rtt_s=point["rtt_s"],
            delta_rtt_ps=(delta_rtt * D("1e12") if not delta_tag else ""),
            same_tag_half_range_change_mm=(
                delta_rtt * D("299792458") * 500 if not delta_tag else ""
            ),
            legacy_photons=int(fields[10]),
            current_photons=point["photons"],
            legacy_sigma_rtt_ps=D(fields[11]) / 10,
            current_provider_sigma_rtt_ps=point["sigma_ps"],
            current_bin_rms_ps=point["bin_rms_ps"],
            legacy_pressure_hpa=D(fields[14]) / 100,
            current_pressure_hpa=point["pressure_hpa"],
            legacy_temperature_k=D(fields[15]) / 10 + D("273.15"),
            current_temperature_k=point["temperature_k"],
            legacy_humidity_percent=D(fields[16]),
            current_humidity_percent=point["humidity_percent"],
            legacy_duration_s=D(fields[19]),
            current_duration_s=point["duration_s"],
            legacy_wavelength_nm=D(fields[17]) / 10,
            current_wavelength_nm=point["wavelength_nm"],
            acs_mode=point["acs_mode"],
            acs_statement=point["acs_statement"],
            current_source=point["source"],
        )
    )
with (OUT / "point_differences.csv").open("w") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
report = dict(
    legacy_source={
        "path": legacy.name,
        "sha256": sha(legacy),
        "url": "https://tmurphy.physics.ucsd.edu/apollo/160912e_apollo_np2.txt",
        "parent_sha256": sha(DATA / "160912e_apollo_np2.txt"),
    },
    current_sources=files,
    matching_rule="Same reflector, same UTC date, unique current acquisition window containing original transmit tag. No range values used for matching.",
    current_points=len(current),
    matched_legacy_points=len(rows),
    same_exact_tags=sum(row["delta_tag_s"] == 0 for row in rows),
    unmatched_current_points=[point for n, point in enumerate(current) if n not in used],
    current_header_records=[point["headers"] for point in current],
    current_calibration_records=[point["calibration_record"] for point in current],
    conventions={
        "station_id": "Legacy provider ID70610; current ILRS APOL7045 system95 occupancy1 time scale4.",
        "flags": "Current H4 after endpoints: release1,trop0,COM0,receiveAmplitude0,stationSystemDelay1,spacecraftSystemDelay0,rangeType2,qualityAlert0; record11 event2/transmit. Legacy Group e has qualityB and two undocumented A fields; CRD correction flags cannot be recovered or declared identical from fixed-width rows.",
        "uncertainty": "Current total uncertainty is the explicit provider comment, not CRD binRMS nor binRMS/sqrtN. Provider caveat: excepting a potential scale factor or RSS value, which will follow.",
        "quarantine": "EDC current historical APOLLO records are served under quarantine/7045. WGL021: No XYZ station coordinates found for calculating orbit fit. No rejection based on residuals is inferred.",
        "comparison": "No observation simulation or residuals. Same-tag half-range changes are only c*deltaRTT/2; For unequal tags, no RTT subtraction is reported.",
    },
    retrieved_utc="2026-10-16",
    metadata_urls=[
        "https://earth.gsfc.nasa.gov/geo/projects/sgp/apache-point/updates/complete-dataset-apache-point-lunar-laser-ranging-station",
        "https://edc.dgfi.tum.de/en/data/npt2/885635/",
        "https://edc.dgfi.tum.de/en/oc/crd/v2/#WGL021",
    ],
)
(OUT / "summary.json").write_text(json.dumps(report, indent=2, default=str) + "\n")
for row in rows:
    print(
        row["point"],
        row["target"],
        "tag",
        row["delta_tag_s"],
        "RTTps",
        row["delta_rtt_ps"],
        "sameTagMM",
        row["same_tag_half_range_change_mm"],
        "sigma",
        row["legacy_sigma_rtt_ps"],
        row["current_provider_sigma_rtt_ps"],
        "P",
        row["legacy_pressure_hpa"],
        row["current_pressure_hpa"],
        "T",
        row["legacy_temperature_k"],
        row["current_temperature_k"],
        "N",
        row["legacy_photons"],
        row["current_photons"],
    )
print("Current", len(current), "matched", len(rows), "same tags", report["same_exact_tags"])
print("Extra", report["unmatched_current_points"])
print("Sources", json.dumps(files, indent=2))
