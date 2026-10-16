"""Compare raw station nightly means at nearby actual epochs without time alignment."""

import csv
from collections import defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

source = Path(sys.argv[1])

output = Path(sys.argv[2])

rows = list(csv.DictReader(source.open(newline="")))

groups = defaultdict(list)

for row in rows:
    groups[(row["station"], row["reflector"], row["utc_launch"][:10])].append(row)
nightly = []

for (station, reflector, night), points in sorted(groups.items()):
    epochs = [
        datetime.fromisoformat(r["utc_launch"]).replace(tzinfo=timezone.utc).timestamp()
        for r in points
    ]

    nightly.append(
        dict(
            station=station,
            reflector=reflector,
            utc_night=night,
            mean_epoch_s=sum(epochs) / len(epochs),
            count=len(points),
            raw_mean_m=sum(float(r["one_way_residual_m"]) for r in points) / len(points),
        )
    )
matched = []

used = set()

# Choose the nearest actual mean epoch, allowing either adjacent UTC date, within twelve hours.
for apollo in [r for r in nightly if r["station"] == "APOLLO"]:
    candidates = [
        r for r in nightly if r["station"] == "GRASSE" and r["reflector"] == apollo["reflector"]
    ]

    grasse = (
        min(candidates, key=lambda r: abs(r["mean_epoch_s"] - apollo["mean_epoch_s"]))
        if candidates
        else None
    )

    if grasse is None or abs(grasse["mean_epoch_s"] - apollo["mean_epoch_s"]) > 12 * 3600:
        continue
    used.add((grasse["reflector"], grasse["utc_night"]))
    matched.append(
        dict(
            reflector=apollo["reflector"],
            apollo_night=apollo["utc_night"],
            grasse_night=grasse["utc_night"],
            apollo_mean_utc=datetime.fromtimestamp(
                apollo["mean_epoch_s"], timezone.utc
            ).isoformat(),
            grasse_mean_utc=datetime.fromtimestamp(
                grasse["mean_epoch_s"], timezone.utc
            ).isoformat(),
            grasse_minus_apollo_hours=(grasse["mean_epoch_s"] - apollo["mean_epoch_s"]) / 3600,
            apollo_count=apollo["count"],
            grasse_count=grasse["count"],
            apollo_raw_mean_m=apollo["raw_mean_m"],
            grasse_raw_mean_m=grasse["raw_mean_m"],
            grasse_minus_apollo_raw_mean_m=grasse["raw_mean_m"] - apollo["raw_mean_m"],
        )
    )
output.mkdir(parents=True, exist_ok=True)
with (output / "nearby_nightly_raw_means.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(matched[0]))

    writer.writeheader()
    writer.writerows(matched)
matched_apollo = {(r["reflector"], r["apollo_night"]) for r in matched}

unmatched = []

for row in nightly:
    paired = matched_apollo if row["station"] == "APOLLO" else used

    if (row["reflector"], row["utc_night"]) in paired:
        continue
    other = [
        r for r in nightly if r["station"] != row["station"] and r["reflector"] == row["reflector"]
    ]

    nearest = (
        min(other, key=lambda r: abs(r["mean_epoch_s"] - row["mean_epoch_s"])) if other else None
    )

    unmatched.append(
        dict(
            station=row["station"],
            reflector=row["reflector"],
            utc_night=row["utc_night"],
            count=row["count"],
            nearest_other_station_night=nearest["utc_night"] if nearest else None,
            nearest_other_minus_this_hours=(
                (nearest["mean_epoch_s"] - row["mean_epoch_s"]) / 3600 if nearest else None
            ),
        )
    )

summary = dict(
    pairing="Nearest Grasse station/reflector nightly mean observation epoch to each APOLLO one, within 12 hours, allowing adjacent UTC dates; no interpolation or fitted alignment",
    matched_pairs=len(matched),
    adjacent_date_pairs=sum(r["apollo_night"] != r["grasse_night"] for r in matched),
    unmatched_station_reflector_nights=unmatched,
    comparison_scope="Raw means at distinct times; differences contain changing lunar/station geometry and are not simultaneous station biases",
)
(output / "nearby_nightly_comparison.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps({k: summary[k] for k in ["matched_pairs", "adjacent_date_pairs"]}))
