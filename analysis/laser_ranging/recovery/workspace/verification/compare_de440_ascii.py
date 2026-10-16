"""Compare original DE440 ASCII coefficients with the pinned SPK and lunar PCK."""

from pathlib import Path
import hashlib
import json
import numpy as np
import mpmath as mp
from jplephem.spk import SPK
from jplephem.pck import PCK

mp.mp.dps = 80
root = Path("/home/dominic/Tudat/laser-ranging-work")
source = root / "research/de440-coefficients/ascp01950.440"
header_path = root / "data/llr/header.440"
header = header_path.read_text()
names = header.split("GROUP   1040")[1].split("GROUP")[0].split()[1:]
values = header.split("GROUP   1041")[1].split("GROUP")[0].split()[1:]
constants = dict(zip(names, [mp.mpf(v.replace("D", "e")) for v in values]))
rows = header.split("GROUP   1050")[1].split("GROUP")[0].splitlines()
rows = [list(map(int, r.split())) for r in rows if r.strip()]
indices, counts, intervals = rows
records = {}
selected = {}
for case, directory in [
    ("apollo", "apollo-20170309-survey-stable-cubic30"),
    ("matera", "matera-20220821-stable-cubic30"),
]:
    events = json.loads((root / "results" / directory / "event_diagnostics.json").read_text())
    selected[case] = [
        (
            e["index"],
            mp.mpf(e["launch_tdb_time_split"][0]) * 3600 + mp.mpf(e["launch_tdb_time_split"][1]),
        )
        for e in events
    ]
required_days = [
    mp.mpf("2451545") + epoch / 86400 for points in selected.values() for _, epoch in points
]
with source.open() as stream:
    while line := stream.readline():
        if not line.strip():
            continue
        number, length = map(int, line.split())
        text = [stream.readline() for _ in range((length + 2) // 3)]
        fields = "".join(text).split()[:length]
        start, stop = [mp.mpf(s.replace("D", "e")) for s in fields[:2]]
        if any(start <= day < stop for day in required_days):
            records[number] = [mp.mpf(s.replace("D", "e")) for s in fields]

spk = SPK.open(str(root / "data/llr/kernels/de440.bsp"))
pck = PCK.open(str(root / "data/llr/kernels/moon_pa_de440_200625.bpc"))


def clenshaw(coefficients, x):
    """Evaluate an ascending Chebyshev coefficient list without rounding intermediate values."""
    previous = previous_previous = mp.mpf(0)
    for coefficient in coefficients[:0:-1]:
        current = coefficient + 2 * x * previous - previous_previous
        previous_previous, previous = previous, current
    return coefficients[0] + x * previous - previous_previous


def ascii_value(quantity, epoch):
    """Select the actual JPL subinterval using the original header's pointer table."""
    day = mp.mpf("2451545") + epoch / 86400
    record = next(row for row in records.values() if row[0] <= day < row[1])
    duration = (record[1] - record[0]) / intervals[quantity]
    subinterval = int(mp.floor((day - record[0]) / duration))
    x = 2 * (day - record[0] - subinterval * duration) / duration - 1
    first = indices[quantity] - 1 + 3 * counts[quantity] * subinterval
    return mp.matrix(
        [
            clenshaw(
                record[first + axis * counts[quantity] : first + (axis + 1) * counts[quantity]], x
            )
            for axis in range(3)
        ]
    )


def binary_value(segment, epoch):
    """Evaluate exact binary coefficients with the same independent high-precision arithmetic."""
    initial, duration, coefficients = (
        segment._data if hasattr(segment, "target") else segment._load()
    )
    initial, duration = mp.mpf(float(initial)), mp.mpf(float(duration))
    subinterval = int(mp.floor((epoch - initial) / duration))
    x = 2 * (epoch - initial - subinterval * duration) / duration - 1
    return mp.matrix(
        [
            clenshaw([mp.mpf(float(v)) for v in coefficients[::-1, axis, subinterval]], x)
            for axis in range(3)
        ]
    )


def norm(vector):
    """Keep precision while differencing positions or orientation angles."""
    return float(mp.sqrt(sum(v * v for v in vector)))


report = {
    "ascii_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    "header_sha256": hashlib.sha256(header_path.read_bytes()).hexdigest(),
    "DENUM": str(constants["DENUM"]),
    "LENUM": str(constants["LENUM"]),
    "TDATEF": str(constants["TDATEF"]),
    "TDATEB": str(constants["TDATEB"]),
    "arithmetic_decimal_digits": 80,
    "selected_record_numbers": list(records),
    "cases": {},
}
for case, points in selected.items():
    result = []
    for index, epoch in points:
        emb = ascii_value(2, epoch)
        moon = ascii_value(9, epoch)
        earth_from_emb = -moon / (1 + constants["EMRAT"])
        moon_from_emb = moon + earth_from_emb
        differences = {
            "emb_m": norm(1000 * (binary_value(spk[0, 3], epoch) - emb)),
            "earth_from_emb_m": norm(1000 * (binary_value(spk[3, 399], epoch) - earth_from_emb)),
            "moon_from_emb_m": norm(1000 * (binary_value(spk[3, 301], epoch) - moon_from_emb)),
            "sun_m": norm(1000 * (binary_value(spk[0, 10], epoch) - ascii_value(10, epoch))),
        }
        angles = binary_value(pck.segments[0], epoch) - ascii_value(12, epoch)
        differences["libration_angle_differences_rad"] = [float(v) for v in angles]
        differences["libration_angle_norm_rad"] = norm(angles)
        result.append({"index": index, **differences})
    report["cases"][case] = result
(root / "logs/llr-de440-ascii-kernel-comparison.json").write_text(
    json.dumps(report, indent=2) + "\n"
)
print(json.dumps({k: v for k, v in report.items() if k != "cases"}, indent=2))
for case, rows in report["cases"].items():
    print(
        case,
        {
            k: max(abs(row[k]) for row in rows)
            for k in [
                "emb_m",
                "earth_from_emb_m",
                "moon_from_emb_m",
                "sun_m",
                "libration_angle_norm_rad",
            ]
        },
        "firstangles",
        rows[0]["libration_angle_differences_rad"],
    )
