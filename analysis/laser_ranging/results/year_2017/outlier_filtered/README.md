This is user-authorized postprocessing of the saved calendar-2017 residuals.
Exactly two common observation identities are omitted from this derived table
in every model. All raw provider files and residual CSVs remain unchanged.
No model runs, fitted offsets, centering or additional cuts are performed.
The earlier user-requested March 9 Grasse exclusion is already absent from
these input CSVs and is not counted again.

| Source identity (Grasse, L1) | Printed UTC | DE440 residual (m) | INPOP17a + JPL glass (m) | INPOP17a without separate glass (m) |
| --- | --- | ---: | ---: | ---: |
| `luna17_201706.npt:58` | 2017-06-13T03:12:09.609942659400 | 0.5395044088363647 | 0.37395578622817993 | 0.4372681975364685 |
| `luna17_201706.npt:154` | 2017-06-22T09:47:06.763983741400 | 110.45823282003403 | 110.2642195224762 | 110.3275335431099 |

The same exact identity mask is applied to all three models. It removes no
APOLLO points and does not apply a residual threshold to the model-dependent
DE440 December pattern. Each derived case contains 311 APOLLO and 1,266 Grasse
points (1,577 total), compared with 1,579 input rows.

| Model | Station | N | Mean (mm) | RMS (mm) | P95 absolute residual (mm) |
| --- | --- | ---: | ---: | ---: | ---: |
| de440 | APOLLO | 311 | 15.626 | 72.395 | 146.130 |
| de440 | GRASSE | 1266 | -10.382 | 75.992 | 150.192 |
| inpop17a_jpl_nh | APOLLO | 311 | -51.842 | 56.255 | 89.357 |
| inpop17a_jpl_nh | GRASSE | 1266 | -62.619 | 67.296 | 107.297 |
| inpop17a_no_separate_glass | APOLLO | 311 | -8.220 | 22.986 | 46.205 |
| inpop17a_no_separate_glass | GRASSE | 1266 | -13.543 | 26.286 | 50.078 |

[statistics.csv](statistics.csv) retains the numerical values. For residuals
r in metres, mean = 1000 × math.fsum(r)/N and RMS =
1000 × sqrt(math.fsum(r²)/N), about zero. P95 is the linear percentile of
absolute residuals: sort |r|, use index 0.95 × (N − 1), and linearly interpolate
between adjacent entries; multiply by 1000 for millimetres. No mean is removed.

The unchanged source CSVs are `../<model>/residuals.csv`, pinned as follows:

| Model | Raw CSV SHA256 |
| --- | --- |
| de440 | `d169c8cda2f04f5a1eef3aa59be26c7b8ae448d2d16d68d2130eeea2186c2c48` |
| inpop17a_jpl_nh | `378bc8d63f111bd387395babd0796b40b3d1a26fee2147e0774f43d31adcca95` |
| inpop17a_no_separate_glass | `7b9701cb97ea980d5cdc467aebbf231b2fc8dc0a5bbe2a3b18834edcde34ad9c` |

Run this stdlib-only snippet from this directory to reproduce all six rows
from the adjacent raw result directories:

```sh
python - <<'PYTHON'
import csv, math
from pathlib import Path

mask = {
    ("GRASSE", "luna17_201706.npt", "58", "2017-06-13T03:12:09.609942659400"),
    ("GRASSE", "luna17_201706.npt", "154", "2017-06-22T09:47:06.763983741400"),
}
for model in ("de440", "inpop17a_jpl_nh", "inpop17a_no_separate_glass"):
    with (Path("..") / model / "residuals.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    kept = [r for r in rows if tuple(r[k] for k in
            ("station", "source_file", "source_line", "utc_launch")) not in mask]
    assert len(rows) == 1579 and len(kept) == 1577
    for station in ("APOLLO", "GRASSE"):
        r = [float(x["one_way_residual_m"]) for x in kept if x["station"] == station]
        a = sorted(map(abs, r))
        index = .95 * (len(a) - 1)
        lo, hi = math.floor(index), math.ceil(index)
        p95 = a[lo] + (index - lo) * (a[hi] - a[lo])
        print(model, station, len(r), 1000 * math.fsum(r) / len(r),
              1000 * math.sqrt(math.fsum(x*x for x in r) / len(r)), 1000 * p95)
PYTHON
```
