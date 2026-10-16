# Calendar 2017 fixed-ephemeris LLR diagnostic

Each case retains **1,579 observations: 311 APOLLO and 1,268 Grasse**, across
all five reflectors. The frozen [manifest](../../llr_year_2017_dataset.json)
contains 1,582 source rows. Exactly one March 9 Grasse point is excluded at
Dominic's request; two February 20 records are separately excluded because
their printed 17/25 ns round-trip times cannot represent lunar observations.
Their exact identities and reasons remain in [summary.json](summary.json).
No other cuts, fitted parameters, offset removal or trend removal are applied.
Both copies of the provider's repeated August 28 source record remain.

| Case | Raw CSV | Full-span PNG | ±200 mm detail PNG |
| --- | --- | --- | --- |
| DE440 | [CSV](de440/residuals.csv) | [Full](de440/figures/residuals_full_span.png) | [Detail](de440/figures/residuals_y_detail.png) |
| INPOP17a, explicit JPL glass term | [CSV](inpop17a_jpl_nh/residuals.csv) | [Full](inpop17a_jpl_nh/figures/residuals_full_span.png) | [Detail](inpop17a_jpl_nh/figures/residuals_y_detail.png) |
| INPOP17a, no separate glass term | [CSV](inpop17a_no_separate_glass/residuals.csv) | [Full](inpop17a_no_separate_glass/figures/residuals_full_span.png) | [Detail](inpop17a_no_separate_glass/figures/residuals_y_detail.png) |

Matching PDFs are beside each PNG. Points are unconnected; Grasse green/IR
wavelengths have different markers. Detail panels label points outside their
display limits; full-span figures and all statistics retain them.

Raw one-way O−C statistics, **mm**, without centering:

| Case | Station | N | Mean | RMS | Median \|r\| | P95 \|r\| |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| DE440 | APOLLO | 311 | +15.626 | 72.395 | 49.993 | 146.130 |
| DE440 | Grasse | 1268 | +77.172 | 3102.942 | 53.456 | 151.711 |
| INPOP17a, JPL glass | APOLLO | 311 | −51.842 | 56.255 | 50.425 | 89.357 |
| INPOP17a, JPL glass | Grasse | 1268 | +24.734 | 3097.276 | 59.140 | 107.421 |
| INPOP17a, no separate glass | APOLLO | 311 | −8.220 | 22.986 | 14.288 | 46.205 |
| INPOP17a, no separate glass | Grasse | 1268 | +73.833 | 3098.441 | 17.401 | 50.259 |

[raw_statistics.csv](raw_statistics.csv) also records population scatter and
full precision; P95 uses linear interpolation between sorted absolute values.
Grasse RMS is dominated by the **retained** June 22 L1 point
`luna17_201706.npt:154` (+110.328 m in the no-separate-glass case). Its cause
is unestablished. It is one of three explicitly marked long endpoint-weather
holds, retained as source-data limitations; those holds do not establish
millimetre atmospheric accuracy. These residuals are not intrinsic noise
estimates. The INPOP optical reference convention remains unresolved, so both
glass cases are sensitivities, not a choice justified by their residual size.
See the [source audit](../../llr_year_2017_sources/README.md) and
[temporal assessment](../year_2017_temporal/README.md) for dated patterns,
coverage gaps and station pairing.

All 915 overlapping March–April computed ranges are bit-identical to the
three previous 305-point cases; independent source/count/finite-output checks
also pass ([verification](../../llr_year_2017_sources/annual_output_verification.md)).
The four first/last control rows match every final CSV field in each case.
[provenance.json](provenance.json) records commands, source blobs, hashes and
timings: source `9082ce261` is blob-equivalent to integration `5e4934584`,
with kernel `cab49e73…` (native source `a366f436c`). Three independent
single-thread interpreters completed together in 236.6 s; per-case setup,
solver and elapsed times are retained. No native rebuild was needed.

Replay with the pinned provider cache and matching native kernel environment:

```sh
python analysis/laser_ranging/llr_long_arc.py \
  --dataset-manifest analysis/laser_ranging/llr_year_2017_dataset.json \
  --data-root /path/to/llr-year-2017 \
  --kernel-dir /path/to/coherent-kernels \
  --ephemeris inpop17a --glass-convention none --output /new/output
python analysis/laser_ranging/plot_llr_residuals.py /new/output/residuals.csv \
  --output /new/figures --title '2017 — INPOP17a — no separate glass term' \
  --zoom-mm -200 200
```

Use `inpop17a/jpl-nh` or `de440/jpl-nh` for the other cases, with their
corresponding coherent kernel directory. Set numerical thread counts to one;
the exact executed environment and local build path are in provenance.
