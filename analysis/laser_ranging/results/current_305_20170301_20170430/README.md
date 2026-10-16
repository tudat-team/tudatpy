These current results exclude exactly one Grasse A15 normal point at Dominic's
explicit request: `2017-03-09T18:34:15.297622200200` UTC,
`apollo15_201703.npt` line171. Its physical cause is not established. Selection
uses station, source filename, line and verbatim UTC together, never residual
size. All other CSV records, including their computed ranges and metadata, are
copied **byte-for-byte** from the immutable original 306-point packages. No
observation is recomputed, reweighted, demeaned or fitted.

| Model / explicit prism delay | Station | N | Raw mean (mm) | RMS (mm) | Scatter (mm) | Median absolute (mm) | 95th percentile absolute (mm) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| DE440 / JPL n h | APOLLO | 69 | +1.830 | 58.488 | 58.460 | 50.281 | 104.197 |
| DE440 / JPL n h | Grasse | 236 | −12.959 | 70.101 | 68.892 | 48.721 | 140.127 |
| INPOP17a / JPL n h | APOLLO | 69 | −47.439 | 49.408 | 13.807 | 45.630 | 71.250 |
| INPOP17a / JPL n h | Grasse | 236 | −60.826 | 63.920 | 19.648 | 58.420 | 95.991 |
| INPOP17a / no explicit glass | APOLLO | 69 | −3.363 | 12.341 | 11.874 | 8.212 | 22.850 |
| INPOP17a / no explicit glass | Grasse | 236 | −12.807 | 22.095 | 18.004 | 12.368 | 37.351 |

RMS is about zero; scatter is descriptive and is never subtracted. Absolute
order statistics use every retained raw value without centering; the 95th
percentile uses NumPy's default linear interpolation. Full precision,
station/reflector statistics and separately descriptive temporal slopes are in
[summary.json](summary.json). Grasse bin RMS remains metadata, not a formal NP
uncertainty. Neither prism convention is established by these residual sizes;
the published INPOP optical reference definition remains unresolved.

New figures contain raw markers and supplied APOLLO error bars only. Residuals
are never connected; the plotting helper's optional nightly means also use
unconnected markers. Zero reference lines remain. Each case includes full-span
and separately labelled ±200mm-detail PNG/PDF:

| Case | Full-span PNG / PDF | Detail PNG / PDF |
| --- | --- | --- |
| DE440 | [PNG](de440/figures/residuals_full_span.png) / [PDF](de440/figures/residuals_full_span.pdf) | [PNG](de440/figures/residuals_y_detail.png) / [PDF](de440/figures/residuals_y_detail.pdf) |
| INPOP17a n h | [PNG](inpop17a_jpl_nh/figures/residuals_full_span.png) / [PDF](inpop17a_jpl_nh/figures/residuals_full_span.pdf) | [PNG](inpop17a_jpl_nh/figures/residuals_y_detail.png) / [PDF](inpop17a_jpl_nh/figures/residuals_y_detail.pdf) |
| INPOP17a no explicit glass | [PNG](inpop17a_zero_explicit_glass/figures/residuals_full_span.png) / [PDF](inpop17a_zero_explicit_glass/figures/residuals_full_span.pdf) | [PNG](inpop17a_zero_explicit_glass/figures/residuals_y_detail.png) / [PDF](inpop17a_zero_explicit_glass/figures/residuals_y_detail.pdf) |

The remaining curved temporal and reflector-dependent signatures are visible
without a fitted adjustment. Excluding this isolated point changes Grasse's
distribution and its March9 A15 nightly mean; it does not establish a physical
explanation or remove other discrepancies. Two months do not identify a period.
The original [DE440306](../paired_long_arc_20170301_20170430/README.md) and
[both INPOP306 controls](../inpop17a_long_arc_20170301_20170430/README.md), including
their original figures and all reduction limitations, remain unchanged.

Reproduce this postprocessing and all new figures from the repository root:

```sh
python analysis/laser_ranging/llr_requested_exclusion.py --output "$NEW_OUTPUT"
```

Use the existing TudatPy analysis environment; this command imports the shared
statistics/plot helpers but never constructs an environment or solves an
observation. It checks the pinned original CSV hashes, validates source 69+237
before retaining 69+236, and preserves original line bytes. Source and output
hashes are in [provenance.json](provenance.json).

For future reductions, `llr_long_arc_dataset.json` explicitly distinguishes
source and retained counts. The native reader still validates every source
record and its input hash. Environment construction receives all original
session points and measured weather, then the observation loop skips the one
requested identity. The affected session remains 14 source NPs and 14 measured
weather records; its retained output has 13 NPs. Weather interpolation, table
padding, model physics, clocks, calibrations and tolerances remain unchanged.
