These are descriptive summaries of the completed calendar-2017 cases: **311
APOLLO and 1268 Grasse NPs**, 1579 retained per model. Source selection excludes
only the explicitly requested March9 point and two independently impossible
nanosecond-RTT lunar records. All remaining points, including the 110 m Grasse
discrepancy, remain in every raw distribution. No offset, drift, period, station
parameter or other correction is fitted or subtracted.

The [raw annual packages and marker-only figures](../year_2017/README.md) preserve
source fields, counts, uncertainties and reduction limitations. This directory
adds monthly medians/N/P10/P90, dated nightly extrema and fixed 12 h cross-station
comparisons without duplicating raw residuals or figures.

Monthly raw medians below are in millimetres; A/G denotes APOLLO/Grasse.
These ALL-reflector pools have unequal target, wavelength and lunar-phase
sampling, so changes are descriptive and must not be called station drift.
Per-reflector monthly values and empty months remain explicit in each case's
`monthly_raw_medians.csv`.

| UTC month | N A/G | DE440 median A/G (mm) | INPOP n h median A/G (mm) | INPOP no separate glass median A/G (mm) |
| --- | ---: | ---: | ---: | ---: |
| 2017-01 | 4 / 52 | +66.4 / -50.2 | -59.2 / -58.3 | -19.2 / -14.9 |
| 2017-02 | 27 / 64 | +23.5 / -64.2 | -48.0 / -70.0 | -1.4 / -26.0 |
| 2017-03 | 41 / 94 | +22.9 / -39.0 | -49.2 / -59.0 | -5.5 / -13.8 |
| 2017-04 | 28 / 142 | -29.0 / -28.6 | -43.7 / -56.8 | +1.1 / -8.4 |
| 2017-05 | 38 / 127 | +4.1 / -16.9 | -57.5 / -47.4 | -12.0 / -1.0 |
| 2017-06 | 29 / 121 | +24.2 / +29.0 | -58.1 / -67.3 | -18.1 / -20.9 |
| 2017-07 | 23 / 63 | -14.6 / +8.1 | -54.1 / -72.2 | -13.6 / -24.8 |
| 2017-08 | 15 / 116 | -6.8 / -11.6 | -49.7 / -59.7 | -8.0 / -19.5 |
| 2017-09 | 19 / 108 | -44.4 / +16.1 | -55.7 / -58.8 | -15.7 / -18.5 |
| 2017-10 | 48 / 200 | +16.1 / -9.7 | -53.4 / -58.1 | -10.0 / -14.8 |
| 2017-11 | 16 / 118 | +18.4 / -60.2 | -49.8 / -57.2 | -9.8 / -9.4 |
| 2017-12 | 23 / 63 | +68.4 / -16.3 | -17.1 / -57.6 | +24.9 / -2.5 |

The pooled monthly median span is −44.4 to +68.4mm for APOLLO/DE440 and−19.2
to +24.9mm for INPOP without separate glass; Grasse spans −64.2 to +29.0mm and
−26.0 to −1.0mm respectively. INPOP n h spans −59.2 to −17.1mm and −72.2 to −47.4mm.
These alternatives change orbit, rotation, coordinates and published surface/map
conventions together; this comparison does not isolate one omitted term or
establish the unresolved optical reference convention.

The annual figures show repeated nonmonotonic rises/falls in DE440, shared by
both stations at nearby epochs, plus persistent reflector separation. INPOP
reduces those large excursions but leaves structure. For example, A11 monthly
medians have a midyear trough at both stations: INPOP without separate glass
gives APOLLO June/July −42.8/−38.4mm and Grasse −48.5/−48.6mm. L1 changes from
negative February medians (APOLLO −19.7, Grasse −38.0mm) to positive May medians
(+18.8,+22.8mm). A15 remains mostly negative through November; APOLLO's December
median becomes +22.9mm while Grasse is−3.9mm, but their December time coverage
differs. These observations do not establish a unique period or cause.

The dated nightly means retain shorter structure that monthly pools obscure.
For INPOP without separate glass, APOLLO reflector peak-to-peak nightly means
span 63.6–92.2mm, versus 177.8–287.3mm for DE440. Extrema dates/counts are retained
per reflector, including sparse nights; they are not trend fits. Late-December
APOLLO rises occur for several reflectors, while Grasse ends on December 21,
so later APOLLO excursions cannot be treated as contemporaneous station
comparisons. No straight-line detrending or periodic correction is applied.

All-point absolute order statistics supplement means/RMS in the core table:

| Station | N | DE440 raw P95 absolute (mm) | INPOP n h (mm) | INPOP no separate glass (mm) |
| --- | ---: | ---: | ---: | ---: |
| APOLLO | 311 | 146.130 | 89.357 | 46.205 |
| Grasse | 1268 | 151.711 | 107.421 | 50.259 |

These percentiles use every retained value, without centering or rejection.
The Grasse L1 point at 2017-06-22T09:47:06.763983741400 UTC,
`luna17_201706.npt:154`, is +110.458233/+110.264220/+110.327534m in the three
cases. It dominates complete Grasse RMS (about 3.1 m); it is not attributed to
weather or source error and is not removed. The INPOP no-separate-glass June22
L1 nightly mean remains +36.778597m across three NPs. The long-weather-hold flag
is preserved, but that coincidence does not establish causation.

Each APOLLO station/reflector nightly mean is paired with the nearest Grasse
one at most 12 h away, allowing adjacent UTC dates. All models have the same 64
pairs (18 across a date boundary), actual gaps 2.134–11.458h. Pairs are not
interpolated or shifted in time, and this descriptive subset does not remove
unpaired observations from any raw statistics.

| Reflector | Pairs | DE440 median G−A (mm) | INPOP n h (mm) | INPOP no separate glass (mm) |
| --- | ---: | ---: | ---: | ---: |
| ALL | 64 | +0.876 | -4.073 | -4.075 |
| A11 | 14 | -2.943 | -4.569 | -4.569 |
| A14 | 12 | -1.083 | -3.741 | -3.747 |
| A15 | 20 | -9.284 | -7.867 | -7.867 |
| L1 | 8 | +1.534 | -2.036 | -2.034 |
| L2 | 10 | +11.489 | +3.394 | +3.397 |

The no-separate-glass pair differences span −30.131 to +14.437mm, with P10/P90
−20.714/+7.955mm; DE440 spans −31.684 to +39.067mm, with −20.937/+27.046mm.
These distinct-epoch differences include changing geometry and are not
simultaneous station biases. The two INPOP pair medians are nearly identical;
cross-station agreement does not resolve a reflector term common to both.
There are 94 unmatched APOLLO and 319 unmatched Grasse reflector/nights, stored
once in `unmatched_station_reflector_nights.csv`, identical for all three cases.
Every matched pair retains actual mean UTC epochs, signed gap and NP counts.

Reproduce each case from its saved raw CSV without loading kernels or computing
observations:

```sh
python analysis/laser_ranging/llr_temporal_summary.py \
  analysis/laser_ranging/results/year_2017/de440/residuals.csv \
  --output "$NEW_OUTPUT"
```

Use the corresponding INPOP case directory for either alternative. The small
postprocessor reuses the existing raw-CSV reader. It writes an identical
unmatched-night CSV for each case; this package stores their common copy once.
Monthly and paired CSVs reproduce the independently prepared external analysis
exactly. Each case pins its raw-input and script SHA256; `provenance.json` also
pins Git blobs and packaged outputs. No model solve, fit or further source
exclusion is part of this analysis.
