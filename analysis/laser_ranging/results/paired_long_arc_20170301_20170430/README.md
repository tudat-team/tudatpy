All 306 frozen March–April 2017 normal points completed with the same
DE440/ICRF model and independently preselected INPOP17a G+IR station pair.
No mean, offset or trend is fitted to the observation model. No normal points
are clipped or removed, no calibration is reapplied, and all source wavelengths
and solver tolerances remain unchanged. The native wavelength guard now allows
the reporting interval around nominal 1064 nm and evaluates the actual 1064.2 nm.
This does not establish a new atmospheric ray-tracing accuracy envelope.

| Station | Points | UTC nights | Raw mean (mm) | RMS about zero (mm) |
| --- | ---: | ---: | ---: | ---: |
| APOLLO | 69 | 9 | +1.829591 | 58.488480 |
| Grasse 7845 | 237 | 24 | +10.006628 | 359.582805 |

Every APOLLO output row, including its residual and geometry, exactly matches
the preceding 69-point control. All frozen source keys, reflector counts and
wavelengths match: Grasse has 224 IR points at 1064.2 nm and 13 green points at
532.1 nm. Green measurements occur only on April 3–4, so wavelength and time
effects are not independently resolved by separate group means. Transmit
elevations span 15.50–68.95 degrees. Grasse measured weather is interpolated
for 167 points and held at the last session measurement for 70, within the
declared 1000 s limit. Weather is never interpolated across nights.

The Grasse RMS is dominated by one retained +5.429946840 m A15 residual at
2017-03-09T18:34:15.297622200200 UTC, source `apollo15_201703.npt` line 171.
It contributes approximately 96.2% of the complete Grasse sum of squared
residuals. Every other Grasse residual has magnitude below 0.2 m. The exact
printed RTT, epoch, target, configuration and native metadata agree with the
frozen inventory. Weather is bracketed, calibration flags say already applied,
and the other 13 points in its block have residuals −71.3 to −36.7 mm. This is
an isolated unexplained observation/model discrepancy; it is not proved to be
a source error and remains in all statistics and figures. Three extreme-record
source audits, including the negative green extreme, are in
`extreme_source_audit.json`. Some provider H4 end times precede later printed
NPs in the same block; each NP retains its own printed epoch without an H4 cut.

A bounded numerical control reconstructs the original full 14-NP/14-weather
session environment and solves only the extreme and its immediate neighbours.
The 30 s tables reproduce all three original computed ranges exactly. With
15 s native `Time` cubic tables and every other input unchanged, computed
half-ranges change by 0, −9.596 and +10.133 micrometres in source order. The
extreme remains +5.429956436 m. This excludes a 30 s interpolation discontinuity
as the explanation for its metre-scale discrepancy, without classifying or
removing the record. Exact values, epochs and hashes are in
`extreme_precision_control.json`.

The [time figures](../paired_long_arc_time_plots/README.rst) show both stations
by reflector with actual UTC gaps and wavelength markers. Their full-range view
contains every point; a separately labeled ±200 mm detail reports the single
point outside its window. Both stations show a nonmonotonic rise toward March 22,
fall toward early April and subsequent recovery where observations exist.
This shared behavior is inconsistent with describing the arc solely by a
constant station offset or one straight-line drift. Sparse sampling over two
months does not identify a robust period or its physical cause. Separate
descriptive slopes in the plotting metadata are never subtracted.

`nearby_nightly_raw_means.csv` pairs each APOLLO reflector/night with the nearest
Grasse nightly mean observation epoch within 12 h, allowing adjacent UTC dates.
There are 17 pairs, including five across a date boundary, with actual gaps
4.675–11.458 h. Raw Grasse-minus-APOLLO mean differences span −31.684 to
+19.386 mm. The March 22 crest appears for all five reflectors at both stations;
their differences are +6.76 to +19.39 mm at gaps 4.68–5.80 h. These are distinct
epochs with changing geometry, not simultaneous measurements or fitted station
biases. The comparison JSON retains the 16 APOLLO and 49 Grasse reflector/nights
not paired, with their nearest-other time gaps. All their normal points remain
in the complete residual CSV.

`residuals.csv` preserves raw O−C, provider identifiers and uncertainties where
supplied, weather handling, elevation, lunar phase, PA Earth/Sun directions and
reconstructed line-of-sight diagnostics. Geometry uses recorded native event
doubles for interpretation; the production solver retains `Time` internally.
`summary.json` contains compact model/input/source/kernel provenance and raw
station/reflector statistics. Repeated group metadata are consolidated into
`group_timing_and_weather.csv` and common flags; the original external summary
is retained and its hash recorded. `verification.json` records complete counts
and control equality. The frozen input URLs and hashes remain in
`../../llr_long_arc_dataset.json`.

Runtime was 26.694 s: 5.813 s setup and 19.359 s solving across 79 groups.
No runtime optimization was introduced. `provenance.json` pins source revision,
native revision, imported kernel hash and the portable command. The preceding
failed guard attempt and old-kernel APOLLO result remain preserved separately.
