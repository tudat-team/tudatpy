# Controlled two-week SLR deformation comparison

All six runs use the same frozen 23 and 131 accepted normal points, combined
postfit SP3 orbits, SLRF2020 histories, measured weather and dated CoM/DHF values.
Native source is `7d52716f6723c704b20ca5f833570acd0f0cc0ba`, kernel SHA256
`80b924287707124a8efea66455b26233b5511b7856983d35c05d21f82734f215`.
All new runs use ICRF axes, Time-valued body/orbit/rotation evaluation and the
rotation-owned observable converter. No fitting, mean removal or residual cuts.

Raw one-way observed-minus-computed statistics, in millimetres:

| Configuration | Mar31–Apr06 mean | RMS | Apr07–13 mean | RMS |
|---|---:|---:|---:|---:|
| Standard IERS ocean EOP / linear C04; solid + pole + GOT | -0.0603 | 7.2210 | +7.6589 | 17.4349 |
| DS159 / four-point C04; same deformation | -0.0604 | 7.2216 | +7.6586 | 17.4349 |
| DS159 / four-point; add CM ocean pole + RP03 S1/S2 | -0.4052 | 7.1416 | +7.5355 | 17.4564 |

The added tides change residuals by RMS 0.4388/0.4716 mm and maximum absolute
0.9916/1.0097 mm in the two weeks. Standard-to-DS maximum changes are
5.33/27.69 µm. Historical geodetic GOT results to the new standard baseline
change by at most 11.72/21.30 µm; that control combines corrected frame labels,
Time rotation and the native converter/kernel update. It does not isolate one
of those changes. Historical outputs remain unchanged in
`../stable_cubic_two_weeks/`.

`summary.json` records raw means, uncentred RMS, centred scatter and statistics
for each station, exact input/kernel/source hashes and Git blob IDs. Reconstruct
executed sources with `git show COMMIT:PATH`; no source snapshots are duplicated.
The six CSVs preserve all numerical strings, reducing only source paths to
basenames. Original external CSV/summary hashes remain in the summary.
`residual_changes.csv` pairs every accepted point, with new-minus-old residual
changes for the three comparisons. Station IDs, event tags, observed intervals
and CoM/DHF fields were checked against the historical frozen CSVs.

The high-tidal profile was declared from supported real provider models before
its residuals were evaluated. The two weeks respond differently; neither a
profile nor an offset is selected from these results. GOT remains the explicit
17-term geodetic HARPOS realization, distinct from exact ASI BLQ/HARDISP.
Non-tidal loading is excluded pending basis and reference-displacement support.
These raw residuals against a fitted, millimetre-quantized orbit are not an
independent orbit validation or a normal-point noise-level claim. Detailed
model choices and reproduction arguments are in the parent analysis README.
