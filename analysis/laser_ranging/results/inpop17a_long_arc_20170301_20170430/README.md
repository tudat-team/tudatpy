Both predeclared INPOP17a optical sensitivities retain all 306 frozen March–April
2017 points, their source fields, wavelengths, measured weather and solver
tolerances. No offset, trend, empirical libration, calibration or point cut is
fitted or applied. DE440 remains the primary/default model. These controls change
the matched orbit/rotation, lunar coordinates, GMs and documented surface/map
conventions together; they do not isolate one cause of the DE440 discrepancy.

| Model / explicit prism delay | Station | N | Raw mean (mm) | RMS (mm) | Scatter (mm) | Median absolute (mm) | 95th percentile absolute (mm) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| DE440 reference / JPL n h | APOLLO | 69 | +1.830 | 58.488 | 58.460 | 50.281 | 104.197 |
| DE440 reference / JPL n h | Grasse | 237 | +10.007 | 359.583 | 359.444 | 48.864 | 142.072 |
| INPOP17a / JPL n h | APOLLO | 69 | −47.439 | 49.408 | 13.807 | 45.630 | 71.250 |
| INPOP17a / JPL n h | Grasse | 237 | −37.637 | 358.749 | 356.770 | 58.571 | 97.228 |
| INPOP17a / no explicit glass | APOLLO | 69 | −3.363 | 12.341 | 11.874 | 8.212 | 22.850 |
| INPOP17a / no explicit glass | Grasse | 237 | +10.348 | 356.316 | 356.166 | 12.529 | 37.777 |

RMS is about zero; scatter is a descriptive standard deviation and is never
subtracted. Absolute order statistics use every raw station residual, without
centering or rejection; the 95th percentile uses NumPy's default linear
order-statistic interpolation. Full precision and station/reflector/wavelength
statistics are in [summary.json](summary.json).

The same unexplained Grasse A15 point remains at
2017-03-09T18:34:15.297622200200 UTC, `apollo15_201703.npt` line171:
+5.434879780 m with n h, +5.474904299 m without, versus +5.429946840 m DE440.
Its source audit remains in the [DE440 package](../paired_long_arc_20170301_20170430/README.md).
It is not classified as erroneous or removed. Grasse formal NP uncertainties
remain unavailable; its bin RMS is metadata, not an inferred standard error.
Neither case establishes intrinsic measurement-noise performance.

The [n h full-span PNG](full/jpl_nh/figures/residuals_full_span.png) and
[PDF](full/jpl_nh/figures/residuals_full_span.pdf), and the
[zero-explicit-glass full-span PNG](full/zero_explicit_glass/figures/residuals_full_span.png)
and [PDF](full/zero_explicit_glass/figures/residuals_full_span.pdf), retain every
point. Separate ±200 mm detail views report points outside their window:
[n h PNG](full/jpl_nh/figures/residuals_y_detail.png) /
[PDF](full/jpl_nh/figures/residuals_y_detail.pdf),
[zero PNG](full/zero_explicit_glass/figures/residuals_y_detail.png) /
[PDF](full/zero_explicit_glass/figures/residuals_y_detail.pdf).
Each full case also contains APOLLO/Grasse time-detail PNG/PDF files. Plot metadata
keeps any descriptive nightly slope separate from unchanged O−C.

The time plots change both the large reflector ordering and the temporal shape.
On March22 the APOLLO five-reflector raw nightly-mean spread is 139.25 mm for
DE440, 32.04 mm for INPOP with n h, and 17.06 mm without explicit glass. A shared
nonmonotonic signature remains: all five APOLLO reflectors have a March22 trough
and recover toward April3. With zero explicit glass, March22 nightly means
A11/A14/A15/L1/L2 are −28.475/−19.863/−21.537/−11.548/−28.609 mm.
Grasse also shows a trough near March22 for the sampled reflectors and subsequent
recovery, with remaining station/reflector differences. The smaller common curve
cannot be reduced to a fixed reflector offset, but two sparsely sampled months
do not establish its period or physical cause. The raw nightly comparison is
in [nightly_model_comparison.csv](nightly_model_comparison.csv); no trend or
constant is fitted to the model or subtracted from these values.

This remains an incomplete INPOP reduction. The [source/table audit](../llr_physical_terms_audit/ephemeris_controls.rst)
records the reviewed Lunokhod baseline repair, with no recovered author erratum.
Manche2011 Eq9.1 and the published coordinate tables do not identify a prism
front-face/vertex/effective-bounce reference. Both n h per leg and zero explicit
glass were therefore declared before evaluating the twelve control points.
Their different RMS values do not establish which convention is correct.

The official archive README labels only INPOP17a, not G+IR. Its header constants
strongly associate it with Table7's G+IR column: CMR2M=.3931482590053706,
TAUM=.08655346141536852 d, KCMB/CMR2M=1.7531636063880315e−8 and
C32M=4.844405991342345e−6 match that column within printed precision. An explicit
archive-to-column identifier was not recovered, so the association remains an
inference independent of residuals. The final [2018 paper](https://doi.org/10.1093/mnras/sty096)
describes IR as part of INPOP17a but does not uniquely identify this archive.

The native reference-radius degree2 tide is rescaled once to each reflector's
actual radius before the LB spatial map. Degree3 surface constants are not
imported from JPL; the documented <0.2 mm spin deformation is omitted. The map
retains only the existing Sun potential: the Earth contribution to a literal
all-external-body Eq9.32 map is about 20 micrometres at lunar distance. Current
Earth EOP/clock/loading conventions are explicit and are not a replay of the
historical INPOP fit. The 2017 report Table10 p23 row15 additionally gives a
Grasse 2016-06-01–2018-01-01 G+IR estimated bias of +2.46±.01 cm in **two-way**
light-time units (12.3 mm half-range magnitude). Its sign/application is not
inferred here. This unapplied solution-level parameter is distinct from the
already-applied source system calibration and cannot explain the shared time
dependence by itself. No fitted radiation/thermal or equipment bias is imported.
The same report Table4 ends APOLLO groupD coverage at 2016-11-25, Grasse green
at 2017-02-07 and IR/total at 2017-02-19: our March–April2017 groupE/Grasse
observations are later than that fit coverage, making this a fixed-product
extrapolation control rather than a reproduction of published fit RMS.

The predeclared [twelve-point controls](twelve) retain exact raw CSVs and their
PNG/PDF views; two sampled nights do not support a trend/period claim.
[Default equivalence](controls/default_equivalence.json) compares all eleven
original March9 DE440 ranges against the preceding committed core, with exact
equality and maximum difference 0 m. The independent
[Eq9.25 radius control](controls/tide_radius_control.json) agrees with the native
rescaled Earth+Sun vector to 5.48e−11 m at one frozen event.

The twelve-point raw vector exports remain in the external run directories
identified by [provenance.json](provenance.json), with exact hashes and explicit
`event_vector_conventions.json` sidecars. Their Earth body-fixed vectors were
requested with origin `Earth` and precede native BCRS station spatial mapping;
they are **not complete observable endpoints**. Moon custom spatial terms are
included. All exported vectors reconstruct diagnostic double epochs, while the
observable solver retains `Time`. Future exports request origin `SSB` and record
that choice; this diagnostic-only change does not alter any range result.

The executed source and manifest are committed as `e26fa7443`; provenance pins
their Git blobs, SHA256, the cleared eight-kernel pool and native binary hash.
Retrieve originals with `git show REV:analysis/laser_ranging/PATH`; full copied
source snapshots are unnecessary. The later manifest addition only documents
the release-to-fit limitation. Frozen input URLs/hashes remain in
`../../llr_long_arc_dataset.json`. Full original summaries and event vectors
remain outside the repository; their hashes are retained. To reproduce either
case, use the portable command in provenance with `CONVENTION=jpl-nh` or `none`,
or append `--control-twelve` for its predeclared subset. Plot using the existing
`plot_llr_residuals.py`, `--nightly-means --zoom-mm -200 200 --station-details`.
