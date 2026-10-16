Fixed-observation station and ephemeris controls
==============================================

All cases use ICRF axes, the same eleven APOLLO or two Matera observations,
30-second native Time cubic tables, no loading, and the original 50-iteration,
``1e-13`` fractional light-time criterion. No offset, mean or trend is fitted;
no point or uncertainty is changed. DE440 controls reproduce the preceding
corrected-frame raw CSVs exactly on the new verified kernel.

===========================  ===========  ==========  ============
Configuration                Mean (mm)   RMS (mm)    Scatter (mm)
===========================  ===========  ==========  ============
DE440 / APOLLO INPOP17a       -25.500      49.900      42.893
DE440 / Matera SLRF2020       -65.200      65.228      1.887
DE440 / APOLLO Hofmann        -3.119       42.083      41.967
DE440 / Matera Hofmann        -24.635      25.318      5.844
DE421 MER / APOLLO INPOP17a   -267.197     282.428     91.493
DE421 MER / Matera SLRF2020   -305.971     306.031     6.041
===========================  ===========  ==========  ============

RMS is about zero; scatter is a descriptive statistic, without subtracting
anything from the raw residuals. Hofmann is the independently preselected
standard station solution. DE421 uses its own translation, PCK, explicit MER
frame, source-header mass parameters and published frame-compatible reflector
coordinates. These are controlled sensitivities, not reproductions of their
original fitted reductions. No default is selected by these residuals.

The APOLLO Hofmann reflector means remain A11 -61.477, A14 +35.805,
A15 -28.315 and L1 +47.247 mm. The short sequence alone cannot separate
reflector-dependent errors from temporal signatures. The two Matera points
have no supplied normal-point uncertainty; their bin RMS is not substituted.

``summary.json`` records the actual source/kernel/input hashes, station and
lunar-frame conventions, per-reflector statistics and original output hashes.
The six CSVs preserve every numerical field. ``computed_half_range_changes.csv``
compares each station or ephemeris change with its declared DE440 control.

``native_interpolation_audit.json`` compares the existing direct and 60/30/15/10 s
Time routes at every frozen launch and recorded bounce/reception epoch with
independent 80-digit evaluation of the stored DE440 coefficients. Absolute
BCRS table errors remain below 48 micrometres; shorter spacing does not
monotonically improve the rounding floor. The existing Moon-relative-Earth
Time tables remain below 0.103 micrometres. These state controls are not
full-range errors, and a simultaneous relative vector cannot replace a
retarded two-epoch BCRS light-time leg.
