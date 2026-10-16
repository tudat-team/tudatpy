Paired APOLLO and Grasse time plots
==================================

The plots retain all 69 APOLLO and 237 Grasse normal points from the fixed
March--April 2017 interval, using the DE440/INPOP17a station-pair reduction.
The raw CSV and physical/source/kernel conventions are in
``../paired_long_arc_20170301_20170430/``. No mean, offset, trend, uncertainty or
point has been changed for these figures.

Each panel shows one reflector in millimetres of raw one-way O-C against UTC
launch time. Station colors remain fixed. Grasse's recorded 1064.2 nm infrared
points are triangles; its 532.1 nm green points are circles. APOLLO error bars
use its provider normal-point uncertainties; Grasse has no supplied NP sigma,
and its return-bin RMS is not substituted. Dotted diamond connections are
unweighted raw nightly means, retained separately from the raw points; these
connections are visual guides, not a fitted or interpolated prediction. Grasse
nightly means pool both wavelengths, as explicitly labeled, so color changes
must not be mistaken for a temporal trend.

``residuals_full_span`` retains the +5.429946840 m Grasse A15 observation on
2017-03-09. ``residuals_y_detail`` is an additional, explicitly labeled
[-200, 200] mm view; its A15 panel states that one raw point lies outside that
window. The full view, nightly means, statistics and source CSV retain that
point. Both figures have PNG and PDF versions, and both were visually reviewed.

Several reflectors show repeated rises and falls visible at both stations.
``plot_summary.json`` records every input hash/count and separate unweighted
nightly-mean slopes by station and reflector. Those slopes are descriptive,
sensitive to sampling, wavelengths and the retained extreme, and are never
subtracted. Two months do not establish a robust monthly period or identify a
physical mechanism.

Reproduce into a new output directory from the repository root::

    python analysis/laser_ranging/plot_llr_residuals.py \
      analysis/laser_ranging/results/paired_long_arc_20170301_20170430/residuals.csv \
      --nightly-means --zoom-mm -200 200 \
      --title 'DE440 / INPOP17a station pair: March–April 2017' \
      --output /path/to/new-output

The optional y-detail argument always saves the full-span figure first and
never changes residual statistics. The preceding short-control figures remain
unchanged in ``../short_de440_time_plots/``, with their original source identity.
