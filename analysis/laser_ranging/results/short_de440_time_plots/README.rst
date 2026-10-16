Short DE440 time views
======================

These figures retain all eleven APOLLO INPOP17a and two Matera SLRF2020 control
points from ``../coherent_station_ephemeris_controls/``. Each panel shows one
reflector in millimetres of raw one-way observed-minus-computed range. No offset,
trend or point is removed. APOLLO error bars use its provider normal-point
uncertainty; Matera has no normal-point uncertainty supplied, so its return-bin
RMS is not substituted.

``residuals_full_span`` retains the actual 2017--2022 gap. The labeled APOLLO and
Matera detail figures resolve their separate intra-night time spans. PNG and PDF
versions are included. Missing reflector panels explicitly say that there are
no observations. ``plot_summary.json`` records source CSV hashes and counts,
raw statistics, and the plot script hash. A single night cannot establish a
multi-night trend; in particular, the two Matera points cannot establish one.

From the repository root, reproduce the figures in a new output directory::

    python analysis/laser_ranging/plot_llr_residuals.py \
      analysis/laser_ranging/results/coherent_station_ephemeris_controls/de440-inpop17a.csv \
      analysis/laser_ranging/results/coherent_station_ephemeris_controls/de440-matera.csv \
      --station APOLLO --station MATERA --station-details \
      --title 'DE440 controls: APOLLO INPOP17a / Matera SLRF2020' \
      --output /path/to/new-output

The same plotter accepts the long-arc CSV columns ``utc_launch``, ``station``,
``reflector``, ``one_way_residual_m`` and optional ``one_way_sigma_m``. It needs
NumPy and Matplotlib, with no native Tudat import or model evaluation. Add
``--nightly-means`` to connect explicitly labeled unweighted raw nightly means.
Their descriptive slopes are reported separately for each station/reflector
only when at least three UTC nights exist; these slopes are never subtracted
from the observations. Even two months do not establish a robust monthly period.
