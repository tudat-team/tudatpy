The frozen March–April 2017 APOLLO arc contains all 69 provider normal points
on nine UTC nights, with the same DE440/ICRF, INPOP17a station, weather and
complete correction model used for the short control. No mean, offset or trend
is fitted to the observation model; no points are removed. The raw one-way
O−C mean is +1.829591 mm, RMS about zero 58.488480 mm and scatter about the
mean 58.459857 mm. The eleven March 9 residuals reproduce the earlier control
exactly. Runtime was 5.332 s, including 0.633 s setup and 4.244 s solving.

The PNG and PDF retain every raw point and provider NP uncertainty. Diamonds
and dotted lines show unweighted raw nightly means, with their actual time
gaps. The plotted source is `analysis/laser_ranging/plot_llr_residuals.py`,
revision `17eacc162`; its exact hash and the CSV hash are in `plot_summary.json`.

The visible variation is not a single monotonic drift. A11, A14, A15 and L1
rise toward March 22 and fall toward early April; later measurements recover
for A11, A14 and A15. Raw peak-to-peak spans are 88–159 mm across reflectors.
The separate descriptive slopes of unweighted nightly means are −0.552,
−0.8158, −0.9854, −0.0244 and −2.5694 mm/day for A11, L1, A14, A15 and L2.
These slopes are not subtracted or used as corrections. Four to nine observed
nights per reflector do not establish a period or distinguish a lunar geometry
signature from station and reduction effects. The same-reflector Grasse
comparison remains necessary.

`residuals.csv` retains source identities, measured weather handling and lunar
geometry diagnostics. `nightly_raw_statistics.csv` contains raw nightly
statistics. `summary.json` records input, model, kernel and source hashes;
`provenance.json` pins the executed revisions and portable command. The frozen
inputs and station conventions are in `../../llr_long_arc_dataset.json`.

The preceding full 306-point attempt completed APOLLO and stopped before the
first Grasse solve: the native wavelength guard rejected the recorded 1064.2 nm
against its 1064.0 nm upper bound. Its source/kernel hashes and error are retained
in the provenance. No wavelength was rounded and no Grasse result is claimed
in this APOLLO-only package.
