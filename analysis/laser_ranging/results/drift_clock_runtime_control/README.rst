Linear station-drift clock control
=================================

Only ``station_kinematic_displacement`` now uses geocentric TT for the
published linear station velocities. The measured topocentric/geocentric
TT difference times velocity bounds the displacement change by
8.62e-16 m at APOLLO and 1.69e-15 m at Grasse across the frozen 306 launches.
Observable time normalization, Earth simultaneity and the generic converter
retain their existing implementations.

One sequential baseline/candidate pair used the same machine, data,
``cab49e738cfcf0ea5f251cb140ab8919d5372dbc3e5ba09314d0d5be4eedfcef``
kernel (native source ``a366f436c6f6ff03aa57f8b1610075784196bdf8``),
and numerical thread limits of one. The baseline was clean revision
``3dc1574d9``; the candidate changed only this function.

=================== ============ =============
Measured time (s)   Baseline     Candidate
=================== ============ =============
Setup               5.8897       5.5959
Solver              19.6312      12.7785
Driver total        27.0663      18.9936
External wall       28.0381      19.8782
=================== ============ =============

All 306 rows and 39 fields (11,934 fields), including all 15 geometry
columns, were byte-identical. Nightly statistics, source identities,
weather and physical settings also matched. The measured solver reduction
was 34.9%; this local pair is not a general performance guarantee. No
repeat was needed to distinguish the observed benefit. No residuals were
removed, fitted or corrected.

``comparison.json`` records the comparisons and independent displacement
bound; ``run_provenance.json`` records exact commands, source/binary/CSV
hashes and UTC execution times. Original raw outputs remain locally in
``/home/dominic/Tudat/laser-ranging-work/results/llr-drift-clock-{baseline,candidate}-306``.
Runner and comparison source remain under that work directory's
``verification/`` as ``run_llr_drift_clock_control.py`` and
``compare_llr_drift_clock_controls.py``. The independent bound source is
``verification/audit_grasse_wavelength_and_drift_clock.py`` with its hashed
result identified in ``comparison.json``. Existing primary raw results
remain in ``../paired_long_arc_20170301_20170430/``.
