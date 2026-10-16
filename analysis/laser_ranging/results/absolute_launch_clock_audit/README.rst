Absolute UTC launch-epoch audit
===============================

Three original APOLLO/Grasse observations spanning March-April 2017 were
solved unchanged. Independent ERFA 2.0.1.5 calendar, UTC→TAI→TT and full
SOFA-derived TDB−TT calculations agree with the native absolute TDB launch
epoch within 2.402 ps. The native 1.410–1.575 ms TDB−TT term is present once.
The selected native UT1 is supplied only for the independent topocentric
term; this does not independently validate the EOP data.

Both native tabulated settings and the returned observation collection
retain the original ``Time`` epoch exactly. The Grasse dependent-variable
transmit timestamp is a double, rounded by +33.447 ns; this is output
serialization, not the solver input. All three diagnostic doubles equal
``float(launch)``. Split-Time differences were evaluated before float
serialization. Raw UTC labels remain verbatim, including the Grasse
12-digit fraction; its CRD double-seconds parse differs from the independent
decimal-label parse by 3.477 ps.

``summary.json`` records all source/binary hashes, split epochs, results and
source routing at revision ``55ce93920``. Kernel ``cab49e738cfcf0ea5f251cb140ab8919d5372dbc3e5ba09314d0d5be4eedfcef``
was loaded from ``/home/dominic/Tudat/laser-ranging-build/src/tudatpy/kernel.so``;
native source was ``a366f436c6f6ff03aa57f8b1610075784196bdf8``.

The portable script changes only path handling and formatting from the
recorded external probe; packaging did not rerun observations. Reproduce
with the same pinned source products and environment::

    python analysis/laser_ranging/results/absolute_launch_clock_audit/audit_absolute_launch_clock.py \
      --data-root /path/to/longarc-20170304 --kernel-dir /path/to/llr/kernels \
      --de440-header /path/to/header.440 --sofa-source /path/to/sofa/src/dtdb.c \
      --output /new/path/absolute-launch-clock-audit.json

The `SOFA manual <https://www.iausofa.org/s/sofa_ts_c.pdf>`_ describes the
independent calendar/time route. SOFA ``dtdb.c`` note 7 gives better than
±3 ns relative to DE405 numerical time ephemerides over 1950–2050, not a
verified bound versus DE440. The cached original DE440 header has zero
coefficients/subintervals for TT−TDB (column 15), so an independent DE440
time-ephemeris comparison was unavailable. This audit excludes a missed
millisecond conversion; it does not claim an independent full light-time
model or exact DE440 clock integration. Hypothetical extra-shift range
changes in the JSON are sensitivities, never residual corrections.
