Laser ranging: first validation milestone
========================================

This checkpoint establishes passing implementation checks and a frozen
fixed-orbit SLR reduction. The lunar diagnostic remains unresolved; these
results do not establish general millimetre accuracy or an exact replay of an
orbit provider's fitted reduction. Model conventions and input selection are
documented in the `SLR analysis notes <../analysis/laser_ranging/README.md>`_ and
the `fixed-ephemeris LLR diagnostic <llr_fixed_ephemeris_diagnostic.rst>`_.

Build and implementation evidence
---------------------------------

The local build is ``/home/dominic/Tudat/laser-ranging-build``, configured against
``/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-ranging``. The latest parser
rebuild used source ``2fe4f1cae73d2f8ac5c778f170f27844b330b256``. Its imported
``src/tudatpy/kernel.so`` is 67,129,648 bytes, SHA-256
``1864ed04530210fa056e78128bdbbfe1ad536cae6a9d3ade4f23295fbeda1744``.
`Kernel provenance <../../../laser-ranging-work/logs/final-kernel-provenance.json>`_
records the interpreter, import origins, environment and successful Python read
of the complete SLRF2020 provider fixture. Numerical-library threads were one;
incremental builds used at most eight jobs, and native suites ran serially.

The `native evidence summary <../../../laser-ranging-work/logs/native-current-summary.json>`_
contains six passing executables, 34 cases and 11,154 assertions. It combines
the affected reader reruns with retained passes; it is not a claim that every
suite was rebuilt and rerun at the final source revision.

.. list-table:: Passing native evidence
   :header-rows: 1

   * - Suite
     - Cases
     - Assertions
     - Evidence origin
   * - LaserRangingData
     - 1
     - 104
     - Latest parser rebuild, 2fe4f1cae
   * - CrdAndSinexReaders
     - 6
     - 59
     - Retained from 9f9f868d7
   * - NWayRangeObservationModel
     - 5
     - 1,019
     - Retained
   * - AtmosphereCorrection
     - 15
     - 9,877
     - Retained
   * - SinexStationHistory
     - 3
     - 74
     - Latest parser rebuild, 2fe4f1cae

The `focused Python report <../../../laser-ranging-work/logs/python-focused.xml>`_
retains 59 passing cases from ``test_data_input_earth_orientation.py``,
``test_station_reference_point_bindings.py`` and ``test_time_conversions.py``.
The `documentation and runtime check <../../../laser-ranging-work/logs/final-python-surface.json>`_
resolved all 87 autosummary names in edited RST files since ``ebe3115ec``.
Its binding smoke exercised the optical correction factory, body-deformation station motion, actual SP3/CRD
records and the selected Earth-rotation time converter. That check used source
``9f9f868d7`` and kernel
``029f3ea61f723093697a8b055e304189d3906ffb541d6fae929008890625deff``;
it was retained after the private SINEX coordinate-parser change. These checks
establish accessible interfaces and their tested behavior, not observational
accuracy.

Frozen SLR result
-----------------

The combined ILRSA fitted LAGEOS-1 SP3 orbit and external SLRF2020 station
solution retain 23 of the 29 preselected normal points over 2024-03-31 through
2024-04-06. No bias, station, orbit or other parameter was fitted; no mean or
trend was removed, and no point was clipped by its residual. Residuals are
equivalent one-way observed minus computed ranges; RMS is about zero.

.. list-table:: Unfitted SLR residuals in millimetres
   :header-rows: 1

   * - Selection
     - Points
     - Mean
     - RMS
   * - Combined
     - 23
     - -0.474
     - 8.843
   * - Graz, 7839 / 11001S002
     - 13
     - -0.402
     - 11.526
   * - Herstmonceux, 7840 / 13212S001
     - 10
     - -0.567
     - 2.675

Six fixed physical/coverage exclusions remain: session/normal-point keys
``53/2``, ``65/3`` and ``65/4`` fall below the provider's 20-degree elevation
cutoff; ``172/4``, ``172/5`` and ``172/6`` lie outside the SP3 interpolation
interior. All six are Graz points. Ocean loading was omitted in this frozen
baseline. A subsequently prepared GOT4.10c helper is not part of these results.
The fitted orbit, station conventions and omitted loading prevent treating the
small combined mean as proof of millimetre model accuracy.

The `SLR script <../analysis/laser_ranging/fixed_orbit_slr.py>`_ and pinned
`dataset manifest <../analysis/laser_ranging/dataset.json>`_ accompany the
`raw residuals <../analysis/laser_ranging/results/slrf2020_combined_2024-03-31_2024-04-06/residuals.csv>`_
and `full summary <../analysis/laser_ranging/results/slrf2020_combined_2024-03-31_2024-04-06/summary.json>`_.
The summary records the executed script and all input hashes and uses the latest
parser kernel identified above. Later precision experiments are separate
diagnostics and do not replace this frozen baseline.

A fresh run reproduced the complete 23-row residual CSV byte for byte, together
with all station/session statistics, six exclusions, interpolation controls and
input hashes. The residual CSV SHA-256 is
``2a17ed360d346c2efa8dc7afd39fc576caff7762e563e62796b9094da46ace33``.
Local repeatability evidence is recorded in
``laser-ranging-work/logs/slr-baseline-repeatability-20261016.json``.

Unresolved LLR diagnostic
-------------------------

All eleven APOLLO points were retained in the 30-second cubic translation-table
experiment: one-way mean **219.550 mm**, RMS about zero **223.164 mm**. No
parameters were fitted and no residual selection was applied. This is an
unresolved diagnostic, with no millimetre-accuracy claim. The approximate station
survey, unmatched fitted station solution/calibration and omitted displacement
terms are detailed in the existing LLR diagnostic, together with the separate
numerical convergence controls.

See the `LLR script <../analysis/laser_ranging/llr_residuals.py>`_,
`raw normal-point results <../../../laser-ranging-work/results/apollo-20170309-cubic30/normal_points.csv>`_
and `full summary <../../../laser-ranging-work/results/apollo-20170309-cubic30/summary.json>`_.
This retained experiment used kernel
``c475a3e90a39eed4332442f58c631398cd608f2195d25bc08cd59e3cfd0bd866``;
its summary preserves its own script, input and resource identities. It was not
rerun merely to update this report or after the unrelated SINEX parser changes.

Native/raw build logs and provider source products remain in the external
``/home/dominic/Tudat/laser-ranging-work`` evidence cache or resource directory.
Selected portable derived results are committed in the repository; the SLR links
above point to those copies. LLR artifact links currently refer to the local cache.

Precise-Time binding and lunar elastic-tide audit
------------------------------------------------

The subsequent binding build used native source ``8506f9188`` and kernel SHA-256
``bd8d893574bd29281188934ff5a62647e4c08494eb4bfb846c17cfc06a062968``
(67,149,448 bytes, same build/import path). The two environment tests passed,
including precise custom-ephemeris and Body dispatch for epochs 20 ns apart
that share one double value. All 57 environment autosummary names resolved.
`Binding provenance <../../../laser-ranging-work/logs/time-callback-kernel-provenance.json>`_
distinguishes the native source from later script-only changes; the six-suite
native evidence above remains retained, without unrelated rebuilds or reruns.

An independent Legendre-polynomial degree-2/3 displacement calculation checked
all eleven saved native reflection events from the INPOP station sensitivity
run. The actual factory supplied ``R=1738000 m``, ``h2=0.0423``, ``l2=0.0107``,
``h3=0.0234`` and ``l3=0.0030`` with loaded DE440 gravitational parameters.
Native intrinsic reflector displacement agreed within ``1.49e-10 m``; the
final reflector agreed with one full elastic tide plus the stated spatial
correction within ``2.61e-10 m``. The native PA rotation and an independent
SPICE ``J2000`` to ``MOON_PA_DE440`` rotation differed by at most ``5.56e-16``
per matrix entry. Runtime tracing confirmed exact TDB ``Time`` routing to the
Moon, Earth and Sun ephemerides. The saved event outputs are doubles; changing
the input over their full ULP changed the calculated tide by at most
``5.21e-13 m``. No elastic-tide implementation bug or duplicate/missing full
equilibrium tide was found by this bounded audit.

The `audit script <../../../laser-ranging-work/verification/audit_lunar_elastic_tides.py>`_
and `detailed evidence <../../../laser-ranging-work/logs/lunar-elastic-tide-audit.json>`_
retain per-event vectors, separate Earth/Sun degree contributions and source
hashes outside the repository. This agreement does not explain the remaining
LLR centimetre discrepancies. JPL2020 gives omitted dissipation of about 2 mm
or less and thermal motion of 1–1.5 mm for Apollo and about 5 mm vertically for
Lunokhod, with horizontal components; these are not uniformly submillimetre
effects. No residual correction, fitting or tide-convention toggle was used.

Atomic interval and cubic arithmetic controls
--------------------------------------------

An independent clock audit used all eleven saved APOLLO INPOP station-profile
transmit/receive events, the selected rotation/converter, DE440 Earth velocity,
external Sun/planet/Moon potentials and the station ``v_E dot r / c^2`` endpoint
term. It used ``LC=(LB-LG)/(1-LG)`` and split ``Time`` offsets, without subtracting
large epoch floats. The geocentric interval agreed within **0.397 micrometres**
of equivalent one-way range; the complete station interval agreed within
**0.094 mm**. Native atomic-minus-TDB interval corrections ranged from -63.487
to -53.815 mm. The sign, normalization and single station clock term agreed
with the independent physical calculation, providing no centimetre-scale clock
explanation for the residuals.

The existing topocentric inverse stops at a 1 ps forward residual per endpoint,
which permits about 0.300 mm of interval error in one-way range. Accounting for
those measured inverse residuals left at most 4.913 micrometres of difference
from the leading-post-Newtonian calculation; this is a bounded diagnostic,
not an exact DE440 time ephemeris or an adjusted production correction.
The `clock evidence <../../../laser-ranging-work/logs/llr-atomic-interval-audit.json>`_
and `script <../../../laser-ranging-work/verification/audit_llr_atomic_interval.py>`_
record truncation allowances and the saved receive-event precision limit.

For the archived final A14 point, the receive-event double advances one ULP
between atmospheric and complete solves. The resulting native clock-offset
change accounts for -96.922 micrometres of the -118.318 micrometre difference
between the sequential glass increment and its fixed-event value. The remaining
-21.396 micrometres involves event geometry/solve arithmetic; the registered
per-leg glass correction remains correct. This
`separate check <../../../laser-ranging-work/logs/llr-glass-clock-check.json>`_
used the preceding ``bd8d8935`` kernel and changed no convergence setting.

The subsequent cubic interpolation change evaluates relative to the nearer
endpoint before adding curvature, avoiding artificial jumps from large absolute
position offsets. At native source ``b1aa772503070347e8ac0b263c8de9278d390c17``,
the existing cubic suite passed **6 cases / 3,033 assertions**, including its new
precise-time regression; NWayRange passed **5 cases / 1,019 assertions**, with
unchanged tolerances. The rebuilt kernel is 67,159,336 bytes, SHA-256
``10f8307ba799847f11d9614d2c3f3da3db23962c46a739d73f536f9dc31ae5d5``.
`Build provenance <../../../laser-ranging-work/logs/stable-cubic-kernel-provenance.json>`_
links the focused reports and preserved old kernel. The earlier six-suite
checkpoint remains historical evidence; this rebuild did not repeat that whole
set. Subsequent observation runs must retain their own kernel identities.

Selected range clock and precise rotation bindings
-------------------------------------------------

Native source ``7d52716f6723c704b20ca5f833570acd0f0cc0ba`` shares the selected
Earth-station rotation clock with non-TDB one-way/n-way range normalization,
and exposes the four existing precise-``Time`` rotation/derivative methods in
Python while retaining the float methods. The incremental build at
``/home/dominic/Tudat/laser-ranging-build`` produced kernel SHA-256
``80b924287707124a8efea66455b26233b5511b7856983d35c05d21f82734f215``
(67,172,672 bytes), imported from ``src/tudatpy/kernel.so`` in that build.
NWayRange passed **5 cases / 1,023 assertions** with the selected Desai--Sibois
and four-point EOP context; all **3 environment Python tests passed**. The native
test-only frame alignment is ``a17e7bb44``; the final Python test is ``26ed42495``.
The preceding unrelated native evidence remains retained.

The new binding test resolves station motion at epochs 20 ns apart that share
one double value, exercises all four Time methods and their inverse identities,
and retains the float route. Its additional 20 ns finite-difference comparison
with the analytic derivative was removed after an independent probe traced its
0.12 micrometre discrepancy to inherited SOFA Earth-rotation-angle double
arithmetic. UT1 retained the 20 ns interval; the ERA increment had a 267 ps
equivalent quantization error at this epoch. No production accuracy setting or
existing numerical tolerance was weakened.

`Kernel provenance <../../../laser-ranging-work/logs/selected-clock-time-rotation-kernel-provenance.json>`_
records the actual interpreter/import origins, source and focused Python result.
The `build log <../../../laser-ranging-work/logs/selected-clock-time-rotation-build.log>`_,
`native result <../../../laser-ranging-work/logs/selected-clock-nway-native-framefix-result.json>`_
and `Python log <../../../laser-ranging-work/logs/selected-clock-environment-python.log>`_
retain the focused checks. The
`SOFA probe <../../../laser-ranging-work/verification/probe_earth_rotation_20ns.py>`_
and `quantization evidence <../../../laser-ranging-work/logs/earth-rotation-20ns-probe.json>`_
retain the measured floor and the original failed assertion separately. The
preceding kernel and its provenance are preserved externally under
``laser-ranging-work/binaries/10f8307ba799847f11d9614d2c3f3da3db23962c46a739d73f536f9dc31ae5d5/``.

The laser target signature correction was later removed: its only use (the lunar reflector glass path)
is applied as a constant range bias in the LLR driver. Annual DE440 residuals agree with the saved results
to within light-time convergence noise (median 2 micrometres, maximum 0.26 mm).
