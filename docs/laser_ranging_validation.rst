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
   * - LaserTargetSignature
     - 4
     - 21
     - Retained
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
Its binding smoke exercised scalar and array target signatures, optical and
target correction factories, body-deformation station motion, actual SP3/CRD
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
