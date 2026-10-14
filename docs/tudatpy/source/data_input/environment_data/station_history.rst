.. _station_history:

``station_history``
===================

This input-only module reads fixed-width SINEX 2.00--2.02 station products,
including IGN SSC catalogues, discontinuity histories, PSD parameters and
seasonal coefficients. It creates neutral records, without propagating stations,
joining independently supplied histories or applying geophysical corrections.

.. code-block:: python

   from tudatpy.data_input.environment_data import station_history

   history = station_history.read_sinex_station_history("ITRF2020-IGS-TRF.SSC")
   key = station_history.SinexSolutionKey("AB07", "A", "1")
   solution = station_history.get_cartesian_solution(history, key)

A Cartesian query requires a complete position triplet. Velocity components are
optional and remain ``None`` when absent; ``require_complete_velocity=True``
requires all three. Positions are metres and velocities are metres per second,
using 365.25 days per source ``m/y`` year. Original values, units, standard
deviations, parameter indices, constraint codes and source lines remain available.
Parameter-only files do not require a Cartesian solution.

Epochs retain their civil UTC year, day of year, seconds of day and original
``YY:DDD:SSSSS`` token. SINEX uses 2050 for year 50 and 1951 for year 51.
A possible seconds-of-day 86400 label is retained; this input layer does not
consult a leap-second table or collapse that label to a scalar UTC epoch.
``00:000:00000`` remains explicitly unspecified, a header-start/end default,
or an ITRF validity open bound according to its field. It is never turned into
a calendar date or infinity. Source observation spans are separate from position
and velocity discontinuity intervals; equal solution identifiers do not imply
that these independent histories should be joined.

``history.sites`` is a list, preserving repeated raw site/point identifiers and
all their DOMES/technique/source records. Combined IGN seasonal products contain
site/point codes for different monuments. Cartesian selection rejects ambiguous
monument metadata; it does not infer a DOMES association for seasonal coefficients.
Identical repeated metadata stays in the raw list and can be recognized as one
identity only during an explicit Cartesian query.

IGN PSD files label relaxation parameters ``TEXP_*`` and ``TLOG_*`` with the
literal unit ``m`` although the IGN model routines use years. These records
retain the original fields and report ``unresolved_extension`` with no SI value.
Repeated same-type/event terms also remain separate parameter-index records.
A declared provider model must interpret their pairing and units; the reader
does not do so. Frame, release, origin and tide-system metadata stay empty when
the source does not explicitly supply them. These four metadata fields can be
assigned explicitly by the caller; the reader never infers them from a filename.

Receiver and antenna histories retain their intervals, model names, serial
numbers, firmware and complete source records. Eccentricities retain the source
``UNE`` or ``XYZ`` component order. Unknown blocks and covariance declarations
remain opaque records with source/line provenance and are not evaluated.
Provider extensions after the standard equipment fields, such as the IGS
antenna azimuth column, remain available in the original source line.

``read_sinex_station_history_text`` applies the file reader's same validation
to an in-memory string and accepts a diagnostic source label. Malformed fields,
nonfinite values, reversed intervals, duplicate Cartesian components, estimate
counts inconsistent with the header, and truncated blocks raise exceptions with
source and line information. The declared parameter count remains source
metadata when a caller explicitly edits C++ preloaded records. This reader uses
only standard C++17 and is compiled in the independently linkable
``tudat_environment_data`` target. The target also provides the other environment
readers and uses Eigen, without dynamics, estimation or SOFA dependencies.

The ``ilrs.read_sinex_station_data`` map represents one solution
per DOMES monument and rejects ambiguous multi-solution histories. Use the full
history reader and an explicit solution key in that case. Both readers require
a SINEX header, fixed-width records, matching block delimiters and ``%ENDSNX``.

.. currentmodule:: tudatpy.data_input.environment_data.station_history

.. autofunction:: read_sinex_station_history

.. autofunction:: read_sinex_station_history_text

.. autofunction:: validate_sinex_station_history

.. autofunction:: get_cartesian_solution_keys

.. autofunction:: get_cartesian_solution

.. autoclass:: SinexStationHistory
   :members:

.. autoclass:: SinexSolutionKey
   :members:

.. autoclass:: SinexParameterRecord
   :members:

.. autoclass:: SinexEpoch
   :members:
