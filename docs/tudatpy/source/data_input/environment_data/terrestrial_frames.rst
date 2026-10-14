.. _terrestrial_frame_data:

``terrestrial_frames``
======================

Load provenance-bearing terrestrial-frame Helmert parameters before constructing
a transformation model. The strict IGN table reader retains reference hub,
source/target realizations and release provenance, parameter decimal-year epoch,
original unit declarations, row/rate source lines, original records and source
notices. Translation, rotation, scale and rates are normalized to SI using the
published Julian-year rate unit (365.25 days). This unit conversion does not
establish a physical-to-decimal-year time mapping.

.. code-block:: python

   from tudatpy.data_input.environment_data import terrestrial_frames as tf

   data = tf.official_itrf2020_parameters()
   print(data.reference_frame.realization)  # ITRF2020
   print(data.get_parameters("ITRF2014").reference_decimal_year)  # 2015.0
   data.validate()

File and text readers are also available. Units, column order, paired rates and
the explicit published position-vector equation are checked; malformed records,
nonfinite data, duplicates or truncation raise with source and line context.
Pinned 2014 and superseding 2020 catalogues retain their distinct original hashes.
Release metadata is provenance rather than an automatically different frame
datum; an aligned ITRF2020 catalogue update does not imply a nonzero Helmert shift.

The dynamics API is :ref:`terrestrial_frame_transformations`. Its structured epoch
contains the published decimal year and dDecimalYear/dInputSecond. The model
accepts preloaded data and an explicitly supplied time callback. No scalar epoch
is guessed to be UTC, TT, TDB or calendar decimal year. Model evaluation performs
no file I/O. State velocity follows the callback's input seconds:

.. math::

   y=T+A x, \qquad \dot y=\dot T+\dot A x+A\dot x.

The inverse uses the actual inverse affine matrix and exact derivative, and
composition follows the declared catalogue hub. The linearized matrix includes
scale and is not an orthogonal rotation. Rotation sign conventions must be
explicit. The native default cannot silently select an older table or datum.

The existing ``convert_ground_station_state_between_itrf_frames`` compatibility
function retains the explicitly pinned 2014 table and its existing uniform
Julian-year seconds since J2000: year=2000+epoch/31557600. This is not a native
TDB-to-calendar-year conversion. Its corrected requested-epoch position and
velocity replace the former propagation through a fixed reference epoch.
Same nonempty realization labels give identity without catalogue lookup,
including unknown identical labels. This behavior does not validate catalogue
membership; distinct labels require the declared hub and known records.

.. automodule:: tudatpy.data_input.environment_data.terrestrial_frames
   :members:
