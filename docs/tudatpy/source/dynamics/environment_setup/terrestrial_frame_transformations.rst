.. _terrestrial_frame_transformations:

``terrestrial_frame_transformations``
====================================

Transform Cartesian positions and moving states between terrestrial-frame
realizations using the published linearized Helmert equations. Load parameters
through :ref:`terrestrial_frame_data`, then construct a model from the neutral
records. Select the original IGN ITRF2014 or superseding ITRF2020 catalogue
explicitly; each retains its published reference hub and source provenance.

The evaluation epoch supplies both the parameter decimal year and its derivative
with respect to the input state's seconds. The following example deliberately
uses a uniform 365.25-day year coordinate. It does not convert a calendar date,
UTC, TT or TDB into a year:

.. code-block:: python

   from tudatpy.data_input.environment_data import terrestrial_frames as data
   from tudatpy.dynamics.environment_setup import terrestrial_frame_transformations as tf

   model = tf.TerrestrialFrameTransformationModel(data.official_itrf2020_parameters())
   epoch = tf.TerrestrialFrameEvaluationEpoch(2025.0, 1.0 / 31557600.0)
   source_state = [4075578.385, 931852.890, 4801570.154, 0.012345, -0.023456, 0.034567]
   target_state = model.transform_state(source_state, "ITRF2020", "ITRF2014", epoch)

For native time evaluation, supply the model's ``clock`` callback. It accepts a
``Time`` and returns a ``TerrestrialFrameEvaluationEpoch`` with the same physical
second convention as the input velocity. The ``*_at_time`` methods require this
callback; no clock mapping is inferred from a scalar epoch or a catalogue label.

The affine matrix includes scale and linearized signed rotations. Position and
velocity are evaluated at the requested epoch:

.. math::

   y=T+A x, \qquad \dot y=\dot T+\dot A x+A\dot x.

``inverse`` uses the actual matrix inverse and its derivative. ``compose_after``
includes all translation and matrix derivative terms. Model composition follows
the declared catalogue hub. Equal nonempty realization labels return identity
before catalogue lookup, including unknown identical labels; this does not
validate catalogue membership. Aligned catalogue-release updates retain their
realization identity and do not create additional datum transformations.

The existing :func:`~tudatpy.dynamics.environment_setup.convert_ground_station_state_between_itrf_frames`
retains the IGN ITRF2014 catalogue and its documented uniform-Julian-year epoch
convention. It delegates to the same requested-epoch affine state evaluation.

.. automodule:: tudatpy.dynamics.environment_setup.terrestrial_frame_transformations
   :members:
   :undoc-members:
