.. _tracking_data:

``tracking_data``
=================

This module contains the common containers used by the individual tracking-data
readers in the submodules below. Those readers load file- or service-specific
tracking data into :class:`~tudatpy.data_input.tracking_data.TrackingData`
objects, optionally with
:class:`~tudatpy.data_input.tracking_data.TrackingSupplementaryData` objects
containing auxiliary information such as station frequency ramps or camera
settings. The tracking data can then be converted to an
:class:`~tudatpy.estimation.observations.ObservationCollection` with
:func:`~tudatpy.estimation.observations.create_observation_collection_from_tracking_data`.
Supplementary data that updates bodies, ground stations, or related environment
objects is applied with
:func:`~tudatpy.estimation.observations.set_tracking_supplementary_data_in_bodies`.

.. toctree::
   :maxdepth: 2
   :caption: Modules

   tracking_data/atdf
   tracking_data/obs_80_cols
   tracking_data/fdets
   tracking_data/generic_text_file
   tracking_data/ifms
   tracking_data/mpc
   tracking_data/odf
   tracking_data/optical_utilities
   tracking_data/psf
   tracking_data/tnf

.. automodule:: tudatpy.data_input.tracking_data
   :members:

Tracking data containers
------------------------
.. currentmodule:: tudatpy.data_input.tracking_data

These classes hold tracking data after it has been read from a source, but
before it is used in estimation. :class:`TrackingData` stores the observations
themselves, while :class:`TrackingSupplementaryData` stores additional data from
the same source that must be applied to the simulation environment.

.. autosummary::

   TrackingData
   TrackingSupplementaryData

.. autoclass:: tudatpy.data_input.tracking_data.TrackingData
   :members:

.. autoclass:: tudatpy.data_input.tracking_data.TrackingSupplementaryData
   :members:

TrackingSupplementaryData constituent classes
------------------------------------------------

The classes in this section store the individual categories of supplementary
data that may be attached to a :class:`TrackingSupplementaryData` object. They
are normally created by the tracking-data readers, but can also be inspected
before the supplementary data is applied to the bodies.

.. autosummary::

   FrequencySupplementaryData
   FrequencyRamp
   RampedFrequencySupplementaryData
   PiecewiseConstantFrequencySupplementaryData
   InstrumentSupplementaryData
   CameraInstrumentSupplementaryData
   TranslationalStateSupplementaryData
   RotationalStateSupplementaryData

.. autoclass:: tudatpy.data_input.tracking_data.FrequencySupplementaryData
   :members:

.. autoclass:: tudatpy.data_input.tracking_data.FrequencyRamp
   :members:

.. autoclass:: tudatpy.data_input.tracking_data.RampedFrequencySupplementaryData
   :members:

.. autoclass:: tudatpy.data_input.tracking_data.PiecewiseConstantFrequencySupplementaryData
   :members:

.. autoclass:: tudatpy.data_input.tracking_data.InstrumentSupplementaryData
   :members:

.. autoclass:: tudatpy.data_input.tracking_data.CameraInstrumentSupplementaryData
   :members:

.. autoclass:: tudatpy.data_input.tracking_data.TranslationalStateSupplementaryData
   :members:

.. autoclass:: tudatpy.data_input.tracking_data.RotationalStateSupplementaryData
   :members:
