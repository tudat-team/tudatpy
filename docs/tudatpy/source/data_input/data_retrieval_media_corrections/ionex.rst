.. _data_retrieval_media_corrections_ionex:

``ionex``
=========

.. automodule:: tudatpy.data_input.data_retrieval.media_corrections.ionex

This submodule contains functionality to download IONEX global ionospheric TEC
maps from NASA CDDIS. The :func:`download_ionex` function is the main interface
for retrieving IONEX files for a UTC date range. Product and temporal-resolution
selection are controlled with :class:`IonexProduct` and :class:`IonexResolution`.

.. currentmodule:: tudatpy.data_input.data_retrieval.media_corrections.ionex

Functions
---------

.. autosummary::

   download_ionex

.. autofunction:: download_ionex

Classes
-------

.. autosummary::

   IonexProduct
   IonexResolution

.. autoclass:: IonexProduct
   :members:

.. autoclass:: IonexResolution
   :members:
