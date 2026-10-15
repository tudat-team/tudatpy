.. _data_retrieval_media_corrections:

``media_corrections``
=====================

.. toctree::
   :maxdepth: 2
   :caption: Modules

   data_retrieval_media_corrections/ionex
   data_retrieval_media_corrections/vmf

.. automodule:: tudatpy.data_input.data_retrieval.media_corrections

This submodule contains download utilities for ancillary media-correction data,
including IONEX ionospheric maps and VMF tropospheric mapping functions. The
:func:`download_ancillary` function is the main convenience interface for
downloading both data sources for a UTC date range. Source-specific interfaces
for IONEX and VMF products are documented in the submodules listed above. The
:class:`DownloadAtmosphericData` class provides the legacy atmospheric-data
downloader interface in this media-correction namespace.

.. currentmodule:: tudatpy.data_input.data_retrieval.media_corrections

Functions
---------

.. autosummary::

   download_ancillary

.. autofunction:: download_ancillary

Classes
-------

.. autosummary::

   DownloadAtmosphericData
   DownloadResult

.. autoclass:: DownloadAtmosphericData
   :members:

.. autoclass:: DownloadResult
   :members:

Exceptions
----------

.. autosummary::

   AncillaryDownloadError
   AuthenticationError

.. autoclass:: AncillaryDownloadError

.. autoclass:: AuthenticationError
