.. _data_retrieval_media_corrections_vmf:

``vmf``
=======

.. automodule:: tudatpy.data_input.data_retrieval.media_corrections.vmf

This submodule contains functionality to download Vienna Mapping Function
tropospheric mapping files. The :func:`download_vmf` function is the main
interface for retrieving VMF products for a UTC date range. Geodetic technique
and processing-mode selection are controlled with :class:`VmfTechnique` and
:class:`VmfProcessing`.

.. currentmodule:: tudatpy.data_input.data_retrieval.media_corrections.vmf

Functions
---------

.. autosummary::

   download_vmf

.. autofunction:: download_vmf

Classes
-------

.. autosummary::

   VmfTechnique
   VmfProcessing

.. autoclass:: VmfTechnique
   :members:

.. autoclass:: VmfProcessing
   :members:
