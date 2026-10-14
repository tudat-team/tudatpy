.. _tracking_data_sum_lmk:

``sum_lmk``
===========

.. automodule:: tudatpy.data_input.tracking_data.sum_lmk

This submodule contains functionality to load optical landmark observations
from SUM and LMK files, as written by the stereophotoclinometry (SPC) software.
Each SUM file describes one image, with its epoch, camera calibration and
pointing, and the pixel coordinates of the landmarks seen in it. Each LMK file
describes one landmark and its position on the target body. The
:func:`read_sum_lmk_data` function is the main interface for loading the data
and converting it to objects that Tudat can process further. See also
:ref:`tracking_data`. All other functionality in this module is reserved for
better understanding what data is being loaded before it is processed into
Tudat-compatible objects.

Repeated measurements of the same landmark in the same image are loaded once
when their pixel coordinates agree. Conflicting coordinates, image metadata,
or landmark definitions raise an error. This also applies when converting raw
data directly with :func:`convert_raw_sum_lmk_data`.

.. currentmodule:: tudatpy.data_input.tracking_data.sum_lmk

.. autofunction:: read_sum_lmk_data

Supporting API
--------------

The functions and classes below expose the raw contents of SUM and LMK files.
They are used internally by :func:`read_sum_lmk_data` and can be used to
inspect the files, or to convert the same files several times with
:func:`convert_raw_sum_lmk_data`.

.. autofunction:: convert_raw_sum_lmk_data

.. autofunction:: read_raw_sum_files

.. autofunction:: read_raw_lmk_files

.. autoclass:: SumImageData
   :members:

.. autoclass:: SumLandmarkObservation
   :members:

.. autoclass:: LmkLandmarkData
   :members:
