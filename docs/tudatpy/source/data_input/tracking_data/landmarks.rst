.. _tracking_data_landmarks:

``landmarks``
=============

.. automodule:: tudatpy.data_input.tracking_data.landmarks

This submodule contains functionality to load optical landmark observations
from SUM and LMK files, as written by the stereophotoclinometry (SPC) software
(see `Gaskell et al., 2008 <https://doi.org/10.1111/j.1945-5100.2008.tb00692.x>`_).
Each SUM file describes one image, with its epoch, camera calibration and
pointing, and the pixel coordinates of the landmarks seen in it. Each LMK file
describes one landmark and its position on the target body. These files can be
converted to pixel-coordinate tracking data, and to supplementary data that
adds the landmarks as ground stations on the target body, and the cameras, with
the calibration and pointing of each image, on the receiver body. The
:func:`read_sum_lmk_data` function is the main interface for loading the data
and converting it to objects that Tudat can process further; see also
:ref:`tracking_data`. All other functionality in this module is reserved for
better understanding what data is being loaded, and in some cases manipulating
it, before it is processed into Tudat-compatible objects.

.. currentmodule:: tudatpy.data_input.tracking_data.landmarks

.. autofunction:: read_sum_lmk_data

Supporting API
--------------

The functions and classes below expose the contents of SUM and LMK files in
dedicated containers. They are used internally by :func:`read_sum_lmk_data` and
can be used to inspect the files, or to convert the same files several times
with :func:`convert_raw_sum_lmk_data`.

.. autofunction:: read_raw_sum_files

.. autofunction:: read_raw_lmk_files

.. autofunction:: convert_raw_sum_lmk_data

.. autoclass:: SumImageData
   :members:

.. autoclass:: SumLandmarkObservation
   :members:

.. autoclass:: LmkLandmarkData
   :members:
