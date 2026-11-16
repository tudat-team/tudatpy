.. _tracking_data_psf:

``psf``
=======

.. automodule:: tudatpy.data_input.tracking_data.psf

This submodule contains functionality to load tracking data from PSF files of
the type distributed by the JPL Solar System Dynamics group as spacecraft
optical observations of planetary satellites
(https://ssd.jpl.nasa.gov/sats/obs_data.html). These files contain optical
image measurements, camera metadata, and image orientation information that can
be converted to pixel-line tracking data. The :func:`read_psf_data` function is
the main interface for loading the data and converting it to objects that Tudat
can process further; see also
:ref:`tracking_data`. All other functionality in this module is
reserved for better understanding what data is being loaded, and in some cases
manipulating it, before it is processed into Tudat-compatible objects.

.. currentmodule:: tudatpy.data_input.tracking_data.psf

.. autofunction:: read_psf_data

Supporting API
--------------

The functions and classes below expose the raw contents of PSF files in
dedicated containers. They are used internally by :func:`read_psf_data` and can
be used to inspect image records, measurements and camera metadata before they
are converted to Tudat pixel-line tracking data.

.. autofunction:: read_raw_psf_file_contents

.. autoclass:: RawPsfFileContents
   :members:

.. autoclass:: RawPsfFileImageContents
   :members:

.. autoclass:: RawPsfMeasurement
   :members:

.. autoclass:: RawPsfStarMeasurement
   :members:

.. autoclass:: OpticalImageType
