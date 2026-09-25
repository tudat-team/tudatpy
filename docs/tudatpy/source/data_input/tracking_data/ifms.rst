.. _tracking_data_ifms:

``ifms``
========

.. automodule:: tudatpy.data_input.tracking_data.ifms

IFMS files contain radio tracking data, typically closed-loop Doppler data, for
a number of ESA deep space missions, optionally with station tropospheric
corrections that can be applied during loading.
This submodule contains functionality to load Doppler tracking data from Level 2 IFMS files.
The file format is described in `IFMS Doppler Processing Software : Level 1a to Level 2, Table 3-2 and 3-3 <https://archives.esac.esa.int/psa/ftp/MARS-EXPRESS/MRS/MEX-M-MRS-1-2-3-EXT9-4441-V1.0/DOCUMENT/MRS_DOC/MEX_MRS_IGM_DS_3035.PDF>`_.

The :func:`read_ifms_data` function is the main interface for loading the data and
converting it to objects that Tudat can process further; see also
:ref:`tracking_data`.

.. currentmodule:: tudatpy.data_input.tracking_data.ifms

.. autofunction:: read_ifms_data
