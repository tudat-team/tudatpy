``ocean_loading``
=================

Read BLQ ocean tidal loading coefficients independently of dynamics. The data
records preserve up/radial, west and south amplitudes in metres, Greenwich phase
lags in degrees, the standard eleven-constituent order, station identifiers and
provider headers. Explicit CMC declarations identify CE or CM; absent declarations
remain unspecified. These records do not convert the supplied reference origin.

.. currentmodule:: tudatpy.data_input.environment_data.ocean_loading

.. autofunction:: read_ocean_loading_data

.. autofunction:: read_ocean_loading_data_text

.. autofunction:: constituent_names

.. autoclass:: OceanLoadingData
   :members:

.. autoclass:: OceanLoadingStation
   :members:

.. autoclass:: OceanLoadingReferenceOrigin
   :members:

To construct a loading model from these records, use
:func:`~tudatpy.dynamics.environment_setup.shape_deformation.ocean_tidal_from_data`.
The default expands the eleven coefficients to 342 constituents using the
corrected IERS HARDISP admittance model. The explicit
``OceanTideSynthesis.legacy_eleven_constituents`` option retains the previous
synthesis for reproducibility. File input occurs before displacement evaluation.
