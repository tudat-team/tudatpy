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

Ocean pole-load grids
--------------------

The IERS conventional ocean pole-load grid is available from the
`IERS Conventions centre <https://iers-conventions.obspm.fr/content/chapter7/additional_info/opoleloadcoefcmcor.txt.gz>`_.
Decompress it locally and pass the plain-text filename to ``read_ocean_pole_tide_grid``.
This input API performs no interpolation, astronomical evaluation or reference-origin conversion.
The neutral record retains source headers and the declared origin; its dimensionless
coefficient matrices have radial/north/east rows and real/imaginary columns. Nodes
are stored longitude-fast, then latitude. Longitude and latitude axes retain the
provider degrees without inferring a geodetic or geocentric coordinate convention.

The official grid already includes the correction to the instantaneous centre of
mass of the Earth and its loading ocean pole tide. The dynamics model therefore
applies no further geocentre correction.

.. autoclass:: OceanPoleTideGrid
   :members:

.. autofunction:: read_ocean_pole_tide_grid

.. autofunction:: read_ocean_pole_tide_grid_text
