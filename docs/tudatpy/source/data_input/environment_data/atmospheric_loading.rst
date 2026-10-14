.. _atmospheric_loading_data:

``atmospheric_loading``
=======================

Read conventional atmospheric S1/S2 loading coefficients before constructing
an Earth deformation model. The Bernese grid format has longitude, latitude,
and twelve signed coefficients: radial/north/east rows, each ordered
cos(S1), sin(S1), cos(S2), sin(S2). ATL station blocks contain the same twelve
coefficients following an exact station identifier (including numeric IDs).
Readers normalize millimetres to metres and retain original units, source
headers, model, Green function and origin declarations. Both longitude seam
endpoints remain in the raw grid. They perform no interpolation or synthesis.

.. code-block:: python

   from tudatpy.data_input.environment_data import atmospheric_loading as atl

   metadata = atl.AtmosphericTideMetadata()
   metadata.expected_longitude_count = 361
   metadata.expected_latitude_count = 181
   grid = atl.read_atmospheric_tide_grid("s1_s2_def_ce.dat", metadata)
   grid.validate()

The supported grid has no count header; expected counts additionally detect
truncation of complete terminal longitude blocks. Malformed, nonfinite,
incomplete, duplicate or inconsistent data raise with source and line context.

The provider and Bernese manual establish the component labels and units, but
do not settle the grid latitude convention or whether radial means spherical
radial versus ellipsoid-normal up for automatic station use. These declarations
remain independent and unspecified by default. A dynamics model requires an
explicit lookup latitude convention and displacement basis. Selecting geodetic
lookup does not establish geodetic displacement basis. Geodetic conversion also
requires Earth's configured shape. Precomputed station coefficients require an
explicit displacement basis, with no automatic geographic lookup.

Native deformation settings are available through
``environment_setup.shape_deformation.atmospheric_tide_from_grid`` (and the
station/file variants). The environment scalar time argument is TDB seconds
since J2000; phase is measured from midnight UT1 using the model-owned shared
Earth time converter. Context position evaluation remains continuous through
inserted leap seconds. State/velocity are analytic loading derivatives with
respect to TDB. They chain dUT1/dTAI from the selected continuous EOP interval
and harmonic corrections with the model's geocentric dTT/dTDB clock. Values at
knots use the interval to the right (left at the final sample). Unsupported or
missing rate stencils raise for rates/states while independently valid position
evaluation remains available.

Source origin is preserved by default. The opt-in RP03 CE-to-CM policy requires
explicit identification of the IERS RP03 no-IB Farrell source profile, declared
model/Green function and CE or CM origin. It adds the IERS Table 7.6 Cartesian
translation after local rotation exactly once for CE inputs; already-CM data
receive none. Unknown products, including the separate Vienna APL grid, cannot
use this conversion.

The pinned RP03 provider has rounded seam and polar coefficient differences.
Strict endpoint consistency is the default. The explicitly named averaged-seam
and polar-Cartesian approximation records their maximum discrepancies and
modifies only model coefficients. The original grid remains available.
The Bernese manual station sample does not agree with this downloaded grid;
that discrepancy remains documented in the test provenance and is not claimed
as independent-software conformance.

.. automodule:: tudatpy.data_input.environment_data.atmospheric_loading
   :members:
