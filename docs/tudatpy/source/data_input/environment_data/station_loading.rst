.. _station_loading_data:

``station_loading``
===================

Load precomputed station-displacement series while retaining their scientific
inventory and source geometry. The module reads
`EPHEDISP 2005.06.30 <https://astrogeo.org/malo/ephedisp_format.txt>`_ and accepts
generic preloaded ``StationLoadingData``. It performs no station matching,
interpolation, clock conversion or reference-origin translation.

.. code-block:: python

   from tudatpy.data_input.environment_data import station_loading

   metadata = station_loading.imls_merra2_atmospheric_loading_metadata(
       release="my documented provider release")
   data = station_loading.read_ephedisp_station_loading_data(
       "provider.eph", metadata=metadata)

EPHEDISP's authoritative clock is TAI. Its first split epoch (integer MJD plus
fractional seconds) and sample interval define the indexed epoch grid. The
informational calendar fields are retained as text and ignored. Provider
Cartesian locations and validity radii use metres; displacement components
use metres in geocentric radial/east/north order, defined at those source
locations. Identifiers are internal labels and do not establish station identity.

Each site's sample dictionary uses zero-based indexes into ``data.epochs``.
Leading/trailing missing samples and sites without any samples are preserved;
internal gaps are rejected. The reader validates fixed fields, counted sections,
unique sites/samples, sampling, chronological order, finite values and the
trailer. Both the specification's T-before-A layout and the pinned IMLS file's
A-before-T layout are supported. A singleton can declare a zero sample interval.
Source lines, source records, headers and informational fields remain available.

Metadata retains product/release, source profile, physical components, source
origin and its definition, removed harmonics, retained seasonal/mean/trend content,
ocean-response assumptions and reference frame. Unknown fields stay unknown.
``imls_merra2_atmospheric_loading_metadata`` explicitly selects the documented
IMLS atmospheric profile and requires a release identifier. Filenames and
labels never select that profile. The profile uses whole-Earth CM including
atmosphere; it establishes no identity between CE and CF.

Other providers can supply a strictly increasing explicit epoch grid, clock,
Cartesian or geocentric radial/east/north basis, metadata, locations/radii and
samples, then call ``data.validate()``. A nominal sampling interval is optional
for an irregular grid. Source epochs remain split; they are never rounded through
a large Julian date. Declared EPHEDISP begin/end epochs remain separate from its
index-derived grid, preserving the original fields' stated precision.

The separately linkable C++ target ``tudat_station_loading_data`` and header
``tudat/io/readStationLoadingData.h`` require only standard C++17.

.. currentmodule:: tudatpy.data_input.environment_data.station_loading

.. autoclass:: StationLoadingData
   :members:

.. autoclass:: StationLoadingMetadata
   :members:

.. autoclass:: StationLoadingEpoch
   :members:

.. autoclass:: StationLoadingSite
   :members:

.. autoclass:: StationLoadingSample
   :members:

.. autofunction:: read_ephedisp_station_loading_data

.. autofunction:: read_ephedisp_station_loading_data_text

.. autofunction:: imls_merra2_atmospheric_loading_metadata
