.. _shape_deformation:

``shape_deformation``
=====================
This module contains a set of factory functions for setting up the
shape deformation models of celestial bodies in an environment.

The main interfaces with Tudat is the :attr:`~tudatpy.dynamics.environment_setup.BodySettings.shape_deformation_settings` list
attribute  (with entries of type :class:`~tudatpy.dynamics.environment_setup.shape_deformation.BodyShapeSettings`) of the body settings, which defines settings for the global shape deformation model of a body.
**The functions in this submodule are used to create these settings objects.** When creating a body (typically using the
:func:`~tudatpy.dynamics.environment_setup.create_system_of_bodies` function), a list of objects of type
:class:`~tudatpy.dynamics.environment.BodyShapeDeformationModel` (or a derived class) is created
and added to the associated :class:`~tudatpy.dynamics.environment.Body` object based on the settings object, which can
be retrieved using the :attr:`~tudatpy.dynamics.environment.Body.shape_deformation_model` attribute.

The shape deformation models are used for global shape variations of a celestial body, which are used (for instance)
for the high-accuracy modelling of the positions of ground stations on a body. Station-specific models for position variations
of a station (such as plate motion) are set in the ground stations themselves (see :ref:`ground_station`).


Station motion receives TDB seconds since J2000. Native Earth deformation models
route that physical epoch explicitly: solid-tide astronomical arguments and
HARDISP loading use TT and the corresponding UTC phase, ephemerides use TDB,
and pole tides use the model's precise TT/EOP path. A GCRS/ITRS rotation uses
its own converter and declared input scale. Generic models and explicitly
supplied custom argument callbacks retain their original scalar forwarding.

Native C++ deformation models expose ``calculateDisplacement(Time, station)``
and ``calculateDisplacementState(Time, station)``. Both receive precise TDB;
the model owns its selected time converter. Outputs are body-fixed Cartesian
metres and, for states, metres per second with respect to TDB. Position evaluation
does not require rates. Velocity callers extract the state tail; ENU callers
apply the station's body-fixed-to-topocentric rotation. TT-only models apply the
geocentric dTT/dTDB factor, and mixed solid tides differentiate the complete
routed displacement. Existing explicit scalar argument callbacks receive TDB
seconds and retain their own argument convention.
A configured TT/TDB interpolator without a supported clock derivative permits
positions but rejects native states/rates. Current scalar UTC tide phases reject
ambiguous inserted-second positions and derivative intervals crossing an offset
step; continuous UT1 pole-tide positions remain available through a leap second.

The solid-tide displacement is equivalent to DEHANTTIDEINEL Steps 1 and 2 and
includes the time-independent degree-2 displacement. Nominal coordinates must
therefore be conventional tide-free coordinates. No Step 3 permanent-tide
"uncorrection" is applied. A mean-tide catalogue cannot receive this full
addition unchanged. This is the displacement convention in
`IERS Conventions 2010, Chapter 7, section 7.1.1.2
<https://iers-conventions.obspm.fr/content/chapter7/icc7.pdf>`_; it does not set a
gravity-field tide system.

Functions
---------
.. currentmodule:: tudatpy.dynamics.environment_setup.shape_deformation

.. autosummary::

   basic_solid_body_tidal

   degree_two_basic_solid_body_tidal

   iers_2010_solid_body_tidal

   ocean_tidal

   ocean_tidal_from_data

   pole_tidal
.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.basic_solid_body_tidal

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.degree_two_basic_solid_body_tidal

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.iers_2010_solid_body_tidal

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.ocean_tidal

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.ocean_tidal_from_data

.. autoclass:: tudatpy.dynamics.environment_setup.shape_deformation.OceanTideSynthesis
   :members:

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.pole_tidal





Conventional ocean pole-tide loading implements IERS Section 7.1.5 with the
2018 secular pole and IERS 2010 Chapter 1 constants. The official grid covers
latitudes -89.75 through +89.75 degrees. ``require_grid_coverage`` raises outside
the supplied extent; ``clamp_to_nearest_grid_latitude_approximation`` explicitly
uses the nearest covered latitude. No polar-cap interpolation is specified by
the IERS coefficient grid. Longitude is periodic for a complete grid. Cartesian
displacement uses the geocentric spherical basis; public station ENU uses its
actual topocentric frame. Standalone epochs are TT; station environment calls convert their TDB epoch
through the rotation model's own clock and preserve precise EOP coverage.

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.ocean_pole_tidal

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.ocean_pole_tidal_from_data

.. autoclass:: tudatpy.dynamics.environment_setup.shape_deformation.OceanPoleTideLatitudePolicy
   :members:

Atmospheric S1/S2 loading accepts preloaded grid or station coefficients and
file-backed settings. Lookup latitude and displacement basis are declared
separately; geodetic choices require the configured body shape. The model uses
continuous UT1 from its own time converter, including through inserted UTC leap
seconds. Loading velocities differentiate the harmonics analytically and apply
the selected dUT1/dTAI and geocentric dTT/dTDB factors to produce metres per TDB
second. At daily EOP knots, rates use the interval to the right, except at the
final sample where the left interval applies. Missing or unsupported rate
stencils reject state/rate requests while independently valid positions remain
available. An explicitly selected RP03 CE-to-CM translation and its derivative
are added once in Cartesian coordinates. See
:doc:`/data_input/environment_data/atmospheric_loading` for source, geometry and
origin policies and the documented limits of the reference comparisons.

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.atmospheric_tide_from_grid

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.atmospheric_tide_from_stations

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.atmospheric_tide_from_grid_file

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.atmospheric_tide_from_station_file

Precomputed station loading
---------------------------

``station_loading(data, evaluation_settings)`` consumes neutral series loaded
with :mod:`tudatpy.data_input.environment_data.station_loading`. It supports
atmospheric, hydrological and ocean loading as explicitly inventoried source
components. It adds no implicit loading model to the environment.

.. code-block:: python

   from tudatpy.data_input.environment_data import station_loading
   from tudatpy.dynamics.environment_setup import shape_deformation

   metadata = station_loading.imls_merra2_atmospheric_loading_metadata(
       release="my documented provider release")
   metadata.reference_frame = "ITRS"  # Explicit caller declaration of Cartesian axes.
   data = station_loading.read_ephedisp_station_loading_data("provider.eph", metadata)
   policy = shape_deformation.StationLoadingEvaluationSettings()
   policy.station_bindings = {"MyStation": "weuc"}
   policy.output_origin_definition = data.metadata.origin_definition
   settings = shape_deformation.station_loading(data, policy)
   body_settings.get("Earth").shape_deformation_settings = [settings]

Bindings map environment station IDs to arbitrary source labels and still
enforce each provider point's validity radius. Alternatively, explicitly select
``within_provider_radius``; zero matches are unavailable and multiple matches
raise a configuration error. Zero-radius points cover only their exact source
location. Source radial/east/north axes are geocentric axes at the provider
Cartesian point, independent of the evaluated station's geodetic topocentric
axes. Displacements and their rates are returned in Cartesian body-fixed axes.
The factory requires the declared Cartesian frame to match the body's fixed
frame; it performs no implicit transformation between realizations.

Source TAI, TT, UTC or TDB epochs are converted with the model-owned clock to
continuous TAI before interpolation. EPHEDISP supplies TAI. The supported
policy is explicit linear interpolation, with a default maximum support gap of
86400 SI seconds. This is not a claim of equivalence to provider B-spline
interpolation. A rate differentiates the selected interpolant and applies
``dTAI/dTDB`` once. At an exact knot it uses the interval to the right, except
at the final endpoint where it uses the interval to the left. No extrapolation,
edge order reduction or silent gap bridging occurs. Exact source-node positions
remain available beside unsupported gaps, and singleton positions remain
available when their rates are unavailable.

Native C++ calls use the common ``Time`` argument in environment TDB and a
``GroundStationState`` for displacement or combined state. Source epochs remain
split and precise. Native queries classify source-node TDB images exactly,
without rounding or tolerance snapping. Preloaded source TDB images remain the
original precise epochs, and uniform TAI sampling need not be uniform in source
TDB. Only genuine scalar station adapters
recognize exact rounded node images; collapsed scalar images raise an error.
Continuous TAI interpolation works through an inserted UTC second.

The requested output-origin definition must match the declared source definition
or have a supplied ``StationLoadingOriginTranslation``. Its Cartesian series
adds ``source origin - target origin`` once after source-basis rotation, including
its derivative. Both origin definitions, Cartesian axes and the complete known
physical component/signal scope must match. CE and CF are never equated
implicitly. Reusing the same global translation identifier across selected
contributions raises an error even when overlaps are explicitly permitted;
separate component-scoped translations may each apply once.

Unknown product/release/component or signal inventory requires explicit
``allow_unknown_inventory`` or resolved source metadata. Axes and origins still
must be declared. Assemblies reject repeated physical components, retained tidal
bands combined with another tidal-loading model, and physical seasonal loading
combined with fitted station seasonal terms. ``allow_declared_overlap`` records
an explicit decision to combine these contributions. Scientific conflicts
propagate independently of optional coverage: optional configuration may skip
an unavailable central displacement, while unavailable derivatives of a selected
displacement still raise.

.. autofunction:: tudatpy.dynamics.environment_setup.shape_deformation.station_loading

.. autoclass:: tudatpy.dynamics.environment_setup.shape_deformation.StationLoadingEvaluationSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.shape_deformation.StationLoadingOriginTranslation
   :members:

Classes
-------
.. currentmodule:: tudatpy.dynamics.environment_setup.shape_deformation

.. autosummary::

   BodyDeformationSettings

   BasicSolidBodyDeformationSettings



.. autoclass:: tudatpy.dynamics.environment_setup.shape_deformation.BodyDeformationSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.shape_deformation.BasicSolidBodyDeformationSettings
   :members:
