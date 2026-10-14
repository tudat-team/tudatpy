.. _ground_station:

``ground_station``
==================
This module contains a set of factory functions for setting up ground
stations and associated models. Note that in Tudat, no distinction is
made between a ground station/lander on Earth or a different body.

The main interfaces with Tudat is the :attr:`~tudatpy.dynamics.environment_setup.BodySettings.ground_station_settings` list
attribute (with entries of type :class:`~tudatpy.dynamics.environment_setup.ground_station.GroundStationSettings`) of the body settings, which defines settings for the ground stations of a body.
The functions in this submodule are used to create these settings objects. When creating a body (typically using the
:func:`~tudatpy.dynamics.environment_setup.create_system_of_bodies` function), a list of objects of type
:class:`~tudatpy.dynamics.environment.GroundStation` (or a derived class) is created
and added to the associated :class:`~tudatpy.dynamics.environment.Body` object based on the settings object, which can
be retrieved using the :attr:`~tudatpy.dynamics.environment.Body.ground_station_list` attribute. The ground station
can also be safely added to the :class:`~tudatpy.dynamics.environment.Body` after its creation using the :func:`~tudatpy.dynamics.environment_setup.add_ground_station` function.

A ground station defines a reference point (and other relevant properties)
on a celestial body. Although ground stations are considered part of the 
environment in Tudat (as properties of a :class:`~tudatpy.dynamics.environment.Body` object), they do not 
influence the numerical propagation (unless a custom model imposing this 
is implemented by the user).

Properties such as a transmitting frequency of the station is added after its creation (to the :class:`~tudatpy.dynamics.environment.GroundStation` object). When creating the settings for a ground stations (using the functions in this module)
the body may be endowed with a list of :class:`~tudatpy.dynamics.environment_setup.ground_station.GroundStationMotionSettings`
settings, which define station dependent models to define the station-specific motion of a ground station (such as plate motion).
Models for the deformation of the full body (such as tidal shape variations) are to be defined through the :ref:`shape_deformation`
module.

For Earth, we provide several options to create default stations, such as the :func:`~tudatpy.dynamics.environment_setup.ground_station.dsn_stations` and
:func:`~tudatpy.dynamics.environment_setup.ground_station.evn_stations`.

Station histories, post-seismic deformation and seasonal motion
--------------------------------------------------------------

``ign_station_history`` builds an absolute terrestrial station trajectory from
preloaded :mod:`tudatpy.data_input.environment_data.station_history` products.
Supply explicit and matching ``frame``, ``release``, ``origin`` (for example CM
or CF) and ``tide_system`` metadata on each product. The reader preserves source
records and does not infer these scientific conventions from filenames.
``ITRF2020`` and its aligned catalogue updates have separate release identities;
mixing solution, PSD or seasonal releases raises an error.

Select a site code and point code explicitly. Every selected Cartesian solution
requires complete position and velocity coefficients and one matching position
discontinuity interval. Physical TAI intervals are half-open: the new solution
applies exactly at its start. Gaps, overlaps, unresolved default endpoints and
extrapolation raise errors. Observational spans do not define model validity.
The independent velocity discontinuity ID is retained in
``trajectory.solutions[i].velocity_solution_id`` and ``velocity_validity``.
Published SSC velocities remain attached to their Cartesian/position solution;
the velocity discontinuity ID is not an SSC solution lookup. AB07 position
solution 6, for example, belongs to velocity interval 4. Each position interval
must fit one velocity interval. Coefficients sharing a velocity interval must
agree within the explicit componentwise
``velocity_consistency_tolerance_metres_per_julian_year`` (default ``1e-8``, or
10 nm/year). This is a numerical equivalence budget, not printed-field precision;
the supplied coefficients are preserved. Set zero for strict equality.

PSD amplitudes retain their original parameter indices. The IGN profile pairs
an amplitude at index ``i`` with the matching relaxation parameter at ``i+1``;
record order is irrelevant. Every applicable event and every repeated EXP/LOG
term accumulates. It interprets the native IGN ``TEXP``/``TLOG`` numbers as Julian
years even when the published unit field says ``m``; raw records remain unchanged.
Relaxation times must be positive and finite. At an event the correction is zero
and the returned velocity is the right-hand rate. At a solution change the new
solution supplies position and right-hand velocity. Neither value asserts a
two-sided derivative at a discontinuity.

ENU coefficients require an explicit, constant, proper ``enu_to_terrestrial``
rotation whose columns are East, North and Up. The model does not recompute that
basis as the station moves. Duplicate SITE/ID entries with distinct DOMES
monuments require a ``domes_selector``; the selector resolves metadata, and the
caller must establish the parameter-to-monument association from companion
products. Metadata-only selection cannot recover information absent from a
parameter-only source.

The named native time convention is ``IGN nominal UTC``. Linear and PSD elapsed
time is the difference of nominal civil MJD divided by 365.25. Annual and
semiannual phase uses the elapsed fraction of the actual calendar year,
including February 29; it is not a fixed 365.25-day sinusoid. ``StationHistoryTime``
keeps physical TAI, split MJD and calendar-year coordinates separate and accepts
optional derivatives with respect to environment TDB seconds. ``position`` does
not need rates; ``state`` returns metres and metres per TDB second and requires
the corresponding source rates. On ordinary UTC intervals these are
``dTAI/dTDB / 86400`` and ``dTAI/dTDB / (days_in_year * 86400)``.

The native conversion uses the installed SOFA leap table and supports UTC from
1972 onward. Dubious dates, inserted-second source/reference/event labels and
state derivatives exactly at a UTC offset step raise errors. Physical interval
conversion can represent a real inserted-second label. At an ordinary January 1
the seasonal phase is continuous and its rate uses the new year's length.
Caller-supplied source coordinates or a ``source_time_function(TDB_Time,
require_rates)`` callback may implement another explicitly named convention;
its name must agree with the trajectory metadata.

Attach the result with ``station_history_motion(trajectory)``. It subtracts the
station's nominal position so that the resulting station position equals the
absolute trajectory. The default callback uses the body's terrestrial time
converter when available. For example, with matching metadata already supplied:

.. code-block:: python

    from tudatpy.data_input.environment_data import station_history as sh
    from tudatpy.dynamics.environment_setup import ground_station as gs

    trajectory = gs.ign_station_history(
        coordinates, discontinuities, sh.SinexSiteKey("AB07", "A"),
        post_seismic=psd, enu_to_terrestrial=enu_basis)
    motion = gs.station_history_motion(trajectory)
    station = gs.basic_station("AB07", nominal_position,
                               station_motion_settings=[motion])

``trajectory.metadata`` retains the release, origin, source inventory and signal
inventory; ``trajectory.has_fitted_seasonal_terms`` reports whether empirical
seasonal terms are active. When also applying physical loading, explicitly assess
the overlap between fitted seasonal signals and the loading correction.

Validation includes offline outputs from the unchanged hash-pinned IGN
``parametric.f`` for selectors 0–5, signed amplitudes, very short and very long
elapsed times, cumulative events, analytic rate checks, leap/year boundaries,
physical interval selection and published repeated/seasonal coefficient inputs.
The fixture provenance is in ``tests/test_tudat/reference_files/station_history``.
No independent complete trajectory comparison is claimed.

The direct C++ ``LinearStationMotionModel`` constructor now defaults its reference
epoch to zero TDB seconds since J2000, matching the existing settings/factory
default. Its previous default was the Julian day number used as seconds.

Transforming a station history with a terrestrial frame catalogue
---------------------------------------------------------------

``transform_station_history`` applies a preloaded conventional or supplied
linearized Helmert catalogue to the complete absolute trajectory. Its calendar
decimal year and derivative come from the same ``StationHistoryTime`` as the
marker coordinates; ``station_history_frame_epoch`` exposes that epoch for
explicit affine model queries. Parameter rates still use their declared Julian
year units. An IGN nominal UTC leap-year query uses the actual calendar-year
fraction and its TDB derivative, including the documented right-hand year rate.

Supply a :class:`~tudatpy.data_input.environment_data.terrestrial_frames.TerrestrialFrameIdentity`
with both target realization and release, plus explicit target origin and tide
system. Origin and tide conventions must match the source. This operation does
not convert CM to CF or change a tide system. A release different from the
catalogue's declared identity requires ``acknowledge_aligned_releases=True``;
that declaration warrants alignment of the release to the named realization,
not the existence of an additional transform. In particular, an aligned
ITRF2020-u2024 catalogue can retain the ITRF2020 realization. Same-frame identity
preserves origin/tide semantics and requires an explicit alignment declaration
when the release label changes. Distinct realizations require actual catalogue
records.

The returned trajectory uses full affine kinematics:
``y = T + A*x`` and ``dy/dTDB = dT/dTDB + dA/dTDB*x + A*dx/dTDB``. Inverses use
the actual matrix inverse, including scale. Position evaluation remains available
when source-coordinate rates are absent; state evaluation requires the rates
needed by the source and, for a nonidentity transformation, the calendar-year
rate. No derivative is invented to obtain a position.

.. code-block:: python

    from tudatpy.data_input.environment_data import terrestrial_frames as tf
    from tudatpy.dynamics.environment_setup import ground_station as gs

    target = tf.TerrestrialFrameIdentity()
    target.realization = target.release = "ITRF2014"
    transformed = gs.transform_station_history(
        trajectory, tf.official_itrf2020_parameters(), target,
        target_origin=trajectory.metadata.origin,
        target_tide_system=trajectory.metadata.tide_system)
    motion = gs.station_history_motion(transformed)

``transformed.metadata`` describes the evaluated output. The input identity,
release, signal inventory and original coefficients remain available through
``transformed.source_trajectory``; inherited solution/PSD/seasonal records also
remain source-coordinate records. ``transformed.catalogue`` retains the pinned
parameter source and hash, and ``acknowledged_release_alignment`` records the
explicit release declaration. Empirical seasonal activity is propagated for
subsequent loading-overlap decisions.













Functions
---------
.. currentmodule:: tudatpy.dynamics.environment_setup.ground_station

.. autosummary::

    basic_station

    dsn_station

    dsn_stations

    radio_telescope_stations

    evn_stations

    default_ilrs_sinex_state_file

    default_ilrs_sinex_eccentricity_file

    ilrs_stations_from_sinex_domes

    linear_station_motion

    piecewise_constant_station_motion

    custom_station_motion

    get_approximate_dsn_ground_station_positions

    get_vlbi_station_positions

    get_vlbi_station_velocities

    get_radio_telescope_positions

    optical_telescope_stations

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.basic_station

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.dsn_station

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.dsn_stations

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.radio_telescope_stations

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.evn_stations

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.default_ilrs_sinex_state_file

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.default_ilrs_sinex_eccentricity_file

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.ilrs_stations_from_sinex_domes

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.linear_station_motion

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.piecewise_constant_station_motion

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.custom_station_motion

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.get_approximate_dsn_ground_station_positions

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.get_vlbi_station_positions

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.get_vlbi_station_velocities

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.get_radio_telescope_positions

.. autofunction:: tudatpy.dynamics.environment_setup.ground_station.optical_telescope_stations






Classes
-------
.. currentmodule:: tudatpy.dynamics.environment_setup.ground_station

.. autosummary::

   GroundStationSettings

   GroundStationMotionSettings

   LinearGroundStationMotionSettings

   PiecewiseConstantGroundStationMotionSettings

   CustomGroundStationMotionSettings



.. autoclass:: tudatpy.dynamics.environment_setup.ground_station.GroundStationSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.ground_station.GroundStationMotionSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.ground_station.LinearGroundStationMotionSettings
   :members:
   :special-members: __init__

.. autoclass:: tudatpy.dynamics.environment_setup.ground_station.PiecewiseConstantGroundStationMotionSettings
   :members:

.. autoclass:: tudatpy.dynamics.environment_setup.ground_station.CustomGroundStationMotionSettings
   :members:

