.. _release_notes:

Release notes
=============

Version 1.1.0
-------------

Highlights of this release
~~~~~~~~~~~~~~~~~~~~~~~~~~


* Full two-body spherical-harmonic gravity, including figure-figure effects,
  provides consistent translational and rotational models for two extended bodies.
  This allows the high-fidelity orbital-rotational dynamics of for instance binary
  asteroid or close-in moons to be modelled. Stay tuned for the consistent incorporation
  of tidal effects into this dynamics for the next release!
* New readers, converters, and download tools support DSN ATDF tracking data and
  optical-navigation pixel/line data (tested with Voyager data). The DSN TNF and
  MPC readers have also been improved. This facilitates analyses using data coming from
  even more sources/missions. See our `EPSC abstract and poster
  <https://meetingorganizer.copernicus.org/EPSC2026/EPSC2026-349.html>`_
* Orbit determination now supports camera pixel/line, azimuth/elevation,
  differenced time-of-arrival, differenced frequency-of-arrival, and measured
  one-way-frequency observables. This has further broadened the data types that
  can be used, permitting even broader application and data fusion applications
* A post-Newtonian space-time framework provides a common basis for translational
  dynamics and coordinate- and proper-time propagation. This allow for high-fidelity
  solar system dynamics and tests of relativity
* The Mars Climate Database interface has been updated for high-fidelity Martian
  atmospheric density calculations. This is essential in high-accuracy orbit
  estimation of Martian orbiters
* Constrained multi-arc estimation can connect consecutive arcs through soft
  constraints. This allows for a flexible way to balance between the advantages
  of a single- and multi-arc estimation
* Supported simulation and estimation objects can be serialized to binary or
  JSON files, or with Python pickle. This can make life much easier by simply
  saving all analysis results to a file, and then delegating all the post-processing to a
  separate script - with no manual file reading/writing!


New features - Description
~~~~~~~~~~~~~~~~~~~~~~~~~~


Extended-body, radiation, and custom force modelling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Full two-body spherical-harmonic dynamics account for figure-figure interactions
between the spherical-harmonic fields of two bodies. Unlike standard spherical-
harmonic gravity, both the acceleration and torque depend on the gravity fields
and orientations of both bodies.

* The `acceleration.full_two_body
  <https://py.api.tudat.space/en/latest/dynamics/propagation_setup/acceleration.html#tudatpy.dynamics.propagation_setup.acceleration.full_two_body_spherical_harmonic_gravity>`_
  and `torque.full_two_body
  <https://py.api.tudat.space/en/latest/dynamics/propagation_setup/torque.html#tudatpy.dynamics.propagation_setup.torque.full_two_body_spherical_harmonic_gravitational>`_
  models implement the full two-body formulation of `Dirkx et al. (2019)
  <https://doi.org/10.1007/s10509-019-3521-4>`_. Users can select truncation
  degrees or specify the coefficient combinations directly.
* A closed-form degree-two by degree-two torque model implements the classical
  formulation of `Schutz et al. (1981)
  <https://doi.org/10.1007/BF01229194>`_.

Other dynamics additions include:

* The `acceleration.three_coefficient
  <https://py.api.tudat.space/en/latest/dynamics/propagation_setup/acceleration.html#tudatpy.dynamics.propagation_setup.acceleration.three_coefficient_radiation_pressure>`_
  model follows `McMahon and Scheeres (2015)
  <https://doi.org/10.2514/1.G000666>`_ and provides an intermediate-fidelity
  radiation-pressure model whose coefficients can be estimated.
* `torque.radiation_pressure_torque
  <https://py.api.tudat.space/en/latest/dynamics/propagation_setup/torque.html#tudatpy.dynamics.propagation_setup.torque.radiation_pressure_torque>`_
  computes radiation-pressure torque for panelled bodies defined by a macromodel.
  Additional dependent variables provide panel geometry, illumination, irradiance,
  and force information.
* The `acceleration.rtg
  <https://py.api.tudat.space/en/latest/dynamics/propagation_setup/acceleration.html#tudatpy.dynamics.propagation_setup.acceleration.rtg>`_
  model represents acceleration due to a spacecraft RTG, with constant or
  arc-wise parameter estimation for single- and multi-arc problems.
* Custom-atmosphere outputs and estimatable panel material properties are also
  available.


Mission, tracking, and catalogue data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* New or improved processing is available for DSN ATDF and TNF, FDETS,
  optical-navigation PSF, and ILRS CRD/SINEX data. The ODF, IFMS, ramp-table,
  and Doppler-compression workflows have also been extended.
* MPC support includes a separate parser for 80-column observations,
  standardized table handling, and custom target names.
* Space-Track queries use more selective API calls and configurable timeouts.
  OMM records without embedded TLE fields are handled correctly, and TLE parsing
  is more robust.
* `environment.fit_tle
  <https://py.api.tudat.space/en/latest/dynamics/environment.html#tudatpy.dynamics.environment.fit_tle_to_cartesian_state_history>`_
  fits a TLE to a Cartesian state history, allowing TLEs to be generated from
  propagated or tabulated states.
* SBDB objects expose minimum-orbit-intersection-distance properties. Failures
  from the remote Horizons, DISCO, and Space-Track services are now reported more
  consistently.


Observation simulation and tracking analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The observation framework now covers additional radio, optical, laser-ranging,
and navigation-camera applications.

* The `model.pixel_coordinates
  <https://py.api.tudat.space/en/latest/estimation/observable_models_setup/model_settings.html#tudatpy.estimation.observable_models_setup.model_settings.pixel_coordinates>`_
  model allows optical-navigation pixel/line data to be used in estimation. It
  has been tested with Voyager 2 data at Neptune.
* Additional observation models are available for `model.azimuth_elevation
  <https://py.api.tudat.space/en/latest/estimation/observable_models_setup/model_settings.html#tudatpy.estimation.observable_models_setup.model_settings.azimuth_elevation>`_,
  `model.dtoa
  <https://py.api.tudat.space/en/latest/estimation/observable_models_setup/model_settings.html#tudatpy.estimation.observable_models_setup.model_settings.differenced_time_of_arrival>`_,
  `model.dfoa
  <https://py.api.tudat.space/en/latest/estimation/observable_models_setup/model_settings.html#tudatpy.estimation.observable_models_setup.model_settings.differenced_frequency_of_arrival>`_,
  and `model.measured_frequency
  <https://py.api.tudat.space/en/latest/estimation/observable_models_setup/model_settings.html#tudatpy.estimation.observable_models_setup.model_settings.doppler_measured_frequency>`_
  observables.
* Frequency-dependent corrections can now be used with all observation models.
  Individual light-time correction contributions can be saved as observation
  dependent variables.
* `light_time.nequick2
  <https://py.api.tudat.space/en/latest/estimation/observable_models_setup/light_time_corrections.html#tudatpy.estimation.observable_models_setup.light_time_corrections.nequick2_ionospheric_light_time_correction>`_
  implements the ionospheric model of `Nava et al. (2008)
  <https://doi.org/10.1016/j.jastp.2008.01.015>`_. IONEX data can rescale the
  model and supply uncertainty information.
* Initial ILRS laser-ranging support covers CRD observation parsing, SINEX
  station states and eccentricities, and VMF3o tropospheric corrections.
* DSN additions include custom seasonal troposphere calibrations, improved
  ramp-table handling, N-way range partials, and transponder delays.
* Standalone geometry functions check whether a target is visible from a station
  and whether the station is dark.


Relativistic dynamics, time propagation, and estimation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The new `space_time
<https://py.api.tudat.space/en/latest/dynamics/environment_setup/space_time.html>`_
module defines a common post-Newtonian model for relativistic dynamics and time
transformations. It is intended primarily for Solar-System dynamics and tests of
relativity.

* Space-time can be configured with a Schwarzschild metric around one central
  body or with the multi-body Solar-System metric of
  `Soffel et al. (2003) <https://doi.org/10.1086/378162>`_. The models include
  first- and second-order post-Newtonian terms, PPN parameters, and optional
  spherical-harmonic, angular-momentum, and body-acceleration contributions.
* Time differences between barycentric coordinate time (such as TCB/TDB),
  body-centred coordinate time (such as TCG/TT), and proper time can be
  propagated with `propagator.first_order_bodycentric_relativistic_time_settings
  <https://py.api.tudat.space/en/latest/dynamics/propagation_setup/propagator.html#tudatpy.dynamics.propagation_setup.propagator.first_order_bodycentric_relativistic_time_settings>`_,
  `propagator.bodycentered_to_topocentric_time_settings
  <https://py.api.tudat.space/en/latest/dynamics/propagation_setup/propagator.html#tudatpy.dynamics.propagation_setup.propagator.bodycentered_to_topocentric_time_settings>`_, and
  `propagator.direct_relativistic_time_settings
  <https://py.api.tudat.space/en/latest/dynamics/propagation_setup/propagator.html#tudatpy.dynamics.propagation_setup.propagator.direct_relativistic_time_settings>`_.
  The equations can be integrated separately or together with the translational
  dynamics. Relativistic time converters can also be attached to bodies, so that
  the dynamics, clocks, and observation models use the same transformations.
* `acceleration.relativistic_from_metric
  <https://py.api.tudat.space/en/latest/dynamics/propagation_setup/acceleration.html#tudatpy.dynamics.propagation_setup.acceleration.relativistic_from_metric>`_
  evaluates the translational acceleration directly from the configured
  space-time metric.



Atmospheres and environmental modelling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Mars Climate Database is available as an optional source of atmospheric
density data for high-fidelity Mars orbiter simulations.
See `atmosphere.mcd
<https://py.api.tudat.space/en/latest/dynamics/environment_setup/atmosphere.html#tudatpy.dynamics.environment_setup.atmosphere.mars_climate_database_climate_model>`_
and `Forget et al. (1999) <https://doi.org/10.1029/1999JE001025>`_. Tudat does not
ship MCD; it must be compiled separately by the user.

A new cometary-coma framework implements the four-dimensional spherical-harmonic
polynomial model of `Reichel et al. (2026)
<https://doi.org/10.5194/epsc2026-370>`_ for density, temperature, and wind fields.

* `atmosphere.coma_density
  <https://py.api.tudat.space/en/latest/dynamics/environment_setup/atmosphere.html#tudatpy.dynamics.environment_setup.atmosphere.coma_model_from_stokes_data>`_
  constructs density and temperature models from polynomial or precomputed Stokes
  spherical-harmonic data.
* `atmosphere.coma_wind
  <https://py.api.tudat.space/en/latest/dynamics/environment_setup/atmosphere.html#tudatpy.dynamics.environment_setup.atmosphere.coma_wind_model>`_
  constructs spherical-harmonic wind fields with optional nucleus co-rotation.

Other environment additions include IAU rotation models with configurable
reference epochs and base frames; estimatable pole and pole-libration parameters;
subtraction of mean tidal forcing from solid-body tides; ocean- and pole-tide
deformation settings for Earth-based ground stations; optical- and radio-telescope
station data; and configurable station motion.


Orbit determination and covariance analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Constrained multi-arc estimation applies optional soft continuity constraints
between neighbouring arcs, following `Lari et al. (2021)
<https://doi.org/10.1007/s10686-021-09823-8>`_. The `estimation.continuity
<https://py.api.tudat.space/en/latest/estimation/estimation_analysis.html#tudatpy.estimation.estimation_analysis.position_only_continuity>`_
functions provide position, velocity, full-state, and general weighted constraints.

Other estimation additions include:

* Observation-only estimation can be run without propagator settings when the
  estimated parameters do not require integrated dynamics.
* Estimation ephemerides accept user-defined interpolator settings, for cases
  where the default sixth-order Lagrange interpolation is not suitable.
* New estimatable parameters cover arc-wise aerodynamic force scaling,
  exponential-atmosphere base density and scale height, rotation-model pole and
  libration properties, RTG force properties, and panel material coefficients.


Serialization of simulation and estimation results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Selected Tudat objects can now be stored in binary or JSON files and serialized
with Python pickle. Supported objects include propagation results, observation
collections and sets, and estimation outputs. Serialized data contains schema
and build metadata, which are checked when the data is loaded. Compatibility
across arbitrary TudatPy versions or build configurations is not guaranteed.


Some smaller additions
^^^^^^^^^^^^^^^^^^^^^^

* Date/time conversion handles Python ``datetime`` edge cases more consistently.
  Historical Earth-orientation conversion also accepts pre-1972 UT1 inputs.
* Python interfaces are available for direct TLE and multi-arc ephemeris settings,
  empty tabulated ephemerides, and custom-state retrieval.
* Ground-station setup exposes VLBI station positions and velocities, as well as
  combined DSN/VLBI position catalogues. ILRS station settings can be loaded from
  SINEX state and eccentricity files.


Main improvements to existing features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Propagation and environment modelling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Spherical-harmonic acceleration is now regularized at the poles.
* NRLMSISE-00 storm-condition ``Ap`` inputs and warnings for user overrides have
  been corrected.
* The state derivative can be modified during propagation, and termination
  callbacks can inspect the current propagated state.
* Multi-arc and hybrid-arc propagation have additional consistency checks and
  corrected result interfaces.


Observations and estimation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Residuals for periodic angular observables are wrapped correctly at
  discontinuities.
* Fixes cover Doppler compression, Doppler bias partials, IFMS cadence-gap
  processing, occultation checks, and several light-time correction paths.
* Estimation checks propagator and parameter consistency earlier and gives clearer
  errors when setup fails.


Data access
^^^^^^^^^^^

* Network-backed data clients use timeouts and recognize common remote-service
  outages, avoiding indefinite waits and misleading local errors.
* The motion of DSN station DSS-65 is represented piecewise to account for its
  relocation in 2005.


Deprecated features
~~~~~~~~~~~~~~~~~~~

* ``ionex_ionospheric_light_time_correction`` remains available for compatibility,
  but is documented as deprecated in favour of
  ``nequick2_ionospheric_light_time_correction``.
* Observation ancillary APIs now use the correct ``ancillary`` spelling. The old
  ``ancilliary`` names remain as deprecated aliases and issue warnings.


Backwards-incompatible changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* The deprecated radiation-pressure interface API has been removed. Radiation
  pressure must now be configured using radiation source and target models in
  ``tudatpy.dynamics.environment_setup.radiation_pressure``. This includes the
  removal of the legacy cannonball interface settings and the
  ``add_radiation_pressure_interface`` helper.
* The deprecated thrust direction and magnitude settings API has been removed.
  Configure thrust using engine and body rotation model settings together with
  ``thrust_from_engine``, ``thrust_from_engines``, or
  ``thrust_from_all_engines``.
* The deprecated aerodynamic guidance API has been removed, including
  ``AerodynamicGuidance``, ``set_aerodynamic_guidance``,
  ``set_aerodynamic_orientation_functions``, and
  ``set_constant_aerodynamic_orientation``.
* Deprecated generic ``custom`` factory aliases in the ephemeris, acceleration,
  dependent-variable, mass-rate, and torque setup modules have been removed.
  Use the corresponding explicitly named ``custom_*`` factories.
* Other deprecated setup and analysis conveniences have been removed, including
  ``spherical_harmonic_triaxial_body``, ``single_observation_set``, and the
  ``*_split_output`` covariance and formal-error propagation functions.
* Deprecated names containing the misspelling ``ancilliary`` have been removed.
  Use the corresponding ``ancillary`` names.
* The deprecated ``tudatpy.io`` module has been removed. Import the corresponding
  functionality from ``tudatpy.data``.
* ``tudatpy.util.split_history`` has been removed.


Build, packaging, and documentation changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* The CMake build now supports compiler caching, unity builds, and precompiled
  headers.
* Several Boost utilities have been replaced by standard-library equivalents,
  reducing Boost to a header-only dependency.
* Generated API signatures and type information have been cleaned up. Module-level
  documentation has been expanded, and formats such as aerodynamic coefficient
  files are now documented.


Contributors
~~~~~~~~~~~~

The following GitHub users contributed to this release.

* `alfonsoSR <https://github.com/alfonsoSR>`_
* `DominicDirkx <https://github.com/DominicDirkx>`_
* `HannoSpreeuw <https://github.com/HannoSpreeuw>`_
* `jo11he <https://github.com/jo11he>`_
* `kimonito98 <https://github.com/kimonito98>`_
* `larshinueber <https://github.com/larshinueber>`_
* `luigigisolfi <https://github.com/luigigisolfi>`_
* `Markus-Reichel <https://github.com/Markus-Reichel>`_
* `mrevellino <https://github.com/mrevellino>`_
* `nikemais <https://github.com/nikemais>`_
* `PascalChristiaanse <https://github.com/PascalChristiaanse>`_
* `Spinola000 <https://github.com/Spinola000>`_
* `Strtus <https://github.com/Strtus>`_
* `suvayu <https://github.com/suvayu>`_
* `valeriof7 <https://github.com/valeriof7>`_
* `yiquintero <https://github.com/yiquintero>`_

