=========
Changelog
=========

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/) and this project adheres to [Semantic Versioning](http://semver.org/).

`0.2.0` (unreleased)
---------------------

* Support for [2. Perturbed Earth-orbiting Satellite](http://tudat.tudelft.nl/tutorials/applicationWalkthroughs/perturbedEarthOrbitingSatellite.html)

`0.1.0`
-------

:Date: Apr 29, 2020

Added
-----

* Support for [1. Unperterbed Earth-orbiting Satellite](http://tudat.tudelft.nl/tutorials/applicationWalkthroughs/unperturbedEarthOrbitingSatellite.html) added. The Python version of the script included as
  ``tutorial_1.py`` under ``examples``.
* Added module ``tudatpy.constants``.

  * Contains all constants as named in ``celestialBodyConstants.h``, ``physicalConstants.h``,
    and ``mathematicalConstants.h`` in Tudat.

* Added module ``tudatpy.interpolators``.

  * Currently only acts as a placeholder for the ``LangrangeInterpolatorSettings`` needed for the default argument of
    ``InterpolatedSpiceEphemerisSettings``.

* Added module ``tudatpy.spice_interface``.

  * Exposed function ``spice_interface.load_standard_spice_kernels()``.
  * Exposed function ``spice_interface.clear_spice_kernels()``.

* Added module ``tudatpy.basic_astrodynamics``.

  * Exposed enumeration ``AvailableAcceleration``.
  * Partial exposure of 3D specialised ``AccelerationModel`` (declaration only).

* Added module ``tudatpy.gravitation``.

  * Partially exposure of class ``GravityFieldModel`` for ``gravitational_parameter`` property.

* Added module ``tudatpy.numerical_integrators``.

  * Exposed enumeration ``AvailableIntegrators`` for ``IntegratorSettings`` argument.
  * Partial exposure of class ``IntegratorSettings`` (only constructor).

* Added module ``tudatpy.propagators``.

  * Exposed function ``propagators.get_single_integration_size()``
  * Exposed function ``propagators.get_single_integration_differential_equation_order()``
  * Exposed function ``propagators.get_generalized_acceleration_size()``
  * Partial exposure of ``SingleStateTypeDerivative`` (only declaration).
  * Partial exposure of ``NBodyStateDerivative`` (only declaration).
  * Partial exposure of ``NBodyCowellStateDerivative`` (only constructor).
  * Exposed class ``SingleArcDynamicsSimulator`` (complete).
  * Exposed enumeration ``TranslationalPropagatorType`` (complete).
  * Partial exposure of class ``DependentVariableSaveSettings`` (only constructor).
  * Partial exposure of class ``PropagatorSettings`` (only declaration).
  * Partial exposure of class ``SingleArcPropagatorSettings`` (only declaration).
  * Partial exposure of class ``TranslationalStatePropagatorSettings`` (only all constructors).


* Added module ``tudatpy.orbital_element_conversions``.

  * Exposed enumeration ``KeplerianElementIndices``.
  * Exposed function ``tudatpy.convert_keplerian_to_cartesian_elements()``.

* Added module ``tudatpy.simulation_setup``.

  * Exposed class ``BodySettings`` (complete).
  * Exposed class ``Body`` (complete).
  * Exposed function ``tudatpy.get_default_body_settings()`` (both overloads).
  * Partially exposed class ``Ephemeris`` (only declaration).
  * Exposed class ``ConstantEphemeris`` (complete).
  * Exposed enumeration ``EphemerisType``.
  * Exposed class ``EphemerisSettings`` (complete).
  * Exposed class ``DirectSpiceEphemerisSettings`` (complete).
  * Exposed class ``InterpolatedSpiceEphemerisSettings`` (complete).
  * Exposed class ``ApproximatePlanetPositionSettings`` (complete).
  * Exposed class ``ConstantEphemerisSettings`` (complete).
  * Exposed class ``CustomEphemerisSettings`` (complete).
  * Exposed class ``KeplerEphemerisSettings`` (complete).
  * Exposed class ``TabulatedEphemerisSettings`` (complete).
  * Exposed function ``create_body_ephemeris()``.
  * Exposed function ``get_safe_interpolation_interval()``.
  * Exposed function ``set_global_frame_body_ephemerides()``.
  * Exposed function ``create_bodies()``.
  * Exposed function ``create_acceleration_models_dict()`` (both overloads).

Fixed
-----

* Fixed hard requirement for Python 3.7
* Modified installation procedure for ``python3-dev`` requirment prior to compilation of the tudatBundle using the
  ``USE_TUDATPY`` option.


Other Changes
--------------