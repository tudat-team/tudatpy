=========
Changelog
=========

All notable changes to this project will be documented in this file.

The format is based on `Keep a Changelog`_ and this project adheres to `Semantic Versioning`_.

.. _`Keep a Changelog` : http://keepachangelog.com/
.. _`Semantic Versioning` : http://semver.org/

0.3.0 (02-05-2020)
------------------

* Support for `3. Un-guided Capsule Entry`_ added.

* Implemented module :ref:`tudatpy.unittests`

    * Exposed function ``unittests.get_apollo_coefficient_interface()``

* Exposed method ``Body.get_flight_conditions()``

* Added module :ref:`tudatpy.aerodynamics`.

    * Exposed class ``FlightConditions()``
    * Exposed method ``FlightConditions.get_aerodynamic_angle_calculator()``

* Added module :ref:`tudatpy.reference_frames`.

    * Exposed class ``AerodynamicAngleCalculator``.

* Exposed function :py:meth:`tudatpy.orbital_element_conversions.convert_spherical_orbital_to_cartesian_state`

* Added module :ref:`tudatpy.ephemerides`.

    * Added function :py:meth:`tudatpy.ephemerides.transform_state_to_global_frame()`

* Added class :py:meth:`tudatpy.propagators.SingleDependentVariableSaveSettings`.
* Added enum :py:meth:`tudatpy.propagators.PropagationDependentVariables`.
* Added class :py:meth:`tudatpy.propagators.PropagationDependentVariableTerminationSettings`.

.. _`3. Un-guided Capsule Entry` : http://tudat.tudelft.nl/tutorials/applicationWalkthroughs/unguidedCapsuleEntry.html

0.2.1 (01-05-2020)
------------------

* Implemented shortened ``AccelerationSettings`` assignment for ``aerodynamic``, ``spherical_harmonic_gravity``
  and ``cannonball_radiation_pressure`` using methods calls of Python class ``simulation_setup.Acceleration``.

0.2.0 (30-04-2020)
------------------

Added
~~~~~

* Support for `2. Perturbed Earth-orbiting Satellite`_ added.

    * Exposed method ``BodySettings::ephemerisSettings()``
    * Exposed method ``BodySettings::rotationModelSettings()``
    * Exposed method ``Body::setConstantBodyMass()``
    * Exposed method ``Body::setAerodynamicCoefficientInterface()``
    * Exposed method ``Body::setRadiationPressureInterface()``
    * Exposed class ``ConstantAerodynamicCoefficientSettings()``
    * Exposed class ``SphericalHarmonicAccelerationSettings()``
    * Exposed function ``createAerodynamicCoefficientInterface()``
    * Exposed function ``createRadiationPressureInterface()``


.. _`2. Perturbed Earth-orbiting Satellite` : http://tudat.tudelft.nl/tutorials/applicationWalkthroughs/perturbedEarthOrbitingSatellite.html



0.1.0 (29-04-2020)
------------------

Added
~~~~~

* Support for `1. Unperturbed Earth-orbiting Satellite`_ added. The Python version of the script included as
  ``tutorial_1.py`` under ``examples``.

.. _`1. Unperturbed Earth-orbiting Satellite` : http://tudat.tudelft.nl/tutorials/applicationWalkthroughs/unperturbedEarthOrbitingSatellite.html

* Added module :ref:`tudatpy.constants`.

  * Contains all constants as named in ``celestialBodyConstants.h``, ``physicalConstants.h``,
    and ``mathematicalConstants.h`` in Tudat.

* Added module :ref:`tudatpy.interpolators`.

  * Currently only acts as a placeholder for the ``LangrangeInterpolatorSettings`` needed for the default argument of
    ``InterpolatedSpiceEphemerisSettings``.

* Added module :ref:`tudatpy.spice_interface`.

  * Exposed function ``spice_interface.load_standard_spice_kernels()``.
  * Exposed function ``spice_interface.clear_spice_kernels()``.

* Added module :ref:`tudatpy.basic_astrodynamics`.

  * Exposed enumeration ``AvailableAcceleration``.
  * Partial exposure of 3D specialised ``AccelerationModel`` (declaration only).

* Added module :ref:`tudatpy.gravitation`.

  * Partially exposure of class ``GravityFieldModel`` for ``gravitational_parameter`` property.

* Added module :ref:`tudatpy.numerical_integrators`.

  * Exposed enumeration ``AvailableIntegrators`` for ``IntegratorSettings`` argument.
  * Partial exposure of class ``IntegratorSettings`` (only constructor).

* Added module :ref:`tudatpy.propagators`.

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


* Added module :ref:`tudatpy.orbital_element_conversions`.

  * Exposed enumeration ``KeplerianElementIndices``.
  * Exposed function ``tudatpy.convert_keplerian_to_cartesian_elements()``.

* Added module :ref:`tudatpy.simulation_setup`.

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
~~~~~

* Fixed hard requirement for Python 3.7
* Modified installation procedure for ``python3-dev`` requirment prior to compilation of the tudatBundle using the
  ``USE_TUDATPY`` option.
