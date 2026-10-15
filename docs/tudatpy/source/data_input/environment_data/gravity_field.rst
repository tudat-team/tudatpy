.. _gravity_field_data:

``gravity_field``
=================

.. automodule:: tudatpy.data_input.environment_data.gravity_field

These readers retain provider coefficients, units, normalization and provenance.
They do not construct gravity models, infer missing GM/radius or permanent-tide
conventions, or add loading and phase corrections. A missing coefficient remains
distinct from a published zero. Coefficient dictionaries use ``(degree, order)``
keys. Input files must already be available locally and decompressed.

Use :func:`~tudatpy.dynamics.environment_setup.gravity_field.from_icgem` for a
static field, and the factories in :ref:`gravity_field_variation` for gravity
tides. Ocean-pole gravity A/B integrals are different from station-displacement
grids in :ref:`ocean_loading_data`.

.. currentmodule:: tudatpy.data_input.environment_data.gravity_field

Functions
---------

.. autosummary::

   read_icgem_gravity_field
   read_ocean_pole_gravity_data
   read_harmonic_potential_data

.. autofunction:: read_icgem_gravity_field

.. autofunction:: read_ocean_pole_gravity_data

.. autofunction:: read_harmonic_potential_data

Enumerations
------------

.. autosummary::

   TideSystem
   IcgemNormalization
   IcgemCoefficientErrors
   HarmonicPotentialProfile

.. autoclass:: TideSystem
   :members:

.. autoclass:: IcgemNormalization
   :members:

.. autoclass:: IcgemCoefficientErrors
   :members:

.. autoclass:: HarmonicPotentialProfile
   :members:

Records
-------

.. autosummary::

   IcgemGravityCoefficient
   IcgemTemporalCoefficient
   IcgemGravityFieldData
   OceanPoleGravityCoefficient
   OceanPoleGravityData
   HarmonicPotentialCoefficient
   HarmonicPotentialConstituent
   HarmonicPotentialData

.. autoclass:: IcgemGravityCoefficient
   :members:

.. autoclass:: IcgemTemporalCoefficient
   :members:

.. autoclass:: IcgemGravityFieldData
   :members:

.. autoclass:: OceanPoleGravityCoefficient
   :members:

.. autoclass:: OceanPoleGravityData
   :members:

.. autoclass:: HarmonicPotentialCoefficient
   :members:

.. autoclass:: HarmonicPotentialConstituent
   :members:

.. autoclass:: HarmonicPotentialData
   :members:

Non-tidal RL07 records
---------------------

RL07 AOD1B fields retain their printed GPS epochs and separate ATM, OCN, GLO and
OBA components. They are potential anomalies relative to the provider's
2007–2014 mean. The reader does not sum components or reinterpret GPS dates as
UTC. Use ``Aod1bComponent.global_mass`` for GLO; OBA has a different physical
domain and is not an alternative decomposition of GLO.

.. autosummary::

   read_aod1b_gravity_field
   read_aod1b_gravity_fields
   Aod1bComponent
   Aod1bGravityCoefficient
   Aod1bGravitySample
   Aod1bGravityFieldData

.. autofunction:: read_aod1b_gravity_field

.. autofunction:: read_aod1b_gravity_fields

.. autoclass:: Aod1bComponent
   :members:

.. autoclass:: Aod1bGravityCoefficient
   :members:

.. autoclass:: Aod1bGravitySample
   :members:

.. autoclass:: Aod1bGravityFieldData
   :members:

ICGEM temporal records
----------------------

``IcgemGravityFieldData.temporal_coefficients`` retains each gfct, trnd, dot,
asin and acos record, including its original values, uncertainties, independent
printed dates and optional period in years. The reader does not infer a physical
time scale or a year convention, and does not evaluate these records. The static
``from_icgem`` adapter rejects a temporal collection instead of discarding it.
