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
   IcgemGravityFieldData
   OceanPoleGravityCoefficient
   OceanPoleGravityData
   HarmonicPotentialCoefficient
   HarmonicPotentialConstituent
   HarmonicPotentialData

.. autoclass:: IcgemGravityCoefficient
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
