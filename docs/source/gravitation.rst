``gravitation``
===============
Utility functions for calculations related to (spherical harmonic) gravity fields


This module contains a list of utility functions for calculations related to (spherical harmonic) gravity fields.
Note that the calculations relating to gravity fields that are relevant for a numerical propagation and estimation are
done in the relevant environment, acceleration, etc. models in the \`\`numerical_simulation\`\``. The functions in
this module are meant to support the user on relevant pre- and post-processing steps.












Functions
---------
.. currentmodule:: tudatpy.astro.gravitation

.. autosummary::

   legendre_normalization_factor

   normalize_spherical_harmonic_coefficients

   unnormalize_spherical_harmonic_coefficients

   spherical_harmonic_coefficients_from_inertia



.. autofunction:: tudatpy.astro.gravitation.legendre_normalization_factor

.. autofunction:: tudatpy.astro.gravitation.normalize_spherical_harmonic_coefficients

.. autofunction:: tudatpy.astro.gravitation.unnormalize_spherical_harmonic_coefficients

.. autofunction:: tudatpy.astro.gravitation.spherical_harmonic_coefficients_from_inertia







