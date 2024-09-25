``element_conversion``
======================
Functions for converting between sets of orbital elements.

This module provide a variety of functions and classes to
convert between different representations of translational and
rotational states (e.g. Cartesian ↔ Keplerian).

.. note:: Rotations between different reference frames are provided in
          the :ref:`\`\`frame_conversion\`\`` module.






Notes
-----

Unless specified otherwise, the Keplerian elements are ordered as:

+-------+---------------------------------------------------------------------------------------+
| Index | Keplerian Element                                                                     |
+-------+---------------------------------------------------------------------------------------+
| ``0`` | Semi-major axis (except if eccentricity = ``1.0``, then represents semi-latus rectum) |
+-------+---------------------------------------------------------------------------------------+
| ``1`` | Eccentricity                                                                          |
+-------+---------------------------------------------------------------------------------------+
| ``2`` | Inclination                                                                           |
+-------+---------------------------------------------------------------------------------------+
| ``3`` | Argument of periapsis                                                                 |
+-------+---------------------------------------------------------------------------------------+
| ``4`` | Longitude of ascending node                                                           |
+-------+---------------------------------------------------------------------------------------+
| ``5`` | True anomaly                                                                          |
+-------+---------------------------------------------------------------------------------------+

Unless specified otherwise, the spherical elements are ordered as:

+-------+--------------------------------------------------------------------------------------+
| Index | Spherical State Element                                                              |
+-------+--------------------------------------------------------------------------------------+
| ``0`` | Radial distance                                                                      |
+-------+--------------------------------------------------------------------------------------+
| ``1`` | Latitude                                                                             |
+-------+--------------------------------------------------------------------------------------+
| ``2`` | Longitude                                                                            |
+-------+--------------------------------------------------------------------------------------+
| ``3`` | Speed                                                                                |
+-------+--------------------------------------------------------------------------------------+
| ``4`` | Flight-path angle                                                                    |
+-------+--------------------------------------------------------------------------------------+
| ``5`` | Heading angle                                                                        |
+-------+--------------------------------------------------------------------------------------+

Unless specified otherwise, the Modified equinoctial elements (MEE; see `here <https://arc.aiaa.org/doi/pdf/10.2514/1.32237>`_, element set *k*) are ordered as follows.
The element :math:`I`, which defines the location of the MEE singularity is not treated as a state element, but is provided/determined separately:

+-------+--------------------------------------------------------------------------------------+
| Index | Modified equinoctial Element                                                         |
+-------+--------------------------------------------------------------------------------------+
| ``0`` | Semi-parameter                                                                       |
+-------+--------------------------------------------------------------------------------------+
| ``1`` | f-element                                                                            |
+-------+--------------------------------------------------------------------------------------+
| ``2`` | g-element                                                                            |
+-------+--------------------------------------------------------------------------------------+
| ``3`` | h-element                                                                            |
+-------+--------------------------------------------------------------------------------------+
| ``4`` | k-element                                                                            |
+-------+--------------------------------------------------------------------------------------+
| ``5`` | True longitude                                                                       |
+-------+--------------------------------------------------------------------------------------+








Functions
---------
.. currentmodule:: tudatpy.astro.element_conversion

.. autosummary::

   cartesian_to_keplerian

   keplerian_to_cartesian

   keplerian_to_cartesian_elementwise

   mean_to_true_anomaly

   true_to_mean_anomaly

   true_to_eccentric_anomaly

   eccentric_to_true_anomaly

   eccentric_to_mean_anomaly

   mean_to_eccentric_anomaly

   elapsed_time_to_delta_mean_anomaly

   delta_mean_anomaly_to_elapsed_time

   mean_motion_to_semi_major_axis

   semi_major_axis_to_mean_motion

   keplerian_to_mee_manual_singularity

   keplerian_to_mee

   flip_mee_singularity

   mee_to_keplerian

   cartesian_to_mee

   cartesian_to_mee_manual_singularity

   mee_to_cartesian

   quaternion_entries_to_rotation_matrix

   rotation_matrix_to_quaternion_entries

   cartesian_to_spherical

   spherical_to_cartesian

   spherical_to_cartesian_elementwise



.. autofunction:: tudatpy.astro.element_conversion.cartesian_to_keplerian

.. autofunction:: tudatpy.astro.element_conversion.keplerian_to_cartesian

.. autofunction:: tudatpy.astro.element_conversion.keplerian_to_cartesian_elementwise

.. autofunction:: tudatpy.astro.element_conversion.mean_to_true_anomaly

.. autofunction:: tudatpy.astro.element_conversion.true_to_mean_anomaly

.. autofunction:: tudatpy.astro.element_conversion.true_to_eccentric_anomaly

.. autofunction:: tudatpy.astro.element_conversion.eccentric_to_true_anomaly

.. autofunction:: tudatpy.astro.element_conversion.eccentric_to_mean_anomaly

.. autofunction:: tudatpy.astro.element_conversion.mean_to_eccentric_anomaly

.. autofunction:: tudatpy.astro.element_conversion.elapsed_time_to_delta_mean_anomaly

.. autofunction:: tudatpy.astro.element_conversion.delta_mean_anomaly_to_elapsed_time

.. autofunction:: tudatpy.astro.element_conversion.mean_motion_to_semi_major_axis

.. autofunction:: tudatpy.astro.element_conversion.semi_major_axis_to_mean_motion

.. autofunction:: tudatpy.astro.element_conversion.keplerian_to_mee_manual_singularity

.. autofunction:: tudatpy.astro.element_conversion.keplerian_to_mee

.. autofunction:: tudatpy.astro.element_conversion.flip_mee_singularity

.. autofunction:: tudatpy.astro.element_conversion.mee_to_keplerian

.. autofunction:: tudatpy.astro.element_conversion.cartesian_to_mee

.. autofunction:: tudatpy.astro.element_conversion.cartesian_to_mee_manual_singularity

.. autofunction:: tudatpy.astro.element_conversion.mee_to_cartesian

.. autofunction:: tudatpy.astro.element_conversion.quaternion_entries_to_rotation_matrix

.. autofunction:: tudatpy.astro.element_conversion.rotation_matrix_to_quaternion_entries

.. autofunction:: tudatpy.astro.element_conversion.cartesian_to_spherical

.. autofunction:: tudatpy.astro.element_conversion.spherical_to_cartesian

.. autofunction:: tudatpy.astro.element_conversion.spherical_to_cartesian_elementwise







