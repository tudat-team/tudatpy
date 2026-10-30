.. _close_encounters:

``close_encounters``
====================
Functions for the analysis of close encounters between a small body and a planet.


This module contains functions to describe a close encounter with the osculating two-body hyperbola
that matches the state of the small body relative to the encountered body: the incoming and outgoing
asymptotes, the B-plane (target plane) coordinates in a user-defined reference frame, the coordinates
of the Öpik theory, and the capture cross section of the encountered body.



Functions
---------
.. currentmodule:: tudatpy.astro.close_encounters

.. autosummary::

   compute_b_plane_coordinates

   compute_opik_coordinates

   compute_opik_reference_direction

   compute_incoming_asymptote_direction

   compute_outgoing_asymptote_direction

   compute_capture_radius



.. autofunction:: tudatpy.astro.close_encounters.compute_b_plane_coordinates

.. autofunction:: tudatpy.astro.close_encounters.compute_opik_coordinates

.. autofunction:: tudatpy.astro.close_encounters.compute_opik_reference_direction

.. autofunction:: tudatpy.astro.close_encounters.compute_incoming_asymptote_direction

.. autofunction:: tudatpy.astro.close_encounters.compute_outgoing_asymptote_direction

.. autofunction:: tudatpy.astro.close_encounters.compute_capture_radius




Classes
-------
.. currentmodule:: tudatpy.astro.close_encounters

.. autosummary::

   BPlaneData



.. autoclass:: tudatpy.astro.close_encounters.BPlaneData
   :members:
