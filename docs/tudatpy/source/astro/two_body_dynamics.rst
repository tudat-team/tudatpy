.. _two_body_dynamics:

``two_body_dynamics``
=====================
Functions for (semi-)analytical calculations in a simple two-body point-mass system.

This module provides functions and classes for analytical solutions to two-body orbital mechanics problems, including:

- Kepler orbit propagation
- Escape and capture delta-v calculations
- Lambert problem solvers (orbit determination between two positions)
- Gravity assist calculations





Classes
-------
.. currentmodule:: tudatpy.astro.two_body_dynamics

.. autosummary::

   LambertTargeter
   LambertTargeterGooding
   LambertTargeterIzzo
   ZeroRevolutionLambertTargeterIzzo
   MultiRevolutionLambertTargeterIzzo
   PericenterFindingFunctions
   EccentricityFindingFunctions


Functions
---------
.. autosummary::

   propagate_kepler_orbit
   compute_escape_or_capture_delta_v


Lambert Problem Solvers
~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: tudatpy.astro.two_body_dynamics.LambertTargeter
   :members:
   :special-members: __init__

.. autoclass:: tudatpy.astro.two_body_dynamics.LambertTargeterGooding
   :members:
   :special-members: __init__
   :show-inheritance:

.. autoclass:: tudatpy.astro.two_body_dynamics.LambertTargeterIzzo
   :members:
   :special-members: __init__
   :show-inheritance:

.. autoclass:: tudatpy.astro.two_body_dynamics.ZeroRevolutionLambertTargeterIzzo
   :members:
   :special-members: __init__
   :show-inheritance:

.. autoclass:: tudatpy.astro.two_body_dynamics.MultiRevolutionLambertTargeterIzzo
   :members:
   :special-members: __init__
   :show-inheritance:


Gravity Assist Utilities
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: tudatpy.astro.two_body_dynamics.PericenterFindingFunctions
   :members:
   :special-members: __init__

.. autoclass:: tudatpy.astro.two_body_dynamics.EccentricityFindingFunctions
   :members:
   :special-members: __init__


Orbital Propagation Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: tudatpy.astro.two_body_dynamics.propagate_kepler_orbit

.. autofunction:: tudatpy.astro.two_body_dynamics.compute_escape_or_capture_delta_v




Notes
-----

**Lambert's Problem:**

Lambert's problem is the determination of an orbit connecting two position vectors in a specified time of flight. 
This is a fundamental problem in:

- Preliminary mission design
- Interplanetary trajectory optimization  
- Rendezvous and intercept calculations

**Multi-Revolution Transfers:**

For long time-of-flight scenarios, multiple solutions may exist corresponding to different numbers of complete 
revolutions around the central body. The :class:`~tudatpy.astro.two_body_dynamics.MultiRevolutionLambertTargeterIzzo` 
class can compute these solutions.

**Solver Selection:**

- Use :class:`~tudatpy.astro.two_body_dynamics.ZeroRevolutionLambertTargeterIzzo` for most practical applications (fastest)
- Use :class:`~tudatpy.astro.two_body_dynamics.MultiRevolutionLambertTargeterIzzo` when considering phasing orbits
- Use :class:`~tudatpy.astro.two_body_dynamics.LambertTargeterGooding` for maximum robustness across edge cases




References
----------

.. bibliography::
   :filter: docname in docnames
   :keyprefix: two_body_dynamics-