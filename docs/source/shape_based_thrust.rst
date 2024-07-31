``shape_based_thrust``
======================
Functionalities for shape-based low-thrust trajectory design.

This module provides the functionality for creating shape-based low-thrust trajectory design legs.
Without using any numerical propagation, this generates a preliminary low-thrust trajectory using
relatively simple (semi-)analytical methods








References
----------
.. [1] Gondelach, David J., and Ron Noomen. "Hodographic-shaping method for 
       low-thrust interplanetary trajectory design." Journal of Spacecraft and Rockets 52.3 (2015): 728-738.






Functions
---------
.. currentmodule:: tudatpy.trajectory_design.shape_based_thrust

.. autosummary::

   recommended_radial_hodograph_functions

   recommended_normal_hodograph_functions

   recommended_axial_hodograph_functions

   hodograph_constant

   hodograph_sine

   hodograph_cosine

   hodograph_exponential

   hodograph_exponential_sine

   hodograph_exponential_cosine

   hodograph_power

   hodograph_power_sine

   hodograph_power_cosine



.. autofunction:: tudatpy.trajectory_design.shape_based_thrust.recommended_radial_hodograph_functions

.. autofunction:: tudatpy.trajectory_design.shape_based_thrust.recommended_normal_hodograph_functions

.. autofunction:: tudatpy.trajectory_design.shape_based_thrust.recommended_axial_hodograph_functions

.. autofunction:: tudatpy.trajectory_design.shape_based_thrust.hodograph_constant

.. autofunction:: tudatpy.trajectory_design.shape_based_thrust.hodograph_sine

.. autofunction:: tudatpy.trajectory_design.shape_based_thrust.hodograph_cosine

.. autofunction:: tudatpy.trajectory_design.shape_based_thrust.hodograph_exponential

.. autofunction:: tudatpy.trajectory_design.shape_based_thrust.hodograph_exponential_sine

.. autofunction:: tudatpy.trajectory_design.shape_based_thrust.hodograph_exponential_cosine

.. autofunction:: tudatpy.trajectory_design.shape_based_thrust.hodograph_power

.. autofunction:: tudatpy.trajectory_design.shape_based_thrust.hodograph_power_sine

.. autofunction:: tudatpy.trajectory_design.shape_based_thrust.hodograph_power_cosine






Classes
-------
.. currentmodule:: tudatpy.trajectory_design.shape_based_thrust

.. autosummary::

   BaseFunctionHodographicShaping



.. autoclass:: tudatpy.trajectory_design.shape_based_thrust.BaseFunctionHodographicShaping
   :members:



