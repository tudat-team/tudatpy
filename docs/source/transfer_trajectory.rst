``transfer_trajectory``
=======================
Functionalities for multiple-leg interplanetary transfer trajectories.

This module provides the functionality for creating transfer trajectories consisting of
multiple transfer legs separated by transfer nodes. In particular, this allows defining multi-gravity assist
transfers with high-thrust legs, low-thrust legs or a mix of the two.








References
----------
.. [1] Izzo, D. 1st ACT global trajectory optimisation competition: Problem
      description and summary of the results. Acta Astronautica 61, 9 (2007).
      https://doi.org/10.1016/j.actaastro.2007.03.003
.. [2] Musegaas, P. Optimization of Space Trajectories Including Multiple
      Gravity Assists and Deep Space Maneuvers. Delft University of Technology,
      MSc thesis (2012).
      http://resolver.tudelft.nl/uuid:02468c77-5c64-4df8-9a24-1ed7ad9d1408
.. [3] Roegiers, T. Application of the Spherical Shaping Method to a Low-Thrust Multiple Asteroid Rendezvous
      Mission: Implementation, Limitations and Solutions. Delft University of Technology,
      MSc thesis (2014).
      http://resolver.tudelft.nl/uuid:9994980c-e1df-47a6-880a-3a226797e34a
.. [4] Gondelach, D. A Hodographic-Shaping Method for Low-Thrust Trajectory Design.
      Delft University of Technology, MSc thesis (2012).
      http://resolver.tudelft.nl/uuid:6a4f1673-88b1-4823-b2ef-9d864c84ab11






Functions
---------
.. currentmodule:: tudatpy.trajectory_design.transfer_trajectory

.. autosummary::

   mga_settings_unpowered_unperturbed_legs

   mga_settings_dsm_velocity_based_legs

   mga_settings_dsm_position_based_legs

   mga_settings_spherical_shaping_legs

   mga_settings_hodographic_shaping_legs

   mga_settings_hodographic_shaping_legs_with_recommended_functions

   unpowered_leg

   dsm_position_based_leg

   dsm_velocity_based_leg

   spherical_shaping_leg

   hodographic_shaping_leg

   swingby_node

   departure_node

   capture_node

   print_parameter_definitions

   create_transfer_trajectory



.. autofunction:: tudatpy.trajectory_design.transfer_trajectory.mga_settings_unpowered_unperturbed_legs

.. autofunction:: tudatpy.trajectory_design.transfer_trajectory.mga_settings_dsm_velocity_based_legs

.. autofunction:: tudatpy.trajectory_design.transfer_trajectory.mga_settings_dsm_position_based_legs

.. autofunction:: tudatpy.trajectory_design.transfer_trajectory.mga_settings_spherical_shaping_legs

.. autofunction:: tudatpy.trajectory_design.transfer_trajectory.mga_settings_hodographic_shaping_legs

.. autofunction:: tudatpy.trajectory_design.transfer_trajectory.mga_settings_hodographic_shaping_legs_with_recommended_functions

.. autofunction:: tudatpy.trajectory_design.transfer_trajectory.unpowered_leg

.. autofunction:: tudatpy.trajectory_design.transfer_trajectory.dsm_position_based_leg

.. autofunction:: tudatpy.trajectory_design.transfer_trajectory.dsm_velocity_based_leg

.. autofunction:: tudatpy.trajectory_design.transfer_trajectory.spherical_shaping_leg

.. autofunction:: tudatpy.trajectory_design.transfer_trajectory.hodographic_shaping_leg

.. autofunction:: tudatpy.trajectory_design.transfer_trajectory.swingby_node

.. autofunction:: tudatpy.trajectory_design.transfer_trajectory.departure_node

.. autofunction:: tudatpy.trajectory_design.transfer_trajectory.capture_node

.. autofunction:: tudatpy.trajectory_design.transfer_trajectory.print_parameter_definitions

.. autofunction:: tudatpy.trajectory_design.transfer_trajectory.create_transfer_trajectory




Enumerations
------------
.. currentmodule:: tudatpy.trajectory_design.transfer_trajectory

.. autosummary::

   TransferLegTypes



.. autoclass:: tudatpy.trajectory_design.transfer_trajectory.TransferLegTypes
   :members:




Classes
-------
.. currentmodule:: tudatpy.trajectory_design.transfer_trajectory

.. autosummary::

   TransferLeg

   SphericalShapingLeg

   HodographicShapingLeg

   TransferNodeSettings

   SwingbyNodeSettings

   EscapeAndDepartureNodeSettings

   CaptureAndInsertionNodeSettings

   TransferLegSettings

   TransferTrajectory



.. autoclass:: tudatpy.trajectory_design.transfer_trajectory.TransferLeg
   :members:

.. autoclass:: tudatpy.trajectory_design.transfer_trajectory.SphericalShapingLeg
   :members:

.. autoclass:: tudatpy.trajectory_design.transfer_trajectory.HodographicShapingLeg
   :members:

.. autoclass:: tudatpy.trajectory_design.transfer_trajectory.TransferNodeSettings
   :members:

.. autoclass:: tudatpy.trajectory_design.transfer_trajectory.SwingbyNodeSettings
   :members:

.. autoclass:: tudatpy.trajectory_design.transfer_trajectory.EscapeAndDepartureNodeSettings
   :members:

.. autoclass:: tudatpy.trajectory_design.transfer_trajectory.CaptureAndInsertionNodeSettings
   :members:

.. autoclass:: tudatpy.trajectory_design.transfer_trajectory.TransferLegSettings
   :members:

.. autoclass:: tudatpy.trajectory_design.transfer_trajectory.TransferTrajectory
   :members:



