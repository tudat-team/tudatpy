.. _mpc:

``mpc``
=======
Interface to the Minor Planet Center (MPC) database for asteroid and comet observations.

This module contains a wrapper for the MPC interface of the astroquery package,
with additional functionality for applying observation weights based on
:cite:t:`veres2017` and star catalog biases based on :cite:t:`eggl2020`.


Classes
-------
.. currentmodule:: tudatpy.data.mpc

.. autosummary::

   BatchMPC


.. autoclass:: tudatpy.data.mpc.BatchMPC
   :members:
   :exclude-members: create_observations_from_astropy_table, to_tudat_observation_collection
   :special-members: __init__


Functions
---------
.. autosummary::

   load_bias_file
   get_biases_EFCC18
   get_weights_VFCC17


.. autofunction:: tudatpy.data.mpc.load_bias_file

.. autofunction:: tudatpy.data.mpc.get_biases_EFCC18

.. autofunction:: tudatpy.data.mpc.get_weights_VFCC17
