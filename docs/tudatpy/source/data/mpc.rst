.. _mpc:

``mpc``
=======
Interface to the Minor Planet Center (MPC) database for asteroid and comet observations.

This module contains a wrapper for the MPC interface of the astroquery package,
with additional functionality for applying observation weights and star catalog biases.


Classes
-------
.. currentmodule:: tudatpy.data.mpc

.. autosummary::

   BatchMPC


.. autoclass:: tudatpy.data.mpc.BatchMPC
   :members:
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

