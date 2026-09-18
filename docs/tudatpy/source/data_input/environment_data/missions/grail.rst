.. _environment_data_missions_grail:

``grail``
=========

This module contains readers for GRAIL mission files that provide spacecraft
mass and antenna-switch information. These files are typically used together
with GRAIL SPICE kernels and radio-science tracking data. The mass files can be
used when defining time-dependent spacecraft mass properties, while antenna
files provide the switch history needed to assign the appropriate spacecraft
reference point when processing GRAIL tracking observations.

.. automodule:: tudatpy.data_input.environment_data.missions.grail

Functions
---------
.. currentmodule:: tudatpy.data_input.environment_data.missions.grail

.. autosummary::

   grail_antenna_file_reader
   grail_mass_level_0_file_reader
   grail_mass_level_1_file_reader

.. autofunction:: tudatpy.data_input.environment_data.missions.grail.grail_antenna_file_reader

.. autofunction:: tudatpy.data_input.environment_data.missions.grail.grail_mass_level_0_file_reader

.. autofunction:: tudatpy.data_input.environment_data.missions.grail.grail_mass_level_1_file_reader
