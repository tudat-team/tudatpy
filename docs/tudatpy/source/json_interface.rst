.. _json_interface:

``json_interface``
==================

This module creates Tudat settings objects from JSON documents governed by the
packaged contracts. A JSON object of the form ``{"factory_name": {...}}`` invokes
the corresponding settings factory after its inputs have been validated.
Optional properties omitted from a document are passed to the Tudat factory using
the default recorded in the contract. A ``null`` contract default delegates to the
factory default, for settings that are absent or defaults that JSON cannot represent.
Inputs listed in a factory's ``unsupported`` array remain part of the API contract
but cannot yet be supplied through JSON. Optional unsupported inputs always use
their default; a factory with a required unsupported input cannot yet be loaded.

The special object ``{"$ref": "relative/path.json"}`` includes another JSON
document. Paths are resolved relative to the file containing the reference.

For translational propagation, the JSON document contains acceleration settings
rather than acceleration models. The system of bodies is supplied as runtime
context, and the acceleration models are created before the ordinary
translational propagator factory is called:

.. code-block:: python

   from tudatpy import dynamics
   from tudatpy.json_interface import (
       load_body_list_settings,
       load_translational_propagator_settings,
   )

   body_list_settings = load_body_list_settings("body_list_settings.json")
   bodies = dynamics.environment_setup.create_system_of_bodies(body_list_settings)
   propagator_settings = load_translational_propagator_settings(
       "propagator_settings.json", bodies
   )
   dynamics_simulator = dynamics.simulator.create_dynamics_simulator(
       bodies, propagator_settings
   )

Functions
---------

.. autofunction:: tudatpy.json_interface.load_body_list_settings

.. autofunction:: tudatpy.json_interface.load_acceleration_settings

.. autofunction:: tudatpy.json_interface.load_translational_propagator_settings

.. autofunction:: tudatpy.json_interface.load_settings

.. autofunction:: tudatpy.json_interface.validate_contract_against_module

.. autofunction:: tudatpy.json_interface.validate_all_contracts
