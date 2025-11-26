.. _links:

``links``
=========


This module contains a set of factory functions and classes for setting link ends and associated functionality in Tudat.
The concept of a link end in Tudat denotes a participant in an observation (e.g. a transmitter, receiver, retransmitter, etc.).

Functionality in this module is used to create objects of type :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`,
which represents a single participant in an observation, enums of type :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType`,
which denotes the role a participant has in an observation (transmitter,receiver, etc.) and :class:`~tudatpy.estimation.observable_models_setup.links.LinkDefinition`,
which compines a set of the former two into a complete set of participants and their roles in Tudat. These types are used
in many places in Tudat related to observations.

On the user guide, we provide a more complete overview on the concept of  `link ends in Tudat <https://docs.tudat.space/en/latest/user-guide/state-estimation/link-ends-setup.html>`_ and their interaction with the
`creation of observation models <https://docs.tudat.space/en/latest/user-guide/state-estimation/observation-model-setup.html>`_.


Functions
---------
.. currentmodule:: tudatpy.estimation.observable_models_setup.links

.. autosummary::

   one_way_downlink_link_ends

   one_way_uplink_link_ends

   get_default_reference_link_end

   body_origin_link_end_id

   body_reference_point_link_end_id

   link_definition

.. autofunction:: tudatpy.estimation.observable_models_setup.links.one_way_downlink_link_ends

.. autofunction:: tudatpy.estimation.observable_models_setup.links.one_way_uplink_link_ends

.. autofunction:: tudatpy.estimation.observable_models_setup.links.get_default_reference_link_end

.. autofunction:: tudatpy.estimation.observable_models_setup.links.body_origin_link_end_id

.. autofunction:: tudatpy.estimation.observable_models_setup.links.body_reference_point_link_end_id

.. autofunction:: tudatpy.estimation.observable_models_setup.links.link_definition

Enumerations
------------
.. currentmodule:: tudatpy.estimation.observable_models_setup.links

.. autosummary::

   LinkEndType


.. autoclass:: tudatpy.estimation.observable_models_setup.links.LinkEndType
   :members:




Classes
-------
.. currentmodule:: tudatpy.estimation.observable_models_setup.links

.. autosummary::

   LinkEndId

   LinkDefinition

.. autoclass:: tudatpy.estimation.observable_models_setup.links.LinkEndId
   :members:

.. autoclass:: tudatpy.estimation.observable_models_setup.links.LinkDefinition
   :members:
   :special-members: __init__
