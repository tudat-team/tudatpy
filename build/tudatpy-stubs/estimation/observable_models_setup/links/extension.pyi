import typing
__all__ = ['LinkDefinition', 'LinkEndId', 'LinkEndType', 'body_origin_link_end_id', 'body_reference_point_link_end_id', 'get_default_reference_link_end', 'link_definition', 'observed_body', 'observer', 'one_way_downlink_link_ends', 'one_way_uplink_link_ends', 'receiver', 'reflector1', 'reflector2', 'reflector3', 'reflector4', 'retransmitter', 'transmitter', 'transmitter2', 'unidentified_link_end']

class LinkDefinition:
    """Base class storing the link ends involved in a given observation.
    Instances of this class are typically created defining a *Link_Ends* dictionary via the :func:`~tudatpy.estimation.observable_models_setup.links.link_definition` function,
    whose output is a *LinkDefinition* object, storing the Link Ends involved in a given observation.
    
    Examples
    --------
    .. code-block:: python
    
        # Code Snippet to produce a LinkDefinition object
        from tudatpy.estimation.observable_models_setup import links
    
        link_ends = dict()
        link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
        link_ends[links.transmitter] = links.body_origin_link_end_id("Delfi-C3")
    
        # Show that what we created is a LinkDefinition object
        Link_Definition_Object = links.link_definition(link_ends)
        print(Link_Definition_Object)
    
        # [Optional]: Print the Link Ends (receiver and transmitter)  names
        receiver_name = links.link_definition(link_ends).link_end_id(links.receiver).body_name
        transmitter_name = links.link_definition(link_ends).link_end_id(links.transmitter).body_name
        print(receiver_name)
        print(transmitter_name)"""

    def __init__(self, link_ends: dict[LinkEndType, LinkEndId]) -> None:
        ...

    def link_end_id(self, link_end_type: LinkEndType) -> LinkEndId:
        """
                 Function to provide a dictionary of link ends.
        
                 Function to provide a dictionary of link ends, with the key denoting the role in the observation, and the associated value the identifier for the link end.
        
                 Parameters
                 ----------
                 "link_end_type" : :type: LinkEndType
        
                 Returns
                 -------
                 :type: dict[LinkEndType,LinkEndId]
                     Dictionary of link ends
        
                 Examples
                 --------
                 .. code-block:: python
        
                     # Code Snippet to retrieve the LinkEnds names from a LinkDefinition object,
                     # using the "link_end_id" property of LinkDefinition (LinkDefinition.link_end_id)
                     from tudatpy.estimation.observable_models_setup import links
        
                     link_ends = dict()
                     link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
                     link_ends[links.transmitter] = links.body_origin_link_end_id("Delfi-C3")
        
                     Link_Definition_Object = links.link_definition(link_ends)
        
                     # [Optional] Show that what we created is a LinkDefinition object
                     print(Link_Definition_Object)
        
                     # Print the Link Ends (receiver and transmitter)  names using the "link_end_id" property
                     print(links.link_definition(link_ends).link_end_id(links.receiver).body_name)
                     print(links.link_definition(link_ends).link_end_id(links.transmitter).body_name)
        """

class LinkEndId:
    """Base class serving as identifier of a specific link end.
    
            Base class serving as identifier of a specific link end.
            Instances of this class are typically created via the :func:`~tudatpy.estimation.observable_models_setup.links.body_origin_link_end_id` function,
            whose output is indeed a *LinkEndId* object, representing the center of mass of a body.
    
    Examples
    --------
    .. code-block:: python
    
        # Code Snippet to produce a LinkEndId object
        from tudatpy.estimation.observable_models_setup import links
    
        link_ends = dict()
        link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
        link_ends[links.transmitter] = links.body_origin_link_end_id("Delfi-C3")
    
        # The keys of this dictionary are LinkEndType objects.
        print(link_ends.keys())
        # The values of this dictionary are LinkEndId objects.
        print(link_ends.values())
    
        # Print out (explicitly) the keys (link types) and values (link names).
        # [Note: To accomplish this, we use the "name" property (link_type.name) of the LinkEndType enumeration,
        # and the "body_name" property (link_name.body_name) of the LinkEndId class]
    
        for link_type, link_name in link_ends.items():
            print(f'LinkEndType: {link_type.name}, LinkEndId: {link_name.body_name}')
    
    
         """

    @property
    def body_name(self) -> str:
        """
                 Name of the body where the reference point is located, str
        
             Examples
             --------
             .. code-block:: python
        
                 # Code Snippet to produce a LinkEndId object
                 from tudatpy.estimation.observable_models_setup import links
        
                 link_ends = dict()
                 link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
                 link_ends[links.transmitter] = links.body_origin_link_end_id("Delfi-C3")
        
                 # The keys of this dictionary are LinkEndType objects.
                 print(link_ends.keys())
                 # The values of this dictionary are LinkEndId objects.
                 print(link_ends.values())
        
                 # Print out the keys (link types) and values (link names)
                 for link_type, link_name in link_ends.items():
                     print(f'LinkEndType: {link_type.name}, LinkEndId: {link_name.body_name}')
        """

    @property
    def reference_point(self) -> str:
        """
                 Function for setting a name for the reference point on a body.
        
                 Function for setting a name for the reference point on a body (tipically, the name of a ground station).
        
             Examples
             --------
             .. code-block:: python
        
                 # Code Snippet to produce a LinkEndId object (e.g. ground station) on the Earth Surface
                 # and retrieve the link reference point using the  "reference_point" property
        
                 from tudatpy.estimation.observable_models_setup import links
        
                 # Set CoolTrackingStation (defined as a Reference Point on Earth) as a receiver
                 link_ends = dict()
                 link_ends[links.receiver] = links.body_reference_point_link_end_id("Earth", "CoolTrackingStation")
        
                 # Verify that CoolTracking Station is associated to the key links.receiver
                 link_end_body = link_ends[links.receiver].body_name # body on which the reference point is located
                 link_end_name = link_ends[links.receiver].reference_point #reference point name
                 print(f'Link End Name: {link_end_name} is found on body: {link_end_body}')
        """

class LinkEndType:
    """Enumeration of available link end types.
    
    Examples
    --------
    .. code-block:: python
    
        # Code snippet to print all available Link End Types
        from tudatpy.estimation import observable_models_setup
    
        # Check how many Link End Types are available in Tudatpy
        num_link_end_types = len(observable_models_setup.links.LinkEndType.__members__)
        print(f'The length of all available Tudatpy Link End Types is: {num_link_end_types}')
    
        # Print all available Link End Types using the "name" property
        for i in range(num_link_end_types):
            print(i, observable_models_setup.links.LinkEndType(i).name)
    
    
    
    
          
    
    Members:
    
      unidentified_link_end
    
      transmitter
    
      reflector1
    
      retransmitter
    
      reflector2
    
      reflector3
    
      reflector4
    
      receiver
    
      transmitter2
    
      observer
    
      observed_body"""
    __members__: typing.ClassVar[dict[str, LinkEndType]]
    observed_body: typing.ClassVar[LinkEndType]
    observer: typing.ClassVar[LinkEndType]
    receiver: typing.ClassVar[LinkEndType]
    reflector1: typing.ClassVar[LinkEndType]
    reflector2: typing.ClassVar[LinkEndType]
    reflector3: typing.ClassVar[LinkEndType]
    reflector4: typing.ClassVar[LinkEndType]
    retransmitter: typing.ClassVar[LinkEndType]
    transmitter: typing.ClassVar[LinkEndType]
    transmitter2: typing.ClassVar[LinkEndType]
    unidentified_link_end: typing.ClassVar[LinkEndType]

    def __eq__(self, other: typing.Any) -> bool:
        ...

    def __getstate__(self) -> int:
        ...

    def __hash__(self) -> int:
        ...

    def __index__(self) -> int:
        ...

    def __init__(self, value: int) -> None:
        ...

    def __int__(self) -> int:
        ...

    def __ne__(self, other: typing.Any) -> bool:
        ...

    def __repr__(self) -> str:
        ...

    def __setstate__(self, state: int) -> None:
        ...

    def __str__(self) -> str:
        ...

    @property
    def name(self) -> str:
        ...

    @property
    def value(self) -> int:
        ...

def body_origin_link_end_id(body_name: str) -> LinkEndId:
    """Function to create a link end identifier for the origin (typically center of mass) of a body.
    
    Function to create a link end identifier for the origin (typically center of mass) of a body.
    Using this option will simulate the origin of a body transmitter, receiving, etc. the observation.
    Although this is typically not physically realistic, it can be a useful approximation, in particular for simulation studies.
    
    
    Parameters
    ----------
    body_name : str
        Name of the body
    
    Returns
    -------
    LinkEndId
        A LinkEndId object representing the center of mass of a body
    
    Examples
    --------
    .. code-block:: python
    
        # Code Snippet to showcase the use of the body_origin_link_end_id
        from tudatpy.estimation.observable_models_setup import links
    
        # Input of body_origin_link_end_id are strings (name of bodies, or satellites, or ground stations, etc...)
        receiver = "Earth"
        transmitter = "Delfi-C3"
    
        # Call and print links.body_origin_link_end_id with the proper inputs (receiver, transmitter)
        # a LinkEndId object is returned for both receiver and transmitter
        print(links.body_origin_link_end_id(receiver))
        print(links.body_origin_link_end_id(transmitter))"""

def body_reference_point_link_end_id(body_name: str, reference_point_id: str) -> LinkEndId:
    """Function to create a link end identifier for a reference point on a body.
    
    Function to create a link end identifier for a reference point on a body, where the reference point
    is typically the identifier of a ground stations.
    
    
    Parameters
    ----------
    body_name : str
        Name of the body on which the reference point is located: :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`, str
    
    reference_point_id : str
        Identifier of a specific link end: :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`, str
    
    Returns
    -------
    LinkEndId
        A LinkEndId object representing a reference point on a body
    
    Examples
    --------
    .. code-block:: python
    
        # Code Snippet to showcase the use of the body_reference_point_link_end_id
        from tudatpy.estimation.observable_models_setup import links
    
        # Input of body_reference_point_link_end_id are strings (name of bodies, or satellites, or ground stations, etc...)
        receiver = "Earth"
        reference_point = "CoolTrackingStation"
    
        # Call and print links.body_reference_point_link_end_id with the proper inputs (receiver, reference_point)
        # a LinkEndId object is returned
        print(links.body_reference_point_link_end_id(receiver, reference_point))"""

def get_default_reference_link_end(observabl_type: ...) -> LinkEndType:
    """Function for automatically retrieving the reference link end associated with a given observable type.
    
    
    Parameters
    ----------
    observable_type : :class:`ObservableType`
        Observable type for which the associated reference link end is to be retrieved.
    Returns
    -------
    :class:`LinkEndType`
        Defines the link end (via the :class:`LinkEndType`) which is typically used as a reference for observation times in *e.g.* :func:`~tudatpy.estimation.observations_setup.tabulated_simulation_settings`."""

def link_definition(link_ends: dict[LinkEndType, LinkEndId]) -> LinkDefinition:
    """Function to create a link definition object.
    
    Function to create a link definition object.
    It returns the ``LinkDefinition`` object storing the link ends of the observation.
    
    
    Parameters
    ----------
    link_ends : dict[LinkEndType,LinkEndId]
        Dictionary of link ends, with the key denoting the role in the observation, and the associated value the identifier for the link end.
    Returns
    -------
    LinkDefinition
        The ``LinkDefinition`` object storing the link ends of the observation
    
    Examples
    --------
    .. code-block:: python
    
        # Code Snippet to showcase the use of the link_definition function to return a LinkDefinition object
        from tudatpy.estimation.observable_models_setup import links
    
        # Create link_ends. These are the input parameters of the link_definition function
        link_ends = dict()
        link_ends[links.receiver] = links.body_origin_link_end_id("Earth")
        link_ends[links.transmitter] = links.body_origin_link_end_id("Delfi-C3")
    
        # Show that, using links.link_definition, a LinkDefinition object is returned
        print(links.link_definition(link_ends))"""

def one_way_downlink_link_ends(transmitter: ..., receivers: list[...]) -> list[dict[LinkEndType, ...]]:
    """Function for defining one-way downlinks via LinkDefinition types.
    
    Function for defining single or multiple one-way downlinks via LinkDefinition types.
    Multiple downlinks share the same transmitters, but may each have a different receiver.
    For each downlink, the returned list will contain an additional `LinkDefinition` type.
    
    
    Parameters
    ----------
    transmitter : Tuple[str, str]
        List of :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId` types (tuple of strings), where, for each tuple, the first entry identifies the body and the second entry reference point of the single transmitter link end(s).
    
    receivers : List[ Tuple[str, str] ]
        List of :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId` types (tuple of strings), where for each tuple the first entrance identifies the body and the second entry the reference point of the receiver link end(s).
    
    Returns
    -------
    List[ LinkDefinition ]
        List of one or more :class:`~tudatpy.estimation.observable_models_setup.links.LinkDefinition` types, each defining the geometry for one one-way downlink.
        A `LinkDefinition` type for a one one-way link is composed a dict with one `receiver` and one `transmitter` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` key, to each of which a :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId` type is mapped.
    
    Examples
    --------
    .. code-block:: python
    
        # Code Snippet to showcase the use of the one_way_downlink_link_ends function to return a LinkDefinition object
    
        from tudatpy.estimation.observable_models_setup import links
    
        # Create a dictionary of LinkEndId objects
        link_ends = {
            links.receiver: links.body_origin_link_end_id("Earth"),
            links.transmitter: links.body_origin_link_end_id("Delfi-C3")
        }
    
        # Print individual LinkEndId objects
        print("Transmitter:", link_ends[links.transmitter])
        print("Receiver:", link_ends[links.receiver])
    
        # Call one_way_downlink_link_ends with properly formatted arguments
        # Note: The function expects a transmitter and a list of receivers
        link_definition = links.one_way_downlink_link_ends(
            link_ends[links.transmitter],
            [link_ends[links.receiver]]  # Receivers must be in a list
        )
    
        # Verify that the one_way_downlink_link_ends function returns a LinkDefinition object
        print(link_definition)"""

def one_way_uplink_link_ends(transmitters: list[...], receiver: ...) -> list[dict[LinkEndType, ...]]:
    """Function for defining single or multiple one-way uplinks via LinkDefinition types.
    
    Function for defining single or multiple one-way uplinks via LinkDefinition types.
    Multiple uplinks share the same receiver, but may each have a different transmitter.
    For each uplink, the returned list will contain an additional `LinkDefinition` type.
    
    
    Parameters
    ----------
    transmitters : List[ Tuple[str, str] ]
        List of :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId` types (tuple of strings), where, for each tuple, the first entry identifies the body and the second entry the reference point of the transmitter link end(s).
    
    receivers : Tuple[str, str]
        List of :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId` types (tuple of strings), where, for each tuple, the first entry identifies the body and the second entry the reference point of the single receiver link end(s).
    
    Returns
    -------
    List[ LinkDefinition ]
        List of one or more :class:`~tudatpy.estimation.observable_models_setup.links.LinkDefinition` types, each defining the geometry for one one-way uplink.
        A :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId` type for a one one-way link is made of a dict with one `receiver` and one `transmitter` :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType` key, to each of which a :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId` type is mapped.
    
    Examples
    --------
    .. code-block:: python
    
        # Code Snippet to showcase the use of the one_way_uplink_link_ends function to return a LinkDefinition object
        from tudatpy.estimation.observable_models_setup import links
    
        # Create a dictionary of LinkEndId objects
        link_ends = {
            links.receiver: links.body_origin_link_end_id("Earth"),
            links.transmitter: links.body_origin_link_end_id("Delfi-C3")
        }
    
        # Print individual LinkEndId objects
        print("Transmitter:", link_ends[links.transmitter])
        print("Receiver:", link_ends[links.receiver])
    
        # Call one_way_uplink_link_ends with properly formatted arguments
        # Note: The function expects a transmitter and a list of receivers
        link_definition = links.one_way_uplink_link_ends(
            [link_ends[links.transmitter]], # Transmitters must be in a list
            link_ends[links.receiver]
        )
    
        # Verify that the one_way_uplink_link_ends function returns a LinkDefinition object
        print(link_definition)"""
observed_body: LinkEndType
observer: LinkEndType
receiver: LinkEndType
reflector1: LinkEndType
reflector2: LinkEndType
reflector3: LinkEndType
reflector4: LinkEndType
retransmitter: LinkEndType
transmitter: LinkEndType
transmitter2: LinkEndType
unidentified_link_end: LinkEndType