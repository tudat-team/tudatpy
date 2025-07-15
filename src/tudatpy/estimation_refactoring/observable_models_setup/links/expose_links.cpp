/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_links.h"
#include <pybind11/functional.h>
#include "scalarTypes.h"

#include "tudat/simulation/estimation_setup/createObservationModel.h"

// namespace tss = tudat::simulation_setup;
namespace tom = tudat::observation_models;
// namespace tuc = tudat::unit_conversions;
// namespace ti = tudat::interpolators;

namespace tudatpy
{
namespace estimation_refactoring
{
namespace observable_models_setup
{

namespace links
{

void expose_links( py::module& m )
{
    
     // ################      Link Definition ################

    py::enum_< tom::LinkEndType >( m, "LinkEndType", R"doc(

Enumeration of available link end types.

Examples
--------
.. code-block:: python

    # Code snippet to print all available Link End Types
    from tudatpy.numerical_simulation import estimation_setup

    # Check how many Link End Types are available in Tudatpy
    num_link_end_types = len(estimation_setup.observation.LinkEndType.__members__)
    print(f'The length of all available Tudatpy Link End Types is: {num_link_end_types}')

    # Print all available Link End Types using the "name" property
    for i in range(num_link_end_types):
        print(i, estimation_setup.observation.LinkEndType(i).name)




      )doc" )
            .value( "unidentified_link_end", tom::LinkEndType::unidentified_link_end )
            .value( "transmitter", tom::LinkEndType::transmitter )
            .value( "reflector1", tom::LinkEndType::reflector1 )
            .value( "retransmitter", tom::LinkEndType::retransmitter )
            .value( "reflector2", tom::LinkEndType::reflector2 )
            .value( "reflector3", tom::LinkEndType::reflector3 )
            .value( "reflector4", tom::LinkEndType::reflector4 )
            .value( "receiver", tom::LinkEndType::receiver )
            .value( "transmitter2", tom::LinkEndType::transmitter2 )
            .value( "observer", tom::LinkEndType::observer )
            .value( "observed_body", tom::LinkEndType::observed_body )
            .export_values( );


    m.def( "one_way_downlink_link_ends",
           &tom::getOneWayDownlinkLinkEndsList,
           py::arg( "transmitter" ),
           py::arg( "receivers" ),
           R"doc(

 Function for defining one-way downlinks via LinkDefinition types.

 Function for defining single or multiple one-way downlinks via LinkDefinition types.
 Multiple downlinks share the same transmitters, but may each have a different receiver.
 For each downlink, the returned list will contain an additional `LinkDefinition` type.


 Parameters
 ----------
 transmitter : Tuple[str, str]
     List of :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndId` types (tuple of strings), where, for each tuple, the first entry identifies the body and the second entry reference point of the single transmitter link end(s).

 receivers : List[ Tuple[str, str] ]
     List of :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndId` types (tuple of strings), where for each tuple the first entrance identifies the body and the second entry the reference point of the receiver link end(s).

 Returns
 -------
 List[ LinkDefinition ]
     List of one or more :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkDefinition` types, each defining the geometry for one one-way downlink.
     A `LinkDefinition` type for a one one-way link is composed a dict with one `receiver` and one `transmitter` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` key, to each of which a :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndId` type is mapped.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the one_way_downlink_link_ends function to return a LinkDefinition object

     from tudatpy.kernel.numerical_simulation.estimation_setup import observation

     # Create a dictionary of LinkEndId objects
     link_ends = {
         observation.receiver: observation.body_origin_link_end_id("Earth"),
         observation.transmitter: observation.body_origin_link_end_id("Delfi-C3")
     }

     # Print individual LinkEndId objects
     print("Transmitter:", link_ends[observation.transmitter])
     print("Receiver:", link_ends[observation.receiver])

     # Call one_way_downlink_link_ends with properly formatted arguments
     # Note: The function expects a transmitter and a list of receivers
     link_definition = observation.one_way_downlink_link_ends(
         link_ends[observation.transmitter],
         [link_ends[observation.receiver]]  # Receivers must be in a list
     )

     # Verify that the one_way_downlink_link_ends function returns a LinkDefinition object
     print(link_definition)




     )doc" );

     m.def( "one_way_uplink_link_ends",
           &tom::getOneWayUplinkLinkEndsList,
           py::arg( "transmitters" ),
           py::arg( "receiver" ),
           R"doc(

 Function for defining single or multiple one-way uplinks via LinkDefinition types.

 Function for defining single or multiple one-way uplinks via LinkDefinition types.
 Multiple uplinks share the same receiver, but may each have a different transmitter.
 For each uplink, the returned list will contain an additional `LinkDefinition` type.


 Parameters
 ----------
 transmitters : List[ Tuple[str, str] ]
     List of :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndId` types (tuple of strings), where, for each tuple, the first entry identifies the body and the second entry the reference point of the transmitter link end(s).

 receivers : Tuple[str, str]
     List of :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndId` types (tuple of strings), where, for each tuple, the first entry identifies the body and the second entry the reference point of the single receiver link end(s).

 Returns
 -------
 List[ LinkDefinition ]
     List of one or more :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkDefinition` types, each defining the geometry for one one-way uplink.
     A :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndId` type for a one one-way link is made of a dict with one `receiver` and one `transmitter` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` key, to each of which a :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndId` type is mapped.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the one_way_uplink_link_ends function to return a LinkDefinition object
     from tudatpy.kernel.numerical_simulation.estimation_setup import observation

     # Create a dictionary of LinkEndId objects
     link_ends = {
         observation.receiver: observation.body_origin_link_end_id("Earth"),
         observation.transmitter: observation.body_origin_link_end_id("Delfi-C3")
     }

     # Print individual LinkEndId objects
     print("Transmitter:", link_ends[observation.transmitter])
     print("Receiver:", link_ends[observation.receiver])

     # Call one_way_uplink_link_ends with properly formatted arguments
     # Note: The function expects a transmitter and a list of receivers
     link_definition = observation.one_way_uplink_link_ends(
         [link_ends[observation.transmitter]], # Transmitters must be in a list
         link_ends[observation.receiver]
     )

     # Verify that the one_way_uplink_link_ends function returns a LinkDefinition object
     print(link_definition)




     )doc" );


     m.def( "get_default_reference_link_end",
           &tom::getDefaultReferenceLinkEndType,
           py::arg( "observabl_type" ),
           R"doc(

 Function for automatically retrieving the reference link end associated with a given observable type.


 Parameters
 ----------
 observable_type : :class:`ObservableType`
     Observable type for which the associated reference link end is to be retrieved.
 Returns
 -------
 :class:`LinkEndType`
     Defines the link end (via the :class:`LinkEndType`) which is typically used as a reference for observation times in *e.g.* :func:`~tudatpy.numerical_simulation.estimation_setup.observation.tabulated_simulation_settings`.






     )doc" );

    // ###########      Observation Model Settings
    // ################

    py::class_< tom::LinkEndId, std::shared_ptr< tom::LinkEndId > >( m,
                                                                     "LinkEndId",
                                                                     R"doc(

         Base class serving as identifier of a specific link end.

         Base class serving as identifier of a specific link end.
         Instances of this class are typically created via the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.body_origin_link_end_id` function,
         whose output is indeed a *LinkEndId* object, representing the center of mass of a body.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to produce a LinkEndId object
     from tudatpy.numerical_simulation.estimation_setup import observation

     link_ends = dict()
     link_ends[observation.receiver] = observation.body_origin_link_end_id("Earth")
     link_ends[observation.transmitter] = observation.body_origin_link_end_id("Delfi-C3")

     # The keys of this dictionary are LinkEndType objects.
     print(link_ends.keys())
     # The values of this dictionary are LinkEndId objects.
     print(link_ends.values())

     # Print out (explicitly) the keys (link types) and values (link names).
     # [Note: To accomplish this, we use the "name" property (link_type.name) of the LinkEndType enumeration,
     # and the "body_name" property (link_name.body_name) of the LinkEndId class]

     for link_type, link_name in link_ends.items():
         print(f'LinkEndType: {link_type.name}, LinkEndId: {link_name.body_name}')


      )doc" )
            .def_property_readonly( "body_name",
                                    &tom::LinkEndId::getBodyName,
                                    R"doc(
         Name of the body where the reference point is located, str

     Examples
     --------
     .. code-block:: python

         # Code Snippet to produce a LinkEndId object
         from tudatpy.numerical_simulation.estimation_setup import observation

         link_ends = dict()
         link_ends[observation.receiver] = observation.body_origin_link_end_id("Earth")
         link_ends[observation.transmitter] = observation.body_origin_link_end_id("Delfi-C3")

         # The keys of this dictionary are LinkEndType objects.
         print(link_ends.keys())
         # The values of this dictionary are LinkEndId objects.
         print(link_ends.values())

         # Print out the keys (link types) and values (link names)
         for link_type, link_name in link_ends.items():
             print(f'LinkEndType: {link_type.name}, LinkEndId: {link_name.body_name}')




      )doc" )
            .def_property_readonly( "reference_point",
                                    &tom::LinkEndId::getStationName,
                                    R"doc(
         Function for setting a name for the reference point on a body.

         Function for setting a name for the reference point on a body (tipically, the name of a ground station).

     Examples
     --------
     .. code-block:: python

         # Code Snippet to produce a LinkEndId object (e.g. ground station) on the Earth Surface
         # and retrieve the link reference point using the  "reference_point" property

         from tudatpy.numerical_simulation.estimation_setup import observation

         # Set CoolTrackingStation (defined as a Reference Point on Earth) as a receiver
         link_ends = dict()
         link_ends[observation.receiver] = observation.body_reference_point_link_end_id("Earth", "CoolTrackingStation")

         # Verify that CoolTracking Station is associated to the key observation.receiver
         link_end_body = link_ends[observation.receiver].body_name # body on which the reference point is located
         link_end_name = link_ends[observation.receiver].reference_point #reference point name
         print(f'Link End Name: {link_end_name} is found on body: {link_end_body}')



      )doc" );

    m.def( "body_origin_link_end_id",
           py::overload_cast< const std::string& >( &tom::linkEndId ),
           py::arg( "body_name" ),
           R"doc(

 Function to create a link end identifier for the origin (typically center of mass) of a body.

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
     from tudatpy.numerical_simulation.estimation_setup import observation

     # Input of body_origin_link_end_id are strings (name of bodies, or satellites, or ground stations, etc...)
     receiver = "Earth"
     transmitter = "Delfi-C3"

     # Call and print observation.body_origin_link_end_id with the proper inputs (receiver, transmitter)
     # a LinkEndId object is returned for both receiver and transmitter
     print(observation.body_origin_link_end_id(receiver))
     print(observation.body_origin_link_end_id(transmitter))





     )doc" );

     m.def( "body_reference_point_link_end_id",
           py::overload_cast< const std::string&, const std::string& >( &tom::linkEndId ),
           py::arg( "body_name" ),
           py::arg( "reference_point_id" ),
           R"doc(

 Function to create a link end identifier for a reference point on a body.

 Function to create a link end identifier for a reference point on a body, where the reference point
 is typically the identifier of a ground stations.


 Parameters
 ----------
 body_name : str
     Name of the body on which the reference point is located: :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndId`, str

 reference_point_id : str
     Identifier of a specific link end: :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndId`, str

 Returns
 -------
 LinkEndId
     A LinkEndId object representing a reference point on a body

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the body_reference_point_link_end_id
     from tudatpy.numerical_simulation.estimation_setup import observation

     # Input of body_reference_point_link_end_id are strings (name of bodies, or satellites, or ground stations, etc...)
     receiver = "Earth"
     reference_point = "CoolTrackingStation"

     # Call and print observation.body_reference_point_link_end_id with the proper inputs (receiver, reference_point)
     # a LinkEndId object is returned
     print(observation.body_reference_point_link_end_id(receiver, reference_point))



     )doc" );

    py::class_< tom::LinkDefinition, std::shared_ptr< tom::LinkDefinition > >( m, "LinkDefinition", R"doc(

         Base class storing the link ends involved in a given observation.
         Instances of this class are typically created defining a *Link_Ends* dictionary via the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.link_definition` function,
         whose output is a *LinkDefinition* object, storing the Link Ends involved in a given observation.

         Examples
         --------
         .. code-block:: python

             # Code Snippet to produce a LinkDefinition object
             from tudatpy.numerical_simulation.estimation_setup import observation

             link_ends = dict()
             link_ends[observation.receiver] = observation.body_origin_link_end_id("Earth")
             link_ends[observation.transmitter] = observation.body_origin_link_end_id("Delfi-C3")

             # Show that what we created is a LinkDefinition object
             Link_Definition_Object = observation.link_definition(link_ends)
             print(Link_Definition_Object)

             # [Optional]: Print the Link Ends (receiver and transmitter)  names
             receiver_name = observation.link_definition(link_ends).link_end_id(observation.receiver).body_name
             transmitter_name = observation.link_definition(link_ends).link_end_id(observation.transmitter).body_name
             print(receiver_name)
             print(transmitter_name)




      )doc" )
            .def( py::init< const std::map< tom::LinkEndType, tom::LinkEndId >& >( ), py::arg( "link_ends" ) )
            .def( "link_end_id",
                  &tom::LinkDefinition::at,
                  py::arg( "link_end_type" ),
                  R"doc(

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
             from tudatpy.numerical_simulation.estimation_setup import observation

             link_ends = dict()
             link_ends[observation.receiver] = observation.body_origin_link_end_id("Earth")
             link_ends[observation.transmitter] = observation.body_origin_link_end_id("Delfi-C3")

             Link_Definition_Object = observation.link_definition(link_ends)

             # [Optional] Show that what we created is a LinkDefinition object
             print(Link_Definition_Object)

             # Print the Link Ends (receiver and transmitter)  names using the "link_end_id" property
             print(observation.link_definition(link_ends).link_end_id(observation.receiver).body_name)
             print(observation.link_definition(link_ends).link_end_id(observation.transmitter).body_name)



      )doc" );
    //            .def_property( "link_ends",
    //            &tom::LinkDefinition::linkEnds_,
    //                           get_docstring("LinkDefinition.link_ends").c_str()
    //                           );

    m.def( "link_definition",
           &tom::linkDefinition,
           py::arg( "link_ends" ),
           R"doc(

 Function to create a link definition object.

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
     from tudatpy.numerical_simulation.estimation_setup import observation

     # Create link_ends. These are the input parameters of the link_definition function
     link_ends = dict()
     link_ends[observation.receiver] = observation.body_origin_link_end_id("Earth")
     link_ends[observation.transmitter] = observation.body_origin_link_end_id("Delfi-C3")

     # Show that, using observation.link_definition, a LinkDefinition object is returned
     print(observation.link_definition(link_ends))


     )doc" );

}

}
}
}
}