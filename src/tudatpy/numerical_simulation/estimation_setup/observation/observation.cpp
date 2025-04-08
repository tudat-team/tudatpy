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
#include "observation.hpp"

#include <pybind11/functional.h>

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"
#include "tudat/simulation/estimation_setup/processOdfFile.h"
#include "tudat/simulation/estimation_setup/processTrackingTxtFile.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"

namespace tss = tudat::simulation_setup;
namespace tom = tudat::observation_models;
namespace tuc = tudat::unit_conversions;
namespace ti = tudat::interpolators;

namespace tudat
{

namespace simulation_setup
{

void addNoiseFunctionToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >&
                observationSimulationSettings,
        const std::function< Eigen::VectorXd( const double ) > observationNoiseFunction )
{
    tss::addNoiseFunctionToObservationSimulationSettings< TIME_TYPE, Eigen::VectorXd >(
            observationSimulationSettings, observationNoiseFunction );
}

void addNoiseFunctionToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >&
                observationSimulationSettings,
        const std::function< Eigen::VectorXd( const double ) > observationNoiseFunction,
        const tom::ObservableType observableType )
{
    tss::addNoiseFunctionToObservationSimulationSettings< TIME_TYPE,
                                                          Eigen::VectorXd,
                                                          const tom::ObservableType >(
            observationSimulationSettings, observationNoiseFunction, observableType );
}

void addNoiseFunctionToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >&
                observationSimulationSettings,
        const std::function< Eigen::VectorXd( const double ) > observationNoiseFunction,
        const tom::ObservableType observableType,
        const tom::LinkDefinition& linkEnds )
{
    tss::addNoiseFunctionToObservationSimulationSettings< TIME_TYPE,
                                                          Eigen::VectorXd,
                                                          const tom::ObservableType,
                                                          const tom::LinkDefinition& >(
            observationSimulationSettings, observationNoiseFunction, observableType, linkEnds );
}

void addGaussianNoiseFunctionToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >&
                observationSimulationSettings,
        const double observationNoiseAmplitude )
{
    tss::addGaussianNoiseFunctionToObservationSimulationSettings< TIME_TYPE >(
            observationSimulationSettings, observationNoiseAmplitude );
}

void addGaussianNoiseFunctionToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >&
                observationSimulationSettings,
        const double observationNoiseAmplitude,
        const tom::ObservableType observableType )
{
    tss::addGaussianNoiseFunctionToObservationSimulationSettings< TIME_TYPE,
                                                                  const tom::ObservableType >(
            observationSimulationSettings, observationNoiseAmplitude, observableType );
}

void addGaussianNoiseFunctionToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >&
                observationSimulationSettings,
        const double observationNoiseAmplitude,
        const tom::ObservableType observableType,
        const tom::LinkDefinition& linkEnds )
{
    tss::addGaussianNoiseFunctionToObservationSimulationSettings< TIME_TYPE,
                                                                  const tom::ObservableType,
                                                                  const tom::LinkDefinition& >(
            observationSimulationSettings, observationNoiseAmplitude, observableType, linkEnds );
}

void addViabilityToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >&
                observationSimulationSettings,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >&
                viabilitySettingsList )
{
    tss::addViabilityToObservationSimulationSettings< TIME_TYPE >( observationSimulationSettings,
                                                                   viabilitySettingsList );
}

void addViabilityToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >&
                observationSimulationSettings,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >&
                viabilitySettingsList,
        const tom::ObservableType observableType )
{
    tss::addViabilityToObservationSimulationSettings< TIME_TYPE, const tom::ObservableType >(
            observationSimulationSettings, viabilitySettingsList, observableType );
}

void addViabilityToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >&
                observationSimulationSettings,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >&
                viabilitySettingsList,
        const tom::ObservableType observableType,
        const tom::LinkDefinition& linkEnds )
{
    tss::addViabilityToObservationSimulationSettings< TIME_TYPE,
                                                      const tom::ObservableType,
                                                      const tom::LinkDefinition& >(
            observationSimulationSettings, viabilitySettingsList, observableType, linkEnds );
}

void addDependentVariablesToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >&
                observationSimulationSettings,
        const std::vector< std::shared_ptr< ObservationDependentVariableSettings > >&
                dependentVariableList,
        const SystemOfBodies& bodies )
{
    tss::addDependentVariablesToObservationSimulationSettings< TIME_TYPE >(
            observationSimulationSettings, dependentVariableList, bodies );
}

void addDependentVariablesToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >&
                observationSimulationSettings,
        const std::vector< std::shared_ptr< ObservationDependentVariableSettings > >&
                dependentVariableList,
        const SystemOfBodies& bodies,
        const tom::ObservableType observableType )
{
    tss::addDependentVariablesToObservationSimulationSettings< TIME_TYPE,
                                                               const tom::ObservableType >(
            observationSimulationSettings, dependentVariableList, bodies, observableType );
}

void addDependentVariablesToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >&
                observationSimulationSettings,
        const std::vector< std::shared_ptr< ObservationDependentVariableSettings > >&
                dependentVariableList,
        const SystemOfBodies& bodies,
        const tom::ObservableType observableType,
        const tom::LinkDefinition& linkEnds )
{
    tss::addDependentVariablesToObservationSimulationSettings< TIME_TYPE,
                                                               const tom::ObservableType,
                                                               const tom::LinkDefinition& >(
            observationSimulationSettings,
            dependentVariableList,
            bodies,
            observableType,
            linkEnds );
}

void addAncilliarySettingsToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >&
                observationSimulationSettings,
        const std::shared_ptr< tom::ObservationAncilliarySimulationSettings >& ancilliarySettings,
        const tom::ObservableType observableType )
{
    tss::addAncilliarySettingsToObservationSimulationSettings< TIME_TYPE,
                                                               const tom::ObservableType >(
            observationSimulationSettings, ancilliarySettings, observableType );
}

void addAncilliarySettingsToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >&
                observationSimulationSettings,
        const std::shared_ptr< tom::ObservationAncilliarySimulationSettings >& ancilliarySettings,
        const tom::ObservableType observableType,
        const tom::LinkDefinition& linkEnds )
{
    tss::addAncilliarySettingsToObservationSimulationSettings< TIME_TYPE,
                                                               const tom::ObservableType,
                                                               const tom::LinkDefinition& >(
            observationSimulationSettings, ancilliarySettings, observableType, linkEnds );
}

}  // namespace simulation_setup

}  // namespace tudat

namespace tudatpy
{
namespace numerical_simulation
{
namespace estimation_setup
{
namespace observation
{

void expose_observation_setup( py::module& m )
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

    py::class_< tom::LinkDefinition, std::shared_ptr< tom::LinkDefinition > >(
            m, "LinkDefinition", R"doc(

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
            .def( py::init< const std::map< tom::LinkEndType, tom::LinkEndId >& >( ),
                  py::arg( "link_ends" ) )
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

    py::enum_< tom::ObservableType >( m, "ObservableType", R"doc(

         Enumeration of available observable types.

 Examples
 --------
 .. code-block:: python

     # Code snippet to print all available Observable Types
     from tudatpy.numerical_simulation import estimation_setup

     num_observable_types = len(estimation_setup.observation.ObservableType.__members__)
     print(f'The length of all available Tudatpy Observable Types is: {num_observable_types}')

     # Print all available Observable Types using the "name" property
     for i in range(num_observable_types):
         print(i, estimation_setup.observation.ObservableType(i).name)




      )doc" )
            .value( "one_way_range_type", tom::ObservableType::one_way_range )
            .value( "n_way_range_type", tom::ObservableType::n_way_range )
            .value( "angular_position_type", tom::ObservableType::angular_position )
            .value( "relative_angular_position_type", tom::ObservableType::angular_position )
            .value( "position_observable_type", tom::ObservableType::position_observable )
            .value( "velocity_observable_type", tom::ObservableType::velocity_observable )
            .value( "one_way_instantaneous_doppler_type", tom::ObservableType::one_way_doppler )
            .value( "one_way_averaged_doppler_type",
                    tom::ObservableType::one_way_differenced_range )
            .value( "two_way_instantaneous_doppler_type", tom::ObservableType::two_way_doppler )
            .value( "n_way_averaged_doppler_type", tom::ObservableType::n_way_differenced_range )
            .value( "euler_angle_313_observable_type",
                    tom::ObservableType::euler_angle_313_observable )
            .value( "dsn_one_way_averaged_doppler",
                    tom::ObservableType::dsn_one_way_averaged_doppler )
            .value( "dsn_n_way_averaged_doppler", tom::ObservableType::dsn_n_way_averaged_doppler )
            .value( "doppler_measured_frequency_type",
                    tom::ObservableType::doppler_measured_frequency )
            .value( "dsn_n_way_range", tom::ObservableType::dsn_n_way_range )
            .export_values( );

    py::class_< tom::DopplerProperTimeRateSettings,
                std::shared_ptr< tom::DopplerProperTimeRateSettings > >(
            m,
            "DopplerProperTimeRateSettings",
            R"doc(

         Base class to define proper time rate settings.

         Base class to define proper time rate settings (at a single link end) for instantaneous Doppler observation model settings.
         Specific proper time rate settings must be defined using an object derived from this class.
         The derived classes are made accessible via dedicated functions.





      )doc" );

    py::class_< tom::ObservationModelSettings, std::shared_ptr< tom::ObservationModelSettings > >(
            m, "ObservationSettings", R"doc(

         Base class to define settings of observation models.

         Base class to define settings of observation models.
         Observation model settings define at least the type and geometry of a given observation.
         They can furthermore set observation biases and/or light-time corrections.
         Simple observation models settings that are fully characterised by these elements can be managed by this base class.
         Instances of this class are typically created via functions, such as
         :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`, :func:`~tudatpy.numerical_simulation.estimation_setup.observation.cartesian_position`,
         :func:`~tudatpy.numerical_simulation.estimation_setup.observation.angular_position`, etc.
         Model settings for specific observation models that require additional information such as integration time, retransmission time, etc. must be defined using an object derived from this class.
         The derived classes are made accessible through further functions.

         Examples
         --------
         .. code-block:: python

             # Code snippet to show the creation of an ObservationSettings object
             from tudatpy.numerical_simulation.estimation_setup import observation

             # Create Link Ends dictionary
             link_ends = dict()
             link_ends[observation.receiver] = observation.body_origin_link_end_id("Earth")
             link_ends[observation.transmitter] = observation.body_origin_link_end_id("Delfi-C3")

             # Create a Link Definition Object from link_ends dictionary
             Link_Definition_Object = observation.LinkDefinition(link_ends)

             # Create minimal ObservationSettings object (only required Link_Definition_Object argument is passed)
             # Other optional parameters (bias_settings, light_time_correction_settings,  light_time_convergence_settings) are set by default
             observation_settings = observation.one_way_range(Link_Definition_Object)

             # Show that it is an ObservationSettings object.
             print(observation_settings)




      )doc" );

    py::class_< tom::OneWayDopplerObservationSettings,
                std::shared_ptr< tom::OneWayDopplerObservationSettings >,
                tom::ObservationModelSettings >( m,
                                                 "OneWayDopplerObservationSettings",
                                                 R"doc(

         Derived Class for defining the settings of one-way instantaneous Doppler observation models.

         Derived Class for defining the settings of one-way instantaneous Doppler observation models.
         Settings object can account for additional observation model aspects such as light time corrections and proper time rate settings.
         Instances of this class can be created via the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_doppler_instantaneous` function.
         Associated base class: :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings`.

         Examples
         --------
         .. code-block:: python

             # Code snippet to show the creation of a OneWayDopplerObservationSettings object
             from tudatpy.numerical_simulation.estimation_setup import observation

             # Create Link Ends dictionary
             link_ends = dict()
             link_ends[observation.receiver] = observation.body_origin_link_end_id("Earth")
             link_ends[observation.transmitter] = observation.body_origin_link_end_id("Delfi-C3")

             # Create a Link Definition Object from link_ends dictionary
             Link_Definition_Object = observation.LinkDefinition(link_ends)

             # Use: observation.one_way_doppler_instantaneous to create a OneWayDopplerObservationSettings object (only required Link_Definition_Object argument is passed)
             # Other optional parameters (bias_settings, light_time_correction_settings,  light_time_convergence_settings, proper time rate) are set by default
             doppler_observation_settings = observation.one_way_doppler_instantaneous(Link_Definition_Object)

             # Show that it is an OneWayDopplerObservationSettings object.
             print(doppler_observation_settings)




      )doc" );

    py::class_< tom::NWayRangeObservationSettings,
                std::shared_ptr< tom::NWayRangeObservationSettings >,
                tom::ObservationModelSettings >(
            m, "NWayRangeObservationSettings", R"doc(No documentation found.)doc" );

    py::class_< tom::LightTimeConvergenceCriteria,
                std::shared_ptr< tom::LightTimeConvergenceCriteria > >(
            m,
            "LightTimeConvergenceCriteria",
            R"doc(

         Base class to define criteria of light time convergence.

         Base class to define criteria of light time convergence.
         This class is not used for calculations of corrections, but is used for the purpose of defining the light time convergence criteria.
         Specific light time convergence criteria must be defined using an object derived from this class.
         Instances of this class are typically created via the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.light_time_convergence_settings` function.

         Examples
         --------
         .. code-block:: python

             # Code snippet to show the creation of a LightTimeConvergenceCriteria object
             from tudatpy.numerical_simulation.estimation_setup import observation

             # Create Default Light Time Convergence Settings (no args specified = setting default arguments)
             light_time_convergence_settings = observation.light_time_convergence_settings()

             # Show that it is an LightTimeConvergenceCriteria object.
             print(light_time_convergence_settings)




      )doc" );

    py::enum_< tom::LightTimeFailureHandling >( m, "LightTimeFailureHandling", R"doc(

         Enumeration of behaviour when failing to converge light-time with required settings.

 Examples
 --------
 .. code-block:: python

     # Code snippet to print all available Light Time Failure Handling Types
     from tudatpy.numerical_simulation import estimation_setup

     num_LightTimeFailureHandling_types = len(estimation_setup.observation.LightTimeFailureHandling.__members__)
     print(f'The length of all available Tudatpy Light Time Failure Handling Types is: {num_LightTimeFailureHandling_types}')

     # Print all available Observation Viability Types using the "name" property
     for i in range(num_LightTimeFailureHandling_types):
         print(i, estimation_setup.observation.LightTimeFailureHandling(i).name)



      )doc" )
            .value( "accept_without_warning",
                    tom::LightTimeFailureHandling::accept_without_warning )
            .value( "print_warning_and_accept",
                    tom::LightTimeFailureHandling::print_warning_and_accept )
            .value( "throw_exception", tom::LightTimeFailureHandling::throw_exception )
            .export_values( );

    m.def( "light_time_convergence_settings",
           &tom::lightTimeConvergenceCriteria,
           py::arg( "iterate_corrections" ) = false,
           py::arg( "maximum_number_of_iterations" ) = 50,
           py::arg( "absolute_tolerance" ) = TUDAT_NAN,
           py::arg( "failure_handling" ) = tom::accept_without_warning,
           R"doc(

 Function for creating convergence settings for solving the light-time equation.

 Function for creating convergence settings for solving the light-time equation. Computing the light time
 :math:`s=t_{R}-t_{T}` between two receiver :math:`R` and transmitter :math:`T` requires the iterative
 solution of the following equation:

 .. math::
     t_{R} - t_{T} = c\left(|\mathbf{r}_{R}(t_{R}) - \mathbf{r}_{T}(t_{T})| + \Delta s(t_{R}, t_{T}, \mathbf{r}_{R}(t_{R}), \mathbf{r}_{T}(t_{T}))\right)

 where either the reception time :math:`t_{R}` or the transmission time :math:`t_{T}` is kept fixed (reference link end time). The term :math:`\Delta s` contains any
 deviations in the light-time from straight-line propagation at speed of light (relativistic corrections, media corrections, etc.). The algorithm starts
 at :math:`t_{R}=t_{T}`, and uses this to evaluate the right-hand side of the above equation. This leads to a new value of :math:`t_{R}` or :math:`t_{T}` (depending on which is kept fixed)
 and the right-hand side is re-evaluated in a new iteration. The input to this function defines the settings for when the iteration will terminate.

 Parameters
 ----------
 iterate_corrections : bool, default = False
     Boolean denoting whether the terms :math:`\Delta s` are recomputed at each iteration or not. If false, the corrections are calculated only on the first iteration. Afterwards, the value
     is kept fixed until convergence. Once preliminarily converged, the algorithm recomputes :math:`\Delta s`, and continues the iteration (until proper convergence) while now recomputing
     :math:`\Delta s` each iteration. Setting this input to false is typically safe, and is computationally more efficient.

 maximum_number_of_iterations : int, default = 50
     Maximum number of iterations taken by the algorithm. If this number of iterations is reached without convergence (as defined by ``absolute_tolerance`` input),
     the behaviour of the algorithm is defined by the ``failure_handling`` input.

 absolute_tolerance : float, default = nan
     Difference in :math:`t_{R}-t_{T}` between two consecutive iterations below which the algorithm is considered to be converged. Default value is nan, which means the default value is taken.
     The default value depends on the time representation used (1 ps for float; 1 fs for Time class)

 failure_handling : LightTimeFailureHandling, default = accept_without_warning
     Input defines behaviour when failing to converge within the required number of iterations. NOTE: the default value should be overridden for high-accuracy applications

 Returns
 -------
 :class:`LightTimeConvergenceCriteria`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LightTimeConvergenceCriteria` with the required settings.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the light_time_convergence_settings function
     from tudatpy.numerical_simulation.estimation_setup import observation

     # The light_time_convergence_settings function can be used with default inputs as just:
     light_time_convergence_settings = observation.light_time_convergence_settings()
     # A LightTimeConvergenceCriteria object is returned
     print(light_time_convergence_settings)

     # Users can also specify the following input arguments:
     # iterate_corrections, maximum_number_of_iterations, absolute_tolerance, failure_handling.
     # Let's set the failure_handling argument to LightTimeFailureHandling.print_warning_and_accept (default was LightTimeFailureHandling.accept_without_warning)
     light_time_convergence_settings = observation.light_time_convergence_settings(
         failure_handling = observation.LightTimeFailureHandling.print_warning_and_accept
     )
     # Again, a LightTimeConvergenceCriteria object is returned
     print(light_time_convergence_settings)



     )doc" );

    m.def( "one_way_range",
           &tom::oneWayRangeSettings,
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) =
                   std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) =
                   std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(

 Function for creating settings for a one-way range observable.

 Function for creating observation model settings of one-way range type observables, for a single link definition. The associated observation model creates
 a single-valued observable :math:`h_{_{\text{1-range}}}` as follows (in the unbiased case):

 .. math::
    h_{_{\text{1-range}}}(t_{R},t_{T})=|\mathbf{r}_{R}(t_{R})-\mathbf{r}_{T}(t_{T})| + \Delta s

 where :math:`\mathbf{r}_{R}`, :math:`\mathbf{r}_{T}`, :math:`t_{R}` and :math:`t_{T}` denote the position function of receiver and transmitter, and evaluation time
 of receiver and transmitter. The term :math:`\Delta s` denotes light-time corrections due to e.g relativistic, atmospheric effects (as defined by the ``light_time_correction_settings`` input).
 The transmission and reception time are related to the light-time :math:`T=t_{R}-t_{T}`, which is in turn related to the one-way range as :math:`T=h/c`
 As a result, the calculation of the one-way range (and light-time) requires the iterative solution of the equation:

 .. math::
    t_{R}-t_{T}=c\left(|\mathbf{r}_{R}(t_{R})-\mathbf{r}(t_{R})| + \Delta s\right)

 The method for the iterative solution is described in the :func:`light_time_convergence_settings` entry


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires the
     `transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

 light_time_correction_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LightTimeCorrectionSettings` ], default = list()
     List of corrections for the light-time that are to be used. Default is none, which will result
     in the signal being modelled as moving in a straight line with the speed of light

 bias_settings : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is None (unbiased observation)

 light_time_convergence_settings : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
     Settings for convergence of the light-time

 Returns
 -------
 :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` class defining the settings for the one-way observable.


 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the one_way_range function
     from tudatpy.numerical_simulation.estimation_setup import observation

     # Create Link Ends dictionary
     link_ends = dict()
     link_ends[observation.receiver] = observation.body_origin_link_end_id("Earth")
     link_ends[observation.transmitter] = observation.body_origin_link_end_id("Delfi-C3")

     # Create a Link Definition Object from link_ends dictionary. This will be the input to the function.
     Link_Definition_Object = observation.LinkDefinition(link_ends)

     # Create minimal ObservationSettings object (only required Link_Definition_Object argument is passed)
     # Note: other optional parameters (bias_settings, light_time_correction_settings,  light_time_convergence_settings) can be set
     observation_settings = observation.one_way_range(Link_Definition_Object)

     # Show that this returns an ObservationSettings object.
     print(observation_settings)



     )doc" );

    m.def( "two_way_range",
           &tom::twoWayRangeSimple,
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) =
                   std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) =
                   std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(

 Function for creating settings for a two-way range observable.

 Same as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_range`, with :math:`n=2`. This function is provided
 for convenience.


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires the
     `transmitter`, `retransmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined

 light_time_correction_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LightTimeCorrectionSettings` ], default = list()
     List of corrections for the light-time that are to be used for each constituent one-way range. Default is none, which will result
     in the signal being modelled as moving in a straight line with the speed of light

 bias_settings : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
     Note that only one bias setting is applied to the n-way observable.

 light_time_convergence_settings : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
     Settings for convergence of the light-time

 Returns
 -------
 :class:`~tudatpy.numerical_simulation.estimation_setup.observation.NWayRangeObservationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.NWayRangeObservationSettings` class.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the two_way_range function
     from tudatpy.numerical_simulation.estimation_setup import observation

     # two_way_range() takes a Link Definition Object as input to the function.
     # Note: as for this case, transmitter, retransmitter and receiver are required to define the Link Ends dictionary
     link_ends = dict()
     link_ends[observation.transmitter] = observation.body_origin_link_end_id("Earth")
     link_ends[observation.retransmitter] = observation.body_origin_link_end_id("Delfi-C3")
     link_ends[observation.receiver] = observation.body_origin_link_end_id("Earth")

     # Create the LinkDefinition object
     Link_Definition_Object = observation.LinkDefinition(link_ends)

     # Create minimal ObservationSettings object (only required Link_Definition_Object argument is passed)
     # Note: other optional parameters (bias_settings, light_time_correction_settings,  light_time_convergence_settings) can be set
     observation_settings = observation.two_way_range(Link_Definition_Object)

     # Show that two_way_range() returns an NWayRangeObservationSettings object.
     print(observation_settings)



     )doc" );

    m.def( "two_way_range_from_one_way_links",
           &tom::twoWayRange,
           py::arg( "one_way_range_settings" ),
           py::arg( "bias_settings" ) = nullptr,
           R"doc(

 Function for creating settings for a two-way range observable.

 Same as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_range_from_one_way_links`, with :math:`n=2`. This function is provided
 for convenience.


 Parameters
 ----------
 one_way_range_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` ]
     List of observation model settings of size two, with the first entry the one-way range settings for the uplink, and the second entry the one-way range settings for the downlink.
     The ``LinkDefinition`` of this two-way range observable is created from this list, with the ``transmitter`` and ``retransmitter`` defined by the
     ``transmitter`` and ``receiver`` of the first entry in this list. The ``retransmitter`` and ``receiver`` are defined by the
     ``transmitter`` and ``receiver`` of the second entry of this list.

 bias_settings : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
     Note that only one bias setting is applied to the n-way observable.

 Returns
 -------
 :class:`~tudatpy.numerical_simulation.estimation_setup.observation.NWayRangeObservationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.NWayRangeObservationSettings` class.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the two_way_range_from_one_way_links function
     from tudatpy.numerical_simulation.estimation_setup import observation

     # two_way_range_from_one_way_links() takes a list of ObservationSettings objects
     # Note: as for this case, transmitter, retransmitter and receiver are required to define the Link Ends dictionary
     link_ends = dict()
     link_ends[observation.transmitter] = observation.body_origin_link_end_id("Earth")
     link_ends[observation.retransmitter] = observation.body_origin_link_end_id("Delfi-C3")
     link_ends[observation.receiver] = observation.body_origin_link_end_id("Earth")

     # Create the LinkDefinition object to be used as input
     Link_Definition_Object = observation.LinkDefinition(link_ends) # define LinkDefinition object
     two_way_range_observation_settings_list = [observation.two_way_range(Link_Definition_Object)] # define (minimal) NWayRangeObservationSettings object
     two_way_range_one_way_link_settings = observation.two_way_range_from_one_way_links(two_way_range_observation_settings_list)

     # Show that two_way_range_from_one_way_links() returns an NWayRangeObservationSettings object.
     print(two_way_range_one_way_link_settings)



     )doc" );

    m.def( "n_way_range",
           &tom::nWayRangeSimple,
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) =
                   std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) =
                   std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(

 Function for creating settings for a n-way range observable.

 Function for creating observation model settings of n-way range type observables, for a single link definition. The associated observation model creates
 a single-valued observable :math:`h_{_{\text{N-range}}}` by combining together a series :math:`n` one-way range observations
 (see :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`). By default, the reception time of the :math:`i^{th}` one-way range is set as the
 transmission time of the :math:`(i+1)^{th}` one-way range. A retransmission delay may be defined by ancilliary settings (see :func:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationAncilliarySimulationSettings`) when creating observation
 simulation setings.

 For this function, the settings for each constituent one-way range (with the exception of the link end identifiers) are equal.


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires the
     `transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined, as well
     as a `retransmitter1`, `retransmitter2`, .... (with the number of retransmitters to be defined by the user).

 light_time_correction_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LightTimeCorrectionSettings` ], default = list()
     List of corrections for the light-time that are to be used for each constituent one-way range. Default is none, which will result
     in the signal being modelled as moving in a straight line with the speed of light

 bias_settings : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
     Note that only one bias setting is applied to the n-way observable.

 light_time_convergence_settings : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
     Settings for convergence of the light-time

 Returns
 -------
 :class:`~tudatpy.numerical_simulation.estimation_setup.observation.NWayRangeObservationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.NWayRangeObservationSettings` class.


 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the n_way_range function
     from tudatpy.numerical_simulation.estimation_setup import observation

     # Create Link Ends dictionary
     link_ends = dict()
     link_ends[observation.receiver] = observation.body_origin_link_end_id("Earth")
     link_ends[observation.transmitter] = observation.body_origin_link_end_id("Delfi-C3")

     # n_way_range() takes a Link Definition Object as input to the function.
     # Let's create it from link_ends
     Link_Definition_Object = observation.LinkDefinition(link_ends)

     # Create minimal ObservationSettings object (only required Link_Definition_Object argument is passed)
     # Note: other optional parameters (bias_settings, light_time_correction_settings,  light_time_convergence_settings) can be set
     observation_settings = observation.n_way_range(Link_Definition_Object)

     # Show that n_way_range() returns an NWayRangeObservationSettings object.
     print(observation_settings)





     )doc" );

    m.def( "n_way_range_from_one_way_links",
           &tom::nWayRange,
           py::arg( "one_way_range_settings" ),
           py::arg( "bias_settings" ) = nullptr,
           R"doc(

 Function for creating settings for a n-way range observable.

 Function for creating observation model settings of n-way range type observables, for a single link definition. The
 implementation is the same as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_range`, with the difference
 that the constituent one-way ranges may have different settings.s


 Parameters
 ----------
 one_way_range_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` ]
     List of observation model settings for each of the :math:`n` constituent one-way ranges of the n-way range observable.
     The ``LinkDefinition`` of this n-way range observable is created from this list, with the ``transmitter`` and ``retransmitter`` defined by the
     ``transmitter`` and ``receiver`` of the first entry in this list. The ``retransmitter`` (n-1) and ``receiver`` are defined by the
     ``transmitter`` and ``receiver`` of the :math:`\text{n}^{th}` entry of this list.

 bias_settings : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation).
     Note that only one bias setting is applied to the n-way observable.

 Returns
 -------
 :class:`~tudatpy.numerical_simulation.estimation_setup.observation.NWayRangeObservationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.NWayRangeObservationSettings` class.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the n_way_range_from_one_way_links function
     from tudatpy.numerical_simulation.estimation_setup import observation

     # Create Link Ends dictionary
     link_ends = dict()
     link_ends[observation.receiver] = observation.body_origin_link_end_id("Earth")
     link_ends[observation.transmitter] = observation.body_origin_link_end_id("Delfi-C3")

     # n_way_range_from_one_way_links() takes 1) a list of ObservationSettings objects and 2) bias as input (default is None)
     # Let's create it.
     Link_Definition_Object = observation.LinkDefinition(link_ends) # define LinkDefinition object
     n_way_observation_settings_list = [observation.n_way_range(Link_Definition_Object)] # define (minimal) ObservationSettings object

     n_way_from_one_link_observation_settings = observation.n_way_range_from_one_way_links(n_way_observation_settings_list, bias_settings = None)

     # Show that n_way_range_from_one_way_links() returns an NWayRangeObservationSettings object.
     print(n_way_from_one_link_observation_settings)



     )doc" );

    m.def( "angular_position",
           &tom::angularPositionSettings,
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) =
                   std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) =
                   std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(

 Function for creating settings for an angular position observable.

 Function for creating observation model settings of angular position type observables (as right ascension :math:`\alpha` and declination :math:`\delta`),
 for a single link definition. The associated observation model creates an observable :math:`\mathbf{h}_{_{\text{ang.pos.}}}` of type two as follows (in the unbiased case):

 .. math::
    \Delta\mathbf{r}=\mathbf{r}_{R}(t_{R})-\mathbf{r}_{T}(t_{T})\\
    \tan\alpha=\frac{\Delta r_{y}}{\Delta r_{x}}\\
    \delta=\frac{\Delta r_{z}}{\sqrt{\Delta r_{x}^{2}+\Delta r_{y}^{2}}}\\
    \mathbf{h}_{_{\text{ang.pos.}}} = [\alpha;\delta]

 The relative position vector :math:`\Delta\mathbf{r}` is computed identically as described for the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`
 The angular position observable can be used for optical astrometry, VLBI, etc. Due to the definition of this observable, the xy-plane is defined by the global frame orientation of the
 environment.


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires the
     `transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

 light_time_correction_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LightTimeCorrectionSettings` ], default = list()
     List of corrections for the light-time that are to be used. Default is none, which will result
     in the signal being modelled as moving in a straight line with the speed of light

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
     Settings for convergence of the light-time

 Returns
 -------
 :class:`ObservationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` class defining the settings for the angular position observable.






     )doc" );

    m.def( "relative_angular_position",
           &tom::relativeAngularPositionSettings,
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) =
                   std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) =
                   std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(

 Function for creating settings for an angular position observable.

 Function for creating observation model settings of angular position type observables (as right ascension :math:`\alpha` and declination :math:`\delta`),
 for a single link definition. The associated observation model creates an observable :math:`\mathbf{h}_{_{\text{ang.pos.}}}` of type two as follows (in the unbiased case):

 .. math::
    \Delta\mathbf{r}=\mathbf{r}_{R}(t_{R})-\mathbf{r}_{T}(t_{T})\\
    \tan\alpha=\frac{\Delta r_{y}}{\Delta r_{x}}\\
    \delta=\frac{\Delta r_{z}}{\sqrt{\Delta r_{x}^{2}+\Delta r_{y}^{2}}}\\
    \mathbf{h}_{_{\text{ang.pos.}}} = [\alpha;\delta]

 The relative position vector :math:`\Delta\mathbf{r}` is computed identically as described for the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`
 The angular position observable can be used for optical astrometry, VLBI, etc. Due to the definition of this observable, the xy-plane is defined by the global frame orientation of the
 environment.


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires the
     `transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

 light_time_correction_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LightTimeCorrectionSettings` ], default = list()
     List of corrections for the light-time that are to be used. Default is none, which will result
     in the signal being modelled as moving in a straight line with the speed of light

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
     Settings for convergence of the light-time

 Returns
 -------
 :class:`ObservationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` class defining the settings for the angular position observable.






     )doc" );

    m.def( "cartesian_position",
           &tom::positionObservableSettings,
           py::arg( "link_ends" ),
           py::arg( "bias_settings" ) = nullptr,
           R"doc(

 Function for creating settings for a Cartesian position observable.

 Function for creating observation model settings of Cartesian position type observables.
 Note that this observable is typically not realized in reality, but can be very useful for verification or analysis purposes.
 This observable provides the inertial (w.r.t. global frame origin) Cartesian position of the `observed_body` defined by the `link_ends` input.
 The observable has size 3, and contains the :math:`x`, :math:`y` and :math:`z` position


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires that the
     `observed_body`` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 Returns
 -------
 :class:`ObservationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` class defining the settings for the cartesian position observable.






     )doc" );

    m.def( "relative_cartesian_position",
           &tom::relativePositionObservableSettings,
           py::arg( "link_ends" ),
           py::arg( "bias_settings" ) = nullptr,
           R"doc(

 Function for creating settings for a Cartesian position observable.

 Function for creating observation model settings of Cartesian position type observables.
 Note that this observable is typically not realized in reality, but can be very useful for verification or analysis purposes.
 This observable provides the inertial (w.r.t. global frame origin) Cartesian position of the `observed_body` defined by the `link_ends` input.
 The observable has size 3, and contains the :math:`x`, :math:`y` and :math:`z` position


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires that the
     `observed_body`` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 Returns
 -------
 :class:`ObservationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` class defining the settings for the cartesian position observable.






     )doc" );

    m.def( "cartesian_velocity",
           &tom::velocityObservableSettings,
           py::arg( "link_ends" ),
           py::arg( "bias_settings" ) = nullptr,
           R"doc(

 Function for creating settings for a Cartesian velocity observable.

 Function for creating observation model settings of Cartesian position type observables.
 Note that this observable is typically not realized in reality, but can be very useful for verification or analysis purposes.
 This observable provides the inertial (w.r.t. global frame origin) Cartesian velocity of the `observed_body` defined by the `link_ends` input.
 The observable has size 3, and contains the :math:`x`, :math:`y` and :math:`z` velocity


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires that the
     `observed_body`` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 Returns
 -------
 :class:`ObservationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` class defining the settings for the cartesian velocity observable.






     )doc" );

    m.def( "euler_angles_313",
           &tom::eulerAngle313ObservableSettings,
           py::arg( "link_ends" ),
           py::arg( "bias_settings" ) = nullptr,
           R"doc(No documentation found.)doc" );

    m.def( "one_way_doppler_instantaneous",
           &tom::oneWayOpenLoopDoppler,
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) =
                   std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "transmitter_proper_time_rate_settings" ) = nullptr,
           py::arg( "receiver_proper_time_rate_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) =
                   std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           py::arg( "normalized_with_speed_of_light" ) = false,
           R"doc(

 Function for creating settings for a one-way instantaneous Doppler observable.

 Function for creating settings for a one-way instantaneous Doppler observable for a single link definition. The associated observation model creates
 a single-valued observable :math:`h_{_{\text{1-Dopp.}}}` as follows (in the unbiased case):

 .. math::
    h_{_{\text{1-Dopp.}}}=c\left(\frac{d\tau_{T}}{dt_{T}}\frac{t_{T}}{dt_{R}}\frac{dt_{R}}{d\tau_{R}}-1\right)

 where :math:`t` and :math:`\tau` denote coordinate and proper time of the transmitter T and receiver R, respectively.
 The receiver and transmitter position and coordinate time are computed identically as described for the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`.
 The detailed mathematical implementation are described in: `Moyer, T.D. (2000) Formulation for Observed and Computed Values of Deep Space Network Data Types for Navigation. Monograph 2, Deep Space Communications and Navigation Series, JPL Publication 00-7 <https://www.scirp.org/reference/referencespapers?referenceid=1827210>`_.

 This observable represents the 'instantaneous (non-integrated)' Doppler observable, as obtained from open-loop observations.
 It should *not* be used for the modelling of the typical closed-loop observations used in deep space tracking (for which the
 :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_doppler_averaged` should be used)

 The coordinate
 time derivative :math:`\frac{t_{A}}{dt_{B}}` is always computed when generating this observable. Settings for the proper time
 rates :math:`\frac{d\tau}{dt}` can be specified by the user through the ``transmitter_proper_time_rate_settings`` and ``receiver_proper_time_rate_settings``
 arguments (inputs, see Parameters). Whenever these are left empty, the proper time rates are omitted (set to 1.0).

 The observable may be non-dimensionalized by the speed of light :math:`c`, which results in the observable being equal to thee received and transmitted signal frequencies :math:`f_{R}/f_{T}-1`.


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires that the
     `transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

 light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
     List of corrections for the light-time that are to be used. Default is none, which will result
     in the signal being modelled as moving in a straight line with the speed of light

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 transmitter_proper_time_rate_settings : :class:`DopplerProperTimeRateSettings`, default = None
     Settings for computing the transmitter proper time rate :math:`\frac{d\tau}{dt}`, default is none (:math:`\frac{d\tau}{dt}=1`)

 receiver_proper_time_rate_settings : :class:`DopplerProperTimeRateSettings`, default = None
     Settings for computing the receiver proper time rate :math:`\frac{d\tau}{dt}`, default is none (:math:`\frac{d\tau}{dt}=1`)

 light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
     Settings for convergence of the light-time

 normalized_with_speed_of_light : bool, default = false
     Option to non-dimensionalize the observable with speed of light :math:`c`

 Returns
 -------
 :class:`OneWayDopplerObservationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` derived :class:`OneWayDopplerObservationSettings` class defining the settings for the one-way open doppler observable observable.






     )doc" );

    m.def( "two_way_doppler_instantaneous_from_one_way_links",
           py::overload_cast< const std::shared_ptr< tom::OneWayDopplerObservationSettings >,
                              const std::shared_ptr< tom::OneWayDopplerObservationSettings >,
                              const std::shared_ptr< tom::ObservationBiasSettings > >(
                   &tom::twoWayOpenLoopDoppler ),
           py::arg( "uplink_doppler_settings" ),
           py::arg( "downlink_doppler_settings" ),
           py::arg( "bias_settings" ) = nullptr,
           R"doc(

 Function for creating settings for a two-way instantaneous Doppler observable.


 Function for creating settings for a two-way instantaneous Doppler observable for a single link definition. The
 implementation is the same as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.two_way_doppler_instantaneous`, with the difference
 that the constituent one-way ranges may have different settings.

 The observable may be non-dimensionalized by the speed of light :math:`c` (in the constituent one-way Doppler observable settings),
 which results in the observable being equal to the received and transmitted signal frequencies :math:`f_{R}/f_{T}-1`.


 Parameters
 ----------
 uplink_doppler_settings : :class:`OneWayDopplerObservationSettings`
     Settings for uplink leg of one-way observable, created using :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_open_loop_doppler`

 downlink_doppler_settings : :class:`OneWayDopplerObservationSettings`
     Settings for downlink leg of one-way observable, created using :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_open_loop_doppler`

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the full observation, default is none (unbiased observation). Note that,
     even if no bias is applied to the two-way observable, the constituent one-way observables may still be biased.

 light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
     Settings for convergence of the light-time

 Returns
 -------
 :class:`TwoWayDopplerObservationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` derived :class:`TwoWayDopplerObservationSettings` class defining the settings for the two-way open doppler observable.






     )doc" );

    m.def( "two_doppler_instantaneous",
           py::overload_cast<
                   const tom::LinkDefinition&,
                   const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >&,
                   const std::shared_ptr< tom::ObservationBiasSettings >,
                   const std::shared_ptr< tom::LightTimeConvergenceCriteria >,
                   const bool >( &tom::twoWayOpenLoopDoppler ),
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) =
                   std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) =
                   std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           py::arg( "normalized_with_speed_of_light" ) = false,
           R"doc(No documentation found.)doc" );

    m.def( "one_way_doppler_averaged",
           py::overload_cast<
                   const tom::LinkDefinition&,
                   const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >,
                   const std::shared_ptr< tom::ObservationBiasSettings >,
                   const std::shared_ptr< tom::LightTimeConvergenceCriteria > >(
                   &tom::oneWayClosedLoopDoppler ),
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) =
                   std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) =
                   std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(

 Function for creating settings for a one-way averaged Doppler observable.

 Function for creating observation model settings for one-way averaged Doppler observables, for a single link definition. The associated observation model creates
 a single-valued observable :math:`h_{_{\text{1-\bar{Dopp}}}}` as follows (in the unbiased case):

 .. math::
    h_{_{\text{1-\bar{Dopp}}}}&=c\int_{t-\Delta t}^{t+\Delta t}\frac{t_{T}}{dt_{R}}d\bar{t}\\
                              &=\frac{h_{_{\text{1-range}}}(t_{R}=t+\Delta t,t_{T})-h_{_{\text{1-range}}}(t_{R}=t,t_{T})}{\Delta t}

 where, in the latter formulation (which is the one that is implemented), the observable is referenced to the receiver time. This averaged Doppler observable
 is computed as the difference of two one-way range observables (see :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_range`),
 with the reference time shifted by :math:`\Delta t`. As such, it is sensitive to numerical errors for small :math:`\Delta t`

 The integration time :math:`\Delta t` is defined in the ancilliary settings when simulating the observations (with 60 s as default).


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires that the
     `transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined.

 light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
     List of corrections for the light-time that are to be used. Default is none, which will result
     in the signal being modelled as moving in a straight line with the speed of light

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
     Settings for convergence of the light-time

 Returns
 -------
 :class:`ObservationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` derived `OneWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.






     )doc" );

    m.def( "two_way_doppler_averaged",
           py::overload_cast<
                   const tom::LinkDefinition&,
                   const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >,
                   const std::shared_ptr< tom::ObservationBiasSettings >,
                   const std::shared_ptr< tom::LightTimeConvergenceCriteria > >(
                   &tom::twoWayDifferencedRangeObservationSettings ),
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) =
                   std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) =
                   std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(

 Function for creating settings for an n-way averaged Doppler observable.

 Function for creating observation model settings for n-way averaged Doppler observables, for a single link definition. The implementation is
 analogous to the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_doppler_averaged` observable. But, in the present case
 the observable is computed from the difference of two n-way range observables, with the reference time shifted by :math:`\Delta t`.

 The integration time :math:`\Delta t` is defined in the ancilliary settings when simulating the observations (with 60 s as default).


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires the
     `transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined, as well
     as a `retransmitter1`, `retransmitter2`, .... (with the number of retransmitters to be defined by the user).

 light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
     List of corrections for the light-time that are to be used. Default is none, which will result
     in the signal being modelled as moving in a straight line with the speed of light

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
     Settings for convergence of the light-time

 Returns
 -------
 :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` derived `~tudatpy.numerical_simulation.estimation_setup.observation.NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.






     )doc" );

    m.def( "two_way_doppler_averaged_from_one_way_links",
           py::overload_cast< const std::vector< std::shared_ptr< tom::ObservationModelSettings > >,
                              const std::shared_ptr< tom::ObservationBiasSettings > >(
                   &tom::twoWayDifferencedRangeObservationSettings ),
           py::arg( "one_way_range_settings" ),
           py::arg( "bias_settings" ) = nullptr,
           R"doc(

 Function for creating settings for an n-way averaged Doppler observable.

 Function for creating observation model settings for n-way averaged Doppler observables, for a single link definition. The implemenation is
 analogous to the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_doppler_averaged` observable. But, in the present case
 the observable is computed from the difference of two n-way range observables, with the reference time shifted by :math:`\Delta t`.

 The integration time :math:`\Delta t` is defined in the ancilliary settings when simulating the observations (with 60 s as default).


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires the
     `transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined, as well
     as a `retransmitter1`, `retransmitter2`, .... (with the number of retransmitters to be defined by the user).

 light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
     List of corrections for the light-time that are to be used. Default is none, which will result
     in the signal being modelled as moving in a straight line with the speed of light

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
     Settings for convergence of the light-time

 Returns
 -------
 :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` derived `~tudatpy.numerical_simulation.estimation_setup.observation.NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.






     )doc" );

    m.def( "n_way_doppler_averaged",
           py::overload_cast<
                   const tom::LinkDefinition&,
                   const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >,
                   const std::shared_ptr< tom::ObservationBiasSettings >,
                   const std::shared_ptr< tom::LightTimeConvergenceCriteria > >(
                   &tom::nWayDifferencedRangeObservationSettings ),
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) =
                   std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) =
                   std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(

 Function for creating settings for an n-way averaged Doppler observable.

 Function for creating observation model settings for n-way averaged Doppler observables, for a single link definition. The implementation is
 analogous to the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.one_way_doppler_averaged` observable. But, in the present case
 the observable is computed from the difference of two n-way range observables, with the reference time shifted by :math:`\Delta t`.

 The integration time :math:`\Delta t` is defined in the ancilliary settings when simulating the observations (with 60 s as default).


 Parameters
 ----------
 link_ends : LinkDefinition
     Set of link ends that define the geometry of the observation. This observable requires the
     `transmitter` and `receiver` :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndType` entries to be defined, as well
     as a `retransmitter1`, `retransmitter2`, .... (with the number of retransmitters to be defined by the user).

 light_time_correction_settings : List[ :class:`LightTimeCorrectionSettings` ], default = list()
     List of corrections for the light-time that are to be used. Default is none, which will result
     in the signal being modelled as moving in a straight line with the speed of light

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 light_time_convergence_settings : :class:`LightTimeConvergenceCriteria`, default = :func:`light_time_convergence_settings`
     Settings for convergence of the light-time

 Returns
 -------
 :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` derived `~tudatpy.numerical_simulation.estimation_setup.observation.NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.






     )doc" );

    m.def( "n_way_doppler_averaged_from_one_way_links",
           py::overload_cast< const std::vector< std::shared_ptr< tom::ObservationModelSettings > >,
                              const std::shared_ptr< tom::ObservationBiasSettings >,
                              const std::shared_ptr< tom::LightTimeConvergenceCriteria > >(
                   &tom::nWayDifferencedRangeObservationSettings ),
           py::arg( "one_way_range_settings" ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) =
                   std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(

 Function for creating settings for an n-way averaged Doppler observable.

 Function for creating observation model settings for n-way averaged Doppler observables, for a single link definition.
 The implementation is the same as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_doppler_averaged`, with the difference
 that the constituent one-way range observables may have different settings.


 Parameters
 ----------
 one_way_range_settings : List[ :class:`ObservationSettings` ]
     List of observation model settings for each of the :math:`n` constituent one-way ranges of the n-way averaged range rate observable.
     The ``LinkDefinition`` of this n-way range observable is created from this list, with the ``transmitter`` and ``retransmitter`` defined by the
     ``transmitter`` and ``receiver`` of the first entry in this list. The ``retransmitter`` (n-1) and ``receiver`` are defined by the
     ``transmitter`` and ``receiver`` of the :math:`\text{n}^{th}` entry of this list.

 bias_settings : :class:`ObservationBiasSettings`, default = None
     Settings for the observation bias that is to be used for the observation, default is none (unbiased observation)

 Returns
 -------
 :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` derived `~tudatpy.numerical_simulation.estimation_setup.observation.NWayDifferencedRangeRateObservationSettings` class defining the settings for the one-way closed-loop doppler observable.






     )doc" );

    m.def( "dsn_n_way_doppler_averaged",
           py::overload_cast<
                   const tom::LinkDefinition&,
                   const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >,
                   const std::shared_ptr< tom::ObservationBiasSettings >,
                   const std::shared_ptr< tom::LightTimeConvergenceCriteria >,
                   const bool >( &tom::dsnNWayAveragedDopplerObservationSettings ),
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) =
                   std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) =
                   std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           py::arg( "subtract_doppler_signature" ) = true,
           R"doc(No documentation found.)doc" );

    m.def( "dsn_n_way_doppler_averaged_from_one_way_links",
           py::overload_cast< const std::vector< std::shared_ptr< tom::ObservationModelSettings > >,
                              const std::shared_ptr< tom::ObservationBiasSettings >,
                              const std::shared_ptr< tom::LightTimeConvergenceCriteria >,
                              const bool >( &tom::dsnNWayAveragedDopplerObservationSettings ),
           py::arg( "one_way_range_settings" ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) =
                   std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           py::arg( "subtract_doppler_signature" ) = true,
           R"doc(No documentation found.)doc" );

    m.def( "dsn_n_way_Range",
           py::overload_cast<
                   const tom::LinkDefinition&,
                   const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >,
                   const std::shared_ptr< tom::ObservationBiasSettings >,
                   const std::shared_ptr< tom::LightTimeConvergenceCriteria > >(
                   &tom::dsnNWayRangeObservationSettings ),
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) =
                   std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) =
                   std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(No documentation found.)doc" );

    m.def( "doppler_measured_frequency",
           py::overload_cast<
                   const tom::LinkDefinition&,
                   const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >&,
                   const std::shared_ptr< tom::ObservationBiasSettings >,
                   const std::shared_ptr< tom::LightTimeConvergenceCriteria > >(
                   &tom::dopplerMeasuredFrequencyObservationSettings ),
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) =
                   std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) =
                   std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           R"doc(No documentation found.)doc" );

    m.def( "observation_settings_from_collection",
           &tss::getObservationSimulationSettingsFromObservations< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observation_collection" ),
           py::arg( "bodies" ),
           R"doc(No documentation found.)doc" );

    py::class_< tom::LightTimeCorrectionSettings,
                std::shared_ptr< tom::LightTimeCorrectionSettings > >(
            m,
            "LightTimeCorrectionSettings",
            R"doc(

         Base class to define light time correction settings.

         Base class to define light time correction settings.
         This class is not used for calculations of corrections, but is used for the purpose of defining the light time correction properties.
         Specific light time correction settings must be defined using an object derived from this class.

         Instances of this class are typically created via the
         :func:`~tudatpy.numerical_simulation.estimation_setup.observation.first_order_relativistic_light_time_correction` function

         Examples
         --------
         .. code-block:: python

             # Code snippet to show the creation of a LightTimeCorrectionSettings object
             from tudatpy.numerical_simulation.estimation_setup import observation

             # Create Link Ends dictionary
             link_ends = dict()
             link_ends[observation.receiver] = observation.body_origin_link_end_id("Earth")
             link_ends[observation.transmitter] = observation.body_origin_link_end_id("Delfi-C3")

             # Create a Link Definition Object from link_ends dictionary
             Link_Definition_Object = observation.LinkDefinition(link_ends)

             # Case 1: perturbing body (Earth) involved in the observations
             # In this case, Earth is a receiver, so the body’s state will be evaluated at the reception time.
             perturbing_body = ['Earth']
             doppler_observation_settings = observation.first_order_relativistic_light_time_correction(perturbing_body)

             # Show that it is a LightTimeCorrectionSettings object.
             print(doppler_observation_settings)

             # Case 2: perturbing body (Sun) not involved in the observations
             # In this case, the body's state will be evaluated at the midpoint time between the transmission and reception events.
             perturbing_body = ['Sun']

             # Use: observation.first_order_relativistic_light_time_correction to create a LightTimeCorrectionSettings object
             # Note: first_order_relativistic_light_time_correction only requires the perturbing list of bodies to be passed as arguments
             doppler_observation_settings = observation.first_order_relativistic_light_time_correction(perturbing_body)

             # Show that it is an LightTimeCorrectionSettings object.
             print(doppler_observation_settings.transmitter_proper_time_rate_settings)
             print(dir(doppler_observation_settings))



      )doc" );

    m.def( "first_order_relativistic_light_time_correction",
           &tom::firstOrderRelativisticLightTimeCorrectionSettings,
           py::arg( "perturbing_bodies" ),
           R"doc(

 Function for creating settings for first-order relativistic light-time corrections.

 Function for creating settings for first-order relativistic light-time corrections:  These corrections account for the delay in light travel time caused by stationary point masses, calculated up to
 :math:`c^{-2}` according to general relativity (e.g., Moyer, 2000). A key consideration in the model is the time at which the states of the perturbing bodies are evaluated. This depends on their involvement in the observation link ends:

 * 1. **Perturbing Body as a Link End:**
 If the perturbing body (e.g., Earth) is directly involved in the observation (e.g., as the location of a transmitter or receiver):

     - The body's state is evaluated at the **transmission time** if it acts as the **transmitter**.

     - The body's state is evaluated at the **reception time** if it acts as the **receiver**.

 * 2. **Perturbing Body Not as a Link End:**
 If the perturbing body is not part of the observation link ends, its state is evaluated at the **midpoint time** between the transmission and reception events.

 Parameters
 ----------
 perturbing_bodies : List[str]
     A list containing the names of the bodies due to which the light-time correction is to be taken into account.

 Returns
 -------
 :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LightTimeCorrectionSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LightTimeCorrectionSettings` configured to include
     first-order relativistic light-time corrections.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the first_order_relativistic_light_time_correction function
     from tudatpy.numerical_simulation.estimation_setup import observation

     # Create Link Ends dictionary
     link_ends = dict()
     link_ends[observation.receiver] = observation.body_origin_link_end_id("Earth")
     link_ends[observation.transmitter] = observation.body_origin_link_end_id("Delfi-C3")

     # Create a Link Definition Object from link_ends dictionary
     Link_Definition_Object = observation.LinkDefinition(link_ends)

     # The function first_order_relativistic_light_time_correction() requires a list of strings (perturbing body/bodies) as input
     perturbing_body = ['Earth']
     doppler_observation_settings = observation.first_order_relativistic_light_time_correction(perturbing_body)

     # Show that it returns a LightTimeCorrectionSettings object.
     print(doppler_observation_settings)

     )doc" );

    py::enum_< tom::TroposphericMappingModel >(
            m, "TroposphericMappingModel", R"doc(No documentation found.)doc" )
            .value( "simplified_chao", tom::TroposphericMappingModel::simplified_chao )
            .value( "niell", tom::TroposphericMappingModel::niell )
            .export_values( );

    py::enum_< tom::WaterVaporPartialPressureModel >(
            m, "WaterVaporPartialPressureModel", R"doc(No documentation found.)doc" )
            .value( "tabulated", tom::WaterVaporPartialPressureModel::tabulated )
            .value( "bean_and_dutton", tom::WaterVaporPartialPressureModel::bean_and_dutton )
            .export_values( );

    m.def( "dsn_tabulated_tropospheric_light_time_correction",
           &tom::tabulatedTroposphericCorrectionSettings,
           py::arg( "file_names" ),
           py::arg( "body_with_atmosphere_name" ) = "Earth",
           py::arg( "mapping_model" ) = tom::TroposphericMappingModel::niell,
           R"doc(No documentation found.)doc" );

    m.def( "saastamoinen_tropospheric_light_time_correction",
           &tom::saastamoinenTroposphericCorrectionSettings,
           py::arg( "body_with_atmosphere_name" ) = "Earth",
           py::arg( "mapping_model" ) = tom::TroposphericMappingModel::niell,
           py::arg( "water_vapor_partial_pressure_model" ) =
                   tom::WaterVaporPartialPressureModel::tabulated,
           R"doc(No documentation found.)doc" );

    m.def( "dsn_tabulated_ionospheric_light_time_correction",
           &tom::tabulatedIonosphericCorrectionSettings,
           py::arg( "file_names" ),
           py::arg( "spacecraft_name_per_id" ) = std::map< int, std::string >( ),
           py::arg( "quasar_name_per_id" ) = std::map< int, std::string >( ),
           py::arg( "reference_frequency" ) = 2295e6,
           py::arg( "body_with_atmosphere_name" ) = "Earth",
           R"doc(No documentation found.)doc" );

    m.def( "jakowski_ionospheric_light_time_correction",
           &tom::jakowskiIonosphericCorrectionSettings,
           py::arg( "ionosphere_height" ) = 400.0e3,
           py::arg( "first_order_delay_coefficient" ) = 40.3,
           py::arg( "solar_activity_data" ) =
                   tudat::input_output::solar_activity::readSolarActivityData(
                           tudat::paths::getSpaceWeatherDataPath( ) + "/sw19571001.txt" ),
           py::arg( "geomagnetic_pole_latitude" ) = tuc::convertDegreesToRadians( 80.9 ),
           py::arg( "geomagnetic_pole_longitude" ) = tuc::convertDegreesToRadians( -72.6 ),
           py::arg( "use_utc_for_local_time_computation" ) = false,
           py::arg( "body_with_atmosphere_name" ) = "Earth",
           R"doc(No documentation found.)doc" );

    m.def( "inverse_power_series_solar_corona_light_time_"
           "correction",
           &tom::inversePowerSeriesSolarCoronaCorrectionSettings,
           py::arg( "coefficients" ) = std::vector< double >{ 1.3e14, 0.5e12 },
           py::arg( "positive_exponents" ) = std::vector< double >{ 6.0, 2.0 },
           py::arg( "delay_coefficient" ) = 40.3,
           py::arg( "sun_body_name" ) = "Sun",
           R"doc(No documentation found.)doc" );

    py::class_< tom::ObservationBiasSettings, std::shared_ptr< tom::ObservationBiasSettings > >(
            m, "ObservationBiasSettings", R"doc(

         Base class to defining observation bias settings.

         Base class to defining observation bias settings.
         Specific observation bias settings must be defined using an object derived from this class.
         Instances of this class are typically created via the
         :func:`~tudatpy.numerical_simulation.estimation_setup.observation.absolute_bias` or :func:`~tudatpy.numerical_simulation.estimation_setup.observation.relative_bias` function.


         Examples
         --------
         .. code-block:: python
             # Code snippet to show the creation of an ObservationBiasSettings object
             # using absolute and relative bias settings
             from tudatpy.numerical_simulation.estimation_setup import observation
             import numpy as np

             bias_array = np.array([1e-2])

             # Use absolute_bias function
             absolute_bias_settings = observation.absolute_bias(bias_array)
             # Show that it is an ObservationBiasSettings object.
             print(absolute_bias_settings)

             # Use relative_bias function
             relative_bias_settings = observation.relative_bias(bias_array)
             # Show that it is an ObservationBiasSettings object.
             print(relative_bias_settings)
      )doc" );

    m.def( "clock_induced_bias",
           &tom::clockInducedBias,
           py::arg( "body_name" ),
           py::arg( "station_name" ),
           R"doc(No documentation found.)doc" );

    m.def( "absolute_bias",
           &tom::constantAbsoluteBias,
           py::arg( "bias_value" ),
           R"doc(

 Function for creating settings for an absolute observation bias.

 Function for creating settings for an absolute observation bias. When calculating the observable value, applying this setting
 will take the physically ideal observation :math:`h`, and modify it to obtain the biased observation :math:`\tilde{h}` as follows:

 .. math::
    \tilde{h}=h+K

 where :math:`K` is the `bias_value`. For an observable with size greater than 1, :math:`K` is a vector and the addition is component-wise.


 Parameters
 ----------
 bias_value : numpy.ndarray
     A vector containing the bias that is to be applied to the observable. This vector should be the same size as the observable to which it is
     applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

 Returns
 -------
 :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ConstantObservationBiasSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` defining the settings for a constant, absolute observation bias.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the absolute_bias function
     from tudatpy.numerical_simulation.estimation_setup import observation
     import numpy as np

     # The function absolute_bias() requires a numpy.array of bias values in input
     bias_array = np.array([1e-2])
     absolute_bias_settings = observation.absolute_bias(bias_array)

     # Show that it returns an ObservationBiasSettings object.
     print(absolute_bias_settings)



     )doc" );

    m.def( "relative_bias",
           &tom::constantRelativeBias,
           py::arg( "bias_value" ),
           R"doc(

 Function for creating settings for a relative observation bias.

 Function for creating settings for a relative observation bias. When calculating the observable value, applying this setting
 will take the physically ideal observation :math:`h`, and modify it to obtain the biased observation :math:`\tilde{h}` as follows:

 .. math::
    \tilde{h}=h(1+K)

 where :math:`K` is the`bias_value`. For an observable with size greater than 1, :math:`K` is a vector and the multiplication is component-wise.


 Parameters
 ----------
 bias_value : numpy.ndarray
     A vector containing the bias that is to be applied to the observable. This vector should be the same size as the observable to which it is
     applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

 Returns
 -------
 :class:`ConstantObservationBiasSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` class,
     defining the settings for a constant, relative observation bias.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the relative_bias function
     from tudatpy.numerical_simulation.estimation_setup import observation
     import numpy as np

     # The function relative_bias() requires a numpy.array of bias values in input
     bias_array = np.array([1e-2])
     relative_bias_settings_settings = observation.relative_bias(bias_array)

     # Show that it returns an ObservationBiasSettings object.
     print(relative_bias_settings_settings)




     )doc" );

    m.def( "arcwise_absolute_bias",
           py::overload_cast< const std::vector< double >&,
                              const std::vector< Eigen::VectorXd >&,
                              const tom::LinkEndType >( &tom::arcWiseAbsoluteBias ),
           py::arg( "arc_start_times" ),
           py::arg( "bias_values" ),
           py::arg( "reference_link_end_type" ),
           R"doc(

 Function for creating settings for arc-wise absolute observation biases.

 Function for creating settings for arc-wise absolute observation biases.
 This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.absolute_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.


 Parameters
 ----------
 arc_start_times : List[ float ]
     List containing starting times for each arc.

 bias_values : List[ numpy.ndarray ]
     List of arc-wise bias vectors that are to be applied to the given observable. The vectors should be the same size as the observable to which it is
     applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

 reference_link_end_type : :class:`LinkEndType`
     Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.

 Returns
 -------
 :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the arcwise_absolute_bias function
     from tudatpy.numerical_simulation.estimation_setup import observation
     import numpy as np

     # The function arcwise_absolute_bias() requires:
     # 1) an arc_start_times list ,2) a numpy.array of bias values in input, 3) reference_link_end_type
     # Let's simulate two arcs
     arc_start_times = [0, 60] # define start time in seconds
     arcwise_bias_array = [np.array([1e-2]), np.array([2e-2])] # set arc bias
     reference_link_end_type = observation.receiver # set bias at receiving link end
     arcwise_absolute_bias_settings = observation.arcwise_absolute_bias(arc_start_times, arcwise_bias_array, observation.receiver)

     # Show that it returns an ObservationBiasSettings object.
     print(arcwise_absolute_bias_settings)




     )doc" );

    m.def( "arcwise_absolute_bias_per_time",
           py::overload_cast< const std::map< double, Eigen::VectorXd >&, const tom::LinkEndType >(
                   &tom::arcWiseAbsoluteBias ),
           py::arg( "bias_values_per_start_time" ),
           py::arg( "reference_link_end_type" ),
           R"doc(

 Function for creating settings for arc-wise absolute observation biases.

 Function for creating settings for arc-wise absolute observation biases.
 This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.absolute_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.


 Parameters
 ----------
 bias_values_per_start_time : Dict[float, numpy.ndarray[numpy.float64[m, 1]]]
     Dictionary, in which the bias value vectors for each arc are directly mapped to the starting times of the respective arc.
     The vectors should be the same size as the observable to which it is applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

 reference_link_end_type : :class:`LinkEndType`
     Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.

 Returns
 -------
 :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.

 Examples
 --------
 .. code-block:: python

 # Code Snippet to showcase the use of the arcwise_absolute_bias function
 from tudatpy.numerical_simulation.estimation_setup import observation
 import numpy as np

     # The function arcwise_absolute_bias_settings_per_time() requires:
     # 1) a dictionary with times as keys and bias values as values ,2) a reference_link_end_type
     # Let's simulate two arcs
     bias_value_per_start_time = dict()
     bias_value_per_start_time[0] = np.array([1e-2])
     bias_value_per_start_time[60] = np.array([2e-2])
     reference_link_end_type = observation.receiver # set bias at receiving link end

     arcwise_absolute_bias_settings_per_time = observation.arcwise_absolute_bias_per_time(bias_value_per_start_time, reference_link_end_type)

     # Show that it returns an ObservationBiasSettings object.
     print(arcwise_absolute_bias_settings_per_time)



     )doc" );

    m.def( "arcwise_relative_bias",
           py::overload_cast< const std::vector< double >&,
                              const std::vector< Eigen::VectorXd >&,
                              const tom::LinkEndType >( &tom::arcWiseRelativeBias ),
           py::arg( "arc_start_times" ),
           py::arg( "bias_values" ),
           py::arg( "reference_link_end_type" ),
           R"doc(

 Function for creating settings for arc-wise relative observation biases.

 Function for creating settings for arc-wise relative observation biases.
 This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.relative_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.


 Parameters
 ----------
 arc_start_times : List[ float ]
     List containing starting times for each arc.

 bias_values : List[ numpy.ndarray ]
     List of arc-wise bias vectors that are to be applied to the given observable. The vectors should be the same size as the observable to which it is
     applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

 reference_link_end_type : :class:`LinkEndType`
     Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.

 Returns
 -------
 :class:`ArcWiseConstantObservationBiasSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the arcwise_relative_bias function
     from tudatpy.numerical_simulation.estimation_setup import observation
     import numpy as np

     # The function arcwise_relative_bias() requires:
     # 1) an arc_start_times list ,2) a numpy.array of bias values in input, 3) reference_link_end_type
     # Let's simulate two arcs
     arc_start_times = [0, 60] # define start time in seconds
     arcwise_bias_array = [np.array([1e-2]), np.array([2e-2])] # set arc bias
     reference_link_end_type = observation.receiver # set bias at receiving link end
     arcwise_relative_bias_settings = observation.arcwise_relative_bias(arc_start_times, arcwise_bias_array, reference_link_end_type)

     # Show that it returns an ObservationBiasSettings object.
     print(arcwise_relative_bias_settings)




     )doc" );

    m.def( "arcwise_relative_bias_per_time",
           py::overload_cast< const std::map< double, Eigen::VectorXd >&, const tom::LinkEndType >(
                   &tom::arcWiseRelativeBias ),
           py::arg( "bias_values_per_start_time" ),
           py::arg( "reference_link_end_type" ),
           R"doc(

 Function for creating settings for arc-wise relative observation biases.

 Function for creating settings for arc-wise relative observation biases.
 This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.relative_bias` setting only through the option of setting the `bias_value` :math:`K` to a different values for each arc.


 Parameters
 ----------
 bias_values_per_start_time : Dict[float, numpy.ndarray[numpy.float64[m, 1]]]
     Dictionary, in which the bias value vectors for each arc are directly mapped to the starting times of the respective arc.
     The vectors should be the same size as the observable to which it is applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

 reference_link_end_type : :class:`LinkEndType`
     Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.

 Returns
 -------
 :class:`ArcWiseConstantObservationBiasSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the arcwise_relative_bias function
     from tudatpy.numerical_simulation.estimation_setup import observation
     import numpy as np

     # The function arcwise_relative_bias_per_time() requires:
     # 1) a dictionary with times as keys and bias values as values ,2) a reference_link_end_type
     # Let's simulate two arcs
     bias_value_per_start_time = dict()
     bias_value_per_start_time[0] = np.array([1e-2])
     bias_value_per_start_time[60] = np.array([2e-2])
     reference_link_end_type = observation.receiver # set bias at receiving link end

     arcwise_relative_bias_settings_per_time = observation.arcwise_relative_bias_per_time(bias_value_per_start_time, reference_link_end_type)

     # Show that it returns an ObservationBiasSettings object.
     print(arcwise_relative_bias_settings_per_time)





     )doc" );

    m.def( "time_drift_bias",
           &tom::constantTimeDriftBias,
           py::arg( "bias_value" ),
           py::arg( "time_link_end" ),
           py::arg( "ref_epoch" ),
           R"doc(

 Function for creating settings for a time-drift bias.

 Function for creating settings for a time-drift bias.
 This bias setting generates the configuration for applying a constant time-drift bias to an observation model.

 Parameters
 ----------
 bias_value : numpy.ndarray
     Constant time drift bias that is to be considered for the observation time. This vector should be the same size as the observable to which it is
     assigned (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

 time_link_end : :class:`LinkEndType`
     Defines the link end (via the :class:`LinkEndType`) which is used the current time.

 ref_epoch : float
     Defines the reference epoch at which the effect of the time drift is initialised.

 Returns
 -------
 :class:`constantTimeDriftBias`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.constantTimeDriftBias` class,
     defining the settings for a constant, relative observation bias.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the time_drift_bias function
     from tudatpy.numerical_simulation.estimation_setup import observation
     import numpy as np

     # The function time_drift_bias() requires a numpy.array of bias value, time_link_end and ref_epoch as inputs
     bias_array = np.array([1e-2])
     time_link_end = observation.receiver
     ref_epoch = 0
     time_drift_bias_settings = observation.time_drift_bias(bias_array, time_link_end, ref_epoch)

     # Show that it returns an ObservationBiasSettings object.
     print(time_drift_bias_settings)






     )doc" );

    m.def( "arc_wise_time_drift_bias",
           py::overload_cast< const std::vector< Eigen::VectorXd >&,
                              const std::vector< double >&,
                              const tom::LinkEndType,
                              const std::vector< double >& >( &tom::arcWiseTimeDriftBias ),
           py::arg( "bias_value" ),
           py::arg( "arc_start_times" ),
           py::arg( "time_link_end" ),
           py::arg( "ref_epochs" ),
           R"doc(

 Function for creating settings for arc-wise time-drift biases.

 Function for creating settings for arc-wise time-drift biases.
 This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.time_drift_bias` setting only through the option of setting the `bias_value` (time drift bias) to a different values for each arc.


 Parameters
 ----------
 bias_value : numpy.ndarray
     Constant time drift bias that is to be considered for the observation time. This vector should be the same size as the observable to which it is
     assigned (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

 arc_start_times : List[ float ]
     List containing starting times for each arc.

 time_link_end : :class:`LinkEndType`
     Defines the link end (via the :class:`LinkEndType`) which is used the current time.

 ref_epochs : List[ float ]
     List containing the arc-wise reference epochs at which the effect of the arc-wise time drift is initialised.

 Returns
 -------
 :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`ArcWiseConstantObservationBiasSettings` class.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the arcwise_time_drift_bias function
     from tudatpy.numerical_simulation.estimation_setup import observation
     import numpy as np

     # The function arcwise_time_drift_bias() requires:
     # 1) an arc_start_times list ,2) a numpy.array of bias values in input, 3) reference_link_end_type, 4) list of reference epochs
     # Let's simulate two arcs
     arc_start_times = [0, 60] # define start time in seconds
     arcwise_bias_array = [np.array([1e-2]), np.array([2e-2])] # set arc bias
     reference_link_end_type = observation.receiver # set bias at receiving link end
     ref_epochs = [0,60]
     arcwise_time_drift_bias_settings = observation.arc_wise_time_drift_bias(arcwise_bias_array, arc_start_times, observation.receiver, ref_epochs)

     # Show that it returns an ObservationBiasSettings object.
     print(arcwise_time_drift_bias_settings)




     )doc" );

    m.def( "arc_wise_time_drift_bias_per_time",
           py::overload_cast< const std::map< double, Eigen::VectorXd >&,
                              const tom::LinkEndType,
                              const std::vector< double > >( &tom::arcWiseTimeDriftBias ),
           py::arg( "bias_value_per_start_time" ),
           py::arg( "time_link_end" ),
           py::arg( "ref_epochs" ),
           R"doc(

 Function for creating settings for arc-wise time-drift biases.

 Function for creating settings for arc-wise time-drift biases.
 This bias setting differs from the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.time_drift_bias` setting only through the option of setting the `bias_value` (time drift bias) to a different values for each arc.

 Parameters
 ----------
 bias_value_per_start_time : Dict[float, numpy.ndarray[numpy.float64[m, 1]]]
     Dictionary, in which the time bias value vectors for each arc are directly mapped to the starting times of the respective arc.
     The vectors should be the same size as the observable to which it is applied (*e.g.* size 1 for a range observable, size 2 for angular position, *etc*.)

 time_link_end : :class:`LinkEndType`
     Defines the link end (via the :class:`LinkEndType`) which is used the current time.

 ref_epochs : List[ float ]
     List containing the arc-wise reference epochs at which the effect of the arc-wise time drift is initialised.

 Returns
 -------
 :class:`ArcWiseConstantObservationBiasSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ArcWiseConstantObservationBiasSettings` class.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the arcwise_time_drift_bias_per_time function
     from tudatpy.numerical_simulation.estimation_setup import observation
     import numpy as np

     # The function arcwise_time_drift_bias_per_time() requires:
     # 1) an arc_start_times list ,2)  a dictionary with times as keys and bias values as values, 2) a reference_link_end_type, 3) reference_link_end_type, 4) list of reference epochs
     # Let's simulate two arcs
     bias_value_per_start_time = dict()
     bias_value_per_start_time[0] = np.array([1e-2])
     bias_value_per_start_time[60] = np.array([2e-2])
     reference_link_end_type = observation.receiver # set bias at receiving link end
     ref_epochs = [0,60]
     arcwise_time_drift_bias_settings = observation.arc_wise_time_drift_bias(bias_value_per_start_time, reference_link_end_type, ref_epochs)

     # Show that it returns an ObservationBiasSettings object.
     print(arcwise_time_drift_bias_settings)




     )doc" );

    m.def( "combined_bias",
           &tom::multipleObservationBiasSettings,
           py::arg( "bias_list" ),
           R"doc(

 Function for creating settings for a combined observation bias.

 Function for creating settings for a combined observation bias, calculating by combining any number of bias types.
 Each contribution of the combined bias is computed from the unbiased observable, so when applying both a relative and absolute bias, we get:

 .. math::
    \tilde{h}=h+K_{a}+hK_{r}

 And, crucially:

 .. math::
    \tilde{h}\neq (h+K_{a})(1+K_{r})

 where :math:`K_{r}` and :math:`K_{a}` is the relative and absolute bias, respectively.


 Parameters
 ----------
 bias_list : List[:class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings]
     A list containing the bias settings that are to be applied to the observable.

 Returns
 -------
 :class:`~tudatpy.numerical_simulation.estimation_setup.observation.multipleObservationBiasSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationBiasSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.multipleObservationBiasSettings` class, combining the settings for multiple observation biases.

 Examples
 --------
 .. code-block:: python

     # Code Snippet to showcase the use of the combined_bias function
     from tudatpy.numerical_simulation.estimation_setup import observation
     import numpy as np

     # The function combined_bias() allows to combine multiple ObservationBiasSettings objects.
     # Let's combine absolute, relative and time_drift biases.
     bias_array = np.array([1e-2])

     # Define absolute and relative bias settings
     absolute_bias_settings = observation.absolute_bias(bias_array)
     relative_bias_settings = observation.absolute_bias(bias_array)

     # Define Time Drift Bias
     time_link_end = observation.receiver
     ref_epoch = 0
     time_drift_bias_settings = observation.time_drift_bias(bias_array, time_link_end, ref_epoch)

     # combined_bias takes a list of ObservationBiasSettings objects as input
     bias_list = [absolute_bias_settings, relative_bias_settings, time_drift_bias_settings]
     combined_bias = observation.combined_bias(bias_list)

     # Show that it returns an ObservationBiasSettings object.
     print(combined_bias)




     )doc" );

    m.def( "two_way_time_scale_range_bias",
           &tom::twoWayTimeScaleRangeBias,
           R"doc(No documentation found.)doc" );

    // ###########    Observation Simulation Settings
    // #############

    py::enum_< tom::ObservationViabilityType >( m, "ObservationViabilityType", R"doc(

         Enumeration of observation viability criterion types.

 Examples
 --------
 .. code-block:: python

     # Code snippet to print all available Observation Viability Types
     from tudatpy.numerical_simulation import estimation_setup

     num_observation_viability_types = len(estimation_setup.observation.ObservationViabilityType.__members__)
     print(f'The length of all available Tudatpy Observation Viability Types is: {num_observation_viability_types}')

     # Print all available Observation Viability Types using the "name" property
     for i in range(num_observation_viability_types):
         print(i, estimation_setup.observation.ObservationViabilityType(i).name)




      )doc" )
            .value( "minimum_elevation_angle",
                    tom::ObservationViabilityType::minimum_elevation_angle )
            .value( "body_avoidance_angle", tom::ObservationViabilityType::body_avoidance_angle )
            .value( "body_occultation", tom::ObservationViabilityType::body_occultation )
            .export_values( );

    py::enum_< tom::ObservationAncilliarySimulationVariable >(
            m,
            "ObservationAncilliarySimulationVariable",
            R"doc(

         Enumeration of observation ancillary variable types.

 Examples
 --------
 .. code-block:: python

     # Code snippet to print all available Observation Ancillary Variable Types
     from tudatpy.numerical_simulation import estimation_setup

     num_observation_ancillary_variable_types = len(estimation_setup.observation.ObservationAncilliarySimulationVariable.__members__)
     print(f'The length of all available Tudatpy  Observation Ancillary Variable Types is: {num_observation_ancillary_variable_types}')

     # Print all Observation Ancillary Variable Types using the "name" property
     for i in range(num_observation_ancillary_variable_types):
         print(i, estimation_setup.observation.ObservationAncilliarySimulationVariable(i).name)



      )doc" )
            .value( "link_ends_delays",
                    tom::ObservationAncilliarySimulationVariable::link_ends_delays )
            .value( "doppler_integration_time",
                    tom::ObservationAncilliarySimulationVariable::doppler_integration_time )
            .value( "doppler_reference_frequency",
                    tom::ObservationAncilliarySimulationVariable::doppler_reference_frequency )
            .value( "frequency_bands",
                    tom::ObservationAncilliarySimulationVariable::frequency_bands )
            .value( "reception_reference_frequency_band",
                    tom::ObservationAncilliarySimulationVariable::
                            reception_reference_frequency_band )
            .value( "sequential_range_reference_frequency",
                    tom::ObservationAncilliarySimulationVariable::
                            sequential_range_reference_frequency )
            .value( "sequential_range_lowest_ranging_component",
                    tom::ObservationAncilliarySimulationVariable::
                            sequential_range_lowest_ranging_component )
            .value( "range_conversion_factor",
                    tom::ObservationAncilliarySimulationVariable::range_conversion_factor )
            .export_values( );

    py::class_< tom::ObservationViabilitySettings,
                std::shared_ptr< tom::ObservationViabilitySettings > >(
            m,
            "ObservationViabilitySettings",
            R"doc(

         Class for defining observation viability calculator settings.

         Class for defining the settings for observation viability calculator creation.
         Instances of this class are typically be created through various dedicated functions,such as :func:`~tudatpy.numerical_simulation.estimation_setup.observation.elevation_angle_viability`, :func:`~tudatpy.numerical_simulation.estimation_setup.observation.body_avoidance_viability` and :func:`~tudatpy.numerical_simulation.estimation_setup.observation.body_occultation_viability`

         Examples
         --------
         .. code-block:: python

             # Code snippet to show the creation of an ObservationViabilitySettings object
             import numpy as np
             from tudatpy.numerical_simulation.estimation_setup import observation

             # Create ObservationViabilitySettings object
             # In this case, we exclude observations for which the local elevation angle at link end is less 15 degrees.
             min_elevation = np.deg2rad(15)
             # We apply these settings to every ground station on Earth using the following link_end_id: [“Earth”, “”]
             viability_settings = observation.elevation_angle_viability(["Earth", ""], min_elevation)

             # Show that this is indeed an ObservationViabilitySettings object
             print(viability_settings)




      )doc" );

    m.def( "elevation_angle_viability",
           py::overload_cast< const std::pair< std::string, std::string >, const double >(
                   &tom::elevationAngleViabilitySettings ),
           py::arg( "link_end_id" ),
           py::arg( "elevation_angle" ),
           R"doc(

 Function for defining single elevation angle viability setting.

 Function for defining elevation angle viability settings for single link end.
 When simulating observations, this setting ensures that any applicable observations, for which the local elevation angle at link end is less than some limit value, will be omitted.


 Parameters
 ----------
 link_end_id : Tuple[str,str]
     Link end (as defined by body/reference point pair, see :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndId`), for which the elevation angle viability setting is to be created.
     To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""].

 elevation_angle : float
     Limit elevation angle, below which no observations are produced when using the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.simulate_observations` function. Note: this
     value must be in radians.

 Returns
 -------
 :class:`ObservationViabilitySettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` class, defining the settings for observation viability






     )doc" );

    m.def( "body_avoidance_viability",
           py::overload_cast< const std::pair< std::string, std::string >,
                              const std::string,
                              const double >( &tom::bodyAvoidanceAngleViabilitySettings ),
           py::arg( "link_end_id" ),
           py::arg( "body_to_avoid" ),
           py::arg( "avoidance_angle" ),
           R"doc(

 Function for defining body avoidance observation viability settings.

 Function for defining body avoidance observation viability settings for single link ends.
 When simulating observations, this settings ensures that any applicable observations, for which the signal path passes 'too close' to a body, will be omitted.
 The definition of 'too close' is computed as the angle between:

 * The line-of-sight vector from a link end to a given third body
 * The line-of-sight between two link ends

 This constraint is typically used to prevent the Sun from being too close to the field-of-view of the telescope(s), as defined by
 a so-called 'SPE' (Sun-Probe-Earth) angle constraint. The present viability setting generalizes this constraint.


 Parameters
 ----------
 link_end_id : Tuple[str,str]
     Link end (as defined by body/reference point pair, see :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndId` ), for which the viability settings are to be created.
     To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""] is entry in this list.
     For each link end included in this list, it will be checked if a signal received by and/or transmitted (or reflected) by this
     link end passes too close to the specified body.

 body_to_avoid : str
     Name of the body which the signal path should not pass 'too close' to.

 avoidance_angle : float
     Limit angle (generalization of SPE angle), below which no observations are produced when using the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.simulate_observations` function. Note: this
     value must be in radians.

 Returns
 -------
 :class:`ObservationViabilitySettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings`, defining the settings for observation viability.






     )doc" );

    m.def( "body_occultation_viability",
           py::overload_cast< const std::pair< std::string, std::string >, const std::string >(
                   &tom::bodyOccultationViabilitySettings ),
           py::arg( "link_end_id" ),
           py::arg( "occulting_body" ),
           R"doc(

 Function for defining body occultation viability settings.

 Function for defining body occultation viability settings for single link ends.
 When simulating observations, this setting ensures that any applicable observations, for which the signal path is occulted by a given body, will be omitted.
 The occultation is computed using the shape model of the specified body.


 Parameters
 ----------
 link_end_id : Tuple[str,str]
     Link end (as defined by body/reference point pair, see :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndId`), for which the viability settings are to be created.
     To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""] is entry in this list.

 body_to_avoid : str
     Name of the body which the signal path should not be occulted by.

 Returns
 -------
 :class:`ObservationViabilitySettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings`, defining the settings for observation viability.






     )doc" );

    m.def( "elevation_angle_viability_list",
           py::overload_cast< const std::vector< std::pair< std::string, std::string > >,
                              const double >( &tom::elevationAngleViabilitySettings ),
           py::arg( "link_end_ids" ),
           py::arg( "elevation_angle" ),
           R"doc(

 Function for defining list of elevation angle viability settings.

 Function for defining elevation angle viability settings for multiple link ends.
 Each entry in the returned list contains the observation viability settings for one link end.
 When simulating observations, these settings ensure that any applicable observations, for which the local elevation angle at a link end is less than some limit value, will be omitted.


 Parameters
 ----------
 link_end_ids : List[ Tuple[str,str] ]
     List of individual link ends (as defined by body/reference point pair, see :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndId`), for which the elevation angle viability setting is to be created.
     To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""].
     For each link end included in this list, it will be checked if a signal received by and/or transmitted (or reflected) by this
     link end violates the minimum elevation angle constraint.

 elevation_angle : float
     Limit elevation angle, below which no observations are produced when using the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.simulate_observations` function. Note: this
     value must be in radians.

 Returns
 -------
 :class:`ObservationViabilitySettings`
     List of :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` objects, each defining the settings for observation viability of one link end.






     )doc" );

    m.def( "body_avoidance_viability_list",
           py::overload_cast< const std::vector< std::pair< std::string, std::string > >,
                              const std::string,
                              const double >( &tom::bodyAvoidanceAngleViabilitySettings ),
           py::arg( "link_end_ids" ),
           py::arg( "body_to_avoid" ),
           py::arg( "avoidance_angle" ),
           R"doc(

 Function for defining list of body avoidance viability settings.

 Function for defining body avoidance viability settings for multiple link ends.
 Each entry in the returned list contains the observation viability settings for one link end.
 When simulating observations, these settings ensure that any applicable observations, for which the signal path passes 'too close' to a body, will be omitted.
 The definition of 'too close' is computed as the angle between:

 * The line-of-sight vector from a link end to a given third body
 * The line-of-sight between two link ends

 This constraint is typically used to prevent the Sun from being too close to the field-of-view of the telescope(s), as defined by
 a so-called 'SPE' (Sun-Probe-Earth) angle constraint. The present viability setting generalizes this constraint.


 Parameters
 ----------
 link_end_ids : List[ Tuple[str,str] ]
     List of individual link ends (as defined by body/reference point pair, see :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndId`), for which the elevation angle viability setting is to be created.
     To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""].

 body_to_avoid : str
     Name of the body which the signal path should not pass 'too close' to.

 avoidance_angle : float
     Limit angle (generalization of SPE angle), below which no observations are produced when using the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.simulate_observations` function. Note: this
     value must be in radians.

 Returns
 -------
 :class:`ObservationViabilitySettings`
     List of :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` objects, each defining the settings for observation viability of one link end.






     )doc" );

    m.def( "body_occultation_viability_list",
           py::overload_cast< const std::pair< std::string, std::string >, const std::string >(
                   &tom::bodyOccultationViabilitySettings ),
           py::arg( "link_end_id" ),
           py::arg( "occulting_body" ),
           R"doc(

 Function for defining body occultation viability settings.

 Function for defining body occultation viability settings for multiple link ends.
 Each entry in the returned list contains the observation viability settings for one link end.
 When simulating observations, these settings ensure that any applicable observations, for which the signal path is occulted by a given body, will be omitted.
 The occultation is computed using the shape model of the specified body.


 Parameters
 ----------
 link_end_ids : List[ Tuple[str,str] ]
     List of individual link ends (as defined by body/reference point pair, see :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkEndId`), for which the viability settings are to be created.
     To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""] is entry in this list.
     For each link end included in this list, it will be checked if a signal received by and/or transmitted (or reflected) by this
     link end is occulted by the specified body.

 body_to_avoid : str
     Name of the body which the signal path should not be occulted by.

 Returns
 -------
 :class:`ObservationViabilitySettings`
     List of :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` objects, each defining the settings for observation viability of one link end.






     )doc" );

    py::class_< tss::ObservationSimulationSettings< TIME_TYPE >,
                std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >(
            m,
            "ObservationSimulationSettings",
            R"doc(
         Base class for defining settings for simulated observations.
      )doc" )
            .def_property(
                    "viability_settings_list",
                    &tss::ObservationSimulationSettings< TIME_TYPE >::getViabilitySettingsList,
                    &tss::ObservationSimulationSettings< TIME_TYPE >::setViabilitySettingsList,
                    R"doc(
         viability_settings_list : List[ ObservationViabilitySettings ], default = [ ]) -
         Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.
      )doc" )
            .def_property(
                    "noise_function",
                    &tss::ObservationSimulationSettings< TIME_TYPE >::getObservationNoiseFunction,
                    py::overload_cast< const std::function< double( const double ) >& >(
                            &tss::ObservationSimulationSettings<
                                    TIME_TYPE >::setObservationNoiseFunction ),
                    R"doc(
         noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ], default = None -
         Function providing the observation noise as a function of observation time (can be constant or time-dependent), default is None.
      )doc" );
    //            .def_property("observable_type",
    //                         &tss::ObservationSimulationSettings<double>::getObservableType,
    //                         &tss::ObservationSimulationSettings<double>::setObservableType,
    //                         get_docstring("ObservationSimulationSettings.observable_type").c_str()
    //                         )
    //            .def_property_readonly("link_ends",
    //                         &tss::ObservationSimulationSettings<double>::getLinkEnds,
    //                         get_docstring("ObservationSimulationSettings.link_ends").c_str()
    //                         );

    py::class_< tss::TabulatedObservationSimulationSettings< TIME_TYPE >,
                std::shared_ptr< tss::TabulatedObservationSimulationSettings< TIME_TYPE > >,
                tss::ObservationSimulationSettings< TIME_TYPE > >(
            m,
            "TabulatedObservationSimulationSettings",
            R"doc(

         Class for defining settings for simulating observations at a predefined set of times.

         Class for defining settings for simulating observations at a predefined set of times.
         This type defines predefined time epochs at which applicable observations are to be simulated, stored in a rigid, "tabulated" form.
         Some observation times may be discarded due to the use of viability settings.
         Instances of this class are typically created via the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.tabulated_simulation_settings`
         and :func:`~tudatpy.numerical_simulation.estimation_setup.observation.tabulated_simulation_settings_list` functions.

         Associated base class: :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings`

         Examples
         --------
         .. code-block:: python

             # Code snippet to show the creation of a TabulatedObservationSimulationSettings object
             import numpy as np
             from tudatpy.astro.time_conversion import DateTime
             from tudatpy.numerical_simulation.estimation_setup import observation

             # Set simulation start and end epochs
             simulation_start_epoch = DateTime(2000, 1, 1).epoch()
             simulation_end_epoch   = DateTime(2000, 1, 4).epoch()

             # Define the uplink link ends for one-way observable
             link_ends = dict()
             link_ends[observation.transmitter] = observation.body_origin_link_end_id("Earth")
             link_ends[observation.receiver] = observation.body_origin_link_end_id("Delfi-C3")

             # Create LinkDefinition Object and set observation settings for each link/observable
             link_definition = observation.LinkDefinition(link_ends)
             observation_settings_list = [observation.one_way_doppler_instantaneous(link_definition)]

             # Define observation simulation times (separated by steps of 1 minute)
             observation_times = np.arange(simulation_start_epoch, simulation_end_epoch, 60.0)

             # Create TabulatedObservationSimulationSettings object
             tabulated_observation_simulation_settings = observation.tabulated_simulation_settings(
                 observation.one_way_instantaneous_doppler_type,
                 link_definition,
                 observation_times
             )

             # Show that this is indeed a TabulatedObservationSimulationSettings object
             print(tabulated_observation_simulation_settings)



      )doc" );

    py::class_< tom::ObservationAncilliarySimulationSettings,
                std::shared_ptr< tom::ObservationAncilliarySimulationSettings > >(
            m,
            "ObservationAncilliarySimulationSettings",
            R"doc(

         Base class for holding ancilliary settings for observation simulation.

         Base class for holding ancilliary settings for observation simulation.
         The user can create instances of this class via the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.elevation_angle_dependent_variable` function.

         Examples
         --------
         .. code-block:: python

             # Code snippet to show the creation of an ObservationAncillarySimulationSettings object
             from tudatpy.numerical_simulation.estimation_setup import observation

             # Example 1: Create ObservationAncillarySimulationSettings object using observation.n_way_range_ancilliary_settings function
             # In this case the frequency bands of the retransmitter - we set it to x band.
             n_way_range_ancillary_settings = observation.n_way_range_ancilliary_settings(frequency_bands=[observation.FrequencyBands.x_band])

             # Show that this is indeed an ObservationAncillarySimulationSettings object
             print(n_way_range_ancillary_settings)

             # Example 2: Create ObservationAncillarySimulationSettings object using observation.one_way_doppler_instantaneous function
             # In this case the integration time (in seconds) has to be given as input - we set it to 60s
             doppler_ancillary_settings = observation.doppler_ancilliary_settings(60)

             # Show that this is indeed an ObservationAncillarySimulationSettings object
             print(doppler_ancillary_settings)

             # [OPTIONAL] Verify that we indeed added Frequency Bands as Ancillary Simulation Variables for the n_way_range_ancillary_settings.
             list_num = n_way_range_ancillary_settings.get_float_list_settings(observation.ObservationAncilliarySimulationVariable.frequency_bands)
             for num in list_num:
                 name = observation.ObservationAncilliarySimulationVariable(int(num)).name
                 print(f'Ancillary Simulation Variable(s): {name}, corresponding to enumeration object n. {int(num)} of the ObservationAncilliarySimulationVariable Enumeration')



      )doc" )
            .def( "get_float_settings",
                  &tom::ObservationAncilliarySimulationSettings::getAncilliaryDoubleData,
                  py::arg( "setting_type" ),
                  py::arg( "throw_exception" ) = true,
                  R"doc(


         Parameters
         ----------
         setting_type : ObservationAncilliarySimulationVariable
             Type of the setting for which the value is to be returned

         throw_exception : bool, default = false
             Boolean defining whether to throw an exception if the requested setting does not exist, or does not exist as a floating point value.

         Returns
         -------
         float
             Value of the requested ancilliary variable (or NaN if it does not exist)

     )doc" )
            .def( "get_float_list_settings",
                  &tom::ObservationAncilliarySimulationSettings::getAncilliaryDoubleVectorData,
                  py::arg( "setting_type" ),
                  py::arg( "throw_exception" ) = true,
                  R"doc(


         Parameters
         ----------
         setting_type : ObservationAncilliarySimulationVariable
             Type of the setting for which the value is to be returned

         throw_exception : bool, default = false
             Boolean defining whether to throw an exception if the requested setting does not exist, or does not exist as list of floating point values.

         Returns
         -------
         list[ float ]
             Value of the requested ancilliary variable (or empty list if it does not exist)

         Examples
         --------
         .. code-block:: python

             # Code snippet to show how to retrieve ObservationAncillarySimulationSettings info
             # using the ObservationAncilliarySimulationSettings.get_float_settings function

             from tudatpy.numerical_simulation.estimation_setup import observation

             # Create Ancillary Settings
             n_way_range_ancillary_settings = observation.n_way_range_ancilliary_settings(frequency_bands=[observation.FrequencyBands.x_band])

             # Verify that we indeed added Frequency Bands as Ancillary Simulation Variables, using n_way_range_ancillary_settings.get_float_list_settings
             list_num = n_way_range_ancillary_settings.get_float_list_settings(observation.ObservationAncilliarySimulationVariable.frequency_bands)
             for num in list_num:
                 name = observation.ObservationAncilliarySimulationVariable(int(num)).name
                 print(f'Ancillary Simulation Variable(s): {name}, corresponding to enumeration object n. {int(num)} of the ObservationAncilliarySimulationVariable Enumeration')




     )doc" );

    m.def( "doppler_ancilliary_settings",
           &tom::getAveragedDopplerAncilliarySettings,
           py::arg( "integration_time" ) = 60.0,
           R"doc(

 Function for creating ancilliary settings for averaged Doppler observable.

 Function for creating ancilliary settings for an averaged Doppler observable. Specifically, this
 function can be used to create settings for the integration time of the observable. Note: in case no retransmission
 delays (or other additional ancilliary settings) are to be defined, this setting may be used for one-, two-, or N-way
 averaged Doppler.


 Parameters
 ----------
 integration_time : float, default = 60.0
     Integration time that is to be used for the averaged Doppler observable
 Returns
 -------
 :class:`ObservationAncilliarySimulationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationAncilliarySimulationSettings` with the required settings.






     )doc" );

    m.def( "two_way_range_ancilliary_settings",
           &tom::getTwoWayRangeAncilliarySettings,
           py::arg( "retransmission_delay" ) = 0.0,
           // py::arg("frequency_band") =
           // tom::FrequencyBands::x_band,
           R"doc(

 Function for creating ancilliary settings for two-way range observable.

 Function for creating ancilliary settings for a two-way range observable. Specifically, this
 function can be used to create settings for the retransmission delay of the observable. NOTE:
 this function is provided for convenience, and is equivalent to calling :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_range_ancilliary_settings`
 with a single retransmission delay.


 Parameters
 ----------
 retransmission_delay : float, default = 0.0
     Retransmission delay that is to be applied to the simulation of the two-way observable
 Returns
 -------
 :class:`ObservationAncilliarySimulationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationAncilliarySimulationSettings` with the required settings.






     )doc" );

    m.def( "two_way_doppler_ancilliary_settings",
           &tom::getTwoWayAveragedDopplerAncilliarySettings,
           py::arg( "integration_time" ) = 60.0,
           py::arg( "retransmission_delay" ) = 0.0,
           R"doc(

 Function for creating ancilliary settings for two-way averaged Doppler observable.

 Function for creating ancilliary settings for a two-way range observable. Specifically, this
 function can be used to create settings for the retransmission delay of the observable.  NOTE:
 this function is provided for convenience, and is equivalent to calling :func:`~tudatpy.numerical_simulation.estimation_setup.observation.n_way_doppler_ancilliary_settings`
 with a single retransmission delay.


 Parameters
 ----------
 integration_time : float, default = 60.0
     Integration time that is to be used for the averaged Doppler observable
 retransmission_delay : float, default = 0.0
     Retransmission delay that is to be applied to the simulation of the two-way observable
 Returns
 -------
 :class:`ObservationAncilliarySimulationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationAncilliarySimulationSettings` with the required settings.






     )doc" );

    m.def( "n_way_range_ancilliary_settings",
           &tom::getNWayRangeAncilliarySettings,
           py::arg( "link_end_delays" ) = std::vector< double >( ),
           py::arg( "frequency_bands" ) = std::vector< tom::FrequencyBands >( ),
           R"doc(

 Function for creating ancilliary settings for n-way range observable.

 Function for creating ancilliary settings for a n-way range observable. Specifically, this
 function can be used to create settings for the retransmission delays of the observable, for each of the retransmitters.


 Parameters
 ----------
 retransmission_delays : list[ float ], default = None
     Retransmission delays that are to be applied to the simulation of the n-way observable. If kept empty, this results in 0 retransmission delay at each retransmitter. If defined, this list must be the same length as the number of retransmitters, and the :math:`i^{th}` entry contains the retransmission delay of the :math:`i^{th}` retrasmitter
 Returns
 -------
 :class:`ObservationAncilliarySimulationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationAncilliarySimulationSettings` with the required settings.






     )doc" );

    m.def( "n_way_doppler_ancilliary_settings",
           &tom::getNWayAveragedDopplerAncilliarySettings,
           py::arg( "integration_time" ) = 60.0,
           py::arg( "link_end_delays" ) = std::vector< double >( ),
           py::arg( "frequency_bands" ) = std::vector< tom::FrequencyBands >( ),
           R"doc(

 Function for creating ancilliary settings for n-way averaged Doppler observable.

 Function for creating ancilliary settings for a n-way averaged Doppler observable. Specifically, this
 function can be used to create settings for the integration time of the observable, and the  retransmission delays for each of the retransmitters.


 Parameters
 ----------
 integration_time : float, default = 60.0
     Integration time that is to be used for the averaged Doppler observable
 retransmission_delays : list[ float ], default = None
     Retransmission delays that are to be applied to the simulation of the n-way observable. If kept empty, this results in 0 retransmission delay at each retransmitter. If defined, this list must be the same length as the number of retransmitters, and the :math:`i^{th}` entry contains the retransmission delay of the :math:`i^{th}` retrasmitter
 Returns
 -------
 :class:`ObservationAncilliarySimulationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationAncilliarySimulationSettings` with the required settings.






     )doc" );

    m.def( "dsn_n_way_doppler_ancilliary_settings",
           &tom::getDsnNWayAveragedDopplerAncillarySettings,
           py::arg( "frequency_bands" ),
           py::arg( "reference_frequency_band" ),
           py::arg( "reference_frequency" ),
           py::arg( "integration_time" ) = 60.0,
           py::arg( "link_end_delays" ) = std::vector< double >( ),
           R"doc(No documentation found.)doc" );

    m.def( "dsn_n_way_range_ancilliary_settings",
           &tom::getDsnNWayRangeAncillarySettings,
           py::arg( "frequency_bands" ),
           py::arg( "reference_frequency" ),
           py::arg( "lowest_ranging_component" ),
           py::arg( "link_end_delays" ) = std::vector< double >( ),
           R"doc(No documentation found.)doc" );

    m.def( "doppler_measured_frequency_ancillary_settings",
           &tom::getDopplerMeasuredFrequencyAncilliarySettings,
           py::arg( "frequency_bands" ),
           R"doc(No documentation found.)doc" );

    m.def( "tabulated_simulation_settings",
           &tss::tabulatedObservationSimulationSettings< TIME_TYPE >,
           py::arg( "observable_type" ),
           py::arg( "link_ends" ),
           py::arg( "simulation_times" ),
           py::arg( "reference_link_end_type" ) = tom::receiver,
           py::arg( "viability_settings" ) =
                   std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >( ),
           py::arg( "noise_function" ) = nullptr,
           py::arg( "ancilliary_settings" ) = nullptr,
           R"doc(

 Function for creating settings object for observation simulation, using a predefined list of observation times.

 Function for creating single simulation settings object, using a predefined list of observation times.
 The list of resulting observations may be reduced compared to the ``simulation_times`` should be provided here, as
 only observations that meet the viability settings are retained during observation simulation (these may be
 provide directly here through the ``viability_settings`` input, or added later to the resulting settings object).


 Parameters
 ----------
 observable_type : :class:`ObservableType`
     Observable type of which observations are to be simulated.
 link_ends : LinkDefinition
     Link ends for which observations are to be simulated.
 simulation_times : List[float]
     List of times at which to perform the observation simulation.
 reference_link_end_type : :class:`LinkEndType`, default = :class:`LinkEndType.receiver`
     Defines the link end (via the :class:`LinkEndType`) which is used as a reference time for the observation.
 viability_settings : List[ :class:`ObservationViabilitySettings` ], default = [ ]
     Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.

 noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ], default = None
     Function providing the observation noise factors as a function of observation time.
 Returns
 -------
 :class:`TabulatedObservationSimulationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.TabulatedObservationSimulationSettings` class.






     )doc" );

    m.def( "tabulated_simulation_settings_list",
           &tss::createTabulatedObservationSimulationSettingsList< TIME_TYPE >,
           py::arg( "link_ends_per_observable" ),
           py::arg( "simulation_times" ),
           py::arg( "reference_link_end_type" ) = tom::receiver,
           py::arg( "viability_settings" ) =
                   std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >( ),
           R"doc(

 Function for creating a list of settings object for observation simulation, using a predefined list of observation times.

 Function for creating multiple tabulated observation simulation settings objects in a list. This function is
 equivalent to calling the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.tabulated_simulation_settings` repeatedly, with the different
 observables and link definition provided here through `link_ends_per_observable`.
 During a single call to this function, one simulation settings object is created for each combination of observable type and link geometry given by the `link_ends_per_observable` parameter.


 Parameters
 ----------
 link_ends_per_observable : Dict[:class:`ObservableType`, List[LinkDefinition]]]
     Link geometry per observable type of which observations are to be simulated.
 simulation_times : List[ float ]
     List of times at which to perform the observation simulation.
 reference_link_end_type : :class:`LinkEndType`, default = :class:`LinkEndType.receiver`
     Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.
     The single link end specified here will be considered as the reference link end for all simulation settings object created in the function call.

 viability_settings : List[ :class:`ObservationViabilitySettings` ], default = [ ]
     Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.
     The single settings list given here will be considered as potential viability settings for all simulation settings object created in the function call.

 Returns
 -------
 List[ TabulatedObservationSimulationSettings ]
     List of :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.TabulatedObservationSimulationSettings` objects.






     )doc" );

    m.def( "continuous_arc_simulation_settings",
           &tss::perArcObservationSimulationSettings< TIME_TYPE >,
           py::arg( "observable_type" ),
           py::arg( "link_ends" ),
           py::arg( "start_time" ),
           py::arg( "end_time" ),
           py::arg( "interval_between_observations" ),
           py::arg( "arc_limiting_constraints" ),
           py::arg( "minimum_arc_duration" ),
           py::arg( "maximum_arc_duration" ),
           py::arg( "minimum_time_between_arcs" ),
           py::arg( "reference_link_end_type" ) = tom::receiver,
           py::arg( "additional_viability_settings" ) =
                   std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >( ),
           py::arg( "noise_function" ) = nullptr,
           R"doc(

 Function for creating settings object for observation simulation, using observation times according to a requirement for a continuous tracking arc.

 Function for creating settings object for observation simulation. Unlike the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.tabulated_simulation_settings`
 function, the resulting settings do not define the observation times explicitly. Instead, this settings object determines the observation times adaptively during the
 simulation of the observation, with the requirement that observations should be simulated over a set of contiguous arcs (if possible). The exact algorithm meets the following conditions:

 * Observations are only simulated within the time span of ``start_time`` and ``end_time``
 * A contiguous tracking arc has simulated observations separated by ``interval_between_observations``
 * Starting from ``start_time``, an observation is simulated each ``interval_between_observations``. Once an observation is unviable, as defined by
   the ``arc_limiting_constraints`` input, it is checked whether the arc up until that point
   is longer in duration than ``minimum_arc_duration``. If it is, the arc is added to the simulated observations. If not, the arc is discarded. In either case, a new arc is started once a
   viable is observation is encountered
 * If the current arc reaches a duration greater than ``maximum_arc_duration``, the arc is added to the existing observations, and a new arc is started
 * If defined (e.g. if not NaN), the current observation time is incremented by ``minimum_time_between_arcs`` when an arc has been added to the observations.

 Nominally, this algorithm ensures that any arc of observations has a minimum and maximum duration. In addition, it ensures that (if desired) there is a minimum time interval
 between two tracking arcs. This behaviour can be modified by adding ``additional_viability_settings``, which are *not* used when computing the tracking arcs, but which are instead only used
 to reduce the set of simulated observations afterwards.


 Parameters
 ----------
 observable_type : :class:`ObservableType`
     Observable type of which observations are to be simulated.
 link_ends : LinkDefinition
     Link ends for which observations are to be simulated.
 start_time : float
     First time at which an observation is to be simulated (and checked for viability).
 end_time : float
     Maximum time at which an observation is to be simulated (and checked for viability).
 interval_between_observations : float
     Cadence (in seconds) of subsequent observations in an arc
 arc_limiting_constraints : List[ :class:`ObservationViabilitySettings` ], default = [ ]
     List of settings for the creation of the viability criteria calculators, which are used to check if an observation is viable, and define
     whether an arc should be terminated.

 minimum_arc_duration : float
     Minimum permissible time for a tracking arc
 maximum_arc_duration : float
     Maximum permissible time for a tracking arc
 minimum_time_between_arc : float, default = NaN
     Minimum time between two tracking arcs. If NaN, this is effectively set to the ``interval_between_observations``
 additional_viability_settings : List[ :class:`ObservationViabilitySettings` ], default = [ ]
     Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.
     These settings are *not* used to determine whether an arc is to be terminated, but are instead applied after the arcs have been computed.

 noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ], default = None
     Function providing the observation noise factors as a function of observation time.
 Returns
 -------
 :class:`TabulatedObservationSimulationSettings`
     Instance of the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.TabulatedObservationSimulationSettings` class.






     )doc" );

    m.def( "continuous_arc_simulation_settings_list",
           &tss::perArcObservationSimulationSettingsList< TIME_TYPE >,
           py::arg( "link_ends_per_observable" ),
           py::arg( "start_time" ),
           py::arg( "end_time" ),
           py::arg( "interval_between_observations" ),
           py::arg( "arc_limiting_constraints" ),
           py::arg( "minimum_arc_duration" ),
           py::arg( "maximum_arc_duration" ),
           py::arg( "minimum_time_between_arcs" ),
           py::arg( "reference_link_end_type" ) = tom::receiver,
           py::arg( "additional_viability_settings" ) =
                   std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >( ),
           R"doc(

 Function for creating a list of settings object for observation simulation, using observation times according to a requirement for a continuous tracking arc.

 Function for creating multiple settings objects for observation simulation in a list. This function is
 equivalent to calling the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.continuous_arc_simulation_settings` repeatedly, with the different
 observables and link definition provided here through `link_ends_per_observable`.
 During a single call to this function, one simulation settings object is created for each combination of observable type and link geometry given by the `link_ends_per_observable` parameter.


 Parameters
 ----------
 link_ends_per_observable : Dict[:class:`ObservableType`, List[LinkDefinition]]]
     Link geometry per observable type of which observations are to be simulated.
 start_time : float
     First time at which an observation is to be simulated (and checked for viability).
 end_time : float
     Maximum time at which an observation is to be simulated (and checked for viability).
 interval_between_observations : float
     Cadence (in seconds) of subsequent observations in an arc
 arc_limiting_constraints : List[ :class:`ObservationViabilitySettings` ], default = [ ]
     List of settings for the creation of the viability criteria calculators, which are used to check if an observation is viable, and define
     whether an arc should be terminated.

 minimum_arc_duration : float
     Minimum permissible time for a tracking arc
 maximum_arc_duration : float
     Maximum permissible time for a tracking arc
 minimum_time_between_arc : float, default = NaN
     Minimum time between two tracking arcs. If NaN, this is effectively set to the ``interval_between_observations``
 additional_viability_settings : List[ :class:`ObservationViabilitySettings` ], default = [ ]
     Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.
     These settings are *not* used to determine whether an arc is to be terminated, but are instead applied after the arcs have been computed.

 Returns
 -------
 List[ :class:`TabulatedObservationSimulationSettings` ]
     List of :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` derived :class:`~tudatpy.numerical_simulation.estimation_setup.observation.TabulatedObservationSimulationSettings` objects.






     )doc" );

    m.def( "add_noise_function_to_all",
           py::overload_cast< const std::vector< std::shared_ptr<
                                      tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                              const std::function< Eigen::VectorXd( const double ) > >(
                   &tss::addNoiseFunctionToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "noise_amplitude" ),
           R"doc(

 Function for adding a custom noise function to all existing observation simulation settings.

 Function for including a custom noise function to the simulation settings of all observables.
 The noise settings are added to all :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) in the `observation_simulation_settings`
 list.

 Note: the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects are modified in-place by this function,
 and thus the function does not return anything.


 Parameters
 ----------
 observation_simulation_settings_list : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.

 noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ]
     Function providing the observation noise factors as a function of observation time.

 Returns
 -------
 None
     The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.







     )doc" );

    m.def( "add_noise_function_to_observable",
           py::overload_cast< const std::vector< std::shared_ptr<
                                      tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                              const std::function< Eigen::VectorXd( const double ) >,
                              const tom::ObservableType >(
                   &tss::addNoiseFunctionToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "noise_amplitude" ),
           py::arg( "observable_type" ),
           R"doc(

 Function for adding a custom noise function to selected existing observation simulation settings of a given observable type.

 As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_noise_function_to_all`, except that the function only adds noise to entries of the
 `observation_simulation_settings` list that matches the specified `observable_type`.


 Parameters
 ----------
 observation_simulation_settings_list : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.

 noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ]
     Function providing the observation noise factors as a function of observation time.

 observable_type : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservableType`
     Identifies the observable type in the observation simulation settings to which the noise is to be added.

 Returns
 -------
 None
     The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.







     )doc" );

    m.def( "add_noise_function_to_observable_for_link_ends",
           py::overload_cast< const std::vector< std::shared_ptr<
                                      tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                              const std::function< Eigen::VectorXd( const double ) >,
                              const tom::ObservableType,
                              const tom::LinkDefinition& >(
                   &tss::addNoiseFunctionToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "noise_amplitude" ),
           py::arg( "observable_type" ),
           py::arg( "link_ends" ),
           R"doc(

 Function for adding a custom noise function to existing observation simulation settings of a given observable type and link definition.

 As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_noise_function_to_all`, except that the function only adds noise to entries of the
 `observation_simulation_settings` list that matches the specified `observable_type` and `link_definition`.


 Parameters
 ----------
 observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.

 noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ]
     Function providing the observation noise factors as a function of observation time.

 observable_type : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservableType`
     Identifies the observable type in the observation simulation settings to which the noise is to be added.

 link_definition : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkDefinition`
     Identifies the link definition in the observation simulation settings for which the noise is to be added.

 Returns
 -------
 None
     The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.







     )doc" );

    m.def( "add_gaussian_noise_to_all",
           py::overload_cast< const std::vector< std::shared_ptr<
                                      tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                              const double >(
                   &tss::addGaussianNoiseFunctionToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "noise_amplitude" ),
           R"doc(

 Function for adding gaussian noise function to all existing observation simulation settings.

 Function for including simple time-independent and time-uncorrelated Gaussian noise function to the simulation settings of one or more observable(s).
 The noise settings are added to all :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) in the `observation_simulation_settings`
 list.

 Note: the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects are modified in-place by this function,
 and thus the function does not return anything.


 Parameters
 ----------
 observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
 noise_amplitude : float
     Standard deviation defining the un-biased Gaussian distribution for the noise.
 Returns
 -------
 None
     The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.







     )doc" );

    m.def( "add_gaussian_noise_to_observable",
           py::overload_cast< const std::vector< std::shared_ptr<
                                      tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                              const double,
                              const tom::ObservableType >(
                   &tss::addGaussianNoiseFunctionToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "noise_amplitude" ),
           py::arg( "observable_type" ),
           R"doc(

 Function for adding gaussian noise function to existing observation simulation settings of a given observable type.

 As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_gaussian_noise_to_all`, except that the function only adds noise to entries of the
 `observation_simulation_settings` list that matches the specified `observable_type`.


 Parameters
 ----------
 observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
 noise_amplitude : float
     Standard deviation defining the un-biased Gaussian distribution for the noise.
 observable_type : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservableType`
     Identifies the observable type in the observation simulation settings to which the noise is to be added.

 Returns
 -------
 None
     The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.







     )doc" );

    m.def( "add_gaussian_noise_to_observable_for_link_ends",
           py::overload_cast< const std::vector< std::shared_ptr<
                                      tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                              const double,
                              const tom::ObservableType,
                              const tom::LinkDefinition& >(
                   &tss::addGaussianNoiseFunctionToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "noise_amplitude" ),
           py::arg( "observable_type" ),
           py::arg( "link_definition" ),
           R"doc(

 Function for adding gaussian noise function to existing observation simulation settings of a given observable type and link definition.

 As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_gaussian_noise_to_all`, except that the function only adds noise to entries of the
 `observation_simulation_settings` list that matches the specified `observable_type` and `link_definition`.


 Parameters
 ----------
 observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
 noise_amplitude : float
     Standard deviation defining the un-biased Gaussian distribution for the noise.
 observable_type : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservableType`
     Identifies the observable type in the observation simulation settings to which the noise is to be added.

 link_definition : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkDefinition`
     Identifies the link definition in the observation simulation settings for which the noise is to be added.

 Returns
 -------
 None
     The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.







     )doc" );

    m.def( "add_viability_check_to_all",
           py::overload_cast<
                   const std::vector<
                           std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                   const std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >& >(
                   &tss::addViabilityToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "viability_settings" ),
           R"doc(

 Function for including viability checks into existing observation simulation settings.

 Function for adding viability checks to the observation simulation settings, such that only observations meeting certain conditions are retained.
 The noise settings are added to all :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) in the `observation_simulation_settings`
 list.
 Note: the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects are modified in-place by this function,
 and thus the function does not return anything.


 Parameters
 ----------
 observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
 viability_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` ]
     List of one or more :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` objects, defining the viability checks to be included.

 Returns
 -------
 None
     The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.







     )doc" );

    m.def( "add_viability_check_to_observable",
           py::overload_cast<
                   const std::vector<
                           std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                   const std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >&,
                   const tom::ObservableType >(
                   &tss::addViabilityToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "viability_settings" ),
           py::arg( "observable_type" ),
           R"doc(

 Function for including viability checks into existing observation simulation settings.

 As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_viability_check_to_all`, except that the function only adds viabilitt settings to entries of the
 `observation_simulation_settings` list that matches the specified `observable_type`.


 Parameters
 ----------
 observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
 viability_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` ]
     List of one or more :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` objects, defining the viability checks to be included.

 observable_type : :class:`ObservableType`
     Identifies the observable type in the observation simulation settings for which the viability checks are to be considered.

 Returns
 -------
 None
     The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.







     )doc" );

    m.def( "add_viability_check_to_observable_for_link_ends",
           py::overload_cast<
                   const std::vector<
                           std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                   const std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >&,
                   const tom::ObservableType,
                   const tom::LinkDefinition& >(
                   &tss::addViabilityToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "viability_settings" ),
           py::arg( "observable_type" ),
           py::arg( "link_ends" ),
           R"doc(

 Function for including viability checks into existing observation simulation settings.

 As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_viability_check_to_all`, except that the function only adds noise to entries of the
 `observation_simulation_settings` list that matches the specified `observable_type` and `link_definition`.


 Parameters
 ----------
 observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
 viability_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` ]
     List of one or more :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationViabilitySettings` objects, defining the viability checks to be included.

 observable_type : :class:`ObservableType`
     Identifies the observable type in the observation simulation settings for which the viability checks are to be considered.

 link_definition : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkDefinition`
     Identifies the link definition in the observation simulation settings for which the viability checks are to be considered.

 Returns
 -------
 None







     )doc" );

    m.def( "add_ancilliary_settings_to_observable",
           py::overload_cast<
                   const std::vector<
                           std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                   const std::shared_ptr< tom::ObservationAncilliarySimulationSettings >&,
                   const tom::ObservableType >(
                   &tss::addAncilliarySettingsToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "ancilliary_settings" ),
           py::arg( "observable_type" ),
           R"doc(No documentation found.)doc" );

    m.def( "add_ancilliary_settings_to_observable_for_link_ends",
           py::overload_cast<
                   const std::vector<
                           std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                   const std::shared_ptr< tom::ObservationAncilliarySimulationSettings >&,
                   const tom::ObservableType,
                   const tom::LinkDefinition& >(
                   &tss::addAncilliarySettingsToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "ancilliary_settings" ),
           py::arg( "observable_type" ),
           py::arg( "link_ends" ),
           R"doc(No documentation found.)doc" );

    py::class_< tss::ObservationDependentVariableSettings,
                std::shared_ptr< tss::ObservationDependentVariableSettings > >(
            m,
            "ObservationDependentVariableSettings",
            R"doc(

         Base class for setting observation dependent variables as part of the observation output.

         Base class for setting observation dependent variables as part of the observation output.
         The user can create instances of this class via the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.elevation_angle_dependent_variable` function.
         Note: The associated functionality is not yet mature enough for the end user. Class is exposed for development purposes only.

         Examples
         --------
         .. code-block:: python

             # Code snippet to show the creation of an ObservationDependentVariableSettings object
             from tudatpy.numerical_simulation.estimation_setup import observation

             # Create ObservationDependentVariableSettings object
             elevation_angle_settings = observation.elevation_angle_dependent_variable(observation.receiver)

             # Show that this is indeed an ObservationDependentVariableSettings object
             print(elevation_angle_settings)



      )doc" );

    m.def( "add_dependent_variables_to_all",
           py::overload_cast<
                   const std::vector<
                           std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                   const std::vector<
                           std::shared_ptr< tss::ObservationDependentVariableSettings > >&,
                   const tss::SystemOfBodies& >(
                   &tss::addDependentVariablesToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings" ),
           py::arg( "dependent_variable_settings" ),
           py::arg( "bodies" ),
           R"doc(

 Function for including dependent variables into all existing observation simulation settings.

 Function for including the computation and reporting of dependent variables into the observation simulation settings of all observables.
 Note: The associated functionality is not yet mature enough for the end user. Function is exposed for development purposes only.

 Modifications are applied to all given :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s),
 matching each :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object with the corresponding :class:`ObservationDependentVariableSettings` entry in the `dependent_variable_settings` parameter.
 Note that the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects are modified in-place and thus the function does not return anything.


 Parameters
 ----------
 observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.

 dependent_variable_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationDependentVariableSettings` ]
     List of one or more :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationDependentVariableSettings` objects, defining the dependent variables to be considered.

 bodies : :class:`~tudatpy.numerical_simulation.environment_setup.SystemOfBodies`
     Object consolidating all bodies and environment models that constitute the physical environment.






     )doc" );

    m.def( "add_dependent_variables_to_observable",
           py::overload_cast<
                   const std::vector<
                           std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                   const std::vector<
                           std::shared_ptr< tss::ObservationDependentVariableSettings > >&,
                   const tss::SystemOfBodies&,
                   const tom::ObservableType >(
                   &tss::addDependentVariablesToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings" ),
           py::arg( "dependent_variable_settings" ),
           py::arg( "bodies" ),
           py::arg( "observable_type" ),
           R"doc(

 Function for including dependent variables into selected existing observation simulation settings.

 As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_dependent_variables_to_all`, except that the function only adds includes the
 computation and reporting of dependent variables to entries of the `observation_simulation_settings` list that matches the specified `observable_type`.


 Parameters
 ----------
 observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.

 dependent_variable_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationDependentVariableSettings` ]
     List of one or more :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationDependentVariableSettings` objects, defining the dependent variables to be considered.

 bodies : :class:`~tudatpy.numerical_simulation.environment_setup.SystemOfBodies`
     Object consolidating all bodies and environment models that constitute the physical environment.

 observable_type : :class:`ObservableType`
     Identifies the observable type in the observation simulation settings for which the dependent variables are to be included.






     )doc" );

    m.def( "add_dependent_variables_to_observable_for_link_ends",
           py::overload_cast<
                   const std::vector<
                           std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                   const std::vector<
                           std::shared_ptr< tss::ObservationDependentVariableSettings > >&,
                   const tss::SystemOfBodies&,
                   const tom::ObservableType,
                   const tom::LinkDefinition& >(
                   &tss::addDependentVariablesToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings" ),
           py::arg( "dependent_variable_settings" ),
           py::arg( "bodies" ),
           py::arg( "observable_type" ),
           py::arg( "link_ends" ),
           R"doc(No documentation found.)doc" );

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // DEPENDENT VARIABLES
    /////////////////////////////////////////////////////////////////////////////////////////////////

    py::enum_< tss::IntegratedObservationPropertyHandling >(
            m, "IntegratedObservationPropertyHandling", R"doc(No documentation found.)doc" )
            .value( "interval_start", tss::IntegratedObservationPropertyHandling::interval_start )
            .value( "interval_end", tss::IntegratedObservationPropertyHandling::interval_end )
            .value( "interval_undefined",
                    tss::IntegratedObservationPropertyHandling::interval_undefined )
            .export_values( );

    m.def( "elevation_angle_dependent_variable",
           &tss::elevationAngleDependentVariable,
           py::arg( "link_end_type" ) = tom::unidentified_link_end,
           py::arg( "link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "originating_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "originating_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "integrated_observation_handling" ) = tss::interval_start,
           R"doc(No documentation found.)doc" );

    m.def( "azimuth_angle_dependent_variable",
           &tss::azimuthAngleDependentVariable,
           py::arg( "link_end_type" ) = tom::unidentified_link_end,
           py::arg( "link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "originating_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "originating_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "integrated_observation_handling" ) = tss::interval_start,
           R"doc(No documentation found.)doc" );

    m.def( "target_range_between_link_ends_dependent_variable",
           &tss::targetRangeBetweenLinkEndsDependentVariable,
           py::arg( "start_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "end_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "start_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "end_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "integrated_observation_handling" ) = tss::interval_start,
           R"doc(No documentation found.)doc" );

    m.def( "avoidance_angle_dependent_variable",
           &tss::bodyAvoidanceAngleDependentVariable,
           py::arg( "body_name" ),
           py::arg( "start_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "end_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "start_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "end_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "integrated_observation_handling" ) = tss::interval_start,
           R"doc(No documentation found.)doc" );

    m.def( "body_center_distance_dependent_variable",
           &tss::linkBodyCenterDistanceDependentVariable,
           py::arg( "body_name" ),
           py::arg( "start_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "end_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "start_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "end_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "integrated_observation_handling" ) = tss::interval_start,
           R"doc(No documentation found.)doc" );

    m.def( "body_limb_distance_dependent_variable",
           &tss::linkLimbDistanceDependentVariable,
           py::arg( "body_name" ),
           py::arg( "start_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "end_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "start_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "end_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "integrated_observation_handling" ) = tss::interval_start,
           R"doc(No documentation found.)doc" );

    m.def( "angle_wrt_orbital_plane_dependent_variable",
           &tss::linkAngleWrtOrbitalPlaneDependentVariable,
           py::arg( "body_name" ),
           py::arg( "start_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "end_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "start_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "end_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "integrated_observation_handling" ) = tss::interval_start,
           R"doc(No documentation found.)doc" );

    m.def( "integration_time_dependent_variable",
           &tss::integrationTimeDependentVariable,
           py::arg( "observable_type" ) = tom::undefined_observation_model,
           R"doc(No documentation found.)doc" );

    m.def( "retransmission_delays_dependent_variable",
           &tss::retransmissionDelaysDependentVariable,
           py::arg( "observable_type" ) = tom::undefined_observation_model,
           R"doc(No documentation found.)doc" );

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // FREQUENCIES
    /////////////////////////////////////////////////////////////////////////////////////////////////

    py::enum_< tom::FrequencyBands >( m, "FrequencyBands", R"doc(No documentation found.)doc" )
            .value( "s_band", tom::FrequencyBands::s_band )
            .value( "x_band", tom::FrequencyBands::x_band )
            .value( "ka_band", tom::FrequencyBands::ka_band )
            .value( "ku_band", tom::FrequencyBands::ku_band );

    m.def( "dsn_default_turnaround_ratios",
           &tom::getDsnDefaultTurnaroundRatios,
           py::arg( "uplink_band" ),
           py::arg( "downlink_band" ),
           R"doc(No documentation found.)doc" );

    m.def( "cassini_turnaround_ratios",
           &tom::getCassiniTurnaroundRatio,
           py::arg( "uplink_band" ),
           py::arg( "downlink_band" ),
           R"doc(No documentation found.)doc" );

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // ODF OBSERVATIONS
    /////////////////////////////////////////////////////////////////////////////////////////////////

    py::class_< tom::ProcessedOdfFileContents< TIME_TYPE >,
                std::shared_ptr< tom::ProcessedOdfFileContents< TIME_TYPE > > >(
            m, "ProcessedOdfFileContents", R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "ground_station_names",
                    &tom::ProcessedOdfFileContents< TIME_TYPE >::getGroundStationsNames,
                    R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "processed_observable_types",
                    &tom::ProcessedOdfFileContents< TIME_TYPE >::getProcessedObservableTypes,
                    R"doc(No documentation found.)doc" )
            .def_property_readonly( "start_and_end_time",
                                    &tom::ProcessedOdfFileContents< TIME_TYPE >::getStartAndEndTime,
                                    R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "ignored_odf_observable_types",
                    &tom::ProcessedOdfFileContents< TIME_TYPE >::getIgnoredRawOdfObservableTypes,
                    R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "ignored_ground_stations",
                    &tom::ProcessedOdfFileContents< TIME_TYPE >::getIgnoredGroundStations,
                    R"doc(No documentation found.)doc" )
            .def_property_readonly( "raw_odf_data",
                                    &tom::ProcessedOdfFileContents< TIME_TYPE >::getRawOdfData,
                                    R"doc(No documentation found.)doc" )
            .def( "define_antenna_id",
                  py::overload_cast< const std::string&, const std::string& >(
                          &tom::ProcessedOdfFileContents< TIME_TYPE >::defineSpacecraftAntennaId ),
                  py::arg( "spacecraft_name" ),
                  py::arg( "antenna_name" ),
                  R"doc(No documentation found.)doc" );

    m.def( "process_odf_data_multiple_files",
           py::overload_cast< const std::vector< std::string >&,
                              const std::string&,
                              const bool,
                              const std::map< std::string, Eigen::Vector3d >& >(
                   &tom::processOdfData< TIME_TYPE > ),
           py::arg( "file_names" ),
           py::arg( "spacecraft_name" ),
           py::arg( "verbose" ) = true,
           py::arg( "earth_fixed_ground_station_positions" ) =
                   tss::getApproximateDsnGroundStationPositions( ),
           R"doc(No documentation found.)doc" );

    m.def( "process_odf_data_single_file",
           py::overload_cast< const std::string&,
                              const std::string&,
                              const bool,
                              const std::map< std::string, Eigen::Vector3d >& >(
                   &tom::processOdfData< TIME_TYPE > ),
           py::arg( "file_name" ),
           py::arg( "spacecraft_name" ),
           py::arg( "verbose" ) = true,
           py::arg( "earth_fixed_ground_station_positions" ) =
                   tss::getApproximateDsnGroundStationPositions( ),
           R"doc(No documentation found.)doc" );

    // Create wrapper function
    py::cpp_function getDsnDefaultTurnaroundRatios_wrapper =
            []( tudat::observation_models::FrequencyBands band1,
                tudat::observation_models::FrequencyBands band2 ) {
                return tom::getDsnDefaultTurnaroundRatios( band1, band2 );
            };

    m.def( "set_odf_information_in_bodies",
           &tom::setOdfInformationInBodies< TIME_TYPE >,
           py::arg( "processed_odf_file" ),
           py::arg( "bodies" ),
           py::arg( "body_with_ground_stations_name" ) = "Earth",
           py::arg( "turnaround_ratio_function" ) = getDsnDefaultTurnaroundRatios_wrapper,
           R"doc(No documentation found.)doc" );

    m.def( "create_odf_observed_observation_collection",
           &tom::createOdfObservedObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "processed_odf_file" ),
           py::arg( "observable_types_to_process" ),
           py::arg( "start_and_end_times_to_process" ),
           R"doc(No documentation found.)doc" );

    m.def( "observations_from_odf_files",
           &tom::createOdfObservedObservationCollectionFromFile< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "bodies" ),
           py::arg( "odf_file_names" ),
           py::arg( "target_name" ),
           py::arg( "verbose_output" ) = true,
           py::arg( "earth_fixed_station_positions" ) =
                   tss::getApproximateDsnGroundStationPositions( ),
           R"doc(No documentation found.)doc" );

    m.def( "observations_from_ifms_files",
           &tom::createIfmsObservedObservationCollectionFromFiles< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "ifms_file_names" ),
           py::arg( "bodies" ),
           py::arg( "target_name" ),
           py::arg( "ground_station_name" ),
           py::arg( "reception_band" ),
           py::arg( "transmission_band" ),
           py::arg( "apply_troposphere_correction" ) = true,
           py::arg( "earth_fixed_station_positions" ) =
                   tss::getCombinedApproximateGroundStationPositions( ),
           R"doc(No documentation found.)doc" );

    m.def( "observations_from_multi_station_ifms_files",
           &tom::createMultiStationIfmsObservedObservationCollectionFromFiles< STATE_SCALAR_TYPE,
                                                                               TIME_TYPE >,
           py::arg( "ifms_file_names" ),
           py::arg( "bodies" ),
           py::arg( "target_name" ),
           py::arg( "ground_station_names" ),
           py::arg( "reception_band" ),
           py::arg( "transmission_band" ),
           py::arg( "apply_troposphere_correction" ) = true,
           py::arg( "earth_fixed_station_positions" ) =
                   tss::getCombinedApproximateGroundStationPositions( ),
           R"doc(No documentation found.)doc" );

    m.def( "observations_from_fdets_files",
           &tom::createFdetsObservedObservationCollectionFromFile< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "ifms_file_name" ),
           py::arg( "base_frequency" ),
           py::arg( "column_types" ),
           py::arg( "target_name" ),
           py::arg( "transmitting_station_name" ),
           py::arg( "receiving_station_name" ),
           py::arg( "reception_band" ),
           py::arg( "transmission_band" ),
           py::arg( "earth_fixed_station_positions" ) =
                   tss::getCombinedApproximateGroundStationPositions( ),
           R"doc(No documentation found.)doc" );

    m.def( "create_compressed_doppler_collection",
           &tom::createCompressedDopplerCollection< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "original_observation_collection" ),
           py::arg( "compression_ratio" ),
           py::arg( "minimum_number_of_observations" ) = 10,
           R"doc(No documentation found.)doc" );

    //    m.def("create_odf_observation_simulation_settings_list",
    //          &tom::createOdfObservationSimulationSettingsList<
    //          STATE_SCALAR_TYPE, TIME_TYPE >,
    //          py::arg("observed_observation_collection"),
    //          get_docstring("create_odf_observation_simulation_settings_list").c_str()
    //          );

    m.def( "change_simulation_settings_observable_types",
           &tom::changeObservableTypesOfObservationSimulationSettings< STATE_SCALAR_TYPE,
                                                                       TIME_TYPE >,
           py::arg( "observation_simulation_settings" ),
           py::arg( "replacement_observable_types" ) =
                   std::map< tom::ObservableType, tom::ObservableType >{
                           { tom::dsn_n_way_averaged_doppler, tom::n_way_differenced_range },
                           { tom::dsn_one_way_averaged_doppler, tom::one_way_differenced_range } },
           R"doc(No documentation found.)doc" );

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // Tracking Txt OBSERVATIONS
    /////////////////////////////////////////////////////////////////////////////////////////////////

    m.def( "create_tracking_txtfile_observation_collection",
           py::overload_cast< const std::shared_ptr< tudat::input_output::TrackingTxtFileContents >,
                              const std::string,
                              const std::vector< tom::ObservableType >,
                              const std::map< std::string, Eigen::Vector3d >,
                              const tom::ObservationAncilliarySimulationSettings& >(
                   &tom::createTrackingTxtFileObservationCollection< double, TIME_TYPE > ),
           py::arg( "raw_tracking_txtfile_contents" ),
           py::arg( "spacecraft_name" ),
           py::arg( "observable_types_to_process" ) = std::vector< tom::ObservableType >( ),
           py::arg( "earth_fixed_ground_station_positions" ) =
                   tss::getApproximateDsnGroundStationPositions( ),
           py::arg( "ancillary_settings" ) = tom::ObservationAncilliarySimulationSettings( ),
           R"doc(No documentation found.)doc" );

    //////////////////////////////////////////// DEPRECATED
    ///////////////////////////////////////////////

    m.def( "one_way_open_loop_doppler",
           &tom::oneWayOpenLoopDoppler,
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) =
                   std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "transmitter_proper_time_rate_settings" ) = nullptr,
           py::arg( "receiver_proper_time_rate_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) =
                   std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           py::arg( "normalized_with_speed_of_light" ) = false );

    m.def( "two_way_open_loop_doppler_from_one_way_links",
           py::overload_cast< const std::shared_ptr< tom::OneWayDopplerObservationSettings >,
                              const std::shared_ptr< tom::OneWayDopplerObservationSettings >,
                              const std::shared_ptr< tom::ObservationBiasSettings > >(
                   &tom::twoWayOpenLoopDoppler ),
           py::arg( "uplink_doppler_settings" ),
           py::arg( "downlink_doppler_settings" ),
           py::arg( "bias_settings" ) = nullptr );

    m.def( "two_way_open_loop_doppler",
           py::overload_cast<
                   const tom::LinkDefinition&,
                   const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >&,
                   const std::shared_ptr< tom::ObservationBiasSettings >,
                   const std::shared_ptr< tom::LightTimeConvergenceCriteria >,
                   const bool >( &tom::twoWayOpenLoopDoppler ),
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) =
                   std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) =
                   std::make_shared< tom::LightTimeConvergenceCriteria >( ),
           py::arg( "normalized_with_speed_of_light" ) = false );

    m.def( "one_way_closed_loop_doppler",
           py::overload_cast<
                   const tom::LinkDefinition&,
                   const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >,
                   const std::shared_ptr< tom::ObservationBiasSettings >,
                   const std::shared_ptr< tom::LightTimeConvergenceCriteria > >(
                   &tom::oneWayClosedLoopDoppler ),
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) =
                   std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) =
                   std::make_shared< tom::LightTimeConvergenceCriteria >( ) );

    m.def( "one_way_closed_loop_doppler",
           py::overload_cast<
                   const tom::LinkDefinition&,
                   const std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >,
                   const std::shared_ptr< tom::ObservationBiasSettings >,
                   const std::shared_ptr< tom::LightTimeConvergenceCriteria > >(
                   &tom::oneWayClosedLoopDoppler ),
           py::arg( "link_ends" ),
           py::arg( "light_time_correction_settings" ) =
                   std::vector< std::shared_ptr< tom::LightTimeCorrectionSettings > >( ),
           py::arg( "bias_settings" ) = nullptr,
           py::arg( "light_time_convergence_settings" ) =
                   std::make_shared< tom::LightTimeConvergenceCriteria >( ) );

    //    m.def("gaussian_noise_function",
    //              &ts::getGaussianDistributionNoiseFunction,
    //          py::arg("standard_deviation"),
    //          py::arg("mean") = 0.0,
    //          py::arg("seed") = time(NULL),
    //          py::arg("observable_size") = 1);

    m.def( "set_vmf_troposphere_data",
           &tom::setVmfTroposphereCorrections,
           py::arg( "data_files" ),
           py::arg( "file_has_meteo" ),
           py::arg( "file_has_gradient" ),
           py::arg( "bodies" ),
           py::arg( "set_troposphere_data" ) = true,
           py::arg( "set_meteo_data" ) = true,
           py::arg( "interpolator_settings" ) = ti::cubicSplineInterpolation( ) );
}

}  // namespace observation
}  // namespace estimation_setup
}  // namespace numerical_simulation
}  // namespace tudatpy
