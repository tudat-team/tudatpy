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
#include "expose_model_settings.h"
#include <pybind11/functional.h>
#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"

namespace tom = tudat::observation_models;

namespace tudatpy
{
namespace estimation_refactoring
{
namespace observable_models_setup
{

namespace model_settings
{

void expose_model_settings( py::module& m )
{
    
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
            .value( "relative_position_observable_type", tom::ObservableType::relative_position_observable )
            .value( "one_way_instantaneous_doppler_type", tom::ObservableType::one_way_doppler )
            .value( "one_way_averaged_doppler_type", tom::ObservableType::one_way_differenced_range )
            .value( "two_way_instantaneous_doppler_type", tom::ObservableType::two_way_doppler )
            .value( "n_way_averaged_doppler_type", tom::ObservableType::n_way_differenced_range )
            .value( "euler_angle_313_observable_type", tom::ObservableType::euler_angle_313_observable )
            .value( "dsn_one_way_averaged_doppler", tom::ObservableType::dsn_one_way_averaged_doppler )
            .value( "dsn_n_way_averaged_doppler", tom::ObservableType::dsn_n_way_averaged_doppler )
            .value( "doppler_measured_frequency_type", tom::ObservableType::doppler_measured_frequency )
            .value( "dsn_n_way_range", tom::ObservableType::dsn_n_way_range )
            .export_values( );


    py::class_< tom::DopplerProperTimeRateSettings, std::shared_ptr< tom::DopplerProperTimeRateSettings > >(
            m,
            "DopplerProperTimeRateSettings",
            R"doc(

         Base class to define proper time rate settings.

         Base class to define proper time rate settings (at a single link end) for instantaneous Doppler observation model settings.
         Specific proper time rate settings must be defined using an object derived from this class.
         The derived classes are made accessible via dedicated functions.





      )doc" );

    py::class_< tom::ObservationModelSettings, std::shared_ptr< tom::ObservationModelSettings > >( m, "ObservationSettings", R"doc(

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

    py::class_< tom::NWayRangeObservationSettings, std::shared_ptr< tom::NWayRangeObservationSettings >, tom::ObservationModelSettings >(
            m, "NWayRangeObservationSettings", R"doc(No documentation found.)doc" );

}

}
}
}
}