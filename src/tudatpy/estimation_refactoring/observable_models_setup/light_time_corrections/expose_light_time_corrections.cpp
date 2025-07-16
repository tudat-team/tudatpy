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
#include "expose_light_time_corrections.h"
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

namespace light_time_corrections
{

void expose_light_time_corrections( py::module& m )
{
    py::class_< tom::LightTimeConvergenceCriteria, std::shared_ptr< tom::LightTimeConvergenceCriteria > >( m,
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

    py::class_< tom::LightTimeCorrectionSettings, std::shared_ptr< tom::LightTimeCorrectionSettings > >( m,
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
}

}
}
}
}