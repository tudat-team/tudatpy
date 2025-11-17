/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_parameters.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"

namespace py = pybind11;
namespace tep = tudat::estimatable_parameters;

namespace tudatpy
{
namespace dynamics
{
namespace parameters
{

void expose_parameters( py::module& m )
{
 
    // ESTIMATABLE PARAMETER SET CLASS

    py::class_< tep::EstimatableParameterSet< STATE_SCALAR_TYPE >,
                std::shared_ptr< tep::EstimatableParameterSet< STATE_SCALAR_TYPE > > >(
            m,
            "EstimatableParameterSet",
            R"doc(

         Class containing a consolidated set of estimatable parameters.

         Class containing a consolidated set of estimatable parameters, linked to the environment and acceleration settings of the simulation.
         The user typically creates instances of this class via the :func:`~tudatpy.dynamics.parameters_setup.create_parameter_set` function.





      )doc" )
            .def_property_readonly( "parameter_set_size",
                                    &tep::EstimatableParameterSet<
                                            STATE_SCALAR_TYPE >::getEstimatedParameterSetSize,
                                    R"doc(

         **read-only**

         Size of the parameter set, i.e. amount of estimatable parameters contained in the set.

         :type: int
      )doc" )
            .def_property_readonly(
                    "initial_states_size",
                    &tep::EstimatableParameterSet<
                            STATE_SCALAR_TYPE >::getInitialDynamicalStateParameterSize,
                    R"doc(

         **read-only**

         Amount of initial state parameters contained in the set.

         :type: int
      )doc" )
            .def_property_readonly(
                    "initial_single_arc_states_size",
                    &tep::EstimatableParameterSet<
                            STATE_SCALAR_TYPE >::getInitialDynamicalSingleArcStateParameterSize,
                    R"doc(

         **read-only**

         Amount of initial state parameters in the set, which are treated in a single-arc fashion.

         :type: int
      )doc" )
            .def_property_readonly(
                    "initial_multi_arc_states_size",
                    &tep::EstimatableParameterSet<
                            STATE_SCALAR_TYPE >::getInitialDynamicalMultiArcStateParameterSize,
                    R"doc(

         **read-only**

         Amount of initial state parameters in the set, which are treated in a multi-arc fashion.

         :type: int
      )doc" )
            .def_property_readonly(
                    "constraints_size",
                    &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::getConstraintSize,
                    R"doc(

         **read-only**

         Total size of linear constraint that is to be applied during estimation.

         :type: int
      )doc" )
            .def_property( "parameter_vector",
                           &tep::EstimatableParameterSet<
                                   STATE_SCALAR_TYPE >::getFullParameterValuesWithoutConsiderParameters< double >,
                           &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::resetParameterValuesWithoutConsiderParameters<
                                   double >,
                           R"doc(

         Vector containing the parameter values of all parameters in the set (excluding consider parameters).

         :type: numpy.ndarray[numpy.float64[m, 1]]
      )doc" )
            .def_property( "estimated_and_consider_parameter_vector",
                           &tep::EstimatableParameterSet<
                                   STATE_SCALAR_TYPE >::getFullParameterValuesWithConsiderParameters< double >,
                           &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::resetParameterValuesWithConsiderParameters<
                                   double >,
                           R"doc(

         Vector containing the parameter values of all parameters in the set (including consider parameters; concatenated after estimated parameters).

         :type: numpy.ndarray[numpy.float64[m, 1]]
      )doc" )
            .def_property_readonly("consider_parameters",
                                   &tep::EstimatableParameterSet<
                                           STATE_SCALAR_TYPE >::getConsiderParameters,
                                   R"doc(
                                   
        **read-only**

        Set of consider parameters that are included in the parameter set.

        :type: :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet`

                                   )doc")
            .def( "indices_for_parameter_type",
                  &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::getIndicesForParameterType,
                  py::arg( "parameter_type" ),
                  R"doc(

         Function to retrieve the indices of a given type of parameter.

         Function to retrieve the index of all parameters of a given type from the parameter set.
         This function can be very useful, since the order of parameters within the parameter set does not necessarily correspond to the order in which the elements were added to the set.


         Parameters
         ----------
         parameter_type : Tuple[ :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterTypes`, Tuple[str, str] ]
             Parameter identifier for which the indices are retrieved. The first element of the tuple is the parameter type, the second element is a tuple containing the body name and, if applicable, a secondary identifier (e.g., station name), else this is an empty string.
         Returns
         -------
         List[ Tuple[int, int] ]
             Indices of the parameters corresponding to the description. The first element of the tuple is the start index, the second element is the size of the parameter.

     )doc" )
        .def("indices_for_parameter_description",
                 &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::getIndicesForParameterDescription,
                 py::arg("parameter_description"),
                 R"doc(

         Function to retrieve the indices of a given type of parameter.

         Function to retrieve the index of all parameters of a given type from the parameter set.
         This function can be very useful, since the order of parameters within the parameter set does not necessarily correspond to the order in which the elements were added to the set.


         Parameters
         ----------
         parameter_description : str
                Parameter description that is to be searched in full list of parameters.
                Parameter descriptions can be retrieved from the :meth:`tudatpy.dynamics.parameters.EstimatableParameterSet.get_parameter_descriptions` method.

         Returns
         -------
         List[ Tuple[int, int] ]
             Indices of the parameters corresponding to the description. The first element of the tuple is the start index, the second element is the size of the parameter.

     )doc")
            .def( "get_parameter_descriptions",
                  &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::getParametersDescriptions,
                  R"doc(

         Function to retrieve the descriptions of all parameters in the set.

         Returns
         -------
         List[ str ]
             List of parameter descriptions formatted as strings.


         )doc")
            .def( "get_parameter_identifiers",
                  &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::getParametersIdentifiers,
                  R"doc(

         Function to retrieve the identifiers of all parameters in the set.

         Returns
         -------
         Tuple[ :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterTypes`, Tuple[str, str] ]
             List of parameter identifiers. The first element of the tuple is the parameter type, the second element is a tuple containing the body name and, if applicable, a secondary identifier (e.g., station name), else this is an empty string.


         )doc");


     
    m.def( "print_parameter_names",
           &tep::printEstimatableParameterEntries< double >,
           py::arg( "parameter_set" ),
           R"doc(

 Function for printing a list of estimatable parameter names.

 Function that allows you to print a verbose list of all parameters that shall be estimated. Consider parameters are listed separately.


 Parameters
 ----------
 parameter_set : :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet`.
     Instance of :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet` class, consolidating all estimatable parameters and simulation models.

 Returns
 -------
 List[]
     Verbose List of all parameters that shall be estimated. Consider parameters are listed separately.

 Examples
 --------
 .. code-block:: python

     import numpy as np
     from tudatpy.interface import spice
     from tudatpy.dynamics import environment_setup, propagation_setup, parameters_setup, parameters
     from tudatpy
     from tudatpy.astro.time_representation import DateTime

     # Load SPICE kernels
     spice.load_standard_kernels()

     # Set simulation epochs
     simulation_start_epoch = DateTime(2000, 1, 1).epoch()
     simulation_end_epoch = DateTime(2000, 1, 4).epoch()

     # Create bodies
     bodies_to_create = ["Sun", "Earth"]
     global_frame_origin = "Earth"
     global_frame_orientation = "J2000"
     body_settings = environment_setup.get_default_body_settings(
       bodies_to_create, global_frame_origin, global_frame_orientation)

     # Create vehicle
     body_settings.add_empty_settings("Delfi-C3")
     body_settings.get("Delfi-C3").constant_mass = 2.2
     bodies = environment_setup.create_system_of_bodies(body_settings)

     # Define propagation settings
     bodies_to_propagate = ["Delfi-C3"]
     central_bodies = ["Earth"]
     accelerations_settings_delfi_c3 = dict(
       Sun=[propagation_setup.acceleration.point_mass_gravity()],
       Earth=[propagation_setup.acceleration.spherical_harmonic_gravity(8, 8)]
     )
     acceleration_settings = {"Delfi-C3": accelerations_settings_delfi_c3}
     acceleration_models = propagation_setup.create_acceleration_models(
       bodies, acceleration_settings, bodies_to_propagate, central_bodies)
     initial_state = np.zeros(6)  # Use a real initial state if needed

     # Integrator settings
     integrator_settings = propagation_setup.integrator.runge_kutta_fixed_step(60.0,
                                                                                    coefficient_set=propagation_setup.integrator.CoefficientSets.rkdp_87)

     # Create propagator
     termination_condition = propagation_setup.propagator.time_termination(simulation_end_epoch)
     propagator_settings = propagation_setup.propagator.translational(
       central_bodies, acceleration_models, bodies_to_propagate, initial_state,
       simulation_start_epoch, integrator_settings, termination_condition)

     # Define parameters to estimate
     parameter_settings = parameters_setup.initial_states(propagator_settings, bodies)
     parameter_settings.append(parameters_setup.gravitational_parameter("Earth"))
     parameters_to_estimate = parameters_setup.create_parameter_set(parameter_settings, bodies)

     # Print parameter names
     print(parameters.print_parameter_names(parameters_to_estimate))




     )doc" );
}

}  // namespace parameters
}  // namespace dynamics
}  // namespace tudatpy
