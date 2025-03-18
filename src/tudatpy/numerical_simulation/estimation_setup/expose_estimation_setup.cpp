/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tom = tudat::observation_models;
namespace tep = tudat::estimatable_parameters;

namespace tudatpy
{
namespace numerical_simulation
{
namespace estimation_setup
{

PYBIND11_MODULE( expose_estimation_setup, m )
{
    // *************** PARAMETER ***************

    m.def( "print_parameter_names",
           &tep::printEstimatableParameterEntries< double >,
           py::arg( "parameter_set" ),
           R"doc(

 Function for printing a list of estimatable parameter names.

 Function that allows you to print a verbose list of all parameters that shall be estimated. Consider parameters are listed separately.


 Parameters
 ----------
 parameter_set : :class:`~tudatpy.numerical_simulation.estimation.EstimatableParameterSet`.
     Instance of :class:`~tudatpy.numerical_simulation.estimation.EstimatableParameterSet` class, consolidating all estimatable parameters and simulation models.

 Returns
 -------
 List[]
     Verbose List of all parameters that shall be estimated. Consider parameters are listed separately.

 Examples
 --------
 .. code-block:: python

     import numpy as np
     from tudatpy.interface import spice
     from tudatpy.numerical_simulation import environment_setup, propagation_setup, estimation_setup
     from tudatpy.astro.time_conversion import DateTime

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
     integrator_settings = propagation_setup.integrator.runge_kutta_fixed_step_size(60.0,
                                                                                    coefficient_set=propagation_setup.integrator.CoefficientSets.rkdp_87)

     # Create propagator
     termination_condition = propagation_setup.propagator.time_termination(simulation_end_epoch)
     propagator_settings = propagation_setup.propagator.translational(
       central_bodies, acceleration_models, bodies_to_propagate, initial_state,
       simulation_start_epoch, integrator_settings, termination_condition)

     # Define parameters to estimate
     parameter_settings = estimation_setup.parameter.initial_states(propagator_settings, bodies)
     parameter_settings.append(estimation_setup.parameter.gravitational_parameter("Earth"))
     parameters_to_estimate = estimation_setup.create_parameter_set(parameter_settings, bodies)

     # Print parameter names
     print(estimation_setup.print_parameter_names(parameters_to_estimate))




     )doc" );

    // # EstimatableParameterSettings --> EstimatableParameterSet #
    m.def( "create_parameter_set",
           &tss::createParametersToEstimate< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "parameter_settings" ),
           py::arg( "bodies" ),
           py::arg( "propagator_settings" ) = nullptr,
           py::arg( "consider_parameters_names" ) =
                   std::vector< std::shared_ptr< tep::EstimatableParameterSettings > >( ),
           R"doc(

 Function for creating a consolidated parameter from the given estimatable parameter settings.

 Function for creating a consolidated parameter from the given estimatable parameter settings.
 The function checks for consistency between the parameter settings and the models contained in the simulation setup (given by the `bodies` & and `propagator_settings` parameters).


 Parameters
 ----------
 parameter_settings : list( :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterSettings` )
     List of objects that define the settings for the parameters that are to be created. Each entry in this list is typically created by a call to a function in the :ref:`\`\`parameter\`\`` module

 bodies : :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies`
     Object consolidating all bodies and environment models, including ground station models, that constitute the physical environment.

 propagator_settings : :class:`~tudatpy.numerical_simulation.propagation_setup.propagator.PropagatorSettings`
     Object containing the consolidated propagation settings of the simulation.

 Returns
 -------
 :class:`~tudatpy.numerical_simulation.estimation.EstimatableParameterSet`
     Instance of :class:`~tudatpy.numerical_simulation.estimation.EstimatableParameterSet` class, consolidating all estimatable parameters and simulation models.

 Examples
 --------
 .. code-block:: python

     # Create bodies
     bodies = ...
     # Define parameters settings
     parameter_settings = ...
     # Create the parameters that will be estimated
     parameters_to_estimate = estimation_setup.create_parameter_set(parameter_settings, bodies)

 This code snippet closely follows what is done in: `Full Estimation Example <https://github.com/tudat-team/tudatpy-examples/blob/master/estimation/full_estimation_example.ipynb>`_.




     )doc" );

    // ************** OBSERVATION ***************

    // #   Observation Model Settings --> Observation Simulator #
    m.def( "create_observation_simulators",
           py::overload_cast<
                   const std::vector< std::shared_ptr< tom::ObservationModelSettings > >&,
                   const tss::SystemOfBodies& >(
                   &tom::createObservationSimulators< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "observation_settings" ),
           py::arg( "bodies" ),
           R"doc(

 Function for creating observation simulator objects.

 Function for creating observation simulator objects from observation settings.
 Note that each observation (i.e. combination of observable and link geometry) requires its own observation simulator object.


 Parameters
 ----------
 observation_settings : List[ ObservationSettings ]
     List of settings objects, each object defining the observation model settings for one combination of observable and link geometry that is to be simulated.

 bodies : :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies`
     Object consolidating all bodies and environment models, including ground station models, that constitute the physical environment.

 Returns
 -------
 List[ :class:`~tudatpy.numerical_simulation.estimation.ObservationSimulator` ]
     List of :class:`~tudatpy.numerical_simulation.estimation.ObservationSimulator` objects, each object hosting the functionality for simulating one combination of observable type and link geometry.

 Examples
 --------
 .. code-block:: python

     # Create bodies
     bodies = ...
     # Define parameters settings
     observation_settings = ...
     # Create observation simulators
     observation_simulators = estimation_setup.create_observation_simulators(observation_settings, bodies)

 This code snippet closely follows what is done in: The following snippet closely follows what is done in: `Galilean Moons State Estimation Example <https://github.com/tudat-team/tudatpy-examples/blob/master/estimation/galilean_moons_state_estimation.ipynb>`_.



     )doc" );

    m.def( "single_type_observation_collection",
           py::overload_cast<
                   const tom::ObservableType,
                   const tom::LinkDefinition&,
                   const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >&,
                   const std::vector< TIME_TYPE >,
                   const tom::LinkEndType,
                   const std::shared_ptr< tom::ObservationAncilliarySimulationSettings > >(
                   &tom::createManualObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "observable_type" ),
           py::arg( "link_ends" ),
           py::arg( "observations_list" ),
           py::arg( "times_list" ),
           py::arg( "reference_link_end" ),
           py::arg( "ancilliary_settings" ) = nullptr,
           R"doc(No documentation found.)doc" );
}

}  // namespace estimation_setup
}  // namespace numerical_simulation
}  // namespace tudatpy
