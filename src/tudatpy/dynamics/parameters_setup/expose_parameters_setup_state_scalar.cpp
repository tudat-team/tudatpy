/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"

namespace py = pybind11;
namespace tep = tudat::estimatable_parameters;
namespace tp = tudat::propagators;
namespace tss = tudat::simulation_setup;

namespace tudatpy
{
namespace dynamics
{
namespace parameters_setup
{
namespace TUDATPY_STATE_SCALAR_BINDING_NAMESPACE
{

void expose_parameters_setup_state_scalar( py::module& m )
{
    m.def( "create_parameter_set",
           &tss::createParametersToEstimate< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "parameter_settings" ),
           py::arg( "bodies" ),
           py::arg_v( "propagator_settings", std::shared_ptr< tp::PropagatorSettings< STATE_SCALAR_TYPE > >( ), "None" ),
           py::arg( "consider_parameters_names" ) = std::vector< std::shared_ptr< tep::EstimatableParameterSettings > >( ),
           py::arg( "print_parameter_order_warning" ) = true,
           R"doc(

 Function for creating a consolidated parameter from the given estimatable parameter settings.

 Function for creating a consolidated parameter from the given estimatable parameter settings.
 The function checks for consistency between the parameter settings and the models contained in the simulation setup (given by the `bodies` & and `propagator_settings` parameters).


 Parameters
 ----------
 parameter_settings : list[:class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`]
     List of objects that define the settings for the parameters that are to be created. Each entry in this list is typically created by a call to a function in the :ref:`parameters_setup` module.

 bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
     Object consolidating all bodies and environment models, including ground station models, that constitute the physical environment.

 propagator_settings : :class:`~tudatpy.dynamics.propagation_setup.propagator.PropagatorSettings`
     Object containing the consolidated propagation settings of the simulation.

 consider_parameters_names : list[:class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`] = []
     List of objects that define the settings for the considered parameters that are to be created. Each entry in this list is typically created by a call to a function in the :ref:`parameters_setup` module.

 print_parameter_order_warning : bool = True
     Flag to determine whether a warning is printed to the console if there are both scalar and vector parameters found.



 Returns
 -------
 :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet`
     Instance of :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet` class, consolidating all estimatable parameters and simulation models.

 Examples
 --------
 .. code-block:: python

     # Create bodies
     bodies = ...
     # Define parameters settings
     parameter_settings = ...
     # Create the parameters that will be estimated
     parameters_to_estimate = dynamics.parameters_setup.create_parameter_set(parameter_settings, bodies)

 This code snippet closely follows what is done in: `Full Estimation Example <https://github.com/tudat-team/tudatpy-examples/blob/master/estimation/full_estimation_example.ipynb>`_.




     )doc" );

    m.def( "initial_states",
           &tss::getInitialStateParameterSettings< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "propagator_settings" ),
           py::arg( "bodies" ),
           py::arg( "arc_initial_times" ) = std::vector< double >( ),
           R"doc(

 Function for creating parameter settings for initial state parameters.

 Function for creating a parameter settings object for initial state parameters.
 The function uses the propagator settings to determine which type of initial state parameter (single/multi/hybrid-arc; translational/rotational/... dynamics) is to be estimated,
 e.g. if a single-arc translational state propagator is defined, the function will automatically create the parameters for the associated initial state parameter

 .. note::

    This function return lists of parameter settings objects.
    This means that the return of this function cannot simply be added to the parameter settings objects of single parameters in a list creation statement.
    Instead, list concatenation is recommended. Please see the following example:

 .. code-block:: python

    # define single estimatable parameters
    single_parameter_1 = ...
    single_parameter_2 = ...
    ...

    # bad: list creation statement --> will result in nested list, undesired!
    list_of_all_parameters = [dynamics.parameters_setup.initial_states(...), single_parameter_1, single_parameter_2, ...]

    # better: list concatenation --> will result in simple list, desired!
    list_of_all_parameters = dynamics.parameters_setup.initial_states(...) + [single_parameter_1, single_parameter_2, ...]


 Parameters
 ----------
 propagator_settings : :class:`~tudatpy.dynamics.propagation_setup.propagator.PropagatorSettings`
     Object containing the consolidated propagation settings of the simulation in the context of which the given model parameters are to be estimated.

 bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
     Object consolidating all bodies and environment models that constitute the physical environment.

 arc_initial_times : List[ float ] = []
     Initial times of arcs, only required if arc-wise propagation settings are passed via the `propagator_settings` object.

 Returns
 -------
 List[ :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` ]
     List of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` objects, one per component of each initial state in the simulation.


     )doc" );
}

}  // namespace TUDATPY_STATE_SCALAR_BINDING_NAMESPACE
}  // namespace parameters_setup
}  // namespace dynamics
}  // namespace tudatpy
