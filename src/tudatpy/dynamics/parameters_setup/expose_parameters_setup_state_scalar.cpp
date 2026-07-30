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

 The function checks for consistency between the parameter settings and the models contained in the simulation setup.

 Parameters
 ----------
 parameter_settings
     List of settings for the parameters to create.
 bodies
     System of bodies containing the simulation models.
 propagator_settings
     Consolidated propagation settings.
 consider_parameters_names
     Settings for consider parameters.
 print_parameter_order_warning
     Whether to warn when both scalar and vector parameters are present.

 Returns
 -------
 dynamics.parameters.EstimatableParameterSet
     Consolidated estimatable parameter set.
     )doc" );

    m.def( "initial_states",
           &tss::getInitialStateParameterSettings< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "propagator_settings" ),
           py::arg( "bodies" ),
           py::arg( "arc_initial_times" ) = std::vector< double >( ),
           R"doc(

 Function for creating parameter settings for the propagated initial states.

 The propagator settings determine whether single-, multi-, or hybrid-arc initial-state
 parameters are created and which dynamical state types are included.
     )doc" );
}

}  // namespace TUDATPY_STATE_SCALAR_BINDING_NAMESPACE
}  // namespace parameters_setup
}  // namespace dynamics
}  // namespace tudatpy
