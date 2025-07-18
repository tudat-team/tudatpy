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
#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "expose_estimation.h"
#include "scalarTypes.h"
#include "tudat/astro/propagators/propagateCovariance.h"
#include "tudat/basics/utilities.h"
#include "tudat/simulation/estimation_setup/fitOrbitToEphemeris.h"

namespace py = pybind11;
namespace tep = tudat::estimatable_parameters;
namespace tp = tudat::propagators;
namespace ts = tudat::statistics;
namespace tss = tudat::simulation_setup;
namespace tni = tudat::numerical_integrators;
namespace tom = tudat::observation_models;
namespace trf = tudat::reference_frames;
namespace te = tudat::ephemerides;



namespace tudatpy
{
namespace numerical_simulation
{
namespace estimation
{

void expose_estimation_propagated_covariance( py::module& m )
{
    /*!
     *************** PARAMETERS ***************
     */

     // SAM NOTE should move to dynamics.parameters
    py::class_< tep::EstimatableParameterSet< STATE_SCALAR_TYPE >,
                std::shared_ptr< tep::EstimatableParameterSet< STATE_SCALAR_TYPE > > >(
            m,
            "EstimatableParameterSet",
            R"doc(

         Class containing a consolidated set of estimatable parameters.

         Class containing a consolidated set of estimatable parameters, linked to the environment and acceleration settings of the simulation.
         The user typically creates instances of this class via the :func:`~tudatpy.numerical_simulation.estimation_setup.create_parameter_set` function.





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
                                   STATE_SCALAR_TYPE >::getFullParameterValues< double >,
                           &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::resetParameterValues<
                                   double >,
                           R"doc(

         Vector containing the parameter values of all parameters in the set.

         :type: numpy.ndarray[numpy.float64[m, 1]]
      )doc" )
            .def( "indices_for_parameter_type",
                  &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::getIndicesForParameterType,
                  py::arg( "parameter_type" ),
                  R"doc(

         Function to retrieve the indices of a given type of parameter.

         Function to retrieve the index of all parameters of a given type from the parameter set.
         This function can be very useful, since the order of parameters within the parameter set does not necessarily correspond to the order in which the elements were added to the set.


         Parameters
         ----------
         parameter_type : Tuple[ :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterTypes`, Tuple[str, str] ]
             help
         Returns
         -------
         List[ Tuple[int, int] ]
             help





     )doc" );

    /*!
     *************** STATE TRANSITION INTERFACE ***************
     */

    // SAM NOTE: should also move to dynamics            
    py::class_< tp::CombinedStateTransitionAndSensitivityMatrixInterface,
                std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface > >(
            m,
            "CombinedStateTransitionAndSensitivityMatrixInterface",
            R"doc(

         Class establishing an interface with the simulation's State Transition and Sensitivity Matrices.

         Class establishing an interface to the State Transition and Sensitivity Matrices.
         Instances of this class are instantiated automatically upon creation of :class:`~tudatpy.numerical_simulation.Estimator` objects,
         using the simulation information in the observation, propagation and integration settings that the :class:`~tudatpy.numerical_simulation.Estimator` instance is linked to.





      )doc" )
            .def( "state_transition_sensitivity_at_epoch",
                  &tp::CombinedStateTransitionAndSensitivityMatrixInterface::
                          getCombinedStateTransitionAndSensitivityMatrix,
                  py::arg( "time" ),
                  py::arg( "add_central_body_dependency" ) = true,
                  py::arg( "arc_defining_bodies" ) = std::vector< std::string >( ),
                  R"doc(

         Function to get the concatenated state transition and sensitivity matrix at a given time.

         Function to get the concatenated state transition and sensitivity matrix at a given time.
         Entries corresponding to parameters which are not active at the current arc are omitted.


         Parameters
         ----------
         time : float
             Time at which concatenated state transition and sensitivity matrix are to be retrieved.
         Returns
         -------
         numpy.ndarray[numpy.float64[m, n]]
             Concatenated state transition and sensitivity matrix at a given time.





     )doc" )
            .def( "full_state_transition_sensitivity_at_epoch",
                  &tp::CombinedStateTransitionAndSensitivityMatrixInterface::
                          getFullCombinedStateTransitionAndSensitivityMatrix,
                  py::arg( "time" ),
                  py::arg( "add_central_body_dependency" ) = true,
                  py::arg( "arc_defining_bodies" ) = std::vector< std::string >( ),
                  R"doc(


         Parameters
         ----------
         time : float
             Time at which full concatenated state transition and sensitivity matrix are to be retrieved.
         Returns
         -------
         numpy.ndarray[numpy.float64[m, n]]
             Full concatenated state transition and sensitivity matrix at a given time.





     )doc" )
            .def_property_readonly( "state_transition_size",
                                    &tp::CombinedStateTransitionAndSensitivityMatrixInterface::
                                            getStateTransitionMatrixSize,
                                    R"doc(

         **read-only**

         Size of the (square) state transition matrix.

         :type: int
      )doc" )
            .def_property_readonly( "sensitivity_size",
                                    &tp::CombinedStateTransitionAndSensitivityMatrixInterface::
                                            getSensitivityMatrixSize,
                                    R"doc(

         **read-only**

         Number of columns in the sensitivity matrix.

         :type: int
      )doc" )
            .def_property_readonly( "full_parameter_size",
                                    &tp::CombinedStateTransitionAndSensitivityMatrixInterface::
                                            getFullParameterVectorSize,
                                    R"doc(

         **read-only**

         Full amount of parameters w.r.t. which partials have been set up via State Transition and Sensitivity Matrices.

         :type: int
      )doc" );

    
}

}  // namespace estimation
}  // namespace numerical_simulation
}  // namespace tudatpy
