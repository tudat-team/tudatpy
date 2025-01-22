/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_estimation.h"

#include "tudat/simulation/estimation_setup/fitOrbitToEphemeris.h"
#include "tudat/astro/propagators/propagateCovariance.h"
#include "tudat/basics/utilities.h"

#include "docstrings.h"
#include "scalarTypes.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

namespace py = pybind11;
namespace tep = tudat::estimatable_parameters;
namespace tp = tudat::propagators;
namespace ts = tudat::statistics;
namespace tss = tudat::simulation_setup;
namespace tni = tudat::numerical_integrators;
namespace tom = tudat::observation_models;
namespace trf = tudat::reference_frames;
namespace te = tudat::ephemerides;


namespace tudatpy {
namespace numerical_simulation {
namespace estimation {



void expose_single_observation_set(py::module &m) {


    py::class_< tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>,
            std::shared_ptr<tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>>>(m, "SingleObservationSet",
                                                          get_docstring("SingleObservationSet").c_str() )
            .def("set_observations", py::overload_cast< const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >&>(
                    &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::setObservations),
                    py::arg("observations"), get_docstring("SingleObservationSet.set_observations").c_str() )
            .def("set_observations", py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >& >(
                    &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::setObservations),
                    py::arg("observations"), get_docstring("SingleObservationSet.set_observations").c_str() )
            .def("set_residuals", py::overload_cast< const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >&>(
                    &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::setResiduals),
                    py::arg("residuals"), get_docstring("SingleObservationSet.set_residuals").c_str() )
            .def("set_residuals", py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >& >(
                    &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::setResiduals),
                    py::arg("residuals"), get_docstring("SingleObservationSet.set_residuals").c_str() )
            .def("set_constant_weight", py::overload_cast< const double >( &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::setConstantWeight ),
                    py::arg("weight"), get_docstring("SingleObservationSet.set_constant_weight").c_str() )
            .def("set_constant_weight", py::overload_cast< const Eigen::Matrix< double, Eigen::Dynamic, 1 >& >( &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::setConstantWeight ),
                 py::arg("weight"), get_docstring("SingleObservationSet.set_constant_weight").c_str() )
            .def("set_tabulated_weights", py::overload_cast< const Eigen::VectorXd& >( &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::setTabulatedWeights ),
                 py::arg("weights"), get_docstring("SingleObservationSet.set_tabulated_weights").c_str() )
            .def("filter_observations", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::filterObservations,
                 py::arg("filter"), py::arg("save_filtered_obs") = true,
                 get_docstring("SingleObservationSet.filter_observations").c_str() )
            .def_property_readonly("observable_type", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getObservableType,
                                   get_docstring("SingleObservationSet.observable_type").c_str() )
            .def_property("link_definition",
                          &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getLinkEnds,
                          &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::setLinkEnds,
                          get_docstring("SingleObservationSet.link_definition").c_str() )
            .def_property_readonly("reference_link_end", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getReferenceLinkEnd,
                                   get_docstring("SingleObservationSet.reference_link_end").c_str() )
            .def_property_readonly("number_of_observables", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getNumberOfObservables,
                                   get_docstring("SingleObservationSet.number_of_observables").c_str())
            .def_property_readonly("single_observable_size", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getSingleObservableSize,
                                   get_docstring("SingleObservationSet.single_observaable_size").c_str())
            .def_property_readonly("total_observation_set_size", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getTotalObservationSetSize,
                                   get_docstring("SingleObservationSet.total_observation_set_size").c_str())
            .def_property_readonly("time_bounds", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getTimeBounds,
                                   get_docstring("SingleObservationSet.time_bounds").c_str())
            .def_property_readonly("list_of_observations", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getObservations,
                                   get_docstring("SingleObservationSet.list_of_observations").c_str() )
            .def_property_readonly("observation_times", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getObservationTimes,
                                   get_docstring("SingleObservationSet.observation_times").c_str() )
            .def_property_readonly("concatenated_observations", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getObservationsVector,
                                   get_docstring("SingleObservationSet.concatenated_observations").c_str() )
            .def_property_readonly("computed_observations", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getComputedObservations,
                                   get_docstring("SingleObservationSet.computed_observations").c_str() )
            .def_property_readonly("concatenated_computed_observations",
                                   &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getComputedObservationsVector,
                                   get_docstring("SingleObservationSet.concatenated_computed_observations").c_str() )
            .def_property_readonly("residuals", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getResiduals,
                                   get_docstring("SingleObservationSet.residuals").c_str() )
            .def_property_readonly("concatenated_residuals", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getResidualsVector,
                                   get_docstring("SingleObservationSet.concatenated_residuals").c_str() )
            .def_property_readonly("rms_residuals", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getRmsResiduals,
                                   get_docstring("SingleObservationSet.rms_residuals").c_str() )
            .def_property_readonly("mean_residuals", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getMeanResiduals,
                                   get_docstring("SingleObservationSet.mean_residuals").c_str() )
            .def_property_readonly("weights", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getWeights,
                                   get_docstring("SingleObservationSet.weights").c_str())
            .def_property_readonly("concatenad_weights", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getWeightsVector,
                                   get_docstring("SingleObservationSet.concatenad_weights").c_str())
            .def_property("dependent_variables", &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsDependentVariables,
                                                 &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setObservationsDependentVariables,
                                   get_docstring("SingleObservationSet.dependent_variables").c_str())
            .def_property_readonly("dependent_variables_history", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getDependentVariableHistory,
                                   get_docstring("SingleObservationSet.dependent_variables_history").c_str())
            .def_property_readonly("observations_history", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getObservationsHistory,
                                   get_docstring("SingleObservationSet.observations_history").c_str() )
            .def_property_readonly("ancilliary_settings", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getAncilliarySettings,
                                   get_docstring("SingleObservationSet.ancilliary_settings").c_str() )
            .def_property("weights_vector", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getWeightsVector,
                                            &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::setTabulatedWeights,
                                   get_docstring("SingleObservationSet.weights_vector").c_str() )
            .def_property_readonly("filtered_observation_set", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getFilteredObservationSet,
                                    get_docstring("SingleObservationSet.filtered_observation_set").c_str() )
            .def_property_readonly("number_filtered_observations", &tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>::getNumberOfFilteredObservations,
                                   get_docstring("SingleObservationSet.number_filtered_observations").c_str() )
            .def( "single_dependent_variable",
                  py::overload_cast< std::shared_ptr< tss::ObservationDependentVariableSettings >, const bool >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleDependentVariable ),
                          py::arg( "dependent_variable_settings" ),
                          py::arg( "return_first_compatible_settings" ) = false,
                          get_docstring("SingleObservationSet.single_dependent_variable").c_str() )
            .def( "compatible_dependent_variable_settings",
                  &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getCompatibleDependentVariablesSettingsList,
                  get_docstring( "SingleObservationSet.compatible_dependent_variable_settings" ).c_str() )
            .def( "compatible_dependent_variables_list",
                  &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getAllCompatibleDependentVariables,
                  get_docstring( "SingleObservationSet.compatible_dependent_variables_list" ).c_str() )
            .def( "single_dependent_variable_history",
                  &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleDependentVariableHistory,
                  get_docstring( "SingleObservationSet.single_dependent_variable_history" ).c_str() )
            .def_property_readonly( "dependent_variables_matrix",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsDependentVariablesMatrix,
                                    get_docstring( "SingleObservationSet.dependent_variables_matrix" ).c_str() );


    m.def("single_observation_set",
          &tss::singleObservationSetWithoutDependentVariables< STATE_SCALAR_TYPE, TIME_TYPE >,
          py::arg("observable_type"),
          py::arg("link_definition" ),
          py::arg("observations" ),
          py::arg("observation_times"),
          py::arg("reference_link_end"),
          py::arg("ancilliary_settings") = nullptr,
          get_docstring("single_observation_set").c_str() );



}

}
}
}// namespace tudatpy
