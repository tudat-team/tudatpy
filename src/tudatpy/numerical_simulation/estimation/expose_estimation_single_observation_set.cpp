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

namespace tudat
{

namespace simulation_setup
{

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< tom::SingleObservationSet< ObservationScalarType, TimeType > >
singleObservationSetWithoutDependentVariables(
        const tom::ObservableType observableType,
        const tom::LinkDefinition& linkEnds,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >&
                observations,
        const std::vector< TimeType > observationTimes,
        const tom::LinkEndType referenceLinkEnd,
        const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings >
                ancilliarySettings = nullptr )
{
    return std::make_shared< tom::SingleObservationSet< ObservationScalarType, TimeType > >(
            observableType,
            linkEnds,
            observations,
            observationTimes,
            referenceLinkEnd,
            std::vector< Eigen::VectorXd >( ),
            nullptr,
            ancilliarySettings );
}

}  // namespace simulation_setup

}  // namespace tudat

namespace tudatpy
{
namespace numerical_simulation
{
namespace estimation
{

void expose_estimation_single_observation_set( py::module& m )
{
    py::class_< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > > >(
            m,
            "SingleObservationSet",
            R"doc(

         Class collecting a single set of observations and associated data, of a given observable type, link ends, and ancilliary data.





      )doc" )
            .def( "set_observations",
                  py::overload_cast< const std::vector<
                          Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >& >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                                      TIME_TYPE >::setObservations ),
                  py::arg( "observations" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_observations",
                  py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >& >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                                      TIME_TYPE >::setObservations ),
                  py::arg( "observations" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_residuals",
                  py::overload_cast< const std::vector<
                          Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >& >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                                      TIME_TYPE >::setResiduals ),
                  py::arg( "residuals" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_residuals",
                  py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >& >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                                      TIME_TYPE >::setResiduals ),
                  py::arg( "residuals" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_constant_weight",
                  py::overload_cast< const double >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                                      TIME_TYPE >::setConstantWeight ),
                  py::arg( "weight" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_constant_weight",
                  py::overload_cast< const Eigen::Matrix< double, Eigen::Dynamic, 1 >& >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                                      TIME_TYPE >::setConstantWeight ),
                  py::arg( "weight" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_tabulated_weights",
                  py::overload_cast< const Eigen::VectorXd& >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                                      TIME_TYPE >::setTabulatedWeights ),
                  py::arg( "weights" ),
                  R"doc(No documentation found.)doc" )
            .def( "filter_observations",
                  &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::filterObservations,
                  py::arg( "filter" ),
                  py::arg( "save_filtered_obs" ) = true,
                  R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "observable_type",
                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservableType,
                    R"doc(

         **read-only**

         Type of observable for which the object stores observations

         :type: ObservableType
      )doc" )
            .def_property( "link_definition",
                           &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkEnds,
                           &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setLinkEnds,
                           R"doc(

         **read-only**

         Definition of the link ends for which the object stores observations

         :type: LinkDefinition
      )doc" )
            .def_property_readonly(
                    "reference_link_end",
                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getReferenceLinkEnd,
                    R"doc(

         **read-only**

         Reference link end for all stored observations

         :type: LinkEndType
      )doc" )
            .def_property_readonly( "number_of_observables",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                                                TIME_TYPE >::getNumberOfObservables,
                                    R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "single_observable_size",
                    &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                                TIME_TYPE >::getSingleObservableSize,
                    R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "total_observation_set_size",
                    &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                                TIME_TYPE >::getTotalObservationSetSize,
                    R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "time_bounds",
                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBounds,
                    R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "list_of_observations",
                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservations,
                    R"doc(

         **read-only**

         List of separate stored observations. Each entry of this list is a vector containing a single observation. In cases where the observation is single-valued (range, Doppler), the vector is size 1, but for multi-valued observations such as angular position, each vector in the list will have size >1

         :type: list[ numpy.ndarray[numpy.float64[m, 1]] ]
      )doc" )
            .def_property_readonly(
                    "observation_times",
                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTimes,
                    R"doc(

         **read-only**

         Reference time for each of the observations in ``list_of_observations``

         :type: list[ float]
      )doc" )
            .def_property_readonly( "concatenated_observations",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                                                TIME_TYPE >::getObservationsVector,
                                    R"doc(

         **read-only**

         Concatenated vector of all stored observations

         :type: numpy.ndarray[numpy.float64[m, 1]]
      )doc" )
            .def_property_readonly(
                    "computed_observations",
                    &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                                TIME_TYPE >::getComputedObservations,
                    R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "concatenated_computed_observations",
                    &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                                TIME_TYPE >::getComputedObservationsVector,
                    R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "residuals",
                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getResiduals,
                    R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "concatenated_residuals",
                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getResidualsVector,
                    R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "rms_residuals",
                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getRmsResiduals,
                    R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "mean_residuals",
                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getMeanResiduals,
                    R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "weights",
                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getWeights,
                    R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "concatenad_weights",
                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightsVector,
                    R"doc(No documentation found.)doc" )
            .def_property(
                    "dependent_variables",
                    &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                                TIME_TYPE >::getObservationsDependentVariables,
                    &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                                TIME_TYPE >::setObservationsDependentVariables,
                    R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "dependent_variables_history",
                    &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                                TIME_TYPE >::getDependentVariableHistory,
                    R"doc(No documentation found.)doc" )
            .def_property_readonly( "observations_history",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                                                TIME_TYPE >::getObservationsHistory,
                                    R"doc(

         **read-only**

         Dictionary of observations sorted by time. Created by making a dictionary with ``observation_times`` as keys and  ``list_of_observations`` as values

         :type: dict[ float, numpy.ndarray[numpy.float64[m, 1]] ]
      )doc" )
            .def_property_readonly( "ancilliary_settings",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                                                TIME_TYPE >::getAncilliarySettings,
                                    R"doc(

         **read-only**

         Ancilliary settings all stored observations

         :type: ObservationAncilliarySimulationSettings
      )doc" )
            .def_property(
                    "weights_vector",
                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightsVector,
                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setTabulatedWeights,
                    R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "filtered_observation_set",
                    &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                                TIME_TYPE >::getFilteredObservationSet,
                    R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "number_filtered_observations",
                    &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                                TIME_TYPE >::getNumberOfFilteredObservations,
                    R"doc(No documentation found.)doc" )
            .def( "single_dependent_variable",
                  py::overload_cast< std::shared_ptr< tss::ObservationDependentVariableSettings >,
                                     const bool >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                                      TIME_TYPE >::getSingleDependentVariable ),
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "return_first_compatible_settings" ) = false,
                  R"doc(No documentation found.)doc" )
            .def( "compatible_dependent_variable_settings",
                  &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::
                          getCompatibleDependentVariablesSettingsList,
                  R"doc(No documentation found.)doc" )
            .def( "compatible_dependent_variables_list",
                  &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                              TIME_TYPE >::getAllCompatibleDependentVariables,
                  R"doc(No documentation found.)doc" )
            .def( "single_dependent_variable_history",
                  &tom::SingleObservationSet< STATE_SCALAR_TYPE,
                                              TIME_TYPE >::getSingleDependentVariableHistory,
                  R"doc(No documentation found.)doc" )
            .def_property_readonly( "dependent_variables_matrix",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::
                                            getObservationsDependentVariablesMatrix,
                                    R"doc(No documentation found.)doc" );

    m.def( "single_observation_set",
           &tss::singleObservationSetWithoutDependentVariables< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observable_type" ),
           py::arg( "link_definition" ),
           py::arg( "observations" ),
           py::arg( "observation_times" ),
           py::arg( "reference_link_end" ),
           py::arg( "ancilliary_settings" ) = nullptr,
           R"doc(No documentation found.)doc" );
}

}  // namespace estimation
}  // namespace numerical_simulation
}  // namespace tudatpy
