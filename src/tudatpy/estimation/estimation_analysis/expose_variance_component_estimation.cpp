/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_variance_component_estimation.h"

#include <memory>
#include <vector>

#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/varianceComponentEstimation.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tom = tudat::observation_models;

namespace tudatpy
{
namespace estimation
{
namespace estimation_analysis
{

void expose_variance_component_estimation( py::module& m )
{
    py::class_< tss::VarianceComponentSettings, std::shared_ptr< tss::VarianceComponentSettings > >( m,
                                                                                                     "VarianceComponentSettings",
                                                                                                     R"doc(

         Settings for estimating the weight of each group of observations with Helmert's method.

         Created with :func:`variance_component_settings`, and used by
         :meth:`~tudatpy.estimation.estimation_analysis.Estimator.perform_estimation_with_variance_components`.

      )doc" )
            .def_readonly( "groups", &tss::VarianceComponentSettings::groups_ )
            .def_readwrite( "tolerance", &tss::VarianceComponentSettings::tolerance_ )
            .def_readwrite( "maximum_iterations", &tss::VarianceComponentSettings::maximumNumberOfIterations_ )
            .def_readwrite( "minimum_redundancy", &tss::VarianceComponentSettings::minimumRedundancy_ )
            .def_readwrite( "warm_start", &tss::VarianceComponentSettings::warmStart_ )
            .def_readwrite( "print_output", &tss::VarianceComponentSettings::printOutput_ );

    m.def( "variance_component_settings",
           &tss::varianceComponentSettings,
           py::arg( "groups" ),
           py::arg( "tolerance" ) = 0.01,
           py::arg( "maximum_iterations" ) = 10,
           py::arg( "minimum_redundancy" ) = 1.0,
           py::arg( "warm_start" ) = true,
           py::arg( "print_output" ) = true,
           R"doc(

         Function to create settings for estimating the weight of each group of observations.

         Each group of observations gets its own variance factor, estimated with Helmert's method as described by
         `Sahin et al. (1992) <https://doi.org/10.1007/BF02033189>`_. After each estimation, the weights of a group are
         divided by its variance factor, and the estimation is repeated until all factors are close to 1. Relative
         weights within a group are kept. Observations that are in no group keep their weights.

         To give each arc of a multi-arc estimation its own weights, define one group per observable type and arc (for
         instance by combining an observable type parser with a time bounds parser). Each observation set must then lie
         within one arc.

         Parameters
         ----------
         groups : list[ObservationCollectionParser]
             One parser per group. Groups must not overlap.
         tolerance : float, default = 0.01
             The iteration stops when all variance factors differ from 1 by less than this value.
         maximum_iterations : int, default = 10
             Maximum number of estimations (each one a full least-squares estimation).
         minimum_redundancy : float, default = 1.0
             Groups with a smaller redundancy (number of observations minus the share taken up by the estimated
             parameters) keep their weights.
         warm_start : bool, default = True
             Whether each estimation starts from the result of the previous one. If False, or if the a priori constraint is
             applied to the total parameter deviation, every estimation starts from the initial parameter values.
         print_output : bool, default = True
             Whether to print the variance factors of each iteration.

         Returns
         -------
         :class:`VarianceComponentSettings`
             Settings for the variance component estimation.

      )doc" );

    py::class_< tss::VarianceComponentGroupResult >( m, "VarianceComponentGroupResult", R"doc(

         Variance component results for one group of observations.

      )doc" )
            .def_readonly( "number_of_observations",
                           &tss::VarianceComponentGroupResult::numberOfObservations,
                           "Number of observations in the group with a non-zero weight." )
            .def_readonly( "redundancy",
                           &tss::VarianceComponentGroupResult::redundancy,
                           "Number of observations minus the share taken up by the estimated parameters." )
            .def_readonly( "weighted_square_sum",
                           &tss::VarianceComponentGroupResult::weightedSquareSum,
                           "Weighted sum of squared residuals of the group." )
            .def_readonly( "variance_factor",
                           &tss::VarianceComponentGroupResult::varianceFactor,
                           "Estimated variance factor. The group weights should be divided by this value." )
            .def_readonly( "simplified_variance_factor",
                           &tss::VarianceComponentGroupResult::simplifiedVarianceFactor,
                           "Variance factor from the simplified estimator (weighted square sum divided by redundancy)." )
            .def_readonly( "rms_residual", &tss::VarianceComponentGroupResult::rmsResidual, "Root mean square of the group residuals." )
            .def_readonly( "is_frozen",
                           &tss::VarianceComponentGroupResult::isFrozen,
                           "Whether the group was kept fixed because its redundancy is too small." );

    py::class_< tss::VarianceComponentEstimate >( m, "VarianceComponentEstimate", R"doc(

         Variance component results for all groups, from a single least-squares solution.

      )doc" )
            .def_readonly(
                    "group_results", &tss::VarianceComponentEstimate::groupResults, "Results per group, in the order of the groups." )
            .def_readonly(
                    "helmert_matrix", &tss::VarianceComponentEstimate::helmertMatrix, "Helmert matrix of the groups that were not frozen." )
            .def_readonly( "simplified_estimator_used",
                           &tss::VarianceComponentEstimate::simplifiedEstimatorUsed,
                           "Whether the simplified estimator was used because Helmert's equations gave a non-positive factor." )
            .def_property_readonly( "variance_factors", &tss::VarianceComponentEstimate::getVarianceFactors, "Variance factor per group." )
            .def_property_readonly( "redundancies", &tss::VarianceComponentEstimate::getRedundancies, "Redundancy per group." );

    py::class_< tss::VarianceComponentEstimationOutput< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tss::VarianceComponentEstimationOutput< STATE_SCALAR_TYPE, TIME_TYPE > > >(
            m, "VarianceComponentEstimationOutput", R"doc(

         Output of an estimation in which the weights of observation groups are estimated.

      )doc" )
            .def_readonly( "estimation_output",
                           &tss::VarianceComponentEstimationOutput< STATE_SCALAR_TYPE, TIME_TYPE >::estimationOutput,
                           "Output of the last least-squares estimation." )
            .def_readonly( "iteration_results",
                           &tss::VarianceComponentEstimationOutput< STATE_SCALAR_TYPE, TIME_TYPE >::iterationResults,
                           "Variance component results of each iteration." )
            .def_readonly( "total_variance_factor_history",
                           &tss::VarianceComponentEstimationOutput< STATE_SCALAR_TYPE, TIME_TYPE >::totalVarianceFactorHistory,
                           "Total weight scaling per group used in each iteration's estimation." )
            .def_readonly( "least_squares_best_iterations",
                           &tss::VarianceComponentEstimationOutput< STATE_SCALAR_TYPE, TIME_TYPE >::leastSquaresBestIterations,
                           "Best iteration of the least-squares estimation, for each iteration." )
            .def_readonly( "is_converged",
                           &tss::VarianceComponentEstimationOutput< STATE_SCALAR_TYPE, TIME_TYPE >::isConverged,
                           "Whether all variance factors converged to 1 within the tolerance." )
            .def_property_readonly( "total_variance_factors",
                                    &tss::VarianceComponentEstimationOutput< STATE_SCALAR_TYPE, TIME_TYPE >::getTotalVarianceFactors,
                                    "Total weight scaling per group. Each group sigma was multiplied by the square root of this value." );

    m.def(
            "compute_variance_components",
            []( const std::shared_ptr< tss::EstimationOutput< STATE_SCALAR_TYPE, TIME_TYPE > > estimationOutput,
                const std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > > observationCollection,
                const std::vector< std::shared_ptr< tom::ObservationCollectionParser > >& groups,
                const double minimumRedundancy ) {
                return tss::computeHelmertVarianceComponents< STATE_SCALAR_TYPE, TIME_TYPE >(
                        estimationOutput, observationCollection, groups, minimumRedundancy );
            },
            py::arg( "estimation_output" ),
            py::arg( "observation_collection" ),
            py::arg( "groups" ),
            py::arg( "minimum_redundancy" ) = 1.0,
            R"doc(

         Function to estimate the variance factor of each group of observations, from one estimation.

         Uses Helmert's method (`Sahin et al., 1992 <https://doi.org/10.1007/BF02033189>`_). This function does not change any
         weights, so it can be used to build a custom iteration.

         Parameters
         ----------
         estimation_output : :class:`EstimationOutput`
             Output of the estimation. The design matrix must have been saved.
         observation_collection : :class:`~tudatpy.estimation.observations.ObservationCollection`
             Observations used in the estimation.
         groups : list[ObservationCollectionParser]
             One parser per group. Groups must not overlap.
         minimum_redundancy : float, default = 1.0
             Groups with a smaller redundancy keep a variance factor of 1.

         Returns
         -------
         :class:`VarianceComponentEstimate`
             Variance factor and related values per group.

      )doc" );

    m.def(
            "compute_variance_components_from_matrices",
            []( const Eigen::MatrixXd& designMatrix,
                const Eigen::VectorXd& weights,
                const Eigen::MatrixXd& normalMatrix,
                const Eigen::VectorXd& residuals,
                const std::vector< std::vector< int > >& groupRowIndices,
                const double minimumRedundancy ) {
                return tss::computeHelmertVarianceComponents(
                        designMatrix, weights, normalMatrix, residuals, groupRowIndices, minimumRedundancy );
            },
            py::arg( "design_matrix" ),
            py::arg( "weights" ),
            py::arg( "normal_matrix" ),
            py::arg( "residuals" ),
            py::arg( "group_row_indices" ),
            py::arg( "minimum_redundancy" ) = 1.0,
            R"doc(

         Function to estimate the variance factor of each group of observations, from the matrices of a least-squares solution.

         Uses Helmert's method (`Sahin et al., 1992 <https://doi.org/10.1007/BF02033189>`_). Parts of the normal matrix that do
         not come from any group (such as a priori information) are treated as known, with a variance factor of 1.

         Parameters
         ----------
         design_matrix : numpy.ndarray
             Design matrix, one row per observation.
         weights : numpy.ndarray
             Weight of each observation.
         normal_matrix : numpy.ndarray
             Full normal matrix (inverse covariance) of the solution, including any a priori information.
         residuals : numpy.ndarray
             Post-fit residuals.
         group_row_indices : list[list[int]]
             Row indices of each group. Groups must not overlap.
         minimum_redundancy : float, default = 1.0
             Groups with a smaller redundancy keep a variance factor of 1.

         Returns
         -------
         :class:`VarianceComponentEstimate`
             Variance factor and related values per group.

      )doc" );
}

}  // namespace estimation_analysis
}  // namespace estimation
}  // namespace tudatpy
