/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_estimation_analysis_outlier_rejection.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/outlierRejectionSettings.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;

using BaseSettings = tss::OutlierRejectionSettings< STATE_SCALAR_TYPE, TIME_TYPE >;
using CarpinoSettings = tss::CarpinoOutlierRejectionSettings< STATE_SCALAR_TYPE, TIME_TYPE >;
using SimpleSettings = tss::SimpleOutlierRejectionSettings< STATE_SCALAR_TYPE, TIME_TYPE >;

namespace tudatpy
{
namespace estimation
{
namespace estimation_analysis
{

void expose_estimation_analysis_outlier_rejection( py::module& m )
{
    py::class_< BaseSettings, std::shared_ptr< BaseSettings > >( m, "OutlierRejectionSettings", R"doc(

         Base class for defining the settings of an outlier rejection algorithm.

         Base class for defining the settings of an outlier rejection algorithm. This object is not instantiated directly but is created
         through one of its derived classes.

         All algorithms use a (recovery, rejection) pair of thresholds per
         :class:`~tudatpy.estimation.observations.ObservationSelectionCondition`. Each condition defines a group of
         observations that share one pair of thresholds. Observations that match no condition are not affected by the outlier
         rejection, and an error is raised at the start of the estimation if an observation matches more than one condition.

      )doc" )
            .def_property_readonly( "first_iteration_with_rejection",
                                    &BaseSettings::getFirstIterationWithRejection,
                                    R"doc(

         **read-only**

         Index of the first estimation iteration (counted from zero) in which observations may be rejected.

         :type: int
      )doc" )
            .def_property_readonly( "thresholds_per_condition",
                                    &BaseSettings::getThresholdsPerCondition,
                                    R"doc(

         **read-only**

         (recovery, rejection) thresholds per observation selection condition, in the order in which they were provided.
         If a single pair of thresholds was provided, this list holds one entry, with a condition that
         selects all observations (:meth:`~tudatpy.estimation.observations.ObservationSelectionCondition.all`).

         :type: list[tuple[ObservationSelectionCondition, tuple[float, float]]]
      )doc" );

    py::class_< CarpinoSettings, std::shared_ptr< CarpinoSettings >, BaseSettings >( m,
                                                                                     "CarpinoOutlierRejectionSettings",
                                                                                     R"doc(

         Class for defining the settings of the outlier rejection algorithm of :cite:t:`carpino2003`.

         Class for defining the settings of the outlier rejection algorithm of :cite:t:`carpino2003`. Objects of this
         class are typically created through the
         :func:`~tudatpy.estimation.estimation_analysis.carpino_outlier_rejection_settings` function.

      )doc" );

    m.def( "carpino_outlier_rejection_settings",
           py::overload_cast< const CarpinoSettings::ThresholdsPerCondition&, const int >(
                   &tss::carpinoOutlierRejectionSettings< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "chi2_thresholds_per_condition" ),
           py::arg( "first_iteration_with_rejection" ) = 1,
           R"doc(

 Function for creating settings for the outlier rejection algorithm of :cite:t:`carpino2003`, with thresholds per observation selection condition.

 As the other overload of this function, but with a (recovery, rejection) pair of :math:`\chi^2` thresholds for each
 :class:`~tudatpy.estimation.observations.ObservationSelectionCondition`. Since the :math:`\chi^2` value of an
 observation grows with its number of components (one for a range observation, two for a right ascension/declination
 pair), a single pair of thresholds is not suitable when observables of different sizes are combined. With this overload,
 for instance, the thresholds for range observations can be based on a :math:`\chi^2` distribution with one degree of
 freedom , and those for angular position observations on one with two degrees of freedom.

 Each condition defines a group of observations. The term with :math:`\chi^2_{\max}` in the adaptive rejection threshold
 and the limit of 20% new rejections per iteration are applied per group, using only the observations of that group.
 Since at least one observation may always be rejected per iteration, small groups can also have outliers rejected. The term
 with :math:`N`, which accounts for how well the solution is determined, uses the number of all observations in the fit,
 including those of other groups and those that match no condition.

 The conditions are evaluated once for each observation, when the estimation starts:

 * observations that match none of the conditions are not affected by the outlier rejection;
 * if an observation matches more than one condition, an error is raised;
 * conditions that depend on values that change during the estimation (such as the residuals, or whether an observation
   is rejected) use the values from before the estimation.


 Parameters
 ----------
 chi2_thresholds_per_condition : list[tuple[ObservationSelectionCondition, tuple[float, float]]]
     List of (condition, (recovery, rejection)) entries. Must not be empty. For each entry, both thresholds must not be
     negative, and the recovery threshold must not exceed the rejection threshold.
 first_iteration_with_rejection : int, default = 1
     Index of the first estimation iteration in which observations may be rejected, counted from zero.

 Returns
 -------
 :class:`~tudatpy.estimation.estimation_analysis.CarpinoOutlierRejectionSettings`
     Instance of the :class:`~tudatpy.estimation.estimation_analysis.OutlierRejectionSettings` derived
     :class:`~tudatpy.estimation.estimation_analysis.CarpinoOutlierRejectionSettings` class.


 Examples
 --------
 .. code-block:: python

     from tudatpy.estimation.observable_models_setup import model_settings
     from tudatpy.estimation.observations import observation_query

     # Range observations have one component, angular position observations two. A rejection threshold with the
     # same false-rejection probability as 5 sigma is therefore 25 for range, and about 29 for angular position.
     settings = estimation_analysis.carpino_outlier_rejection_settings(
         [
             (observation_query.observable_type == model_settings.one_way_range_type, (24.0, 25.0)),
             (observation_query.observable_type == model_settings.angular_position_type, (28.0, 29.0)),
         ]
     )

     )doc" );

    m.def( "carpino_outlier_rejection_settings",
           py::overload_cast< const std::pair< double, double >&, const int >(
                   &tss::carpinoOutlierRejectionSettings< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "chi2_thresholds" ),
           py::arg( "first_iteration_with_rejection" ) = 1,
           R"doc(

 Function for creating settings for the outlier rejection algorithm of :cite:t:`carpino2003`.

 Function for creating settings for the outlier rejection algorithm of :cite:t:`carpino2003`. This algorithm rejects and
 recovers outliers based on the :math:`\chi^2` value of the residual:

 .. math::
        \chi^2 = \xi_i \gamma^{-1}_{\xi_i} \xi_i^T

 where :math:`\xi_i` is the residual and :math:`\gamma_{\xi_i}` is the covariance matrix of the residuals. For a scalar
 observable, this is the square of the ratio of the residual value to its uncertainty. For a vector observable, it is the
 squared Mahalanobis distance of the residual. Note that :math:`\gamma_{\xi_i}` is **not** equal to the observation
 covariance/uncertainty. The residual covariance also takes into account whether an observation was included in the
 least-squares inversion or not.

 When the :math:`\chi^2` value of an observation exceeds the rejection threshold, the observation is marked as an outlier
 and excluded from the next iteration. Only 20% of the available observations (rounded down, but at least one) may be rejected per iteration
 for stability: when more observations exceed the threshold, only the worst ones are rejected. The rejection threshold used by the algorithm is increased when the largest
 :math:`\chi^2` value of the observations that are not rejected is large, or when few of these observations remain:
 it is :math:`\max(\chi^2_{\text{rej}} + 400\cdot 1.2^{-N}, \chi^2_{\max}/4)`, with :math:`\chi^2_{\text{rej}}` the
 rejection threshold from the settings, :math:`N` the number of observations that are not rejected, and
 :math:`\chi^2_{\max}` the largest :math:`\chi^2` value of these observations. If a
 rejected observation has a :math:`\chi^2` value below the recovery threshold, it is no longer an outlier, and is included
 in the next iteration of the least-squares fit again. A recovery threshold below the rejection threshold prevents
 observations from being alternately rejected and recovered in the later iterations. The two may also be equal.
 Observations that are rejected before the estimation remain rejected, and are not considered by the algorithm.

 This overload uses one pair of thresholds for all observations. This may not be suitable when using different observable types
 as :math:`\chi^2` grows with the number of components in the observable. In this case, provide thresholds per observation
 selection condition instead.

 Note that this algorithm requires that the observation covariance (or the inverse weight matrix) is compatible with the
 true errors. If all weights are equal to 1 (the default), a warning is printed, since the results are then only meaningful
 if the observation uncertainty is 1. If
 inaccurate weights are set, it may produce unexpected results. The observation covariance is the inverse of the
 complete weight matrix of the observation dataset, including the observations that are rejected, whereas the
 least-squares inversion uses the rows and columns of the weight matrix of the observations that are not rejected. The
 two are only consistent if the weights do not correlate rejected observations with observations that are not rejected.
 This algorithm cannot be used with weights that correlate observations in different observation sets (see
 :meth:`~tudatpy.estimation.observations.ObservationDataset.set_weight_block`).

 The resulting settings object is provided to the :class:`~tudatpy.estimation.estimation_analysis.EstimationInput`
 class, through its ``outlier_rejection_settings`` input.


 Parameters
 ----------
 chi2_thresholds : tuple[float, float]
     (recovery, rejection) :math:`\chi^2` thresholds, used for all observations. An observation that is used in the
     estimation is rejected when its :math:`\chi^2` value exceeds the rejection threshold, and an observation that was
     rejected in an earlier iteration is recovered when its :math:`\chi^2` value drops below the recovery threshold. Both
     must not be negative, and the recovery threshold must not exceed the rejection threshold.
 first_iteration_with_rejection : int, default = 1
     Index of the first estimation iteration in which observations may be rejected, counted from zero. The default
     value of 1 leaves the first iteration untouched, since the residuals of that iteration are computed with the
     a priori parameter values and can be large for all observations.

 Returns
 -------
 :class:`~tudatpy.estimation.estimation_analysis.CarpinoOutlierRejectionSettings`
     Instance of the :class:`~tudatpy.estimation.estimation_analysis.OutlierRejectionSettings` derived
     :class:`~tudatpy.estimation.estimation_analysis.CarpinoOutlierRejectionSettings` class.


 Examples
 --------
 .. code-block:: python

     # Reject observations with a chi-squared above 25 (roughly 5 sigma for a scalar observable), and recover them
     # once their chi-squared drops below 24
     settings = estimation_analysis.carpino_outlier_rejection_settings((24.0, 25.0))

     )doc" );

    py::class_< SimpleSettings, std::shared_ptr< SimpleSettings >, BaseSettings >( m,
                                                                                   "SimpleOutlierRejectionSettings",
                                                                                   R"doc(

         Class for defining the settings of the simple outlier rejection algorithm.

         Class for defining the settings of an outlier rejection algorithm that compares the residual of an observation
         against a (recovery, rejection) pair of thresholds. Objects of this class are typically created through the
         :func:`~tudatpy.estimation.estimation_analysis.simple_outlier_rejection_settings` function.

      )doc" );

    m.def( "simple_outlier_rejection_settings",
           py::overload_cast< const SimpleSettings::ThresholdsPerCondition&, const int >(
                   &tss::simpleOutlierRejectionSettings< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "residual_thresholds_per_condition" ),
           py::arg( "first_iteration_with_rejection" ) = 1,
           R"doc(

 Function for creating settings for the simple outlier rejection algorithm, with thresholds per observation selection condition.

 As the other overload of this function, but with a (recovery, rejection) pair of residual thresholds for each
 :class:`~tudatpy.estimation.observations.ObservationSelectionCondition`. This allows, for instance, different thresholds
 for observations of the same observable type that come from sources with a different accuracy (such as space-based and
 ground-based astrometry), by selecting observations on their link ends.

 The conditions are evaluated once for each observation, when the estimation starts:

 * observations that match none of the conditions are not affected by the outlier rejection;
 * if an observation matches more than one condition, an error is raised;
 * conditions that depend on values that change during the estimation (such as the residuals, or whether an observation
   is rejected) use the values from before the estimation.


 Parameters
 ----------
 residual_thresholds_per_condition : list[tuple[ObservationSelectionCondition, tuple[float, float]]]
     List of (condition, (recovery, rejection)) entries, in SI units. Must not be empty. For each entry, both thresholds
     must not be negative, and the recovery threshold must not exceed the rejection threshold.
 first_iteration_with_rejection : int, default = 1
     Index of the first estimation iteration in which observations may be rejected, counted from zero.

 Returns
 -------
 :class:`~tudatpy.estimation.estimation_analysis.SimpleOutlierRejectionSettings`
     Instance of the :class:`~tudatpy.estimation.estimation_analysis.OutlierRejectionSettings` derived
     :class:`~tudatpy.estimation.estimation_analysis.SimpleOutlierRejectionSettings` class.


 Examples
 --------
 .. code-block:: python

     from tudatpy.estimation.observable_models_setup import links
     from tudatpy.estimation.observations import observation_query

     # Astrometry from Gaia is far more accurate than ground-based astrometry, so it gets tighter thresholds
     # (in radians)
     from_gaia = observation_query.receiver == links.body_origin_link_end_id("Gaia")
     settings = estimation_analysis.simple_outlier_rejection_settings(
         [
             (from_gaia, (5.0e-9, 1.0e-8)),
             (~from_gaia, (5.0e-6, 1.0e-5)),
         ]
     )

     )doc" );

    m.def( "simple_outlier_rejection_settings",
           py::overload_cast< const std::pair< double, double >&, const int >(
                   &tss::simpleOutlierRejectionSettings< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "residual_thresholds" ),
           py::arg( "first_iteration_with_rejection" ) = 1,
           R"doc(

 Function for creating settings for the simple outlier rejection algorithm.

 Function for creating settings for the simple outlier rejection algorithm. This algorithm decides the rejection status of
 an observation based on its absolute residual (the largest element for vector observations), using a (recovery, rejection)
 pair of thresholds:

 * an observation in the fit is rejected when its residual exceeds the rejection threshold;
 * a rejected observation is recovered when its residual drops below the recovery threshold;
 * otherwise, the status is unchanged.

 A recovery threshold below the rejection threshold prevents observations from being alternately rejected and recovered. The
 two may be equal. A recovery threshold of 0 means rejected observations are never recovered. This overload uses one pair
 for all observations (in SI units); use the overload with thresholds per observation selection condition to use different
 thresholds for different groups of observations. Observations that are rejected before the estimation are not considered.

 The rejection threshold should be larger than the residuals of the observations that are not outliers at the iteration
 ``first_iteration_with_rejection``: otherwise many valid observations are rejected, and the estimation may select an
 iteration that used only part of the data as its best iteration (check
 :attr:`~tudatpy.estimation.estimation_analysis.EstimationOutput.active_flags_per_iteration`).

 The resulting settings object is provided to the :class:`~tudatpy.estimation.estimation_analysis.EstimationInput`
 class, through its ``outlier_rejection_settings`` input.


 Parameters
 ----------
 residual_thresholds : tuple[float, float]
     (recovery, rejection) residual thresholds, used for all observations. Both must not be negative, and the recovery
     threshold must not exceed the rejection threshold.
 first_iteration_with_rejection : int, default = 1
     Index of the first estimation iteration in which observations may be rejected, counted from zero. The default
     value of 1 leaves the first iteration untouched, since the residuals of that iteration are computed with the
     a priori parameter values and can be large for all observations.

 Returns
 -------
 :class:`~tudatpy.estimation.estimation_analysis.SimpleOutlierRejectionSettings`
     Instance of the :class:`~tudatpy.estimation.estimation_analysis.OutlierRejectionSettings` derived
     :class:`~tudatpy.estimation.estimation_analysis.SimpleOutlierRejectionSettings` class.


 Examples
 --------
 .. code-block:: python

     # Reject residuals above 10 m, recover them once below 5 m
     settings = estimation_analysis.simple_outlier_rejection_settings((5.0, 10.0))

     )doc" );
}

}  // namespace estimation_analysis
}  // namespace estimation
}  // namespace tudatpy
