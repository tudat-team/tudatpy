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

#include "tudat/simulation/estimation_setup/outlierRejectionSettings.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;

namespace tudatpy
{
namespace estimation
{
namespace estimation_analysis
{

void expose_estimation_analysis_outlier_rejection( py::module& m )
{
    py::class_< tss::OutlierRejectionSettings, std::shared_ptr< tss::OutlierRejectionSettings > >( m,
                                                                                                   "OutlierRejectionSettings",
                                                                                                   R"doc(

         Base class for defining the settings of an outlier rejection algorithm.

         Base class for defining the settings of an outlier rejection algorithm. This object is not instantiated directly but is created
         through one of its derived classes.

      )doc" );

    py::class_< tss::CarpinoOutlierRejectionSettings,
                std::shared_ptr< tss::CarpinoOutlierRejectionSettings >,
                tss::OutlierRejectionSettings >( m,
                                                 "CarpinoOutlierRejectionSettings",
                                                 R"doc(

         Class for defining the settings of the outlier rejection algorithm of Carpino et al. (2003).

         Class for defining the settings of the outlier rejection algorithm of Carpino et al. (2003). Objects of this
         class are typically created through the
         :func:`~tudatpy.estimation.estimation_analysis.carpino_outlier_rejection_settings` function.

      )doc" )
            .def_property_readonly( "chi2_rejection_threshold",
                                    &tss::CarpinoOutlierRejectionSettings::getChi2RejectionThreshold,
                                    R"doc(

         **read-only**

         Chi-squared value above which an observation that is used in the estimation is rejected.

         :type: float
      )doc" )
            .def_property_readonly( "chi2_recovery_threshold",
                                    &tss::CarpinoOutlierRejectionSettings::getChi2RecoveryThreshold,
                                    R"doc(

         **read-only**

         Chi-squared value below which an observation that was rejected in an earlier iteration is recovered.

         :type: float
      )doc" )
            .def_property_readonly( "maximum_rejected_fraction",
                                    &tss::CarpinoOutlierRejectionSettings::getMaximumRejectedFraction,
                                    R"doc(

         **read-only**

         Maximum fraction of all observations that may be in the rejected state.

         :type: float
      )doc" )
            .def_property_readonly( "first_iteration_with_rejection",
                                    &tss::CarpinoOutlierRejectionSettings::getFirstIterationWithRejection,
                                    R"doc(

         **read-only**

         Index of the first estimation iteration in which observations may be rejected.

         :type: int
      )doc" );

    m.def( "carpino_outlier_rejection_settings",
           &tss::carpinoOutlierRejectionSettings,
           py::arg( "chi2_rejection_threshold" ) = 9.0,
           py::arg( "chi2_recovery_threshold" ) = 8.0,
           py::arg( "maximum_rejected_fraction" ) = 0.25,
           py::arg( "first_iteration_with_rejection" ) = 1,
           R"doc(

 Function for creating settings for the outlier rejection algorithm of Carpino et al. (2003).

 Function for creating settings for the outlier rejection algorithm of Carpino et al. (2003). This algorithm rejects and recovers
outliers based on the :math:`\chi^2` value of the residual:

 .. math::
        \chi^2 = \xi_i \gamma^{-1}_{\xi_i} \xi_i^T

where :math:`\xi_i` is the residual and, :math:`\gamma_{\xi_i}` is the covariance matrix of the residuals. For a scalar observable, this
simply becomes the ratio of the residual value to its uncertainty. For a vector observable, it is the Mahalanobis distance of the residual.
Note that :math:`\gamma_{\xi_i}` is **not** equal to the observation covariance/uncertainty. The residual covariance also takes into account
whether an observation was included in the least-squares inversion or not.

When the :math:`\chi^2` value of the observation exceeds the rejection threshold, the observation is marked as an outlier and excluded
from the next iteration. If a rejected observation has a :math:`\chi^2` below the recovery threshold, it is no longer an outlier, and is
included in the next iteration of the least-squares fit again. Seperate recovery and rejection threshold are adopted to avoid observations
from jumping between rejected and recovered in the later iterations.

Note that this algorithm requires that the observation covariance (or the inverse weight matrix) is compatible with the true errors. If no
weights are set on the observations, the algorithm will throw an error. If inaccurate weights are set, it may produce unexpected results.

 The resulting settings object is provided to the :class:`~tudatpy.estimation.estimation_analysis.EstimationInput`
 class, through its ``outlier_rejection_settings`` input.


 Parameters
 ----------
 chi2_rejection_threshold : float, default = 9.0
     Chi-squared value above which an observation that is used in the estimation is rejected.
 chi2_recovery_threshold : float, default = 8.0
     Chi-squared value below which an observation that was rejected in an earlier iteration is recovered. Must be
     smaller than ``chi2_rejection_threshold``.
 maximum_rejected_fraction : float, default = 0.25
     Maximum fraction (between 0 and 1) of all observations that may be in the rejected state.
 first_iteration_with_rejection : int, default = 1
     Index of the first estimation iteration in which observations may be rejected, counted from zero.

 Returns
 -------
 :class:`~tudatpy.estimation.estimation_analysis.CarpinoOutlierRejectionSettings`
     Instance of the :class:`~tudatpy.estimation.estimation_analysis.OutlierRejectionSettings` derived
     :class:`~tudatpy.estimation.estimation_analysis.CarpinoOutlierRejectionSettings` class.

     )doc" );

    py::class_< tss::SimpleOutlierRejectionSettings,
                std::shared_ptr< tss::SimpleOutlierRejectionSettings >,
                tss::OutlierRejectionSettings >( m,
                                                 "SimpleOutlierRejectionSettings",
                                                 R"doc(

         Class for defining the settings of the simple outlier rejection algorithm.

         Class for defining the settings of an outlier rejection algorithm that compares the residual of an observation
         against a maximum allowed value. Objects of this class are typically created through the
         :func:`~tudatpy.estimation.estimation_analysis.simple_outlier_rejection_settings` function.

      )doc" )
            .def_property_readonly( "maximum_allowed_residual_value",
                                    &tss::SimpleOutlierRejectionSettings::getMaximumAllowedResidualValue,
                                    R"doc(

         **read-only**

         Maximum allowed size of a residual, used for every observable type. This value is -1.0 if a separate value was
         provided for each observable type, in which case ``maximum_allowed_residual_value_per_observable_type`` is
         used instead.

         :type: float
      )doc" )
            .def_property_readonly( "maximum_allowed_residual_value_per_observable_type",
                                    &tss::SimpleOutlierRejectionSettings::getMaximumAllowedResidualValueMap,
                                    R"doc(

         **read-only**

         Maximum allowed size of a residual for each observable type. This dictionary is empty if a single value was
         provided for all observable types, in which case ``maximum_allowed_residual_value`` is used instead.

         :type: dict[ObservableType, float]
      )doc" )
            .def_property_readonly( "first_iteration_with_rejection",
                                    &tss::SimpleOutlierRejectionSettings::getFirstIterationWithRejection,
                                    R"doc(

         **read-only**

         Index of the first estimation iteration in which observations may be rejected.

         :type: int
      )doc" )
            .def_property_readonly( "allow_restore",
                                    &tss::SimpleOutlierRejectionSettings::getAllowRestore,
                                    R"doc(

         **read-only**

         Whether an observation that was rejected in an earlier iteration may be recovered.

         :type: bool
      )doc" );

    m.def( "simple_outlier_rejection_settings",
           static_cast< std::shared_ptr< tss::OutlierRejectionSettings > ( * )( const double, const int, const bool ) >(
                   &tss::simpleOutlierRejectionSettings ),
           py::arg( "maximum_allowed_residual_value" ),
           py::arg( "first_iteration_with_rejection" ) = 1,
           py::arg( "allow_restore" ) = true,
           R"doc(

 Function for creating settings for the simple outlier rejection algorithm.

Function for creating settings for the simple outlier rejection algorithm. This algorithm decides the rejection status of an observation
based on the absolute residual value of the observation (or the largest element in the residual vector in case of vector observations).
If it exceeds the threshold, it will be excluded from the next iteration. If ``allow_restore`` is true, it can be recovered again in a later
iteration. One threshold is used for all observable types in the estimation (in SI units).

 The resulting settings object is provided to the :class:`~tudatpy.estimation.estimation_analysis.EstimationInput`
 class, through its ``outlier_rejection_settings`` input.


 Parameters
 ----------
 maximum_allowed_residual_value : float
     Maximum allowed size of a residual, used for every observable type. Must be positive.
 first_iteration_with_rejection : int, default = 1
     Index of the first estimation iteration in which observations may be rejected, counted from zero.
 allow_restore : bool, default = True
     Whether an observation that was rejected in an earlier iteration may be recovered in a later iteration.

 Returns
 -------
 :class:`~tudatpy.estimation.estimation_analysis.SimpleOutlierRejectionSettings`
     Instance of the :class:`~tudatpy.estimation.estimation_analysis.OutlierRejectionSettings` derived
     :class:`~tudatpy.estimation.estimation_analysis.SimpleOutlierRejectionSettings` class.

     )doc" );

    m.def( "simple_outlier_rejection_settings",
           static_cast< std::shared_ptr< tss::OutlierRejectionSettings > ( * )(
                   const std::map< tss::ObservableType, double >&, const int, const bool ) >( &tss::simpleOutlierRejectionSettings ),
           py::arg( "maximum_allowed_residual_value_per_observable_type" ),
           py::arg( "first_iteration_with_rejection" ) = 1,
           py::arg( "allow_restore" ) = true,
           R"doc(

 Function for creating settings for the simple outlier rejection algorithm.

Function for creating settings for the simple outlier rejection algorithm. This algorithm decides the rejection status of an observation
based on the absolute residual value of the observation (or the largest element in the residual vector in case of vector observations).
If it exceeds the threshold, it will be excluded from the next iteration. If ``allow_restore`` is true, it can be recovered again in a later
iteration. A dictionary containing the residual threshold for each observable type used in the estimation must be provided. In case an
observable type is used in the estimation, but missing from the dict, an error will be thrown.

 The resulting settings object is provided to the :class:`~tudatpy.estimation.estimation_analysis.EstimationInput`
 class, through its ``outlier_rejection_settings`` input.


 Parameters
 ----------
 maximum_allowed_residual_value_per_observable_type : dict[ObservableType, float]
     Maximum allowed size of a residual for each observable type. Must not be empty, and each value must be
     positive.
 first_iteration_with_rejection : int, default = 1
     Index of the first estimation iteration in which observations may be rejected, counted from zero. The default
     value of 1 leaves the first iteration untouched, since the residuals of that iteration are computed with the
     a priori parameter values and can be large for all observations.
 allow_restore : bool, default = True
     Whether an observation that was rejected in an earlier iteration may be recovered in a later iteration.

 Returns
 -------
 :class:`~tudatpy.estimation.estimation_analysis.SimpleOutlierRejectionSettings`
     Instance of the :class:`~tudatpy.estimation.estimation_analysis.OutlierRejectionSettings` derived
     :class:`~tudatpy.estimation.estimation_analysis.SimpleOutlierRejectionSettings` class.

     )doc" );
}

}  // namespace estimation_analysis
}  // namespace estimation
}  // namespace tudatpy
