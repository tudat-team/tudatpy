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

         Base class for defining the settings of an algorithm that rejects (and recovers) outlying observations during
         an estimation. Settings objects of this type are not created directly, but through the factory function of a
         specific algorithm, such as :func:`~tudatpy.estimation.estimation_analysis.carpino_outlier_rejection_settings`.
         The resulting object is provided to the :class:`~tudatpy.estimation.estimation_analysis.EstimationInput`
         class, after which the algorithm is applied once per iteration of the estimation.

         Observations that are rejected during the estimation are marked as rejected in the
         :class:`~tudatpy.estimation.observations.ObservationDataset`, and are excluded from the subsequent iterations.
         Which observations were rejected can therefore be inspected on the observation dataset after the estimation
         has finished.

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

 Function for creating settings for the outlier rejection algorithm of Carpino et al. (2003). The algorithm computes,
 for each observation, a chi-squared value from the residual of that observation and the covariance of that residual,
 and compares this value against two thresholds. An observation that is used in the estimation is rejected when its
 chi-squared value exceeds ``chi2_rejection_threshold``. An observation that was rejected in an earlier iteration is
 recovered when its chi-squared value drops below ``chi2_recovery_threshold``. Using a recovery threshold that is
 lower than the rejection threshold prevents observations from oscillating between the rejected and the accepted
 state in successive iterations.

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
     Index of the first estimation iteration in which observations may be rejected, counted from zero. The default
     value of 1 leaves the first iteration untouched, since the residuals of that iteration are computed with the
     a priori parameter values and can be large for all observations.

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

 Function for creating settings for an outlier rejection algorithm that compares the residual of each observation
 against a maximum allowed value. An observation is rejected when the size of its residual exceeds
 ``maximum_allowed_residual_value``. If ``allow_restore`` is true, an observation that was rejected in an earlier
 iteration is recovered as soon as the size of its residual drops back below that value. Observations that were
 excluded from the estimation by the user are never modified by this algorithm.

 The residual is compared against the value component by component, so an observation of an observable type with
 more than one component is rejected if any of its components exceeds the maximum allowed value.

 The resulting settings object is provided to the :class:`~tudatpy.estimation.estimation_analysis.EstimationInput`
 class, through its ``outlier_rejection_settings`` input.


 Parameters
 ----------
 maximum_allowed_residual_value : float
     Maximum allowed size of a residual, used for every observable type. Must be positive.
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

    m.def( "simple_outlier_rejection_settings",
           static_cast< std::shared_ptr< tss::OutlierRejectionSettings > ( * )(
                   const std::map< tss::ObservableType, double >&, const int, const bool ) >( &tss::simpleOutlierRejectionSettings ),
           py::arg( "maximum_allowed_residual_value_per_observable_type" ),
           py::arg( "first_iteration_with_rejection" ) = 1,
           py::arg( "allow_restore" ) = true,
           R"doc(

 Function for creating settings for the simple outlier rejection algorithm, with a separate maximum allowed residual
 value for each observable type.

 Function for creating settings for an outlier rejection algorithm that compares the residual of each observation
 against a maximum allowed value, where this value is defined separately for each observable type. An observation is
 rejected when the size of its residual exceeds the value that is provided for its observable type. If
 ``allow_restore`` is true, an observation that was rejected in an earlier iteration is recovered as soon as the size
 of its residual drops back below that value. Observations that were excluded from the estimation by the user are
 never modified by this algorithm.

 The residual is compared against the value component by component, so an observation of an observable type with
 more than one component is rejected if any of its components exceeds the maximum allowed value.

 An exception is thrown during the estimation if the observations contain an observable type that is not present in
 ``maximum_allowed_residual_value_per_observable_type``.

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
