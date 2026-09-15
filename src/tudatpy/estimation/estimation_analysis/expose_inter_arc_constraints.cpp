/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_inter_arc_constraints.h"

#include <Eigen/Core>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "tudat/simulation/estimation_setup/interArcStateContinuityConstraintSettings.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;

namespace tudatpy
{
namespace estimation
{
namespace estimation_analysis
{

void expose_inter_arc_constraints( py::module& m )
{
    py::class_< tss::InterArcStateContinuityConstraintSettings, std::shared_ptr< tss::InterArcStateContinuityConstraintSettings > >(
            m,
            "InterArcStateContinuityConstraintSettings",
            R"doc(

         Soft inter-arc translational state continuity-prior settings for one or more multi-arc bodies.

         Created via the factory functions :func:`full_state_continuity`, :func:`position_only_continuity`,
         :func:`velocity_only_continuity`, or :func:`general_continuity`. Pass a list of these objects to
         :meth:`CovarianceAnalysisInput.set_inter_arc_continuity_constraints` or the inherited
         :meth:`EstimationInput.set_inter_arc_continuity_constraints` to attach the feature. The feature is
         currently supported for pure multi-arc translational estimators only.

         For body :math:`b` and constrained arc pair :math:`k=(\ell,r)` at connection epoch :math:`t_{c,bk}`,
         the state discrepancy and its parameter partials are

         .. math::

            \mathbf{d}_{bk} &= \mathbf{x}_{r,b}(t_{c,bk})-\mathbf{x}_{\ell,b}(t_{c,bk}),\\
            \mathbf{D}_{bk} &= \frac{\partial\mathbf{d}_{bk}}{\partial\mathbf{p}}
            =\mathbf{M}_{r,b}(t_{c,bk})-\mathbf{M}_{\ell,b}(t_{c,bk}),

         where :math:`\mathbf{M}` contains the applicable rows of the full state-transition and sensitivity
         matrix. Let :math:`\mathbf{C}_{bk}` be the user-supplied positive semi-definite weight matrix and
         :math:`\overline{\mathbf{C}}_{bk}=\tfrac{1}{2}(\mathbf{C}_{bk}+\mathbf{C}_{bk}^{T})` its symmetric
         part. This distinction only affects matrices whose asymmetry is within the numerical validation tolerance.
         With scaling factor :math:`\mu_s` for settings object :math:`s`, :math:`m` scalar observations, and

         .. math::

            m_d=\sum_{b,k}\operatorname{rank}(\overline{\mathbf{C}}_{bk}),\qquad
            \mathbf{W}_{d,bk}=\frac{m}{\mu_s m_d}\overline{\mathbf{C}}_{bk},

         the implemented continuity cost is

         .. math::

            J_d=\frac{1}{2}\sum_{b,k}\mathbf{d}_{bk}^{T}\mathbf{W}_{d,bk}\mathbf{d}_{bk}.

         This is Eqs. (4) and (28) of
         `Lari et al. (2021) <https://doi.org/10.1007/s10686-021-09823-8>`_, after multiplying the averaged total
         objective by :math:`m/2` to match Tudat's :math:`\frac{1}{2}\mathbf{r}^{T}\mathbf{W}\mathbf{r}`
         convention. The rank sum generalizes Lari et al.'s :math:`m_d=6(n-1)` to component masks and general
         rank-deficient weights. If :math:`\mathbf{N}` contains the observation-design-matrix column
         normalization factors, :math:`\widetilde{\mathbf{D}}_{bk}=\mathbf{D}_{bk}\operatorname{diag}(\mathbf{N})^{-1}`.
         The normal matrix and estimation right-hand side receive

         .. math::

            \Delta\widetilde{\mathbf{P}}^{-1} &=
            \sum_{b,k}\widetilde{\mathbf{D}}_{bk}^{T}\mathbf{W}_{d,bk}\widetilde{\mathbf{D}}_{bk},\\
            \Delta\widetilde{\mathbf{g}} &=-
            \sum_{b,k}\widetilde{\mathbf{D}}_{bk}^{T}\mathbf{W}_{d,bk}\mathbf{d}_{bk}.

         ``bodies`` selects :math:`b`; ``connection_epochs`` supplies :math:`t_{c,bk}`; ``arc_pairs`` supplies
         :math:`(\ell,r)` (or consecutive pairs are inferred); the factory weight arguments construct
         :math:`\mathbf{C}_{bk}`; and ``constraint_scaling_factor`` is :math:`\mu_s`. Larger values of
         :math:`\mu_s` weaken the constraint.

         Attributes
         ----------
         bodies : list[str]
             Bodies whose multi-arc translational states are constrained.
         connection_epochs : dict[str, list[float]]
             Body-specific connection epochs :math:`t_{c,bk}`, in seconds since J2000.
         constraint_scaling_factor : float
             Positive penalty scaling factor :math:`\mu_s`.
         arc_pairs : dict[str, list[tuple[int, int]]]
             Optional body-specific zero-based ``(left_arc, right_arc)`` pairs.
      )doc" )
            .def_property_readonly( "bodies", &tss::InterArcStateContinuityConstraintSettings::bodies )
            .def_property_readonly( "connection_epochs", &tss::InterArcStateContinuityConstraintSettings::connectionEpochsByBody )
            .def_property_readonly( "constraint_scaling_factor", &tss::InterArcStateContinuityConstraintSettings::constraintScalingFactor )
            .def_property_readonly( "arc_pairs", &tss::InterArcStateContinuityConstraintSettings::arcPairsByBody );

    m.def( "full_state_continuity",
           static_cast< std::shared_ptr< tss::InterArcStateContinuityConstraintSettings > ( * )(
                   std::vector< std::string >,
                   std::map< std::string, std::vector< double > >,
                   std::map< std::string, std::variant< double, Eigen::VectorXd > >,
                   std::map< std::string, std::variant< double, Eigen::VectorXd > >,
                   double,
                   std::map< std::string, std::vector< std::pair< int, int > > > ) >( &tss::fullStateContinuity ),
           py::arg( "bodies" ),
           py::arg( "epochs" ),
           py::arg( "position_weights" ),
           py::arg( "velocity_weights" ),
           py::arg( "constraint_scaling_factor" ) = 1.0,
           py::arg( "arc_pairs" ) = std::map< std::string, std::vector< std::pair< int, int > > >( ),
           R"doc(

         Build full translational-state continuity settings.

         The position and velocity entries form the diagonal blocks of :math:`\mathbf{C}_{bk}` defined by
         :class:`InterArcStateContinuityConstraintSettings`. One matrix is broadcast to every connection epoch
         of a body.

         Parameters
         ----------
         bodies : list[str]
             Bodies whose position and velocity jumps are constrained.
         epochs : dict[str, list[float]]
             Connection epochs :math:`t_{c,bk}` for each body, in seconds since J2000.
         position_weights : dict[str, float | numpy.ndarray]
             Isotropic scalar or three Cartesian diagonal position weights for each body.
         velocity_weights : dict[str, float | numpy.ndarray]
             Isotropic scalar or three Cartesian diagonal velocity weights for each body.
         constraint_scaling_factor : float, default = 1.0
             Positive factor :math:`\mu_s`; larger values weaken the constraint.
         arc_pairs : dict[str, list[tuple[int, int]]], optional
             Zero-based consecutive arc pairs. If omitted, epoch entry :math:`k` connects arcs :math:`k` and :math:`k+1`.

         Returns
         -------
         InterArcStateContinuityConstraintSettings
             Settings containing the resulting full-state weight matrices.
      )doc" );

    m.def( "position_only_continuity",
           static_cast< std::shared_ptr< tss::InterArcStateContinuityConstraintSettings > ( * )(
                   std::vector< std::string >,
                   std::map< std::string, std::vector< double > >,
                   std::map< std::string, std::variant< double, Eigen::VectorXd > >,
                   double,
                   std::map< std::string, std::vector< std::pair< int, int > > > ) >( &tss::positionOnlyContinuity ),
           py::arg( "bodies" ),
           py::arg( "epochs" ),
           py::arg( "position_weights" ),
           py::arg( "constraint_scaling_factor" ) = 1.0,
           py::arg( "arc_pairs" ) = std::map< std::string, std::vector< std::pair< int, int > > >( ),
           R"doc(

         Build position-only continuity settings.

         This constructs :math:`\mathbf{C}_{bk}` with the velocity block set to zero, leaving inter-arc velocity
         jumps unconstrained. See :class:`InterArcStateContinuityConstraintSettings` for the mathematical model.

         Parameters
         ----------
         bodies : list[str]
             Bodies whose position jumps are constrained.
         epochs : dict[str, list[float]]
             Connection epochs :math:`t_{c,bk}` for each body, in seconds since J2000.
         position_weights : dict[str, float | numpy.ndarray]
             Isotropic scalar or three Cartesian diagonal position weights for each body.
         constraint_scaling_factor : float, default = 1.0
             Positive factor :math:`\mu_s`; larger values weaken the constraint.
         arc_pairs : dict[str, list[tuple[int, int]]], optional
             Zero-based consecutive arc pairs; inferred from epoch order when omitted.

         Returns
         -------
         InterArcStateContinuityConstraintSettings
             Settings containing position-only weight matrices.
      )doc" );

    m.def( "velocity_only_continuity",
           static_cast< std::shared_ptr< tss::InterArcStateContinuityConstraintSettings > ( * )(
                   std::vector< std::string >,
                   std::map< std::string, std::vector< double > >,
                   std::map< std::string, std::variant< double, Eigen::VectorXd > >,
                   double,
                   std::map< std::string, std::vector< std::pair< int, int > > > ) >( &tss::velocityOnlyContinuity ),
           py::arg( "bodies" ),
           py::arg( "epochs" ),
           py::arg( "velocity_weights" ),
           py::arg( "constraint_scaling_factor" ) = 1.0,
           py::arg( "arc_pairs" ) = std::map< std::string, std::vector< std::pair< int, int > > >( ),
           R"doc(

         Build velocity-only continuity settings.

         This constructs :math:`\mathbf{C}_{bk}` with the position block set to zero, leaving inter-arc position
         jumps unconstrained. See :class:`InterArcStateContinuityConstraintSettings` for the mathematical model.

         Parameters
         ----------
         bodies : list[str]
             Bodies whose velocity jumps are constrained.
         epochs : dict[str, list[float]]
             Connection epochs :math:`t_{c,bk}` for each body, in seconds since J2000.
         velocity_weights : dict[str, float | numpy.ndarray]
             Isotropic scalar or three Cartesian diagonal velocity weights for each body.
         constraint_scaling_factor : float, default = 1.0
             Positive factor :math:`\mu_s`; larger values weaken the constraint.
         arc_pairs : dict[str, list[tuple[int, int]]], optional
             Zero-based consecutive arc pairs; inferred from epoch order when omitted.

         Returns
         -------
         InterArcStateContinuityConstraintSettings
             Settings containing velocity-only weight matrices.
      )doc" );

    m.def( "general_continuity",
           static_cast< std::shared_ptr< tss::InterArcStateContinuityConstraintSettings > ( * )(
                   std::vector< std::string >,
                   std::map< std::string, std::vector< double > >,
                   std::map< std::string, std::variant< Eigen::MatrixXd, std::vector< Eigen::MatrixXd > > >,
                   double,
                   std::map< std::string, std::vector< std::pair< int, int > > > ) >( &tss::generalContinuity ),
           py::arg( "bodies" ),
           py::arg( "epochs" ),
           py::arg( "weight_matrices" ),
           py::arg( "constraint_scaling_factor" ) = 1.0,
           py::arg( "arc_pairs" ) = std::map< std::string, std::vector< std::pair< int, int > > >( ),
           R"doc(

         Build continuity settings from general Cartesian-state weight matrices.

         See :class:`InterArcStateContinuityConstraintSettings` for the definition of
         :math:`\mathbf{C}_{bk}` and the complete mathematical model.

         Parameters
         ----------
         bodies : list[str]
             Bodies whose Cartesian-state jumps are constrained.
         epochs : dict[str, list[float]]
             Connection epochs :math:`t_{c,bk}` for each body, in seconds since J2000.
         weight_matrices : dict[str, numpy.ndarray[float, 6, 6] | list[numpy.ndarray[float, 6, 6]]]
             Symmetric positive semi-definite matrices :math:`\mathbf{C}_{bk}`. One matrix is broadcast to all
             epochs of its body; otherwise supply one matrix per epoch.
         constraint_scaling_factor : float, default = 1.0
             Positive factor :math:`\mu_s`; larger values weaken the constraint.
         arc_pairs : dict[str, list[tuple[int, int]]], optional
             Zero-based consecutive arc pairs; inferred from epoch order when omitted.

         Returns
         -------
         InterArcStateContinuityConstraintSettings
             Settings containing the supplied general weight matrices.
      )doc" );
}

}  // namespace estimation_analysis
}  // namespace estimation
}  // namespace tudatpy
