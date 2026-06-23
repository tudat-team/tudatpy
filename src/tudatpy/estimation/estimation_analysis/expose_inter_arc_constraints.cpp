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
#include <memory>
#include <utility>
#include <vector>

#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "tudat/simulation/estimation_setup/interArcStateContinuityConstraintSettings.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;

namespace
{

//! Build a 6x6 diagonal weight matrix from either a scalar or a length-3 sequence of per-component weights.
Eigen::Matrix< double, 6, 6 > buildDiagonalWeight( const py::object& positionWeight, const py::object& velocityWeight )
{
    Eigen::Matrix< double, 6, 6 > C = Eigen::Matrix< double, 6, 6 >::Zero( );
    auto fillBlock = [ &C ]( const py::object& value, int startIdx ) {
        if( value.is_none( ) )
        {
            return;
        }
        if( py::isinstance< py::float_ >( value ) || py::isinstance< py::int_ >( value ) )
        {
            const double scalar = value.cast< double >( );
            C( startIdx, startIdx ) = scalar;
            C( startIdx + 1, startIdx + 1 ) = scalar;
            C( startIdx + 2, startIdx + 2 ) = scalar;
        }
        else
        {
            const std::vector< double > weights = value.cast< std::vector< double > >( );
            if( weights.size( ) != 3 )
            {
                throw std::runtime_error( "Inter-arc continuity weight must be a scalar or a length-3 sequence." );
            }
            for( int i = 0; i < 3; ++i )
            {
                C( startIdx + i, startIdx + i ) = weights[ static_cast< std::size_t >( i ) ];
            }
        }
    };
    fillBlock( positionWeight, 0 );
    fillBlock( velocityWeight, 3 );
    return C;
}

std::vector< double > buildMuValues( const py::object& mu )
{
    if( py::isinstance< py::float_ >( mu ) || py::isinstance< py::int_ >( mu ) )
    {
        return { mu.cast< double >( ) };
    }
    return mu.cast< std::vector< double > >( );
}

}  // namespace

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

         Soft inter-arc translational state continuity-prior settings for a single multi-arc body.

         Created via the factory functions :func:`full_state_continuity`, :func:`position_only_continuity`,
         :func:`velocity_only_continuity`, or :func:`general_continuity`. Pass a list of these objects to
         :meth:`CovarianceAnalysisInput.set_inter_arc_continuity_constraints` or the inherited
         :meth:`EstimationInput.set_inter_arc_continuity_constraints` to attach the feature. The feature is
         currently supported for pure multi-arc translational estimators only.

         The cost contribution per boundary is ``q = d^T W_d d`` with
         ``W_d = (1 / (mu * m_d_total)) * C`` and ``d = x_right(t_c) - x_left(t_c)``. The linearized
         normal-equation terms are ``D_norm.T @ W_d @ D_norm`` and ``-D_norm.T @ W_d @ d``, where ``D_norm``
         is the right-minus-left state-transition/sensitivity row block after applying the estimator's column
         normalization. ``m_d_total`` is the global rank sum across every settings entry. Larger ``mu`` weakens
         the prior.
      )doc" )
            .def_property_readonly( "body", &tss::InterArcStateContinuityConstraintSettings::body )
            .def_property_readonly( "connection_epochs", &tss::InterArcStateContinuityConstraintSettings::connectionEpochs )
            .def_property_readonly( "mu_values", &tss::InterArcStateContinuityConstraintSettings::muValues )
            .def_property_readonly( "arc_pairs", &tss::InterArcStateContinuityConstraintSettings::arcPairs );

    m.def(
            "full_state_continuity",
            []( std::string body,
                std::vector< double > epochs,
                py::object positionWeight,
                py::object velocityWeight,
                py::object mu,
                std::vector< std::pair< int, int > > arcPairs ) {
                auto C = buildDiagonalWeight( positionWeight, velocityWeight );
                return std::make_shared< tss::InterArcStateContinuityConstraintSettings >(
                        std::move( body ),
                        std::move( epochs ),
                        std::vector< Eigen::Matrix< double, 6, 6 > >{ C },
                        buildMuValues( mu ),
                        std::move( arcPairs ) );
            },
            py::arg( "body" ),
            py::arg( "epochs" ),
            py::arg( "position_weight" ) = py::float_( 1.0 ),
            py::arg( "velocity_weight" ) = py::float_( 1.0 ),
            py::arg( "mu" ) = py::float_( 1.0 ),
            py::arg( "arc_pairs" ) = std::vector< std::pair< int, int > >( ),
            R"doc(

         Build a full-state (position + velocity) soft-prior settings object. ``position_weight`` and
         ``velocity_weight`` may each be a scalar (isotropic) or a length-3 sequence (anisotropic). When
         ``arc_pairs`` is empty the prior is applied to every consecutive arc pair as ``(0, 1)``,
         ``(1, 2)``, ... in the order of ``epochs``.
      )doc" );

    m.def(
            "position_only_continuity",
            []( std::string body,
                std::vector< double > epochs,
                py::object positionWeight,
                py::object mu,
                std::vector< std::pair< int, int > > arcPairs ) {
                auto C = buildDiagonalWeight( positionWeight, py::float_( 0.0 ) );
                return std::make_shared< tss::InterArcStateContinuityConstraintSettings >(
                        std::move( body ),
                        std::move( epochs ),
                        std::vector< Eigen::Matrix< double, 6, 6 > >{ C },
                        buildMuValues( mu ),
                        std::move( arcPairs ) );
            },
            py::arg( "body" ),
            py::arg( "epochs" ),
            py::arg( "position_weight" ) = py::float_( 1.0 ),
            py::arg( "mu" ) = py::float_( 1.0 ),
            py::arg( "arc_pairs" ) = std::vector< std::pair< int, int > >( ),
            R"doc(

         Build a position-only soft-prior settings object (rank-3 ``C``). Velocity-row entries of ``C`` are
         zero so the prior leaves the inter-arc Delta-v free. This is the Rosetta OCM-boundary
         configuration.
      )doc" );

    m.def(
            "velocity_only_continuity",
            []( std::string body,
                std::vector< double > epochs,
                py::object velocityWeight,
                py::object mu,
                std::vector< std::pair< int, int > > arcPairs ) {
                auto C = buildDiagonalWeight( py::float_( 0.0 ), velocityWeight );
                return std::make_shared< tss::InterArcStateContinuityConstraintSettings >(
                        std::move( body ),
                        std::move( epochs ),
                        std::vector< Eigen::Matrix< double, 6, 6 > >{ C },
                        buildMuValues( mu ),
                        std::move( arcPairs ) );
            },
            py::arg( "body" ),
            py::arg( "epochs" ),
            py::arg( "velocity_weight" ) = py::float_( 1.0 ),
            py::arg( "mu" ) = py::float_( 1.0 ),
            py::arg( "arc_pairs" ) = std::vector< std::pair< int, int > >( ),
            R"doc(

         Build a velocity-only soft-prior settings object (rank-3 ``C``). Position-row entries of ``C`` are
         zero so the prior leaves inter-arc position jumps free.
      )doc" );

    m.def(
            "general_continuity",
            []( std::string body,
                std::vector< double > epochs,
                std::vector< Eigen::Matrix< double, 6, 6 > > weightMatrices,
                py::object mu,
                std::vector< std::pair< int, int > > arcPairs ) {
                return std::make_shared< tss::InterArcStateContinuityConstraintSettings >(
                        std::move( body ), std::move( epochs ), std::move( weightMatrices ), buildMuValues( mu ), std::move( arcPairs ) );
            },
            py::arg( "body" ),
            py::arg( "epochs" ),
            py::arg( "weight_matrices" ),
            py::arg( "mu" ) = py::float_( 1.0 ),
            py::arg( "arc_pairs" ) = std::vector< std::pair< int, int > >( ),
            R"doc(

         Build a generic soft-prior settings object from a list of 6x6 PSD weight matrices. The list must have
         either length 1 (applied to every pair) or length equal to ``len(epochs)``. Use this when you need a
         dense or per-boundary heterogeneous weight; otherwise prefer the preset factories.
      )doc" );
}

}  // namespace estimation_analysis
}  // namespace estimation
}  // namespace tudatpy
