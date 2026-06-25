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
    Eigen::Matrix< double, 6, 6 > constraintWeightMatrix = Eigen::Matrix< double, 6, 6 >::Zero( );
    auto fillBlock = [ &constraintWeightMatrix ]( const py::object& value, int startIdx ) {
        if( value.is_none( ) )
        {
            return;
        }
        if( py::isinstance< py::float_ >( value ) || py::isinstance< py::int_ >( value ) )
        {
            const double scalar = value.cast< double >( );
            constraintWeightMatrix( startIdx, startIdx ) = scalar;
            constraintWeightMatrix( startIdx + 1, startIdx + 1 ) = scalar;
            constraintWeightMatrix( startIdx + 2, startIdx + 2 ) = scalar;
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
                constraintWeightMatrix( startIdx + i, startIdx + i ) = weights[ static_cast< std::size_t >( i ) ];
            }
        }
    };
    fillBlock( positionWeight, 0 );
    fillBlock( velocityWeight, 3 );
    return constraintWeightMatrix;
}

std::vector< double > buildConstraintScalingFactors( const py::object& constraintScalingFactor )
{
    if( py::isinstance< py::float_ >( constraintScalingFactor ) || py::isinstance< py::int_ >( constraintScalingFactor ) )
    {
        return { constraintScalingFactor.cast< double >( ) };
    }
    return constraintScalingFactor.cast< std::vector< double > >( );
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

         See Lari et al. (2021), Eq. (28) for the underlying mathematics. The cost contribution per boundary
         uses the state discrepancy, the constraint weight matrix, the constraint scaling factor, and the total
         constrained dimension across every settings entry. Larger constraint scaling factors weaken the penalty.
      )doc" )
            .def_property_readonly( "body", &tss::InterArcStateContinuityConstraintSettings::body )
            .def_property_readonly( "connection_epochs", &tss::InterArcStateContinuityConstraintSettings::connectionEpochs )
            .def_property_readonly( "constraint_scaling_factors",
                                    &tss::InterArcStateContinuityConstraintSettings::constraintScalingFactors )
            .def_property_readonly( "arc_pairs", &tss::InterArcStateContinuityConstraintSettings::arcPairs );

    m.def(
            "full_state_continuity",
            []( std::string body,
                std::vector< double > epochs,
                py::object positionWeight,
                py::object velocityWeight,
                py::object constraintScalingFactor,
                std::vector< std::pair< int, int > > arcPairs ) {
                auto constraintWeightMatrix = buildDiagonalWeight( positionWeight, velocityWeight );
                return std::make_shared< tss::InterArcStateContinuityConstraintSettings >(
                        std::move( body ),
                        std::move( epochs ),
                        std::vector< Eigen::Matrix< double, 6, 6 > >{ constraintWeightMatrix },
                        buildConstraintScalingFactors( constraintScalingFactor ),
                        std::move( arcPairs ) );
            },
            py::arg( "body" ),
            py::arg( "epochs" ),
            py::arg( "position_weight" ) = py::float_( 1.0 ),
            py::arg( "velocity_weight" ) = py::float_( 1.0 ),
            py::arg( "constraint_scaling_factor" ) = py::float_( 1.0 ),
            py::arg( "arc_pairs" ) = std::vector< std::pair< int, int > >( ),
            R"doc(

         Build a full-state (position + velocity) continuity settings object. ``position_weight`` and
         ``velocity_weight`` may each be a scalar (isotropic) or a length-3 sequence (anisotropic).
         ``constraint_scaling_factor`` may be a scalar or a sequence matching the number of connection epochs.
         When ``arc_pairs`` is empty the constraint is applied to every consecutive arc pair as ``(0, 1)``,
         ``(1, 2)``, ... in the order of ``epochs``.
      )doc" );

    m.def(
            "position_only_continuity",
            []( std::string body,
                std::vector< double > epochs,
                py::object positionWeight,
                py::object constraintScalingFactor,
                std::vector< std::pair< int, int > > arcPairs ) {
                auto constraintWeightMatrix = buildDiagonalWeight( positionWeight, py::float_( 0.0 ) );
                return std::make_shared< tss::InterArcStateContinuityConstraintSettings >(
                        std::move( body ),
                        std::move( epochs ),
                        std::vector< Eigen::Matrix< double, 6, 6 > >{ constraintWeightMatrix },
                        buildConstraintScalingFactors( constraintScalingFactor ),
                        std::move( arcPairs ) );
            },
            py::arg( "body" ),
            py::arg( "epochs" ),
            py::arg( "position_weight" ) = py::float_( 1.0 ),
            py::arg( "constraint_scaling_factor" ) = py::float_( 1.0 ),
            py::arg( "arc_pairs" ) = std::vector< std::pair< int, int > >( ),
            R"doc(

         Build a position-only continuity settings object. Velocity-row entries of the constraint weight matrix
         are zero so the constraint leaves the inter-arc Delta-v free. ``constraint_scaling_factor`` may be a
         scalar or a sequence matching the number of connection epochs.
      )doc" );

    m.def(
            "velocity_only_continuity",
            []( std::string body,
                std::vector< double > epochs,
                py::object velocityWeight,
                py::object constraintScalingFactor,
                std::vector< std::pair< int, int > > arcPairs ) {
                auto constraintWeightMatrix = buildDiagonalWeight( py::float_( 0.0 ), velocityWeight );
                return std::make_shared< tss::InterArcStateContinuityConstraintSettings >(
                        std::move( body ),
                        std::move( epochs ),
                        std::vector< Eigen::Matrix< double, 6, 6 > >{ constraintWeightMatrix },
                        buildConstraintScalingFactors( constraintScalingFactor ),
                        std::move( arcPairs ) );
            },
            py::arg( "body" ),
            py::arg( "epochs" ),
            py::arg( "velocity_weight" ) = py::float_( 1.0 ),
            py::arg( "constraint_scaling_factor" ) = py::float_( 1.0 ),
            py::arg( "arc_pairs" ) = std::vector< std::pair< int, int > >( ),
            R"doc(

         Build a velocity-only continuity settings object. Position-row entries of the constraint weight matrix
         are zero so the constraint leaves inter-arc position jumps free. ``constraint_scaling_factor`` may be a
         scalar or a sequence matching the number of connection epochs.
      )doc" );

    m.def(
            "general_continuity",
            []( std::string body,
                std::vector< double > epochs,
                std::vector< Eigen::Matrix< double, 6, 6 > > weightMatrices,
                py::object constraintScalingFactor,
                std::vector< std::pair< int, int > > arcPairs ) {
                return std::make_shared< tss::InterArcStateContinuityConstraintSettings >(
                        std::move( body ),
                        std::move( epochs ),
                        std::move( weightMatrices ),
                        buildConstraintScalingFactors( constraintScalingFactor ),
                        std::move( arcPairs ) );
            },
            py::arg( "body" ),
            py::arg( "epochs" ),
            py::arg( "weight_matrices" ),
            py::arg( "constraint_scaling_factor" ) = py::float_( 1.0 ),
            py::arg( "arc_pairs" ) = std::vector< std::pair< int, int > >( ),
            R"doc(

         Build a generic continuity settings object from a list of 6x6 PSD weight matrices. The list must have
         either length 1 (applied to every pair) or length equal to ``len(epochs)``. ``constraint_scaling_factor``
         may use the same scalar-or-per-epoch broadcasting. Use this when you need a dense or per-boundary
         heterogeneous weight; otherwise prefer the preset factories.
      )doc" );
}

}  // namespace estimation_analysis
}  // namespace estimation
}  // namespace tudatpy
