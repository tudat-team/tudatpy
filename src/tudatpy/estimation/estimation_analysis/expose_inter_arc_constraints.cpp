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
#include <set>
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
Eigen::MatrixXd buildDiagonalWeight( const py::object& positionWeight, const py::object& velocityWeight )
{
    Eigen::MatrixXd constraintWeightMatrix = Eigen::MatrixXd::Zero( 6, 6 );
    auto fillBlock = [ &constraintWeightMatrix ]( const py::object& value, int startIdx ) {
        if( value.is_none( ) )
        {
            return;
        }
        if( PyNumber_Check( value.ptr( ) ) != 0 && PySequence_Check( value.ptr( ) ) == 0 )
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

void validateDictionaryKeys( const py::dict& values, const std::vector< std::string >& bodies, const std::string& argumentName )
{
    std::set< std::string > expectedBodies( bodies.begin( ), bodies.end( ) );
    for( const auto& item : values )
    {
        const std::string body = py::cast< std::string >( item.first );
        if( expectedBodies.count( body ) == 0 )
        {
            throw std::runtime_error( "Inter-arc continuity " + argumentName + " contains unknown body \"" + body + "\"." );
        }
    }
    for( const auto& body : bodies )
    {
        if( !values.contains( py::str( body ) ) )
        {
            throw std::runtime_error( "Inter-arc continuity " + argumentName + " is missing body \"" + body + "\"." );
        }
    }
}

py::object getBodyValue( const py::dict& values, const std::string& body )
{
    return py::reinterpret_borrow< py::object >( values[ py::str( body ) ] );
}

py::dict buildUniformWeightDictionary( const std::vector< std::string >& bodies, const double value )
{
    py::dict weights;
    for( const auto& body : bodies )
    {
        weights[ py::str( body ) ] = py::float_( value );
    }
    return weights;
}

tss::InterArcStateContinuityConstraintSettings::WeightMatrixMap buildDiagonalWeightMatricesByBody( const std::vector< std::string >& bodies,
                                                                                                   const py::dict& positionWeights,
                                                                                                   const py::dict& velocityWeights )
{
    validateDictionaryKeys( positionWeights, bodies, "position_weights" );
    validateDictionaryKeys( velocityWeights, bodies, "velocity_weights" );
    tss::InterArcStateContinuityConstraintSettings::WeightMatrixMap weightMatricesByBody;
    for( const auto& body : bodies )
    {
        weightMatricesByBody[ body ] = { buildDiagonalWeight( getBodyValue( positionWeights, body ),
                                                              getBodyValue( velocityWeights, body ) ) };
    }
    return weightMatricesByBody;
}

tss::InterArcStateContinuityConstraintSettings::WeightMatrixMap buildGeneralWeightMatricesByBody( const std::vector< std::string >& bodies,
                                                                                                  const py::dict& weightMatrices )
{
    validateDictionaryKeys( weightMatrices, bodies, "weight_matrices" );
    tss::InterArcStateContinuityConstraintSettings::WeightMatrixMap weightMatricesByBody;
    for( const auto& body : bodies )
    {
        const py::object bodyWeightMatrices = getBodyValue( weightMatrices, body );
        try
        {
            weightMatricesByBody[ body ] = { bodyWeightMatrices.cast< Eigen::MatrixXd >( ) };
        }
        catch( const py::cast_error& )
        {
            weightMatricesByBody[ body ] = bodyWeightMatrices.cast< std::vector< Eigen::MatrixXd > >( );
        }
    }
    return weightMatricesByBody;
}

tss::InterArcStateContinuityConstraintSettings::ArcPairMap buildArcPairsByBody( const py::object& arcPairs )
{
    if( arcPairs.is_none( ) )
    {
        return {};
    }
    return arcPairs.cast< tss::InterArcStateContinuityConstraintSettings::ArcPairMap >( );
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

         Soft inter-arc translational state continuity-prior settings for one or more multi-arc bodies.

         Created via the factory functions :func:`full_state_continuity`, :func:`position_only_continuity`,
         :func:`velocity_only_continuity`, or :func:`general_continuity`. Pass a list of these objects to
         :meth:`CovarianceAnalysisInput.set_inter_arc_continuity_constraints` or the inherited
         :meth:`EstimationInput.set_inter_arc_continuity_constraints` to attach the feature. The feature is
         currently supported for pure multi-arc translational estimators only.

         See Lari et al. (2021), Eq. (28) for the underlying mathematics. The cost contribution per boundary
         uses the state discrepancy, the constraint weight matrix, the constraint scaling factor, and the total
         constrained dimension across every settings entry. Larger constraint scaling factors weaken the penalty.
      )doc" )
            .def_property_readonly( "bodies", &tss::InterArcStateContinuityConstraintSettings::bodies )
            .def_property_readonly( "connection_epochs", &tss::InterArcStateContinuityConstraintSettings::connectionEpochsByBody )
            .def_property_readonly( "constraint_scaling_factor", &tss::InterArcStateContinuityConstraintSettings::constraintScalingFactor )
            .def_property_readonly( "arc_pairs", &tss::InterArcStateContinuityConstraintSettings::arcPairsByBody );

    m.def(
            "full_state_continuity",
            []( std::vector< std::string > bodies,
                tss::InterArcStateContinuityConstraintSettings::EpochMap epochs,
                py::dict positionWeights,
                py::dict velocityWeights,
                double constraintScalingFactor,
                py::object arcPairs ) {
                auto weightMatricesByBody = buildDiagonalWeightMatricesByBody( bodies, positionWeights, velocityWeights );
                auto arcPairsByBody = buildArcPairsByBody( arcPairs );
                return std::make_shared< tss::InterArcStateContinuityConstraintSettings >( std::move( bodies ),
                                                                                           std::move( epochs ),
                                                                                           std::move( weightMatricesByBody ),
                                                                                           constraintScalingFactor,
                                                                                           std::move( arcPairsByBody ) );
            },
            py::arg( "bodies" ),
            py::arg( "epochs" ),
            py::arg( "position_weights" ),
            py::arg( "velocity_weights" ),
            py::arg( "constraint_scaling_factor" ) = 1.0,
            py::arg( "arc_pairs" ) = py::none( ),
            R"doc(

         Build a full translational-state (position + velocity) continuity settings object. ``epochs``, ``position_weights``,
         and ``velocity_weights`` must be dictionaries keyed by body name. Weight values may each be a scalar
         (isotropic) or a length-3 sequence (anisotropic). ``constraint_scaling_factor`` is one global scaling
         factor for all bodies and connection epochs.
      )doc" );

    m.def(
            "position_only_continuity",
            []( std::vector< std::string > bodies,
                tss::InterArcStateContinuityConstraintSettings::EpochMap epochs,
                py::dict positionWeights,
                double constraintScalingFactor,
                py::object arcPairs ) {
                py::dict zeroVelocityWeights = buildUniformWeightDictionary( bodies, 0.0 );
                auto weightMatricesByBody = buildDiagonalWeightMatricesByBody( bodies, positionWeights, zeroVelocityWeights );
                auto arcPairsByBody = buildArcPairsByBody( arcPairs );
                return std::make_shared< tss::InterArcStateContinuityConstraintSettings >( std::move( bodies ),
                                                                                           std::move( epochs ),
                                                                                           std::move( weightMatricesByBody ),
                                                                                           constraintScalingFactor,
                                                                                           std::move( arcPairsByBody ) );
            },
            py::arg( "bodies" ),
            py::arg( "epochs" ),
            py::arg( "position_weights" ),
            py::arg( "constraint_scaling_factor" ) = 1.0,
            py::arg( "arc_pairs" ) = py::none( ),
            R"doc(

         Build a position-only continuity settings object. Velocity-row entries of the constraint weight matrix
         are zero so the constraint leaves the inter-arc Delta-v free. ``epochs`` and ``position_weights`` must
         be dictionaries keyed by body name. ``constraint_scaling_factor`` is one global scaling factor.
      )doc" );

    m.def(
            "velocity_only_continuity",
            []( std::vector< std::string > bodies,
                tss::InterArcStateContinuityConstraintSettings::EpochMap epochs,
                py::dict velocityWeights,
                double constraintScalingFactor,
                py::object arcPairs ) {
                py::dict zeroPositionWeights = buildUniformWeightDictionary( bodies, 0.0 );
                auto weightMatricesByBody = buildDiagonalWeightMatricesByBody( bodies, zeroPositionWeights, velocityWeights );
                auto arcPairsByBody = buildArcPairsByBody( arcPairs );
                return std::make_shared< tss::InterArcStateContinuityConstraintSettings >( std::move( bodies ),
                                                                                           std::move( epochs ),
                                                                                           std::move( weightMatricesByBody ),
                                                                                           constraintScalingFactor,
                                                                                           std::move( arcPairsByBody ) );
            },
            py::arg( "bodies" ),
            py::arg( "epochs" ),
            py::arg( "velocity_weights" ),
            py::arg( "constraint_scaling_factor" ) = 1.0,
            py::arg( "arc_pairs" ) = py::none( ),
            R"doc(

         Build a velocity-only continuity settings object. Position-row entries of the constraint weight matrix
         are zero so the constraint leaves inter-arc position jumps free. ``epochs`` and ``velocity_weights`` must
         be dictionaries keyed by body name. ``constraint_scaling_factor`` is one global scaling factor.
      )doc" );

    m.def(
            "general_continuity",
            []( std::vector< std::string > bodies,
                tss::InterArcStateContinuityConstraintSettings::EpochMap epochs,
                py::dict weightMatrices,
                double constraintScalingFactor,
                py::object arcPairs ) {
                auto weightMatricesByBody = buildGeneralWeightMatricesByBody( bodies, weightMatrices );
                auto arcPairsByBody = buildArcPairsByBody( arcPairs );
                return std::make_shared< tss::InterArcStateContinuityConstraintSettings >( std::move( bodies ),
                                                                                           std::move( epochs ),
                                                                                           std::move( weightMatricesByBody ),
                                                                                           constraintScalingFactor,
                                                                                           std::move( arcPairsByBody ) );
            },
            py::arg( "bodies" ),
            py::arg( "epochs" ),
            py::arg( "weight_matrices" ),
            py::arg( "constraint_scaling_factor" ) = 1.0,
            py::arg( "arc_pairs" ) = py::none( ),
            R"doc(

         Build a generic continuity settings object from body-specific 6x6 PSD weight matrices. ``epochs`` and
         ``weight_matrices`` must be dictionaries keyed by body name. Each body may provide one matrix, broadcast
         to every connection epoch for that body, or one matrix per body-specific connection epoch.
      )doc" );
}

}  // namespace estimation_analysis
}  // namespace estimation
}  // namespace tudatpy
