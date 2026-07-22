/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_INTERARCSTATECONTINUITYCONSTRAINTSETTINGS_H
#define TUDAT_INTERARCSTATECONTINUITYCONSTRAINTSETTINGS_H

#include <memory>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include <Eigen/Core>

namespace tudat
{
namespace simulation_setup
{

//! User-facing soft-continuity constraint between consecutive multi-arc translational state arcs of one or more bodies.
//! The cost added to Tudat's LSQ target is, per pair (Lari et al. 2021, Eqs. 4 and 28):
//!   pairCost = 0.5 * stateDiscrepancy^T * scaledConstraintWeight * stateDiscrepancy
//! where scaledConstraintWeight is the constraint weight matrix multiplied by the number of scalar observations
//! and divided by the product of the constraint scaling factor and total constrained dimension. Larger constraint
//! scaling factors weaken the penalty.
class InterArcStateContinuityConstraintSettings
{
public:
    using EpochMap = std::map< std::string, std::vector< double > >;
    using WeightMatrixMap = std::map< std::string, std::vector< Eigen::MatrixXd > >;
    using CartesianStateWeight = std::variant< double, Eigen::VectorXd >;
    using CartesianStateWeightMap = std::map< std::string, CartesianStateWeight >;
    using WeightMatrixInput = std::variant< Eigen::MatrixXd, std::vector< Eigen::MatrixXd > >;
    using WeightMatrixInputMap = std::map< std::string, WeightMatrixInput >;
    using ArcPairMap = std::map< std::string, std::vector< std::pair< int, int > > >;

    InterArcStateContinuityConstraintSettings( std::vector< std::string > bodies,
                                               EpochMap connectionEpochsByBody,
                                               WeightMatrixMap weightMatricesByBody,
                                               double constraintScalingFactor = 1.0,
                                               ArcPairMap arcPairsByBody = {} );

    const std::vector< std::string >& bodies( ) const
    {
        return bodies_;
    }
    const EpochMap& connectionEpochsByBody( ) const
    {
        return connectionEpochsByBody_;
    }
    const WeightMatrixMap& weightMatricesByBody( ) const
    {
        return weightMatricesByBody_;
    }
    double constraintScalingFactor( ) const
    {
        return constraintScalingFactor_;
    }
    const ArcPairMap& arcPairsByBody( ) const
    {
        return arcPairsByBody_;
    }

    const std::vector< double >& connectionEpochsForBody( const std::string& body ) const;

    const std::vector< std::pair< int, int > >& arcPairsForBody( const std::string& body ) const;

    //! Resolve the weight matrix for the i-th pair of a body (handles 1-or-n broadcasting).
    const Eigen::MatrixXd& weightMatrixForBodyAndPair( const std::string& body, std::size_t pairIndex ) const;

    //! Number of regularized boundaries for a body.
    std::size_t numberOfPairsForBody( const std::string& body ) const
    {
        return connectionEpochsForBody( body ).size( );
    }

    //! Total number of regularized boundaries across all configured bodies.
    std::size_t totalNumberOfPairs( ) const;

private:
    void validate( ) const;
    void validateWeightMatrix( const std::string& body, const Eigen::MatrixXd& constraintWeightMatrix, std::size_t entryIndex ) const;

    std::vector< std::string > bodies_;
    EpochMap connectionEpochsByBody_;
    WeightMatrixMap weightMatricesByBody_;
    double constraintScalingFactor_;
    ArcPairMap arcPairsByBody_;
};

//! Create the 6x6 diagonal weight matrix for Cartesian-state continuity.
Eigen::MatrixXd createCartesianStateWeightMatrix( const Eigen::Vector3d& positionWeights, const Eigen::Vector3d& velocityWeights );

//! Build a soft-prior settings object with both position and velocity continuity.
std::shared_ptr< InterArcStateContinuityConstraintSettings > fullStateContinuity(
        std::vector< std::string > bodies,
        InterArcStateContinuityConstraintSettings::EpochMap connectionEpochsByBody,
        double positionWeight = 1.0,
        double velocityWeight = 1.0,
        double constraintScalingFactor = 1.0,
        InterArcStateContinuityConstraintSettings::ArcPairMap arcPairsByBody = {} );

//! Build full-state settings with body-specific isotropic or component-wise position and velocity weights.
std::shared_ptr< InterArcStateContinuityConstraintSettings > fullStateContinuity(
        std::vector< std::string > bodies,
        InterArcStateContinuityConstraintSettings::EpochMap connectionEpochsByBody,
        InterArcStateContinuityConstraintSettings::CartesianStateWeightMap positionWeightsByBody,
        InterArcStateContinuityConstraintSettings::CartesianStateWeightMap velocityWeightsByBody,
        double constraintScalingFactor = 1.0,
        InterArcStateContinuityConstraintSettings::ArcPairMap arcPairsByBody = {} );

//! Build a settings object with position-only continuity (velocity rows/columns of the weight matrix zeroed).
std::shared_ptr< InterArcStateContinuityConstraintSettings > positionOnlyContinuity(
        std::vector< std::string > bodies,
        InterArcStateContinuityConstraintSettings::EpochMap connectionEpochsByBody,
        double positionWeight = 1.0,
        double constraintScalingFactor = 1.0,
        InterArcStateContinuityConstraintSettings::ArcPairMap arcPairsByBody = {} );

//! Build position-only settings with body-specific isotropic or component-wise position weights.
std::shared_ptr< InterArcStateContinuityConstraintSettings > positionOnlyContinuity(
        std::vector< std::string > bodies,
        InterArcStateContinuityConstraintSettings::EpochMap connectionEpochsByBody,
        InterArcStateContinuityConstraintSettings::CartesianStateWeightMap positionWeightsByBody,
        double constraintScalingFactor = 1.0,
        InterArcStateContinuityConstraintSettings::ArcPairMap arcPairsByBody = {} );

//! Build a settings object with velocity-only continuity (position rows/columns of the weight matrix zeroed).
std::shared_ptr< InterArcStateContinuityConstraintSettings > velocityOnlyContinuity(
        std::vector< std::string > bodies,
        InterArcStateContinuityConstraintSettings::EpochMap connectionEpochsByBody,
        double velocityWeight = 1.0,
        double constraintScalingFactor = 1.0,
        InterArcStateContinuityConstraintSettings::ArcPairMap arcPairsByBody = {} );

//! Build velocity-only settings with body-specific isotropic or component-wise velocity weights.
std::shared_ptr< InterArcStateContinuityConstraintSettings > velocityOnlyContinuity(
        std::vector< std::string > bodies,
        InterArcStateContinuityConstraintSettings::EpochMap connectionEpochsByBody,
        InterArcStateContinuityConstraintSettings::CartesianStateWeightMap velocityWeightsByBody,
        double constraintScalingFactor = 1.0,
        InterArcStateContinuityConstraintSettings::ArcPairMap arcPairsByBody = {} );

//! Build a soft-prior settings object with arbitrary body-specific, possibly per-boundary, 6x6 PSD weight matrices.
std::shared_ptr< InterArcStateContinuityConstraintSettings > generalContinuity(
        std::vector< std::string > bodies,
        InterArcStateContinuityConstraintSettings::EpochMap connectionEpochsByBody,
        InterArcStateContinuityConstraintSettings::WeightMatrixMap weightMatricesByBody,
        double constraintScalingFactor = 1.0,
        InterArcStateContinuityConstraintSettings::ArcPairMap arcPairsByBody = {} );

//! Build general settings from one or more body-specific 6x6 PSD weight matrices.
std::shared_ptr< InterArcStateContinuityConstraintSettings > generalContinuity(
        std::vector< std::string > bodies,
        InterArcStateContinuityConstraintSettings::EpochMap connectionEpochsByBody,
        InterArcStateContinuityConstraintSettings::WeightMatrixInputMap weightMatricesByBody,
        double constraintScalingFactor = 1.0,
        InterArcStateContinuityConstraintSettings::ArcPairMap arcPairsByBody = {} );

}  // namespace simulation_setup
}  // namespace tudat

#endif  // TUDAT_INTERARCSTATECONTINUITYCONSTRAINTSETTINGS_H
