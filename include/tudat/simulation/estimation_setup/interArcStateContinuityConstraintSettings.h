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
#include <vector>

#include <Eigen/Core>

namespace tudat
{
namespace simulation_setup
{

//! User-facing soft-continuity constraint between consecutive multi-arc translational state arcs of one or more bodies.
//! The cost added to the LSQ target is, per pair (Lari et al. 2021 Eq. 28):
//!   pairCost = stateDiscrepancy^T * scaledConstraintWeight * stateDiscrepancy
//! where scaledConstraintWeight is the constraint weight matrix divided by the product of the constraint scaling
//! factor and the total constrained dimension. Larger constraint scaling factors weaken the penalty.
class InterArcStateContinuityConstraintSettings
{
public:
    using EpochMap = std::map< std::string, std::vector< double > >;
    using WeightMatrixMap = std::map< std::string, std::vector< Eigen::MatrixXd > >;
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

namespace detail
{

inline Eigen::MatrixXd diagonalWeight( double position, double velocity )
{
    Eigen::MatrixXd constraintWeightMatrix = Eigen::MatrixXd::Zero( 6, 6 );
    constraintWeightMatrix( 0, 0 ) = position;
    constraintWeightMatrix( 1, 1 ) = position;
    constraintWeightMatrix( 2, 2 ) = position;
    constraintWeightMatrix( 3, 3 ) = velocity;
    constraintWeightMatrix( 4, 4 ) = velocity;
    constraintWeightMatrix( 5, 5 ) = velocity;
    return constraintWeightMatrix;
}

inline InterArcStateContinuityConstraintSettings::WeightMatrixMap createUniformBodyWeightMap( const std::vector< std::string >& bodies,
                                                                                              const Eigen::MatrixXd& weightMatrix )
{
    InterArcStateContinuityConstraintSettings::WeightMatrixMap weightMatricesByBody;
    for( const auto& body : bodies )
    {
        weightMatricesByBody[ body ] = { weightMatrix };
    }
    return weightMatricesByBody;
}

}  // namespace detail

//! Build a soft-prior settings object with both position and velocity continuity.
inline std::shared_ptr< InterArcStateContinuityConstraintSettings > fullStateContinuity(
        std::vector< std::string > bodies,
        InterArcStateContinuityConstraintSettings::EpochMap connectionEpochsByBody,
        double positionWeight = 1.0,
        double velocityWeight = 1.0,
        double constraintScalingFactor = 1.0,
        InterArcStateContinuityConstraintSettings::ArcPairMap arcPairsByBody = {} )
{
    auto weightMatricesByBody = detail::createUniformBodyWeightMap( bodies, detail::diagonalWeight( positionWeight, velocityWeight ) );
    return std::make_shared< InterArcStateContinuityConstraintSettings >( std::move( bodies ),
                                                                          std::move( connectionEpochsByBody ),
                                                                          std::move( weightMatricesByBody ),
                                                                          constraintScalingFactor,
                                                                          std::move( arcPairsByBody ) );
}

//! Build a settings object with position-only continuity (velocity rows/columns of the weight matrix zeroed).
inline std::shared_ptr< InterArcStateContinuityConstraintSettings > positionOnlyContinuity(
        std::vector< std::string > bodies,
        InterArcStateContinuityConstraintSettings::EpochMap connectionEpochsByBody,
        double positionWeight = 1.0,
        double constraintScalingFactor = 1.0,
        InterArcStateContinuityConstraintSettings::ArcPairMap arcPairsByBody = {} )
{
    auto weightMatricesByBody = detail::createUniformBodyWeightMap( bodies, detail::diagonalWeight( positionWeight, 0.0 ) );
    return std::make_shared< InterArcStateContinuityConstraintSettings >( std::move( bodies ),
                                                                          std::move( connectionEpochsByBody ),
                                                                          std::move( weightMatricesByBody ),
                                                                          constraintScalingFactor,
                                                                          std::move( arcPairsByBody ) );
}

//! Build a settings object with velocity-only continuity (position rows/columns of the weight matrix zeroed).
inline std::shared_ptr< InterArcStateContinuityConstraintSettings > velocityOnlyContinuity(
        std::vector< std::string > bodies,
        InterArcStateContinuityConstraintSettings::EpochMap connectionEpochsByBody,
        double velocityWeight = 1.0,
        double constraintScalingFactor = 1.0,
        InterArcStateContinuityConstraintSettings::ArcPairMap arcPairsByBody = {} )
{
    auto weightMatricesByBody = detail::createUniformBodyWeightMap( bodies, detail::diagonalWeight( 0.0, velocityWeight ) );
    return std::make_shared< InterArcStateContinuityConstraintSettings >( std::move( bodies ),
                                                                          std::move( connectionEpochsByBody ),
                                                                          std::move( weightMatricesByBody ),
                                                                          constraintScalingFactor,
                                                                          std::move( arcPairsByBody ) );
}

//! Build a soft-prior settings object with arbitrary body-specific, possibly per-boundary, 6x6 PSD weight matrices.
inline std::shared_ptr< InterArcStateContinuityConstraintSettings > generalContinuity(
        std::vector< std::string > bodies,
        InterArcStateContinuityConstraintSettings::EpochMap connectionEpochsByBody,
        InterArcStateContinuityConstraintSettings::WeightMatrixMap weightMatricesByBody,
        double constraintScalingFactor = 1.0,
        InterArcStateContinuityConstraintSettings::ArcPairMap arcPairsByBody = {} )
{
    return std::make_shared< InterArcStateContinuityConstraintSettings >( std::move( bodies ),
                                                                          std::move( connectionEpochsByBody ),
                                                                          std::move( weightMatricesByBody ),
                                                                          constraintScalingFactor,
                                                                          std::move( arcPairsByBody ) );
}

}  // namespace simulation_setup
}  // namespace tudat

#endif  // TUDAT_INTERARCSTATECONTINUITYCONSTRAINTSETTINGS_H
