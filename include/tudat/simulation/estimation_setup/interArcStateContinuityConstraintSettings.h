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
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>

namespace tudat
{
namespace simulation_setup
{

//! User-facing soft-continuity constraint between consecutive multi-arc translational arcs of a single body.
//! The cost added to the LSQ target is, per pair (Lari et al. 2021 Eq. 28):
//!   pairCost = stateDiscrepancy^T * scaledConstraintWeight * stateDiscrepancy
//! where scaledConstraintWeight is the constraint weight matrix divided by the product of the constraint scaling
//! factor and the total constrained dimension. Larger constraint scaling factors weaken the penalty.
class InterArcStateContinuityConstraintSettings
{
public:
    InterArcStateContinuityConstraintSettings( std::string body,
                                               std::vector< double > connectionEpochs,
                                               std::vector< Eigen::Matrix< double, 6, 6 > > weightMatrices,
                                               std::vector< double > constraintScalingFactors,
                                               std::vector< std::pair< int, int > > arcPairs = {} );

    const std::string& body( ) const
    {
        return body_;
    }
    const std::vector< double >& connectionEpochs( ) const
    {
        return connectionEpochs_;
    }
    const std::vector< Eigen::Matrix< double, 6, 6 > >& weightMatrices( ) const
    {
        return weightMatrices_;
    }
    const std::vector< double >& constraintScalingFactors( ) const
    {
        return constraintScalingFactors_;
    }
    const std::vector< std::pair< int, int > >& arcPairs( ) const
    {
        return arcPairs_;
    }

    //! Resolve the weight matrix for the i-th pair (handles 1-or-n broadcasting).
    const Eigen::Matrix< double, 6, 6 >& weightMatrixForPair( std::size_t pairIndex ) const
    {
        return weightMatrices_.size( ) == 1 ? weightMatrices_.front( ) : weightMatrices_.at( pairIndex );
    }

    //! Resolve the constraint scaling factor for the i-th pair (handles 1-or-n broadcasting).
    double constraintScalingFactorForPair( std::size_t pairIndex ) const
    {
        return constraintScalingFactors_.size( ) == 1 ? constraintScalingFactors_.front( ) : constraintScalingFactors_.at( pairIndex );
    }

    //! Number of regularized boundaries (either the explicit arcPairs count or the connection-epoch count).
    std::size_t numberOfPairs( ) const
    {
        return connectionEpochs_.size( );
    }

private:
    void validate( ) const;
    void validateWeightMatrix( const Eigen::Matrix< double, 6, 6 >& constraintWeightMatrix, std::size_t entryIndex ) const;

    std::string body_;
    std::vector< double > connectionEpochs_;
    std::vector< Eigen::Matrix< double, 6, 6 > > weightMatrices_;
    std::vector< double > constraintScalingFactors_;
    std::vector< std::pair< int, int > > arcPairs_;
};

namespace detail
{

inline Eigen::Matrix< double, 6, 6 > diagonalWeight( double position, double velocity )
{
    Eigen::Matrix< double, 6, 6 > constraintWeightMatrix = Eigen::Matrix< double, 6, 6 >::Zero( );
    constraintWeightMatrix( 0, 0 ) = position;
    constraintWeightMatrix( 1, 1 ) = position;
    constraintWeightMatrix( 2, 2 ) = position;
    constraintWeightMatrix( 3, 3 ) = velocity;
    constraintWeightMatrix( 4, 4 ) = velocity;
    constraintWeightMatrix( 5, 5 ) = velocity;
    return constraintWeightMatrix;
}

}  // namespace detail

//! Build a soft-prior settings object with both position and velocity continuity.
inline std::shared_ptr< InterArcStateContinuityConstraintSettings > fullStateContinuity(
        std::string body,
        std::vector< double > connectionEpochs,
        double positionWeight = 1.0,
        double velocityWeight = 1.0,
        double constraintScalingFactor = 1.0,
        std::vector< std::pair< int, int > > arcPairs = {} )
{
    return std::make_shared< InterArcStateContinuityConstraintSettings >(
            std::move( body ),
            std::move( connectionEpochs ),
            std::vector< Eigen::Matrix< double, 6, 6 > >{ detail::diagonalWeight( positionWeight, velocityWeight ) },
            std::vector< double >{ constraintScalingFactor },
            std::move( arcPairs ) );
}

//! Build a settings object with position-only continuity (velocity rows/columns of the weight matrix zeroed).
inline std::shared_ptr< InterArcStateContinuityConstraintSettings > positionOnlyContinuity(
        std::string body,
        std::vector< double > connectionEpochs,
        double positionWeight = 1.0,
        double constraintScalingFactor = 1.0,
        std::vector< std::pair< int, int > > arcPairs = {} )
{
    return std::make_shared< InterArcStateContinuityConstraintSettings >(
            std::move( body ),
            std::move( connectionEpochs ),
            std::vector< Eigen::Matrix< double, 6, 6 > >{ detail::diagonalWeight( positionWeight, 0.0 ) },
            std::vector< double >{ constraintScalingFactor },
            std::move( arcPairs ) );
}

//! Build a settings object with velocity-only continuity (position rows/columns of the weight matrix zeroed).
inline std::shared_ptr< InterArcStateContinuityConstraintSettings > velocityOnlyContinuity(
        std::string body,
        std::vector< double > connectionEpochs,
        double velocityWeight = 1.0,
        double constraintScalingFactor = 1.0,
        std::vector< std::pair< int, int > > arcPairs = {} )
{
    return std::make_shared< InterArcStateContinuityConstraintSettings >(
            std::move( body ),
            std::move( connectionEpochs ),
            std::vector< Eigen::Matrix< double, 6, 6 > >{ detail::diagonalWeight( 0.0, velocityWeight ) },
            std::vector< double >{ constraintScalingFactor },
            std::move( arcPairs ) );
}

//! Build a soft-prior settings object with arbitrary (possibly per-boundary, possibly dense) 6x6 PSD weight matrices.
inline std::shared_ptr< InterArcStateContinuityConstraintSettings > generalContinuity(
        std::string body,
        std::vector< double > connectionEpochs,
        std::vector< Eigen::Matrix< double, 6, 6 > > weightMatrices,
        double constraintScalingFactor = 1.0,
        std::vector< std::pair< int, int > > arcPairs = {} )
{
    return std::make_shared< InterArcStateContinuityConstraintSettings >( std::move( body ),
                                                                          std::move( connectionEpochs ),
                                                                          std::move( weightMatrices ),
                                                                          std::vector< double >{ constraintScalingFactor },
                                                                          std::move( arcPairs ) );
}

}  // namespace simulation_setup
}  // namespace tudat

#endif  // TUDAT_INTERARCSTATECONTINUITYCONSTRAINTSETTINGS_H
