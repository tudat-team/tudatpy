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
//!   q_pair = (1 / (mu * m_d)) * d^T C d
//! where d = x_right(t_c) - x_left(t_c). The weight matrix C selects which components are constrained
//! (e.g. position-only, velocity-only, full state) and how tightly. Larger mu weakens the penalty.
class InterArcStateContinuityConstraintSettings
{
public:
    InterArcStateContinuityConstraintSettings( std::string body,
                                               std::vector< double > connectionEpochs,
                                               std::vector< Eigen::Matrix< double, 6, 6 > > weightMatrices,
                                               std::vector< double > muValues,
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
    const std::vector< double >& muValues( ) const
    {
        return muValues_;
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

    //! Resolve mu for the i-th pair (handles 1-or-n broadcasting).
    double muForPair( std::size_t pairIndex ) const
    {
        return muValues_.size( ) == 1 ? muValues_.front( ) : muValues_.at( pairIndex );
    }

    //! Number of constrained boundaries (either the explicit arcPairs count or the connection-epoch count).
    std::size_t numberOfPairs( ) const
    {
        return connectionEpochs_.size( );
    }

private:
    void validate( ) const;
    void validateWeightMatrix( const Eigen::Matrix< double, 6, 6 >& C, std::size_t entryIndex ) const;

    std::string body_;
    std::vector< double > connectionEpochs_;
    std::vector< Eigen::Matrix< double, 6, 6 > > weightMatrices_;
    std::vector< double > muValues_;
    std::vector< std::pair< int, int > > arcPairs_;
};

namespace detail
{

inline Eigen::Matrix< double, 6, 6 > diagonalWeight( double position, double velocity )
{
    Eigen::Matrix< double, 6, 6 > C = Eigen::Matrix< double, 6, 6 >::Zero( );
    C( 0, 0 ) = position;
    C( 1, 1 ) = position;
    C( 2, 2 ) = position;
    C( 3, 3 ) = velocity;
    C( 4, 4 ) = velocity;
    C( 5, 5 ) = velocity;
    return C;
}

}  // namespace detail

//! Build a settings object with both position and velocity continuity (anisotropic weights via per-component values).
inline std::shared_ptr< InterArcStateContinuityConstraintSettings > fullStateContinuity(
        std::string body,
        std::vector< double > connectionEpochs,
        double positionWeight = 1.0,
        double velocityWeight = 1.0,
        double mu = 1.0,
        std::vector< std::pair< int, int > > arcPairs = {} )
{
    return std::make_shared< InterArcStateContinuityConstraintSettings >(
            std::move( body ),
            std::move( connectionEpochs ),
            std::vector< Eigen::Matrix< double, 6, 6 > >{ detail::diagonalWeight( positionWeight, velocityWeight ) },
            std::vector< double >{ mu },
            std::move( arcPairs ) );
}

//! Build a settings object with position-only continuity (velocity rows/columns of C zeroed → rank-deficient C).
inline std::shared_ptr< InterArcStateContinuityConstraintSettings > positionOnlyContinuity(
        std::string body,
        std::vector< double > connectionEpochs,
        double positionWeight = 1.0,
        double mu = 1.0,
        std::vector< std::pair< int, int > > arcPairs = {} )
{
    return std::make_shared< InterArcStateContinuityConstraintSettings >(
            std::move( body ),
            std::move( connectionEpochs ),
            std::vector< Eigen::Matrix< double, 6, 6 > >{ detail::diagonalWeight( positionWeight, 0.0 ) },
            std::vector< double >{ mu },
            std::move( arcPairs ) );
}

//! Build a settings object with velocity-only continuity (position rows/columns of C zeroed → rank-deficient C).
inline std::shared_ptr< InterArcStateContinuityConstraintSettings > velocityOnlyContinuity(
        std::string body,
        std::vector< double > connectionEpochs,
        double velocityWeight = 1.0,
        double mu = 1.0,
        std::vector< std::pair< int, int > > arcPairs = {} )
{
    return std::make_shared< InterArcStateContinuityConstraintSettings >(
            std::move( body ),
            std::move( connectionEpochs ),
            std::vector< Eigen::Matrix< double, 6, 6 > >{ detail::diagonalWeight( 0.0, velocityWeight ) },
            std::vector< double >{ mu },
            std::move( arcPairs ) );
}

//! Build a settings object with arbitrary (possibly per-boundary, possibly dense) 6x6 PSD weight matrices.
inline std::shared_ptr< InterArcStateContinuityConstraintSettings > generalContinuity(
        std::string body,
        std::vector< double > connectionEpochs,
        std::vector< Eigen::Matrix< double, 6, 6 > > weightMatrices,
        double mu = 1.0,
        std::vector< std::pair< int, int > > arcPairs = {} )
{
    return std::make_shared< InterArcStateContinuityConstraintSettings >( std::move( body ),
                                                                          std::move( connectionEpochs ),
                                                                          std::move( weightMatrices ),
                                                                          std::vector< double >{ mu },
                                                                          std::move( arcPairs ) );
}

}  // namespace simulation_setup
}  // namespace tudat

#endif  // TUDAT_INTERARCSTATECONTINUITYCONSTRAINTSETTINGS_H
