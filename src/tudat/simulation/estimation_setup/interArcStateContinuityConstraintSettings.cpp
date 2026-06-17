/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/estimation_setup/interArcStateContinuityConstraintSettings.h"

#include <algorithm>

#include <Eigen/Eigenvalues>

namespace tudat
{
namespace simulation_setup
{

InterArcStateContinuityConstraintSettings::InterArcStateContinuityConstraintSettings(
        std::string body,
        std::vector< double > connectionEpochs,
        std::vector< Eigen::Matrix< double, 6, 6 > > weightMatrices,
        std::vector< double > muValues,
        std::vector< std::pair< int, int > > arcPairs ):
    body_( std::move( body ) ), connectionEpochs_( std::move( connectionEpochs ) ), weightMatrices_( std::move( weightMatrices ) ),
    muValues_( std::move( muValues ) ), arcPairs_( std::move( arcPairs ) )
{
    validate( );
}

void InterArcStateContinuityConstraintSettings::validate( ) const
{
    if( body_.empty( ) )
    {
        throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings: body name is empty." );
    }
    if( connectionEpochs_.empty( ) )
    {
        throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings for body " + body_ + ": connectionEpochs is empty." );
    }
    if( !arcPairs_.empty( ) && arcPairs_.size( ) != connectionEpochs_.size( ) )
    {
        throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings for body " + body_ + ": arcPairs size (" +
                                  std::to_string( arcPairs_.size( ) ) + ") does not match connectionEpochs size (" +
                                  std::to_string( connectionEpochs_.size( ) ) + ")." );
    }
    for( const auto& pair : arcPairs_ )
    {
        if( pair.second != pair.first + 1 )
        {
            throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings for body " + body_ +
                                      ": only consecutive arc pairs are supported (got left=" + std::to_string( pair.first ) +
                                      ", right=" + std::to_string( pair.second ) + ")." );
        }
        if( pair.first < 0 )
        {
            throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings for body " + body_ + ": negative arc index (" +
                                      std::to_string( pair.first ) + ")." );
        }
    }
    if( weightMatrices_.size( ) != 1 && weightMatrices_.size( ) != connectionEpochs_.size( ) )
    {
        throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings for body " + body_ + ": weightMatrices size (" +
                                  std::to_string( weightMatrices_.size( ) ) + ") must be either 1 or " +
                                  std::to_string( connectionEpochs_.size( ) ) + "." );
    }
    if( muValues_.size( ) != 1 && muValues_.size( ) != connectionEpochs_.size( ) )
    {
        throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings for body " + body_ + ": muValues size (" +
                                  std::to_string( muValues_.size( ) ) + ") must be either 1 or " +
                                  std::to_string( connectionEpochs_.size( ) ) + "." );
    }
    for( std::size_t i = 0; i < muValues_.size( ); ++i )
    {
        if( !( muValues_[ i ] > 0.0 ) )
        {
            throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings for body " + body_ +
                                      ": mu must be strictly positive (entry " + std::to_string( i ) + " = " +
                                      std::to_string( muValues_[ i ] ) + ")." );
        }
    }
    for( std::size_t i = 0; i < weightMatrices_.size( ); ++i )
    {
        validateWeightMatrix( weightMatrices_[ i ], i );
    }
}

void InterArcStateContinuityConstraintSettings::validateWeightMatrix( const Eigen::Matrix< double, 6, 6 >& C, std::size_t entryIndex ) const
{
    const double cNorm = C.norm( );
    const double symmetryNorm = ( C - C.transpose( ) ).norm( );
    if( symmetryNorm > 1.0E-12 * std::max( cNorm, 1.0 ) )
    {
        throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings for body " + body_ + ": weight matrix entry " +
                                  std::to_string( entryIndex ) + " is not symmetric (asymmetry norm = " + std::to_string( symmetryNorm ) +
                                  ")." );
    }
    Eigen::SelfAdjointEigenSolver< Eigen::Matrix< double, 6, 6 > > eigenSolver( 0.5 * ( C + C.transpose( ) ) );
    if( eigenSolver.info( ) != Eigen::Success )
    {
        throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings for body " + body_ +
                                  ": eigen-decomposition of weight matrix entry " + std::to_string( entryIndex ) + " failed." );
    }
    const double largestAbsoluteEigenValue = eigenSolver.eigenvalues( ).cwiseAbs( ).maxCoeff( );
    const double psdTolerance = -1.0E-12 * largestAbsoluteEigenValue;
    const double smallestEigenValue = eigenSolver.eigenvalues( ).minCoeff( );
    if( smallestEigenValue < psdTolerance )
    {
        throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings for body " + body_ + ": weight matrix entry " +
                                  std::to_string( entryIndex ) +
                                  " is not positive semi-definite (smallest eigenvalue = " + std::to_string( smallestEigenValue ) + ")." );
    }
}

}  // namespace simulation_setup
}  // namespace tudat
