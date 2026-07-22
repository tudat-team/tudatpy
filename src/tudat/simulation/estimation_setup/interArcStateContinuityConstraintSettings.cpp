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
#include <cmath>
#include <set>

#include <Eigen/Eigenvalues>

namespace tudat
{
namespace simulation_setup
{

static Eigen::Vector3d resolveCartesianStateWeights( const std::variant< double, Eigen::VectorXd >& weights )
{
    if( const auto* scalarWeight = std::get_if< double >( &weights ) )
    {
        return Eigen::Vector3d::Constant( *scalarWeight );
    }
    const Eigen::VectorXd& componentWeights = std::get< Eigen::VectorXd >( weights );
    if( componentWeights.size( ) != 3 )
    {
        throw std::runtime_error( "Inter-arc continuity Cartesian component weights must be a scalar or a length-3 sequence." );
    }
    return componentWeights;
}

static std::map< std::string, std::vector< Eigen::MatrixXd > > createUniformBodyWeightMap( const std::vector< std::string >& bodies,
                                                                                           const Eigen::MatrixXd& weightMatrix )
{
    std::map< std::string, std::vector< Eigen::MatrixXd > > weightMatricesByBody;
    for( const auto& body : bodies )
    {
        weightMatricesByBody[ body ] = { weightMatrix };
    }
    return weightMatricesByBody;
}

static std::map< std::string, std::vector< Eigen::MatrixXd > > createCartesianStateWeightMap(
        const std::vector< std::string >& bodies,
        const std::map< std::string, std::variant< double, Eigen::VectorXd > >& positionWeightsByBody,
        const std::map< std::string, std::variant< double, Eigen::VectorXd > >& velocityWeightsByBody )
{
    std::map< std::string, std::vector< Eigen::MatrixXd > > weightMatricesByBody;
    for( const auto& body : bodies )
    {
        weightMatricesByBody[ body ] = { createCartesianStateWeightMatrix(
                resolveCartesianStateWeights( positionWeightsByBody.at( body ) ),
                resolveCartesianStateWeights( velocityWeightsByBody.at( body ) ) ) };
    }
    return weightMatricesByBody;
}

Eigen::MatrixXd createCartesianStateWeightMatrix( const Eigen::Vector3d& positionWeights, const Eigen::Vector3d& velocityWeights )
{
    Eigen::MatrixXd constraintWeightMatrix = Eigen::MatrixXd::Zero( 6, 6 );
    constraintWeightMatrix.diagonal( ).head< 3 >( ) = positionWeights;
    constraintWeightMatrix.diagonal( ).tail< 3 >( ) = velocityWeights;
    return constraintWeightMatrix;
}

std::shared_ptr< InterArcStateContinuityConstraintSettings > fullStateContinuity(
        std::vector< std::string > bodies,
        std::map< std::string, std::vector< double > > connectionEpochsByBody,
        const double positionWeight,
        const double velocityWeight,
        const double constraintScalingFactor,
        std::map< std::string, std::vector< std::pair< int, int > > > arcPairsByBody )
{
    auto weightMatricesByBody = createUniformBodyWeightMap(
            bodies,
            createCartesianStateWeightMatrix( Eigen::Vector3d::Constant( positionWeight ), Eigen::Vector3d::Constant( velocityWeight ) ) );
    return generalContinuity( std::move( bodies ),
                              std::move( connectionEpochsByBody ),
                              std::move( weightMatricesByBody ),
                              constraintScalingFactor,
                              std::move( arcPairsByBody ) );
}

std::shared_ptr< InterArcStateContinuityConstraintSettings > fullStateContinuity(
        std::vector< std::string > bodies,
        std::map< std::string, std::vector< double > > connectionEpochsByBody,
        std::map< std::string, std::variant< double, Eigen::VectorXd > > positionWeightsByBody,
        std::map< std::string, std::variant< double, Eigen::VectorXd > > velocityWeightsByBody,
        const double constraintScalingFactor,
        std::map< std::string, std::vector< std::pair< int, int > > > arcPairsByBody )
{
    auto weightMatricesByBody = createCartesianStateWeightMap( bodies, positionWeightsByBody, velocityWeightsByBody );
    return generalContinuity( std::move( bodies ),
                              std::move( connectionEpochsByBody ),
                              std::move( weightMatricesByBody ),
                              constraintScalingFactor,
                              std::move( arcPairsByBody ) );
}

std::shared_ptr< InterArcStateContinuityConstraintSettings > positionOnlyContinuity(
        std::vector< std::string > bodies,
        std::map< std::string, std::vector< double > > connectionEpochsByBody,
        const double positionWeight,
        const double constraintScalingFactor,
        std::map< std::string, std::vector< std::pair< int, int > > > arcPairsByBody )
{
    return fullStateContinuity( std::move( bodies ),
                                std::move( connectionEpochsByBody ),
                                positionWeight,
                                0.0,
                                constraintScalingFactor,
                                std::move( arcPairsByBody ) );
}

std::shared_ptr< InterArcStateContinuityConstraintSettings > positionOnlyContinuity(
        std::vector< std::string > bodies,
        std::map< std::string, std::vector< double > > connectionEpochsByBody,
        std::map< std::string, std::variant< double, Eigen::VectorXd > > positionWeightsByBody,
        const double constraintScalingFactor,
        std::map< std::string, std::vector< std::pair< int, int > > > arcPairsByBody )
{
    std::map< std::string, std::variant< double, Eigen::VectorXd > > zeroVelocityWeightsByBody;
    for( const auto& body : bodies )
    {
        zeroVelocityWeightsByBody[ body ] = 0.0;
    }
    return fullStateContinuity( std::move( bodies ),
                                std::move( connectionEpochsByBody ),
                                std::move( positionWeightsByBody ),
                                std::move( zeroVelocityWeightsByBody ),
                                constraintScalingFactor,
                                std::move( arcPairsByBody ) );
}

std::shared_ptr< InterArcStateContinuityConstraintSettings > velocityOnlyContinuity(
        std::vector< std::string > bodies,
        std::map< std::string, std::vector< double > > connectionEpochsByBody,
        const double velocityWeight,
        const double constraintScalingFactor,
        std::map< std::string, std::vector< std::pair< int, int > > > arcPairsByBody )
{
    return fullStateContinuity( std::move( bodies ),
                                std::move( connectionEpochsByBody ),
                                0.0,
                                velocityWeight,
                                constraintScalingFactor,
                                std::move( arcPairsByBody ) );
}

std::shared_ptr< InterArcStateContinuityConstraintSettings > velocityOnlyContinuity(
        std::vector< std::string > bodies,
        std::map< std::string, std::vector< double > > connectionEpochsByBody,
        std::map< std::string, std::variant< double, Eigen::VectorXd > > velocityWeightsByBody,
        const double constraintScalingFactor,
        std::map< std::string, std::vector< std::pair< int, int > > > arcPairsByBody )
{
    std::map< std::string, std::variant< double, Eigen::VectorXd > > zeroPositionWeightsByBody;
    for( const auto& body : bodies )
    {
        zeroPositionWeightsByBody[ body ] = 0.0;
    }
    return fullStateContinuity( std::move( bodies ),
                                std::move( connectionEpochsByBody ),
                                std::move( zeroPositionWeightsByBody ),
                                std::move( velocityWeightsByBody ),
                                constraintScalingFactor,
                                std::move( arcPairsByBody ) );
}

std::shared_ptr< InterArcStateContinuityConstraintSettings > generalContinuity(
        std::vector< std::string > bodies,
        std::map< std::string, std::vector< double > > connectionEpochsByBody,
        std::map< std::string, std::vector< Eigen::MatrixXd > > weightMatricesByBody,
        const double constraintScalingFactor,
        std::map< std::string, std::vector< std::pair< int, int > > > arcPairsByBody )
{
    return std::make_shared< InterArcStateContinuityConstraintSettings >( std::move( bodies ),
                                                                          std::move( connectionEpochsByBody ),
                                                                          std::move( weightMatricesByBody ),
                                                                          constraintScalingFactor,
                                                                          std::move( arcPairsByBody ) );
}

std::shared_ptr< InterArcStateContinuityConstraintSettings > generalContinuity(
        std::vector< std::string > bodies,
        std::map< std::string, std::vector< double > > connectionEpochsByBody,
        std::map< std::string, std::variant< Eigen::MatrixXd, std::vector< Eigen::MatrixXd > > > weightMatricesByBody,
        const double constraintScalingFactor,
        std::map< std::string, std::vector< std::pair< int, int > > > arcPairsByBody )
{
    std::map< std::string, std::vector< Eigen::MatrixXd > > resolvedWeightMatricesByBody;
    for( const auto& body : bodies )
    {
        const auto& weightMatrixInput = weightMatricesByBody.at( body );
        if( const auto* weightMatrix = std::get_if< Eigen::MatrixXd >( &weightMatrixInput ) )
        {
            resolvedWeightMatricesByBody[ body ] = { *weightMatrix };
        }
        else
        {
            resolvedWeightMatricesByBody[ body ] = std::get< std::vector< Eigen::MatrixXd > >( weightMatrixInput );
        }
    }
    return generalContinuity( std::move( bodies ),
                              std::move( connectionEpochsByBody ),
                              std::move( resolvedWeightMatricesByBody ),
                              constraintScalingFactor,
                              std::move( arcPairsByBody ) );
}

InterArcStateContinuityConstraintSettings::InterArcStateContinuityConstraintSettings(
        std::vector< std::string > bodies,
        std::map< std::string, std::vector< double > > connectionEpochsByBody,
        std::map< std::string, std::vector< Eigen::MatrixXd > > weightMatricesByBody,
        double constraintScalingFactor,
        std::map< std::string, std::vector< std::pair< int, int > > > arcPairsByBody ):
    bodies_( std::move( bodies ) ), connectionEpochsByBody_( std::move( connectionEpochsByBody ) ),
    weightMatricesByBody_( std::move( weightMatricesByBody ) ), constraintScalingFactor_( constraintScalingFactor ),
    arcPairsByBody_( std::move( arcPairsByBody ) )
{
    validate( );
}

const std::vector< std::string >& InterArcStateContinuityConstraintSettings::bodies( ) const
{
    return bodies_;
}

const std::map< std::string, std::vector< double > >& InterArcStateContinuityConstraintSettings::connectionEpochsByBody( ) const
{
    return connectionEpochsByBody_;
}

const std::map< std::string, std::vector< Eigen::MatrixXd > >& InterArcStateContinuityConstraintSettings::weightMatricesByBody( ) const
{
    return weightMatricesByBody_;
}

double InterArcStateContinuityConstraintSettings::constraintScalingFactor( ) const
{
    return constraintScalingFactor_;
}

const std::map< std::string, std::vector< std::pair< int, int > > >& InterArcStateContinuityConstraintSettings::arcPairsByBody( ) const
{
    return arcPairsByBody_;
}

const std::vector< double >& InterArcStateContinuityConstraintSettings::connectionEpochsForBody( const std::string& body ) const
{
    const auto iterator = connectionEpochsByBody_.find( body );
    if( iterator == connectionEpochsByBody_.end( ) )
    {
        throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings: no connection epochs found for body " + body + "." );
    }
    return iterator->second;
}

const std::vector< std::pair< int, int > >& InterArcStateContinuityConstraintSettings::arcPairsForBody( const std::string& body ) const
{
    static const std::vector< std::pair< int, int > > emptyArcPairs;
    const auto iterator = arcPairsByBody_.find( body );
    return iterator == arcPairsByBody_.end( ) ? emptyArcPairs : iterator->second;
}

const Eigen::MatrixXd& InterArcStateContinuityConstraintSettings::weightMatrixForBodyAndPair( const std::string& body,
                                                                                              std::size_t pairIndex ) const
{
    const auto iterator = weightMatricesByBody_.find( body );
    if( iterator == weightMatricesByBody_.end( ) )
    {
        throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings: no weight matrices found for body " + body + "." );
    }
    const auto& bodyWeightMatrices = iterator->second;
    if( bodyWeightMatrices.empty( ) )
    {
        throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings: empty weight matrix list for body " + body + "." );
    }
    return bodyWeightMatrices.size( ) == 1 ? bodyWeightMatrices.front( ) : bodyWeightMatrices.at( pairIndex );
}

std::size_t InterArcStateContinuityConstraintSettings::numberOfPairsForBody( const std::string& body ) const
{
    return connectionEpochsForBody( body ).size( );
}

std::size_t InterArcStateContinuityConstraintSettings::totalNumberOfPairs( ) const
{
    std::size_t totalNumberOfPairs = 0;
    for( const auto& body : bodies_ )
    {
        totalNumberOfPairs += numberOfPairsForBody( body );
    }
    return totalNumberOfPairs;
}

void InterArcStateContinuityConstraintSettings::validate( ) const
{
    if( bodies_.empty( ) )
    {
        throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings: body list is empty." );
    }
    std::set< std::string > bodySet;
    for( const auto& body : bodies_ )
    {
        if( body.empty( ) )
        {
            throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings: body name is empty." );
        }
        if( !bodySet.insert( body ).second )
        {
            throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings: duplicate body \"" + body + "\"." );
        }
    }
    if( !std::isfinite( constraintScalingFactor_ ) || !( constraintScalingFactor_ > 0.0 ) )
    {
        throw std::runtime_error(
                "Error in InterArcStateContinuityConstraintSettings: constraint scaling factor must be strictly positive (" +
                std::to_string( constraintScalingFactor_ ) + ")." );
    }
    for( const auto& body : bodies_ )
    {
        const auto& connectionEpochs = connectionEpochsByBody_.at( body );
        if( connectionEpochs.empty( ) )
        {
            throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings for body " + body +
                                      ": connection epochs are empty." );
        }
        for( const double connectionEpoch : connectionEpochs )
        {
            if( !std::isfinite( connectionEpoch ) )
            {
                throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings for body " + body +
                                          ": connection epochs must be finite." );
            }
        }
        const auto arcPairsIterator = arcPairsByBody_.find( body );
        if( arcPairsIterator != arcPairsByBody_.end( ) && arcPairsIterator->second.size( ) != connectionEpochs.size( ) )
        {
            throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings for body " + body + ": arcPairs size (" +
                                      std::to_string( arcPairsIterator->second.size( ) ) + ") does not match connection epochs size (" +
                                      std::to_string( connectionEpochs.size( ) ) + ")." );
        }
        if( arcPairsIterator != arcPairsByBody_.end( ) )
        {
            for( const auto& pair : arcPairsIterator->second )
            {
                if( pair.second != pair.first + 1 )
                {
                    throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings for body " + body +
                                              ": only consecutive arc pairs are supported (got left=" + std::to_string( pair.first ) +
                                              ", right=" + std::to_string( pair.second ) + ")." );
                }
                if( pair.first < 0 )
                {
                    throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings for body " + body +
                                              ": negative arc index (" + std::to_string( pair.first ) + ")." );
                }
            }
        }
        const auto& weightMatrices = weightMatricesByBody_.at( body );
        if( weightMatrices.size( ) != 1 && weightMatrices.size( ) != connectionEpochs.size( ) )
        {
            throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings for body " + body + ": weightMatrices size (" +
                                      std::to_string( weightMatrices.size( ) ) + ") must be either 1 or " +
                                      std::to_string( connectionEpochs.size( ) ) + "." );
        }
        for( std::size_t i = 0; i < weightMatrices.size( ); ++i )
        {
            validateWeightMatrix( body, weightMatrices[ i ], i );
        }
    }
}

void InterArcStateContinuityConstraintSettings::validateWeightMatrix( const std::string& body,
                                                                      const Eigen::MatrixXd& constraintWeightMatrix,
                                                                      std::size_t entryIndex ) const
{
    if( constraintWeightMatrix.rows( ) != 6 || constraintWeightMatrix.cols( ) != 6 )
    {
        throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings for body " + body + ": weight matrix entry " +
                                  std::to_string( entryIndex ) + " must be a 6x6 matrix for translational position/velocity continuity." );
    }
    if( !constraintWeightMatrix.allFinite( ) )
    {
        throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings for body " + body + ": weight matrix entry " +
                                  std::to_string( entryIndex ) + " contains a non-finite value." );
    }
    const double weightMatrixNorm = constraintWeightMatrix.norm( );
    const double symmetryNorm = ( constraintWeightMatrix - constraintWeightMatrix.transpose( ) ).norm( );
    if( symmetryNorm > 1.0E-12 * std::max( weightMatrixNorm, 1.0 ) )
    {
        throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings for body " + body + ": weight matrix entry " +
                                  std::to_string( entryIndex ) + " is not symmetric (asymmetry norm = " + std::to_string( symmetryNorm ) +
                                  ")." );
    }
    Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > eigenSolver( 0.5 * ( constraintWeightMatrix + constraintWeightMatrix.transpose( ) ) );
    if( eigenSolver.info( ) != Eigen::Success )
    {
        throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings for body " + body +
                                  ": eigen-decomposition of weight matrix entry " + std::to_string( entryIndex ) + " failed." );
    }
    const double largestAbsoluteEigenValue = eigenSolver.eigenvalues( ).cwiseAbs( ).maxCoeff( );
    const double psdTolerance = -1.0E-12 * largestAbsoluteEigenValue;
    const double smallestEigenValue = eigenSolver.eigenvalues( ).minCoeff( );
    if( smallestEigenValue < psdTolerance )
    {
        throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings for body " + body + ": weight matrix entry " +
                                  std::to_string( entryIndex ) +
                                  " is not positive semi-definite (smallest eigenvalue = " + std::to_string( smallestEigenValue ) + ")." );
    }
}

}  // namespace simulation_setup
}  // namespace tudat
