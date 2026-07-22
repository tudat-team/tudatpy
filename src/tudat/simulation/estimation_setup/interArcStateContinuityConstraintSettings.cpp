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
#include <set>

#include <Eigen/Eigenvalues>

namespace tudat
{
namespace simulation_setup
{

namespace
{

template< typename BodyMapType >
void validateBodyMapKeys( const BodyMapType& values, const std::vector< std::string >& bodies, const std::string& argumentName )
{
    const std::set< std::string > expectedBodies( bodies.begin( ), bodies.end( ) );
    for( const auto& bodyEntry : values )
    {
        if( expectedBodies.count( bodyEntry.first ) == 0 )
        {
            throw std::runtime_error( "Inter-arc continuity " + argumentName + " contains unknown body \"" + bodyEntry.first + "\"." );
        }
    }
    for( const auto& body : expectedBodies )
    {
        if( values.count( body ) == 0 )
        {
            throw std::runtime_error( "Inter-arc continuity " + argumentName + " is missing body \"" + body + "\"." );
        }
    }
}

Eigen::Vector3d resolveCartesianStateWeights( const InterArcStateContinuityConstraintSettings::CartesianStateWeight& weights )
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

InterArcStateContinuityConstraintSettings::WeightMatrixMap createUniformBodyWeightMap( const std::vector< std::string >& bodies,
                                                                                       const Eigen::MatrixXd& weightMatrix )
{
    InterArcStateContinuityConstraintSettings::WeightMatrixMap weightMatricesByBody;
    for( const auto& body : bodies )
    {
        weightMatricesByBody[ body ] = { weightMatrix };
    }
    return weightMatricesByBody;
}

InterArcStateContinuityConstraintSettings::WeightMatrixMap createCartesianStateWeightMap(
        const std::vector< std::string >& bodies,
        const InterArcStateContinuityConstraintSettings::CartesianStateWeightMap& positionWeightsByBody,
        const InterArcStateContinuityConstraintSettings::CartesianStateWeightMap& velocityWeightsByBody )
{
    validateBodyMapKeys( positionWeightsByBody, bodies, "position_weights" );
    validateBodyMapKeys( velocityWeightsByBody, bodies, "velocity_weights" );

    InterArcStateContinuityConstraintSettings::WeightMatrixMap weightMatricesByBody;
    for( const auto& body : bodies )
    {
        weightMatricesByBody[ body ] = { createCartesianStateWeightMatrix(
                resolveCartesianStateWeights( positionWeightsByBody.at( body ) ),
                resolveCartesianStateWeights( velocityWeightsByBody.at( body ) ) ) };
    }
    return weightMatricesByBody;
}

}  // namespace

Eigen::MatrixXd createCartesianStateWeightMatrix( const Eigen::Vector3d& positionWeights, const Eigen::Vector3d& velocityWeights )
{
    Eigen::MatrixXd constraintWeightMatrix = Eigen::MatrixXd::Zero( 6, 6 );
    constraintWeightMatrix.diagonal( ).head< 3 >( ) = positionWeights;
    constraintWeightMatrix.diagonal( ).tail< 3 >( ) = velocityWeights;
    return constraintWeightMatrix;
}

std::shared_ptr< InterArcStateContinuityConstraintSettings > fullStateContinuity(
        std::vector< std::string > bodies,
        InterArcStateContinuityConstraintSettings::EpochMap connectionEpochsByBody,
        const double positionWeight,
        const double velocityWeight,
        const double constraintScalingFactor,
        InterArcStateContinuityConstraintSettings::ArcPairMap arcPairsByBody )
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
        InterArcStateContinuityConstraintSettings::EpochMap connectionEpochsByBody,
        InterArcStateContinuityConstraintSettings::CartesianStateWeightMap positionWeightsByBody,
        InterArcStateContinuityConstraintSettings::CartesianStateWeightMap velocityWeightsByBody,
        const double constraintScalingFactor,
        InterArcStateContinuityConstraintSettings::ArcPairMap arcPairsByBody )
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
        InterArcStateContinuityConstraintSettings::EpochMap connectionEpochsByBody,
        const double positionWeight,
        const double constraintScalingFactor,
        InterArcStateContinuityConstraintSettings::ArcPairMap arcPairsByBody )
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
        InterArcStateContinuityConstraintSettings::EpochMap connectionEpochsByBody,
        InterArcStateContinuityConstraintSettings::CartesianStateWeightMap positionWeightsByBody,
        const double constraintScalingFactor,
        InterArcStateContinuityConstraintSettings::ArcPairMap arcPairsByBody )
{
    InterArcStateContinuityConstraintSettings::CartesianStateWeightMap zeroVelocityWeightsByBody;
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
        InterArcStateContinuityConstraintSettings::EpochMap connectionEpochsByBody,
        const double velocityWeight,
        const double constraintScalingFactor,
        InterArcStateContinuityConstraintSettings::ArcPairMap arcPairsByBody )
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
        InterArcStateContinuityConstraintSettings::EpochMap connectionEpochsByBody,
        InterArcStateContinuityConstraintSettings::CartesianStateWeightMap velocityWeightsByBody,
        const double constraintScalingFactor,
        InterArcStateContinuityConstraintSettings::ArcPairMap arcPairsByBody )
{
    InterArcStateContinuityConstraintSettings::CartesianStateWeightMap zeroPositionWeightsByBody;
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
        InterArcStateContinuityConstraintSettings::EpochMap connectionEpochsByBody,
        InterArcStateContinuityConstraintSettings::WeightMatrixMap weightMatricesByBody,
        const double constraintScalingFactor,
        InterArcStateContinuityConstraintSettings::ArcPairMap arcPairsByBody )
{
    return std::make_shared< InterArcStateContinuityConstraintSettings >( std::move( bodies ),
                                                                          std::move( connectionEpochsByBody ),
                                                                          std::move( weightMatricesByBody ),
                                                                          constraintScalingFactor,
                                                                          std::move( arcPairsByBody ) );
}

std::shared_ptr< InterArcStateContinuityConstraintSettings > generalContinuity(
        std::vector< std::string > bodies,
        InterArcStateContinuityConstraintSettings::EpochMap connectionEpochsByBody,
        InterArcStateContinuityConstraintSettings::WeightMatrixInputMap weightMatricesByBody,
        const double constraintScalingFactor,
        InterArcStateContinuityConstraintSettings::ArcPairMap arcPairsByBody )
{
    validateBodyMapKeys( weightMatricesByBody, bodies, "weight_matrices" );
    InterArcStateContinuityConstraintSettings::WeightMatrixMap resolvedWeightMatricesByBody;
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

InterArcStateContinuityConstraintSettings::InterArcStateContinuityConstraintSettings( std::vector< std::string > bodies,
                                                                                      EpochMap connectionEpochsByBody,
                                                                                      WeightMatrixMap weightMatricesByBody,
                                                                                      double constraintScalingFactor,
                                                                                      ArcPairMap arcPairsByBody ):
    bodies_( std::move( bodies ) ), connectionEpochsByBody_( std::move( connectionEpochsByBody ) ),
    weightMatricesByBody_( std::move( weightMatricesByBody ) ), constraintScalingFactor_( constraintScalingFactor ),
    arcPairsByBody_( std::move( arcPairsByBody ) )
{
    validate( );
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
    auto validateBodyMapKeys = [ &bodySet ]( const auto& bodyMap, const std::string& mapName ) {
        for( const auto& bodyEntry : bodyMap )
        {
            if( bodySet.count( bodyEntry.first ) == 0 )
            {
                throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings: " + mapName + " contains unknown body \"" +
                                          bodyEntry.first + "\"." );
            }
        }
        for( const auto& body : bodySet )
        {
            if( bodyMap.count( body ) == 0 )
            {
                throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings: " + mapName + " is missing body \"" + body +
                                          "\"." );
            }
        }
    };
    validateBodyMapKeys( connectionEpochsByBody_, "connectionEpochsByBody" );
    validateBodyMapKeys( weightMatricesByBody_, "weightMatricesByBody" );
    for( const auto& bodyEntry : arcPairsByBody_ )
    {
        if( bodySet.count( bodyEntry.first ) == 0 )
        {
            throw std::runtime_error( "Error in InterArcStateContinuityConstraintSettings: arcPairsByBody contains unknown body \"" +
                                      bodyEntry.first + "\"." );
        }
    }
    if( !( constraintScalingFactor_ > 0.0 ) )
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
