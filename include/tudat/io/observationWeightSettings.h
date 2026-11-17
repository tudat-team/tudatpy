/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATION_WEIGHT_SETTINGS_H
#define TUDAT_OBSERVATION_WEIGHT_SETTINGS_H

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>

namespace tudat
{
namespace observation_models
{

//! Function to check that weight values are finite and non-negative
inline void checkWeightValues( const Eigen::VectorXd& weights, const std::string& weightDescription )
{
    if( !weights.allFinite( ) || ( weights.array( ) < 0.0 ).any( ) )
    {
        throw std::runtime_error( "Error in observation weight settings, " + weightDescription + " must be finite and non-negative." );
    }
}

//! Function to check that a weight matrix is square, finite and symmetric, with non-negative diagonal
inline void checkWeightMatrix( const Eigen::MatrixXd& weightMatrix, const std::string& weightDescription )
{
    if( weightMatrix.rows( ) != weightMatrix.cols( ) || !weightMatrix.allFinite( ) )
    {
        throw std::runtime_error( "Error in observation weight settings, " + weightDescription + " must be square and finite." );
    }
    checkWeightValues( weightMatrix.diagonal( ), "the diagonal of " + weightDescription );
    const double matrixScale = weightMatrix.size( ) > 0 ? weightMatrix.cwiseAbs( ).maxCoeff( ) : 0.0;
    for( Eigen::Index i = 0; i < weightMatrix.rows( ); ++i )
    {
        for( Eigen::Index j = i + 1; j < weightMatrix.cols( ); ++j )
        {
            if( std::abs( weightMatrix( i, j ) - weightMatrix( j, i ) ) > 1.0E-12 * matrixScale )
            {
                throw std::runtime_error( "Error in observation weight settings, " + weightDescription + " must be symmetric." );
            }
        }
    }
}

//! Base class for the weights of the observations in a single observation set
class ObservationWeightSettings
{
public:
    virtual ~ObservationWeightSettings( ) {}

    //! Check that the weights are consistent with the number of observations and the size of a single observation
    virtual void validateDimensions( const std::size_t numberOfObservations, const unsigned int singleObservationSize ) const = 0;

    //! Return the weight block of a single observation. Not available for set-level weights.
    virtual Eigen::MatrixXd getWeightBlock( const std::size_t observationIndex, const unsigned int singleObservationSize ) const = 0;

    //! Remove the weights of a single observation
    virtual void removeObservation( const std::size_t observationIndex, const unsigned int singleObservationSize ) = 0;

protected:
    //! Function to throw an error for inconsistent dimensions
    void throwDimensionError( ) const
    {
        throw std::runtime_error( "Error in observation weight settings, dimensions are inconsistent with the observations." );
    }
};

//! One scalar weight for every observation and component
class ConstantScalarWeightSettings : public ObservationWeightSettings
{
public:
    ConstantScalarWeightSettings( const double weight ): weight_( weight )
    {
        checkWeightValues( Eigen::VectorXd::Constant( 1, weight ), "the weight" );
    }

    Eigen::MatrixXd getWeightBlock( const std::size_t, const unsigned int singleObservationSize ) const override
    {
        return weight_ * Eigen::MatrixXd::Identity( singleObservationSize, singleObservationSize );
    }

    void validateDimensions( const std::size_t, const unsigned int ) const override {}

    void removeObservation( const std::size_t, const unsigned int ) override {}

private:
    double weight_;
};

//! One scalar weight per observation, for all its components
class ScalarPerObservationWeightSettings : public ObservationWeightSettings
{
public:
    ScalarPerObservationWeightSettings( const std::vector< double >& weights ): weights_( weights )
    {
        checkWeightValues( Eigen::Map< const Eigen::VectorXd >( weights.data( ), weights.size( ) ), "the weights" );
    }

    Eigen::MatrixXd getWeightBlock( const std::size_t observationIndex, const unsigned int singleObservationSize ) const override
    {
        return weights_.at( observationIndex ) * Eigen::MatrixXd::Identity( singleObservationSize, singleObservationSize );
    }

    void validateDimensions( const std::size_t numberOfObservations, const unsigned int ) const override
    {
        if( weights_.size( ) != numberOfObservations )
        {
            throwDimensionError( );
        }
    }

    void removeObservation( const std::size_t observationIndex, const unsigned int ) override
    {
        weights_.erase( weights_.begin( ) + observationIndex );
    }

private:
    std::vector< double > weights_;
};

//! One diagonal (one weight per component) per observation
class DiagonalPerObservationWeightSettings : public ObservationWeightSettings
{
public:
    DiagonalPerObservationWeightSettings( const std::vector< Eigen::VectorXd >& weights ): weights_( weights )
    {
        for( const Eigen::VectorXd& weight : weights_ )
        {
            checkWeightValues( weight, "the weights" );
            if( weight.size( ) != weights_.front( ).size( ) )
            {
                throw std::runtime_error( "Error in observation weight settings, all diagonals must have the same size." );
            }
        }
    }

    const std::vector< Eigen::VectorXd >& getWeights( ) const
    {
        return weights_;
    }

    //! Replace the diagonal weights of a single observation
    void setWeight( const std::size_t observationIndex, const Eigen::VectorXd& weight )
    {
        checkWeightValues( weight, "the weights" );
        if( weight.size( ) != weights_.at( observationIndex ).size( ) )
        {
            throwDimensionError( );
        }
        weights_.at( observationIndex ) = weight;
    }

    Eigen::MatrixXd getWeightBlock( const std::size_t observationIndex, const unsigned int ) const override
    {
        return weights_.at( observationIndex ).asDiagonal( );
    }

    void validateDimensions( const std::size_t numberOfObservations, const unsigned int singleObservationSize ) const override
    {
        if( weights_.size( ) != numberOfObservations ||
            ( !weights_.empty( ) && weights_.front( ).size( ) != static_cast< Eigen::Index >( singleObservationSize ) ) )
        {
            throwDimensionError( );
        }
    }

    void removeObservation( const std::size_t observationIndex, const unsigned int ) override
    {
        weights_.erase( weights_.begin( ) + observationIndex );
    }

private:
    std::vector< Eigen::VectorXd > weights_;
};

//! One weight matrix (over the components of a single observation), the same for every observation
class ConstantBlockWeightSettings : public ObservationWeightSettings
{
public:
    ConstantBlockWeightSettings( const Eigen::MatrixXd& weightBlock ): weightBlock_( weightBlock )
    {
        checkWeightMatrix( weightBlock_, "the weight block" );
    }

    Eigen::MatrixXd getWeightBlock( const std::size_t, const unsigned int ) const override
    {
        return weightBlock_;
    }

    void validateDimensions( const std::size_t, const unsigned int singleObservationSize ) const override
    {
        if( weightBlock_.rows( ) != static_cast< Eigen::Index >( singleObservationSize ) )
        {
            throwDimensionError( );
        }
    }

    void removeObservation( const std::size_t, const unsigned int ) override {}

private:
    Eigen::MatrixXd weightBlock_;
};

//! One weight matrix (over the components of a single observation) per observation
class BlockPerObservationWeightSettings : public ObservationWeightSettings
{
public:
    BlockPerObservationWeightSettings( const std::vector< Eigen::MatrixXd >& weightBlocks ): weightBlocks_( weightBlocks )
    {
        for( const Eigen::MatrixXd& weightBlock : weightBlocks_ )
        {
            checkWeightMatrix( weightBlock, "the weight blocks" );
            if( weightBlock.rows( ) != weightBlocks_.front( ).rows( ) )
            {
                throw std::runtime_error( "Error in observation weight settings, all weight blocks must have the same size." );
            }
        }
    }

    Eigen::MatrixXd getWeightBlock( const std::size_t observationIndex, const unsigned int ) const override
    {
        return weightBlocks_.at( observationIndex );
    }

    void validateDimensions( const std::size_t numberOfObservations, const unsigned int singleObservationSize ) const override
    {
        if( weightBlocks_.size( ) != numberOfObservations ||
            ( !weightBlocks_.empty( ) && weightBlocks_.front( ).rows( ) != static_cast< Eigen::Index >( singleObservationSize ) ) )
        {
            throwDimensionError( );
        }
    }

    void removeObservation( const std::size_t observationIndex, const unsigned int ) override
    {
        weightBlocks_.erase( weightBlocks_.begin( ) + observationIndex );
    }

private:
    std::vector< Eigen::MatrixXd > weightBlocks_;
};

//! One weight matrix for all components of all observations in the set (including correlations between observations)
class SetBlockWeightSettings : public ObservationWeightSettings
{
public:
    SetBlockWeightSettings( const Eigen::MatrixXd& weightMatrix ): weightMatrix_( weightMatrix )
    {
        checkWeightMatrix( weightMatrix_, "the weight matrix" );
    }

    const Eigen::MatrixXd& getWeightMatrix( ) const
    {
        return weightMatrix_;
    }

    Eigen::MatrixXd getWeightBlock( const std::size_t, const unsigned int ) const override
    {
        throw std::runtime_error( "Error, single-observation weight blocks cannot be retrieved for set-level observation weights." );
    }

    void validateDimensions( const std::size_t numberOfObservations, const unsigned int singleObservationSize ) const override
    {
        if( weightMatrix_.rows( ) != static_cast< Eigen::Index >( numberOfObservations * singleObservationSize ) )
        {
            throwDimensionError( );
        }
    }

    //! Remove the rows and columns of a single observation (the matrix is ordered by observation and then by component)
    void removeObservation( const std::size_t observationIndex, const unsigned int singleObservationSize ) override
    {
        const Eigen::Index removedStart = static_cast< Eigen::Index >( observationIndex * singleObservationSize );
        const Eigen::Index remainingSize = weightMatrix_.rows( ) - singleObservationSize;
        const Eigen::Index trailingSize = remainingSize - removedStart;
        Eigen::MatrixXd reducedMatrix( remainingSize, remainingSize );
        reducedMatrix.topLeftCorner( removedStart, removedStart ) = weightMatrix_.topLeftCorner( removedStart, removedStart );
        reducedMatrix.topRightCorner( removedStart, trailingSize ) = weightMatrix_.topRightCorner( removedStart, trailingSize );
        reducedMatrix.bottomLeftCorner( trailingSize, removedStart ) = weightMatrix_.bottomLeftCorner( trailingSize, removedStart );
        reducedMatrix.bottomRightCorner( trailingSize, trailingSize ) = weightMatrix_.bottomRightCorner( trailingSize, trailingSize );
        weightMatrix_ = reducedMatrix;
    }

private:
    Eigen::MatrixXd weightMatrix_;
};

//! Create settings for one scalar weight for every observation and component
inline std::shared_ptr< ObservationWeightSettings > constantScalarWeights( const double weight )
{
    return std::make_shared< ConstantScalarWeightSettings >( weight );
}

//! Create settings for one scalar weight per observation
inline std::shared_ptr< ObservationWeightSettings > scalarPerObservationWeights( const std::vector< double >& weights )
{
    return std::make_shared< ScalarPerObservationWeightSettings >( weights );
}

//! Create settings for one diagonal (one weight per component) per observation
inline std::shared_ptr< ObservationWeightSettings > diagonalPerObservationWeights( const std::vector< Eigen::VectorXd >& weights )
{
    return std::make_shared< DiagonalPerObservationWeightSettings >( weights );
}

//! Create settings for one weight matrix (over the components of a single observation), the same for every observation
inline std::shared_ptr< ObservationWeightSettings > constantBlockWeights( const Eigen::MatrixXd& weightBlock )
{
    return std::make_shared< ConstantBlockWeightSettings >( weightBlock );
}

//! Create settings for one weight matrix (over the components of a single observation) per observation
inline std::shared_ptr< ObservationWeightSettings > blockPerObservationWeights( const std::vector< Eigen::MatrixXd >& weightBlocks )
{
    return std::make_shared< BlockPerObservationWeightSettings >( weightBlocks );
}

//! Create settings for one weight matrix for all components of all observations in the set
inline std::shared_ptr< ObservationWeightSettings > setBlockWeights( const Eigen::MatrixXd& weightMatrix )
{
    return std::make_shared< SetBlockWeightSettings >( weightMatrix );
}

}  // namespace observation_models
}  // namespace tudat

#endif  // TUDAT_OBSERVATION_WEIGHT_SETTINGS_H
