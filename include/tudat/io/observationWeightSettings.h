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
#include <iterator>
#include <memory>
#include <set>
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
    for( Eigen::Index i = 0; i < weightMatrix.rows( ); ++i )
    {
        for( Eigen::Index j = i + 1; j < weightMatrix.cols( ); ++j )
        {
            if( std::abs( weightMatrix( i, j ) - weightMatrix( j, i ) ) >
                1.0E-12 * std::max( { 1.0, std::abs( weightMatrix( i, j ) ), std::abs( weightMatrix( j, i ) ) } ) )
            {
                throw std::runtime_error( "Error in observation weight settings, " + weightDescription + " must be symmetric." );
            }
        }
    }
}

//! Base class for the weights of the observations in a single observation set (or of a subset of these observations)
class ObservationWeightSettings
{
public:
    //! Constructor
    /*!
     * \param observationIndices Indices of the observations for which weights are given (empty for all observations). Only
     * possible for weights of a TrackingData object; the weights of an observation set in a dataset are given for all observations.
     */
    ObservationWeightSettings( const std::vector< unsigned int >& observationIndices = std::vector< unsigned int >( ) ):
        observationIndices_( observationIndices )
    {
        if( std::set< unsigned int >( observationIndices.begin( ), observationIndices.end( ) ).size( ) != observationIndices.size( ) )
        {
            throw std::runtime_error( "Error in observation weight settings, observation indices must be unique." );
        }
    }

    virtual ~ObservationWeightSettings( ) {}

    //! Return the indices of the observations for which weights are given (empty for all observations)
    const std::vector< unsigned int >& getObservationIndices( ) const
    {
        return observationIndices_;
    }

    //! Check that the weights are consistent with the observations to which they are applied
    void validateDimensions( const std::size_t numberOfObservations, const unsigned int singleObservationSize ) const
    {
        if( !observationIndices_.empty( ) &&
            *std::max_element( observationIndices_.begin( ), observationIndices_.end( ) ) >= numberOfObservations )
        {
            throw std::runtime_error(
                    "Error in observation weight settings, observation indices must be smaller than the number of "
                    "observations (" +
                    std::to_string( numberOfObservations ) + ")." );
        }
        validateWeightDimensions( observationIndices_.empty( ) ? numberOfObservations : observationIndices_.size( ),
                                  singleObservationSize );
    }

    //! Return the weight block of a single weight entry (entry i belongs to observation observationIndices_[ i ], or to observation i
    //! if no indices are given). Not available for set-level weights.
    virtual Eigen::MatrixXd getWeightBlock( const std::size_t weightIndex, const unsigned int singleObservationSize ) const = 0;

    //! Remove a single observation, return whether weights for any observation remain
    bool removeObservation( const unsigned int observationIndex, const unsigned int singleObservationSize )
    {
        if( observationIndices_.empty( ) )
        {
            removeWeightEntry( observationIndex, singleObservationSize );
            return true;
        }

        auto indexIterator = std::find( observationIndices_.begin( ), observationIndices_.end( ), observationIndex );
        if( indexIterator != observationIndices_.end( ) )
        {
            removeWeightEntry( std::distance( observationIndices_.begin( ), indexIterator ), singleObservationSize );
            observationIndices_.erase( indexIterator );
        }
        for( unsigned int& index : observationIndices_ )
        {
            if( index > observationIndex )
            {
                index--;
            }
        }
        return !observationIndices_.empty( );
    }

protected:
    //! Check that the weights are consistent with the given number of weight entries and size of a single observation
    virtual void validateWeightDimensions( const std::size_t numberOfWeights, const unsigned int singleObservationSize ) const = 0;

    //! Remove a single weight entry
    virtual void removeWeightEntry( const std::size_t weightIndex, const unsigned int singleObservationSize ) = 0;

    //! Function to throw an error for inconsistent dimensions
    void throwDimensionError( ) const
    {
        throw std::runtime_error( "Error in observation weight settings, dimensions are inconsistent with the observations." );
    }

private:
    //! Indices of the observations for which weights are given (empty for all observations)
    std::vector< unsigned int > observationIndices_;
};

//! One scalar weight for every observation and component
class ConstantScalarWeightSettings : public ObservationWeightSettings
{
public:
    ConstantScalarWeightSettings( const double weight,
                                  const std::vector< unsigned int >& observationIndices = std::vector< unsigned int >( ) ):
        ObservationWeightSettings( observationIndices ), weight_( weight )
    {
        checkWeightValues( Eigen::VectorXd::Constant( 1, weight ), "the weight" );
    }

    double getWeight( ) const
    {
        return weight_;
    }

    Eigen::MatrixXd getWeightBlock( const std::size_t, const unsigned int singleObservationSize ) const override
    {
        return weight_ * Eigen::MatrixXd::Identity( singleObservationSize, singleObservationSize );
    }

protected:
    void validateWeightDimensions( const std::size_t, const unsigned int ) const override {}

    void removeWeightEntry( const std::size_t, const unsigned int ) override {}

private:
    double weight_;
};

//! One scalar weight per observation, for all its components
class ScalarPerObservationWeightSettings : public ObservationWeightSettings
{
public:
    ScalarPerObservationWeightSettings( const std::vector< double >& weights,
                                        const std::vector< unsigned int >& observationIndices = std::vector< unsigned int >( ) ):
        ObservationWeightSettings( observationIndices ), weights_( weights )
    {
        checkWeightValues( Eigen::Map< const Eigen::VectorXd >( weights.data( ), weights.size( ) ), "the weights" );
    }

    const std::vector< double >& getWeights( ) const
    {
        return weights_;
    }

    Eigen::MatrixXd getWeightBlock( const std::size_t weightIndex, const unsigned int singleObservationSize ) const override
    {
        return weights_.at( weightIndex ) * Eigen::MatrixXd::Identity( singleObservationSize, singleObservationSize );
    }

protected:
    void validateWeightDimensions( const std::size_t numberOfWeights, const unsigned int ) const override
    {
        if( weights_.size( ) != numberOfWeights )
        {
            throwDimensionError( );
        }
    }

    void removeWeightEntry( const std::size_t weightIndex, const unsigned int ) override
    {
        weights_.erase( weights_.begin( ) + weightIndex );
    }

private:
    std::vector< double > weights_;
};

//! One diagonal (one weight per component) per observation
class DiagonalPerObservationWeightSettings : public ObservationWeightSettings
{
public:
    DiagonalPerObservationWeightSettings( const std::vector< Eigen::VectorXd >& weights,
                                          const std::vector< unsigned int >& observationIndices = std::vector< unsigned int >( ) ):
        ObservationWeightSettings( observationIndices ), weights_( weights )
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

    //! Replace the diagonal of a single weight entry
    void setWeight( const std::size_t weightIndex, const Eigen::VectorXd& weight )
    {
        checkWeightValues( weight, "the weights" );
        if( weight.size( ) != weights_.at( weightIndex ).size( ) )
        {
            throwDimensionError( );
        }
        weights_.at( weightIndex ) = weight;
    }

    Eigen::MatrixXd getWeightBlock( const std::size_t weightIndex, const unsigned int ) const override
    {
        return weights_.at( weightIndex ).asDiagonal( );
    }

protected:
    void validateWeightDimensions( const std::size_t numberOfWeights, const unsigned int singleObservationSize ) const override
    {
        if( weights_.size( ) != numberOfWeights ||
            ( !weights_.empty( ) && weights_.front( ).size( ) != static_cast< Eigen::Index >( singleObservationSize ) ) )
        {
            throwDimensionError( );
        }
    }

    void removeWeightEntry( const std::size_t weightIndex, const unsigned int ) override
    {
        weights_.erase( weights_.begin( ) + weightIndex );
    }

private:
    std::vector< Eigen::VectorXd > weights_;
};

//! One weight matrix (over the components of a single observation), the same for every observation
class ConstantBlockWeightSettings : public ObservationWeightSettings
{
public:
    ConstantBlockWeightSettings( const Eigen::MatrixXd& weightBlock,
                                 const std::vector< unsigned int >& observationIndices = std::vector< unsigned int >( ) ):
        ObservationWeightSettings( observationIndices ), weightBlock_( weightBlock )
    {
        checkWeightMatrix( weightBlock_, "the weight block" );
    }

    Eigen::MatrixXd getWeightBlock( const std::size_t, const unsigned int ) const override
    {
        return weightBlock_;
    }

protected:
    void validateWeightDimensions( const std::size_t, const unsigned int singleObservationSize ) const override
    {
        if( weightBlock_.rows( ) != static_cast< Eigen::Index >( singleObservationSize ) )
        {
            throwDimensionError( );
        }
    }

    void removeWeightEntry( const std::size_t, const unsigned int ) override {}

private:
    Eigen::MatrixXd weightBlock_;
};

//! One weight matrix (over the components of a single observation) per observation
class BlockPerObservationWeightSettings : public ObservationWeightSettings
{
public:
    BlockPerObservationWeightSettings( const std::vector< Eigen::MatrixXd >& weightBlocks,
                                       const std::vector< unsigned int >& observationIndices = std::vector< unsigned int >( ) ):
        ObservationWeightSettings( observationIndices ), weightBlocks_( weightBlocks )
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

    const std::vector< Eigen::MatrixXd >& getWeightBlocks( ) const
    {
        return weightBlocks_;
    }

    Eigen::MatrixXd getWeightBlock( const std::size_t weightIndex, const unsigned int ) const override
    {
        return weightBlocks_.at( weightIndex );
    }

protected:
    void validateWeightDimensions( const std::size_t numberOfWeights, const unsigned int singleObservationSize ) const override
    {
        if( weightBlocks_.size( ) != numberOfWeights ||
            ( !weightBlocks_.empty( ) && weightBlocks_.front( ).rows( ) != static_cast< Eigen::Index >( singleObservationSize ) ) )
        {
            throwDimensionError( );
        }
    }

    void removeWeightEntry( const std::size_t weightIndex, const unsigned int ) override
    {
        weightBlocks_.erase( weightBlocks_.begin( ) + weightIndex );
    }

private:
    std::vector< Eigen::MatrixXd > weightBlocks_;
};

//! One weight matrix for all components of all observations in the set (including correlations between observations). Always given
//! for all observations.
class SetBlockWeightSettings : public ObservationWeightSettings
{
public:
    SetBlockWeightSettings( const Eigen::MatrixXd& weightMatrix ): ObservationWeightSettings( ), weightMatrix_( weightMatrix )
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

protected:
    void validateWeightDimensions( const std::size_t numberOfWeights, const unsigned int singleObservationSize ) const override
    {
        if( weightMatrix_.rows( ) != static_cast< Eigen::Index >( numberOfWeights * singleObservationSize ) )
        {
            throwDimensionError( );
        }
    }

    //! Remove the rows and columns of a single observation (the matrix is ordered by observation and then by component)
    void removeWeightEntry( const std::size_t weightIndex, const unsigned int singleObservationSize ) override
    {
        const Eigen::Index removedStart = static_cast< Eigen::Index >( weightIndex * singleObservationSize );
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
inline std::shared_ptr< ObservationWeightSettings > constantScalarWeights(
        const double weight,
        const std::vector< unsigned int >& observationIndices = std::vector< unsigned int >( ) )
{
    return std::make_shared< ConstantScalarWeightSettings >( weight, observationIndices );
}

//! Create settings for one scalar weight per observation
inline std::shared_ptr< ObservationWeightSettings > scalarPerObservationWeights(
        const std::vector< double >& weights,
        const std::vector< unsigned int >& observationIndices = std::vector< unsigned int >( ) )
{
    return std::make_shared< ScalarPerObservationWeightSettings >( weights, observationIndices );
}

//! Create settings for one diagonal (one weight per component) per observation
inline std::shared_ptr< ObservationWeightSettings > diagonalPerObservationWeights(
        const std::vector< Eigen::VectorXd >& weights,
        const std::vector< unsigned int >& observationIndices = std::vector< unsigned int >( ) )
{
    return std::make_shared< DiagonalPerObservationWeightSettings >( weights, observationIndices );
}

//! Create settings for one weight matrix (over the components of a single observation), the same for every observation
inline std::shared_ptr< ObservationWeightSettings > constantBlockWeights(
        const Eigen::MatrixXd& weightBlock,
        const std::vector< unsigned int >& observationIndices = std::vector< unsigned int >( ) )
{
    return std::make_shared< ConstantBlockWeightSettings >( weightBlock, observationIndices );
}

//! Create settings for one weight matrix (over the components of a single observation) per observation
inline std::shared_ptr< ObservationWeightSettings > blockPerObservationWeights(
        const std::vector< Eigen::MatrixXd >& weightBlocks,
        const std::vector< unsigned int >& observationIndices = std::vector< unsigned int >( ) )
{
    return std::make_shared< BlockPerObservationWeightSettings >( weightBlocks, observationIndices );
}

//! Create settings for one weight matrix for all components of all observations in the set
inline std::shared_ptr< ObservationWeightSettings > setBlockWeights( const Eigen::MatrixXd& weightMatrix )
{
    return std::make_shared< SetBlockWeightSettings >( weightMatrix );
}

}  // namespace observation_models
}  // namespace tudat

#endif  // TUDAT_OBSERVATION_WEIGHT_SETTINGS_H
