/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATION_WEIGHTS_H
#define TUDAT_OBSERVATION_WEIGHTS_H

#include <cstddef>
#include <optional>
#include <stdexcept>
#include <vector>

#include <Eigen/Core>

namespace tudat
{

namespace observation_models
{

//! Compact weight representation for one observation event.
/*!
 * Scalar and diagonal-vector weights are stored compactly and are expanded only
 * when flattened observation data with a matrix is materialized. Block weights store the full
 * observable-size matrix for one observation event.
 */
struct PerObservationWeight {
    enum class Type { scalar, diagonal_vector, block };

    PerObservationWeight( ) = default;

    explicit PerObservationWeight( const double scalarWeight ): type_( Type::scalar ), scalarWeight_( scalarWeight ) {}

    explicit PerObservationWeight( const Eigen::VectorXd& diagonalWeight ):
        type_( Type::diagonal_vector ), scalarWeight_( 1.0 ), diagonalWeight_( diagonalWeight )
    {}

    explicit PerObservationWeight( const Eigen::MatrixXd& blockWeight ):
        type_( Type::block ), scalarWeight_( 1.0 ), blockWeight_( blockWeight )
    {}

    //! Representation type used for this observation.
    Type type_ = Type::scalar;
    //! Compact scalar weight used when type_ is scalar.
    double scalarWeight_ = 1.0;
    //! Compact component-wise diagonal weights used when type_ is diagonal_vector.
    Eigen::VectorXd diagonalWeight_;
    //! Full observable-size block used when type_ is block.
    Eigen::MatrixXd blockWeight_;

    Eigen::MatrixXd toMatrix( const int observableSize ) const
    {
        if( type_ == Type::scalar )
        {
            return scalarWeight_ * Eigen::MatrixXd::Identity( observableSize, observableSize );
        }
        if( type_ == Type::diagonal_vector )
        {
            if( diagonalWeight_.size( ) != observableSize )
            {
                throw std::runtime_error( "Error when materializing observation weight, diagonal-vector size is inconsistent." );
            }
            return diagonalWeight_.asDiagonal( );
        }
        if( blockWeight_.rows( ) != observableSize || blockWeight_.cols( ) != observableSize )
        {
            throw std::runtime_error( "Error when materializing observation weight, block size is inconsistent." );
        }
        return blockWeight_;
    }

    Eigen::VectorXd toDiagonalVector( const int observableSize ) const
    {
        if( type_ == Type::scalar )
        {
            return scalarWeight_ * Eigen::VectorXd::Ones( observableSize );
        }
        if( type_ == Type::diagonal_vector )
        {
            if( diagonalWeight_.size( ) != observableSize )
            {
                throw std::runtime_error( "Error when retrieving observation weight vector, diagonal-vector size is inconsistent." );
            }
            return diagonalWeight_;
        }
        return toMatrix( observableSize ).diagonal( );
    }

    bool isDiagonalOnly( const int observableSize ) const
    {
        if( type_ == Type::scalar )
        {
            return true;
        }
        if( type_ == Type::diagonal_vector )
        {
            if( diagonalWeight_.size( ) != observableSize )
            {
                throw std::runtime_error( "Error when checking observation weight, diagonal-vector size is inconsistent." );
            }
            return true;
        }
        return isMatrixDiagonal( toMatrix( observableSize ) );
    }

    static bool isMatrixDiagonal( const Eigen::MatrixXd& matrix )
    {
        for( int row = 0; row < matrix.rows( ); ++row )
        {
            for( int column = 0; column < matrix.cols( ); ++column )
            {
                if( row != column && matrix( row, column ) != 0.0 )
                {
                    return false;
                }
            }
        }
        return true;
    }
};

//! Optional larger weight block tied to selected scalar components.
/*!
 * This is the internal advanced layer for rare correlations that do not fit
 * per-observation or set-level diagonal blocks. The initial public API can be
 * kept small while preserving a representation for arbitrary future blocks.
 */
struct ObservationWeightBlock {
    //! Scalar component ids covered by the block rows.
    std::vector< unsigned int > rowScalarComponentIds_;
    //! Scalar component ids covered by the block columns.
    std::vector< unsigned int > columnScalarComponentIds_;
    //! Dense block value for the selected scalar components.
    Eigen::MatrixXd weightBlock_;
};

//! Storage for all observation weights in an ObservationDataset.
/*!
 * The common path stores one compact PerObservationWeight per observation row.
 * A set-level block stores the full M x M block for a newly added batch/set and
 * is materialized only while flattened observation data is assembled. Extra off-diagonal blocks are kept as
 * an internal extension point for larger cross-observation correlations.
 */
class ObservationWeights
{
public:
    //! Append one scalar weight for a newly inserted observation row.
    void appendScalarWeight( const double scalarWeight )
    {
        perObservationWeights_.push_back( PerObservationWeight( scalarWeight ) );
    }

    //! Append one diagonal component-wise weight vector for a newly inserted observation row.
    void appendDiagonalWeightVector( const Eigen::VectorXd& diagonalWeight )
    {
        perObservationWeights_.push_back( PerObservationWeight( diagonalWeight ) );
    }

    //! Append one observable-size dense weight block for a newly inserted observation row.
    void appendWeightBlock( const Eigen::MatrixXd& blockWeight )
    {
        if( blockWeight.rows( ) != blockWeight.cols( ) )
        {
            throw std::runtime_error( "Error when adding observation weight block, block is not square." );
        }
        perObservationWeights_.push_back( PerObservationWeight( blockWeight ) );
    }

    //! Replace the weight of an existing observation row by one scalar value.
    void setScalarWeight( const std::size_t observationId, const double scalarWeight )
    {
        perObservationWeights_.at( observationId ) = PerObservationWeight( scalarWeight );
    }

    //! Replace the weight of an existing observation row by a diagonal component-wise vector.
    void setDiagonalWeightVector( const std::size_t observationId, const Eigen::VectorXd& diagonalWeight )
    {
        perObservationWeights_.at( observationId ) = PerObservationWeight( diagonalWeight );
    }

    //! Replace the weight of an existing observation row by a dense observable-size block.
    void setWeightBlock( const std::size_t observationId, const Eigen::MatrixXd& blockWeight )
    {
        if( blockWeight.rows( ) != blockWeight.cols( ) )
        {
            throw std::runtime_error( "Error when setting observation weight block, block is not square." );
        }
        perObservationWeights_.at( observationId ) = PerObservationWeight( blockWeight );
    }

    const PerObservationWeight& getObservationWeight( const std::size_t observationId ) const
    {
        return perObservationWeights_.at( observationId );
    }

    //! Append one existing compact per-observation representation without materializing it.
    void appendObservationWeight( const PerObservationWeight& observationWeight )
    {
        perObservationWeights_.push_back( observationWeight );
    }

    //! Replace one existing compact per-observation representation without materializing it.
    void setObservationWeight( const std::size_t observationId, const PerObservationWeight& observationWeight )
    {
        perObservationWeights_.at( observationId ) = observationWeight;
    }

    //! Return whether the row stores a dense block instead of compact scalar/diagonal weights.
    bool hasObservationWeightBlock( const std::size_t observationId ) const
    {
        return perObservationWeights_.at( observationId ).type_ == PerObservationWeight::Type::block;
    }

    //! Return the row weight as a dense observable-size matrix, expanding compact diagonal storage if needed.
    Eigen::MatrixXd getObservationWeightMatrix( const std::size_t observationId, const int observableSize ) const
    {
        return perObservationWeights_.at( observationId ).toMatrix( observableSize );
    }

    //! Return the diagonal entries of the row weight, regardless of compact or dense storage.
    Eigen::VectorXd getObservationWeightVector( const std::size_t observationId, const int observableSize ) const
    {
        return perObservationWeights_.at( observationId ).toDiagonalVector( observableSize );
    }

    //! Return whether the row weight has no off-diagonal entries after expansion.
    bool isObservationWeightDiagonalOnly( const std::size_t observationId, const int observableSize ) const
    {
        return perObservationWeights_.at( observationId ).isDiagonalOnly( observableSize );
    }

    std::size_t getNumberOfObservationWeights( ) const
    {
        return perObservationWeights_.size( );
    }

    //! Store a full set-level block that replaces per-row weights in flattened data for this set.
    void setSetWeightBlock( const std::size_t setId, const Eigen::MatrixXd& setWeightBlock )
    {
        if( setWeightBlock.rows( ) != setWeightBlock.cols( ) )
        {
            throw std::runtime_error( "Error when setting observation set weight block, block is not square." );
        }
        if( setWeightBlocks_.size( ) <= setId )
        {
            setWeightBlocks_.resize( setId + 1 );
        }
        setWeightBlocks_.at( setId ) = setWeightBlock;
    }

    //! Return whether a full set-level block is stored for the requested set.
    bool hasSetWeightBlock( const std::size_t setId ) const
    {
        return setId < setWeightBlocks_.size( ) && setWeightBlocks_.at( setId ).has_value( );
    }

    //! Return the full set-level block for the requested set.
    const Eigen::MatrixXd& getSetWeightBlock( const std::size_t setId ) const
    {
        if( !hasSetWeightBlock( setId ) )
        {
            throw std::runtime_error( "Error when retrieving observation set weight block, no block is stored for the requested set." );
        }
        return setWeightBlocks_.at( setId ).value( );
    }

    //! Return whether the stored set-level block has no off-diagonal entries.
    bool isSetWeightBlockDiagonalOnly( const std::size_t setId ) const
    {
        return PerObservationWeight::isMatrixDiagonal( getSetWeightBlock( setId ) );
    }

    //! Add an arbitrary block over scalar-component ids for cross-observation correlations.
    void addExtraWeightBlock( const ObservationWeightBlock& weightBlock )
    {
        if( weightBlock.weightBlock_.rows( ) != static_cast< int >( weightBlock.rowScalarComponentIds_.size( ) ) ||
            weightBlock.weightBlock_.cols( ) != static_cast< int >( weightBlock.columnScalarComponentIds_.size( ) ) )
        {
            throw std::runtime_error( "Error when adding extra observation weight block, block dimensions are inconsistent." );
        }
        extraWeightBlocks_.push_back( weightBlock );
    }

    const std::vector< ObservationWeightBlock >& getExtraWeightBlocks( ) const
    {
        return extraWeightBlocks_;
    }

    bool hasExtraWeightBlocks( ) const
    {
        return !extraWeightBlocks_.empty( );
    }

private:
    //! Per-observation compact weights, aligned one-to-one with observation rows.
    std::vector< PerObservationWeight > perObservationWeights_;
    //! Optional full M x M blocks for complete observation sets/batches.
    std::vector< std::optional< Eigen::MatrixXd > > setWeightBlocks_;
    //! Optional arbitrary blocks for rare off-diagonal correlations.
    std::vector< ObservationWeightBlock > extraWeightBlocks_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATION_WEIGHTS_H
