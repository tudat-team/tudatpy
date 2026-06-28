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
#include <map>
#include <stdexcept>
#include <vector>

#include <Eigen/Core>

namespace tudat
{

namespace observation_models
{

//! Compact weight representation for one observation event.
/*!
 * Scalar weights are stored as scalar values and are expanded only when a
 * vector or matrix projection is materialized. Matrix weights store the full
 * observable-size block for one observation event.
 */
struct PerObservationWeight {
    enum class Type { scalar, block };

    PerObservationWeight( ) = default;

    explicit PerObservationWeight( const double scalarWeight ): type_( Type::scalar ), scalarWeight_( scalarWeight ) {}

    explicit PerObservationWeight( const Eigen::MatrixXd& blockWeight ):
        type_( Type::block ), scalarWeight_( 1.0 ), blockWeight_( blockWeight )
    {}

    //! Representation type used for this observation.
    Type type_ = Type::scalar;
    //! Compact scalar weight used when type_ is scalar.
    double scalarWeight_ = 1.0;
    //! Full observable-size block used when type_ is block.
    Eigen::MatrixXd blockWeight_;

    Eigen::MatrixXd toMatrix( const int observableSize ) const
    {
        if( type_ == Type::scalar )
        {
            return scalarWeight_ * Eigen::MatrixXd::Identity( observableSize, observableSize );
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
        return toMatrix( observableSize ).diagonal( );
    }

    bool isDiagonalOnly( const int observableSize ) const
    {
        if( type_ == Type::scalar )
        {
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
    std::vector< std::size_t > rowScalarComponentIds_;
    //! Scalar component ids covered by the block columns.
    std::vector< std::size_t > columnScalarComponentIds_;
    //! Dense block value for the selected scalar components.
    Eigen::MatrixXd weightBlock_;
};

//! Storage for all observation weights in an ObservationDataset.
/*!
 * The common path stores one compact PerObservationWeight per observation row.
 * A set-level block stores the full M x M block for a newly added batch/set and
 * is materialized only during projection. Extra off-diagonal blocks are kept as
 * an internal extension point for larger cross-observation correlations.
 */
class ObservationWeights
{
public:
    void appendScalarWeight( const double scalarWeight )
    {
        perObservationWeights_.push_back( PerObservationWeight( scalarWeight ) );
    }

    void appendDiagonalWeightVector( const Eigen::VectorXd& diagonalWeight )
    {
        perObservationWeights_.push_back( PerObservationWeight( diagonalWeight.asDiagonal( ).toDenseMatrix( ) ) );
    }

    void appendWeightBlock( const Eigen::MatrixXd& blockWeight )
    {
        if( blockWeight.rows( ) != blockWeight.cols( ) )
        {
            throw std::runtime_error( "Error when adding observation weight block, block is not square." );
        }
        perObservationWeights_.push_back( PerObservationWeight( blockWeight ) );
    }

    void setScalarWeight( const std::size_t observationId, const double scalarWeight )
    {
        perObservationWeights_.at( observationId ) = PerObservationWeight( scalarWeight );
    }

    void setDiagonalWeightVector( const std::size_t observationId, const Eigen::VectorXd& diagonalWeight )
    {
        perObservationWeights_.at( observationId ) = PerObservationWeight( diagonalWeight.asDiagonal( ).toDenseMatrix( ) );
    }

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

    bool hasObservationWeightBlock( const std::size_t observationId ) const
    {
        return perObservationWeights_.at( observationId ).type_ == PerObservationWeight::Type::block;
    }

    Eigen::MatrixXd getObservationWeightMatrix( const std::size_t observationId, const int observableSize ) const
    {
        return perObservationWeights_.at( observationId ).toMatrix( observableSize );
    }

    Eigen::VectorXd getObservationWeightVector( const std::size_t observationId, const int observableSize ) const
    {
        return perObservationWeights_.at( observationId ).toDiagonalVector( observableSize );
    }

    bool isObservationWeightDiagonalOnly( const std::size_t observationId, const int observableSize ) const
    {
        return perObservationWeights_.at( observationId ).isDiagonalOnly( observableSize );
    }

    std::size_t getNumberOfObservationWeights( ) const
    {
        return perObservationWeights_.size( );
    }

    void setSetWeightBlock( const std::size_t setId, const Eigen::MatrixXd& setWeightBlock )
    {
        if( setWeightBlock.rows( ) != setWeightBlock.cols( ) )
        {
            throw std::runtime_error( "Error when setting observation set weight block, block is not square." );
        }
        setWeightBlocks_[ setId ] = setWeightBlock;
    }

    bool hasSetWeightBlock( const std::size_t setId ) const
    {
        return setWeightBlocks_.count( setId ) > 0;
    }

    const Eigen::MatrixXd& getSetWeightBlock( const std::size_t setId ) const
    {
        return setWeightBlocks_.at( setId );
    }

    bool isSetWeightBlockDiagonalOnly( const std::size_t setId ) const
    {
        return PerObservationWeight::isMatrixDiagonal( setWeightBlocks_.at( setId ) );
    }

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
    std::map< std::size_t, Eigen::MatrixXd > setWeightBlocks_;
    //! Optional arbitrary blocks for rare off-diagonal correlations.
    std::vector< ObservationWeightBlock > extraWeightBlocks_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATION_WEIGHTS_H
