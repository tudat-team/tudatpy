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

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <cereal/types/map.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/vector.hpp>

#include "tudat/io/observationWeightSettings.h"

namespace tudat
{
namespace observation_models
{

//! One effective symmetric weight matrix, indexed by dataset scalar storage.
/*!
 * Diagonal entries occupy one vector. Only nonzero off-diagonal entries are
 * stored, once each in the upper triangle. Small observation blocks and sparse
 * cross-observation blocks both update these same coefficients. There are no
 * precedence layers: an assignment replaces precisely the addressed entries.
 * Restriction is a principal submatrix, with the requested index order.
 */
class ObservationWeights
{
public:
    using Index = unsigned int;
    using Entry = std::pair< Index, Index >;

    //! Return the number of scalar rows represented by this weight matrix.
    std::size_t size( ) const
    {
        return diagonal_.size( );
    }

    //! Return whether at least one nonzero off-diagonal coefficient is stored.
    bool hasOffDiagonalWeights( ) const
    {
        return !offDiagonal_.empty( );
    }

    //! Nonzero upper-triangular coefficients, keyed by their two dataset scalar-storage indices.
    const std::map< Entry, double >& getOffDiagonalEntries( ) const
    {
        return offDiagonal_;
    }

    //! Validate that diagonal weights are finite and nonnegative.
    static void validateDiagonal( const Eigen::VectorXd& diagonal )
    {
        if( !diagonal.allFinite( ) || ( diagonal.array( ) < 0.0 ).any( ) )
        {
            throw std::runtime_error( "Observation weight diagonals must be finite and nonnegative." );
        }
    }

    //! Append validated diagonal coefficients to scalar storage.
    void appendDiagonal( const Eigen::VectorXd& diagonal )
    {
        validateDiagonal( diagonal );
        if( diagonal.size( ) > 0 )
        {
            diagonal_.insert( diagonal_.end( ), diagonal.data( ), diagonal.data( ) + diagonal.size( ) );
        }
    }

    //! Gather diagonal coefficients in the requested scalar-index order.
    Eigen::VectorXd getDiagonal( const std::vector< Index >& indices ) const
    {
        Eigen::VectorXd diagonal( indices.size( ) );
        for( std::size_t i = 0; i < indices.size( ); ++i )
        {
            diagonal( i ) = diagonal_.at( indices.at( i ) );
        }
        return diagonal;
    }

    //! Replace the selected principal block by a diagonal, without allocating a dense block.
    void setDiagonal( const std::vector< Index >& indices, const Eigen::VectorXd& diagonal )
    {
        const auto selected = indexMap( indices );
        if( indices.size( ) != static_cast< std::size_t >( diagonal.size( ) ) )
        {
            throw std::runtime_error( "Observation weight diagonal has an inconsistent size." );
        }
        validateDiagonal( diagonal );
        for( auto entry = offDiagonal_.begin( ); entry != offDiagonal_.end( ); )
        {
            if( selected.count( entry->first.first ) && selected.count( entry->first.second ) )
            {
                entry = offDiagonal_.erase( entry );
            }
            else
            {
                ++entry;
            }
        }
        for( std::size_t i = 0; i < indices.size( ); ++i )
        {
            diagonal_.at( indices.at( i ) ) = diagonal( i );
        }
    }

    //! Replace one contiguous observation block by its diagonal.
    void setObservationDiagonal( const Index first, const Eigen::VectorXd& diagonal )
    {
        validateDiagonal( diagonal );
        if( first > size( ) || static_cast< std::size_t >( diagonal.size( ) ) > size( ) - first )
        {
            throw std::runtime_error( "Observation weight diagonal is outside scalar storage." );
        }
        const Index end = first + diagonal.size( );
        // Ordered sparse keys let a batch of per-observation updates visit each
        // relevant coefficient once, instead of rescanning the whole matrix per row.
        auto entry = offDiagonal_.lower_bound( { first, first } );
        while( entry != offDiagonal_.end( ) && entry->first.first < end )
        {
            if( entry->first.second < end )
            {
                entry = offDiagonal_.erase( entry );
            }
            else
            {
                ++entry;
            }
        }
        if( diagonal.size( ) > 0 )
        {
            std::copy( diagonal.data( ), diagonal.data( ) + diagonal.size( ), diagonal_.begin( ) + first );
        }
    }

    //! Assign a block and its transpose; overlapping selections must agree.
    void setBlock( const std::vector< Index >& rows, const std::vector< Index >& columns, const Eigen::MatrixXd& block )
    {
        indexMap( rows );
        indexMap( columns );
        if( block.rows( ) != static_cast< Eigen::Index >( rows.size( ) ) ||
            block.cols( ) != static_cast< Eigen::Index >( columns.size( ) ) || !block.allFinite( ) )
        {
            throw std::runtime_error( "Observation weight block must have matching dimensions and finite entries." );
        }
        std::map< Entry, double > assignments;
        for( std::size_t i = 0; i < rows.size( ); ++i )
        {
            for( std::size_t j = 0; j < columns.size( ); ++j )
            {
                const Entry key = std::minmax( rows.at( i ), columns.at( j ) );
                const double value = block( i, j );
                if( key.first == key.second && value < 0.0 )
                {
                    throw std::runtime_error( "Observation weight diagonals must be nonnegative." );
                }
                const auto inserted = assignments.emplace( key, value );
                if( !inserted.second )
                {
                    const double previous = inserted.first->second;
                    if( std::abs( previous - value ) > 1.0E-12 * std::max( { 1.0, std::abs( previous ), std::abs( value ) } ) )
                    {
                        throw std::runtime_error( "Overlapping observation weight entries must be symmetric." );
                    }
                    inserted.first->second = 0.5 * previous + 0.5 * value;
                }
            }
        }
        // All validation precedes mutation, including duplicate/permuted selectors.
        for( const auto& assignment : assignments )
        {
            setEntry( assignment.first, assignment.second );
        }
    }

    //! Select, reorder or remove scalar rows and columns together.
    ObservationWeights restricted( const std::vector< Index >& indices ) const
    {
        const auto selected = indexMap( indices );
        ObservationWeights result;
        result.diagonal_.reserve( indices.size( ) );
        for( const Index index : indices )
        {
            result.diagonal_.push_back( diagonal_.at( index ) );
        }
        for( const auto& entry : offDiagonal_ )
        {
            const auto row = selected.find( entry.first.first );
            const auto column = selected.find( entry.first.second );
            if( row != selected.end( ) && column != selected.end( ) )
            {
                result.offDiagonal_.emplace( std::minmax( row->second, column->second ), entry.second );
            }
        }
        return result;
    }

    //! Copy an already restricted matrix into selected target scalar entries.
    void copyBlock( const ObservationWeights& source, const std::vector< Index >& targetIndices )
    {
        if( source.size( ) != targetIndices.size( ) )
        {
            throw std::runtime_error( "Copied observation weights have an inconsistent scalar size." );
        }
        // A value snapshot also makes copying within the same object well defined.
        const ObservationWeights snapshot( source );
        setDiagonal( targetIndices, snapshot.diagonalVector( ) );
        for( const auto& entry : snapshot.offDiagonal_ )
        {
            setEntry( std::minmax( targetIndices.at( entry.first.first ), targetIndices.at( entry.first.second ) ), entry.second );
        }
    }

    //! Return all diagonal coefficients as an Eigen vector.
    Eigen::VectorXd diagonalVector( ) const
    {
        Eigen::VectorXd result( diagonal_.size( ) );
        std::copy( diagonal_.begin( ), diagonal_.end( ), result.data( ) );
        return result;
    }

    //! Materialize the complete symmetric sparse weight matrix.
    Eigen::SparseMatrix< double > sparseMatrix( ) const
    {
        Eigen::SparseMatrix< double > result( size( ), size( ) );
        std::vector< Eigen::Triplet< double > > entries;
        entries.reserve( diagonal_.size( ) + 2 * offDiagonal_.size( ) );
        for( std::size_t i = 0; i < diagonal_.size( ); ++i )
        {
            if( diagonal_.at( i ) != 0.0 )
            {
                entries.emplace_back( i, i, diagonal_.at( i ) );
            }
        }
        for( const auto& entry : offDiagonal_ )
        {
            entries.emplace_back( entry.first.first, entry.first.second, entry.second );
            entries.emplace_back( entry.first.second, entry.first.first, entry.second );
        }
        result.setFromTriplets( entries.begin( ), entries.end( ) );
        return result;
    }

    //! Compare the effective compact weight storage for exact equality.
    bool operator==( const ObservationWeights& other ) const
    {
        return diagonal_ == other.diagonal_ && offDiagonal_ == other.offDiagonal_;
    }

    //! Serialize compact diagonal and upper-triangular weight storage.
    template< class Archive >
    void save( Archive& archive ) const
    {
        archive( diagonal_, offDiagonal_ );
    }

    //! Deserialize and validate compact observation-weight storage.
    template< class Archive >
    void load( Archive& archive )
    {
        ObservationWeights loaded;
        archive( loaded.diagonal_, loaded.offDiagonal_ );
        validateDiagonal( loaded.diagonalVector( ) );
        for( const auto& entry : loaded.offDiagonal_ )
        {
            if( entry.first.first >= entry.first.second || entry.first.second >= loaded.size( ) || !std::isfinite( entry.second ) )
            {
                throw std::runtime_error( "Serialized observation weights have invalid sparse entries." );
            }
        }
        *this = std::move( loaded );
    }

private:
    //! Validate scalar indices and map each selected source index to its output position.
    std::unordered_map< Index, Index > indexMap( const std::vector< Index >& indices ) const
    {
        std::unordered_map< Index, Index > result;
        result.reserve( indices.size( ) );
        for( std::size_t i = 0; i < indices.size( ); ++i )
        {
            if( indices.at( i ) >= size( ) || !result.emplace( indices.at( i ), i ).second )
            {
                throw std::runtime_error( "Observation weight selectors must contain unique, existing scalar indices." );
            }
        }
        return result;
    }

    //! Assign one canonical upper-triangular entry in compact storage.
    void setEntry( const Entry& key, const double value )
    {
        if( key.first == key.second )
        {
            diagonal_.at( key.first ) = value;
        }
        else if( value == 0.0 )
        {
            offDiagonal_.erase( key );
        }
        else
        {
            offDiagonal_[ key ] = value;
        }
    }

    //! Diagonal coefficients in dataset scalar-storage order.
    std::vector< double > diagonal_;
    //! Nonzero off-diagonal coefficients stored once in the upper triangle.
    std::map< Entry, double > offDiagonal_;
};

//! Create the effective weights for a new observation set.
/*!
 * \param numberOfObservations Number of observation events in the set.
 * \param singleObservationSize Number of scalar components in each observation event. For example, this is two for
 * angular-position observations.
 * \param weightSettings Representation and numerical weight values to apply. For scalar_per_observation,
 * scalarWeights_ must contain numberOfObservations entries. For diagonal_per_observation, diagonalWeights_ must contain
 * numberOfObservations vectors of length singleObservationSize. For constant_block, weightBlock_ must be a
 * singleObservationSize-by-singleObservationSize matrix. For block_per_observation, weightBlocks_ must contain
 * numberOfObservations matrices of that size. For set_block, weightBlock_ must have
 * numberOfObservations * singleObservationSize rows and columns.
 * \return Validated weights ordered by observation event and then by scalar component within each event.
 */
inline ObservationWeights createObservationWeightsForSet( const std::size_t numberOfObservations,
                                                          const unsigned int singleObservationSize,
                                                          const ObservationWeightSettings& weightSettings )
{
    using WeightsBlockType = ObservationWeightSettings::WeightsBlockType;
    if( singleObservationSize == 0 ||
        numberOfObservations > std::numeric_limits< ObservationWeights::Index >::max( ) / singleObservationSize )
    {
        throw std::runtime_error( "Observation weight dimensions exceed scalar storage capacity." );
    }

    ObservationWeights result;
    Eigen::VectorXd diagonal = Eigen::VectorXd::Ones( numberOfObservations * singleObservationSize );
    switch( weightSettings.type_ )
    {
        case WeightsBlockType::constant_scalar:
            ObservationWeights::validateDiagonal( Eigen::VectorXd::Constant( 1, weightSettings.scalarWeight_ ) );
            diagonal.setConstant( weightSettings.scalarWeight_ );
            break;
        case WeightsBlockType::scalar_per_observation:
            if( weightSettings.scalarWeights_.size( ) != numberOfObservations )
            {
                throw std::runtime_error( "Observation scalar weight count is inconsistent." );
            }
            for( std::size_t i = 0; i < numberOfObservations; ++i )
            {
                diagonal.segment( i * singleObservationSize, singleObservationSize ).setConstant( weightSettings.scalarWeights_.at( i ) );
            }
            break;
        case WeightsBlockType::diagonal_per_observation:
            if( weightSettings.diagonalWeights_.size( ) != numberOfObservations )
            {
                throw std::runtime_error( "Observation diagonal weight count is inconsistent." );
            }
            for( std::size_t i = 0; i < numberOfObservations; ++i )
            {
                const Eigen::VectorXd& observationDiagonal = weightSettings.diagonalWeights_.at( i );
                if( observationDiagonal.size( ) != singleObservationSize )
                {
                    throw std::runtime_error( "Observation diagonal weight size is inconsistent." );
                }
                ObservationWeights::validateDiagonal( observationDiagonal );
                diagonal.segment( i * singleObservationSize, singleObservationSize ) = observationDiagonal;
            }
            break;
        case WeightsBlockType::default_weights:
        case WeightsBlockType::constant_block:
        case WeightsBlockType::block_per_observation:
        case WeightsBlockType::set_block:
            break;
        default:
            throw std::runtime_error( "Unknown observation weight policy." );
    }

    result.appendDiagonal( diagonal );
    if( weightSettings.type_ == WeightsBlockType::set_block )
    {
        std::vector< ObservationWeights::Index > indices( numberOfObservations * singleObservationSize );
        std::iota( indices.begin( ), indices.end( ), 0 );
        result.setBlock( indices, indices, weightSettings.weightBlock_ );
    }
    else if( weightSettings.type_ == WeightsBlockType::constant_block || weightSettings.type_ == WeightsBlockType::block_per_observation )
    {
        if( weightSettings.type_ == WeightsBlockType::block_per_observation &&
            weightSettings.weightBlocks_.size( ) != numberOfObservations )
        {
            throw std::runtime_error( "Observation weight block count is inconsistent." );
        }

        std::vector< ObservationWeights::Index > indices( singleObservationSize );
        if( weightSettings.type_ == WeightsBlockType::constant_block )
        {
            // Validate a constant block even for an empty observation set.
            ObservationWeights singleObservationWeights;
            singleObservationWeights.appendDiagonal( Eigen::VectorXd::Ones( singleObservationSize ) );
            std::iota( indices.begin( ), indices.end( ), 0 );
            singleObservationWeights.setBlock( indices, indices, weightSettings.weightBlock_ );
        }
        for( std::size_t i = 0; i < numberOfObservations; ++i )
        {
            std::iota( indices.begin( ), indices.end( ), i * singleObservationSize );
            result.setBlock( indices,
                             indices,
                             weightSettings.type_ == WeightsBlockType::constant_block ? weightSettings.weightBlock_
                                                                                      : weightSettings.weightBlocks_.at( i ) );
        }
    }
    return result;
}

}  // namespace observation_models
}  // namespace tudat

#endif  // TUDAT_OBSERVATION_WEIGHTS_H
