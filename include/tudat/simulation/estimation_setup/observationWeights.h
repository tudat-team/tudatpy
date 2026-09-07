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

namespace tudat
{
namespace observation_models
{

//! Weight policy used while adding a new observation set.
/*!
 * This object keeps the add-observation-set interface small while still
 * supporting compact scalar weights, per-observation scalar weights,
 * observable-size per-observation blocks and full set-level blocks.
 */
struct ObservationWeightSettings {
    enum class Type { default_weights, constant_scalar, scalar_per_observation, constant_block, block_per_observation, set_block };

    static ObservationWeightSettings defaultWeights( )
    {
        return ObservationWeightSettings( );
    }

    static ObservationWeightSettings constantScalar( const double weight )
    {
        ObservationWeightSettings settings;
        settings.type_ = Type::constant_scalar;
        settings.scalarWeight_ = weight;
        return settings;
    }

    static ObservationWeightSettings scalarPerObservation( const std::vector< double >& weights )
    {
        ObservationWeightSettings settings;
        settings.type_ = Type::scalar_per_observation;
        settings.scalarWeights_ = weights;
        return settings;
    }

    static ObservationWeightSettings constantBlock( const Eigen::MatrixXd& weightBlock )
    {
        ObservationWeightSettings settings;
        settings.type_ = Type::constant_block;
        settings.weightBlock_ = weightBlock;
        return settings;
    }

    static ObservationWeightSettings blockPerObservation( const std::vector< Eigen::MatrixXd >& weightBlocks )
    {
        ObservationWeightSettings settings;
        settings.type_ = Type::block_per_observation;
        settings.weightBlocks_ = weightBlocks;
        return settings;
    }

    static ObservationWeightSettings setBlock( const Eigen::MatrixXd& weightBlock )
    {
        ObservationWeightSettings settings;
        settings.type_ = Type::set_block;
        settings.weightBlock_ = weightBlock;
        return settings;
    }

    Type type_ = Type::default_weights;
    double scalarWeight_ = 1.0;
    std::vector< double > scalarWeights_;
    Eigen::MatrixXd weightBlock_;
    std::vector< Eigen::MatrixXd > weightBlocks_;
};

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

    //! Validate and normalize an addition policy before changing the dataset.
    static ObservationWeights forSet( const std::size_t count, const unsigned int dimension,
                                      const ObservationWeightSettings& settings )
    {
        if( dimension == 0 || count > std::numeric_limits< Index >::max( ) / dimension )
        {
            throw std::runtime_error( "Observation weight dimensions exceed scalar storage capacity." );
        }
        using Type = ObservationWeightSettings::Type;
        ObservationWeights result;
        Eigen::VectorXd diagonal = Eigen::VectorXd::Ones( count * dimension );
        switch( settings.type_ )
        {
        case Type::constant_scalar:
            validateDiagonal( Eigen::VectorXd::Constant( 1, settings.scalarWeight_ ) );
            diagonal.setConstant( settings.scalarWeight_ );
            break;
        case Type::scalar_per_observation:
            if( settings.scalarWeights_.size( ) != count )
            {
                throw std::runtime_error( "Observation scalar weight count is inconsistent." );
            }
            for( std::size_t i = 0; i < count; ++i )
            {
                diagonal.segment( i * dimension, dimension ).setConstant( settings.scalarWeights_.at( i ) );
            }
            break;
        case Type::default_weights:
        case Type::constant_block:
        case Type::block_per_observation:
        case Type::set_block:
            break;
        default:
            throw std::runtime_error( "Unknown observation weight policy." );
        }
        result.appendDiagonal( diagonal );
        if( settings.type_ == Type::set_block )
        {
            std::vector< Index > indices( count * dimension );
            std::iota( indices.begin( ), indices.end( ), 0 );
            result.setBlock( indices, indices, settings.weightBlock_ );
        }
        else if( settings.type_ == Type::constant_block || settings.type_ == Type::block_per_observation )
        {
            if( settings.type_ == Type::block_per_observation && settings.weightBlocks_.size( ) != count )
            {
                throw std::runtime_error( "Observation weight block count is inconsistent." );
            }
            std::vector< Index > indices( dimension );
            if( settings.type_ == Type::constant_block )
            {
                // Validate a constant block even for an empty observation set.
                ObservationWeights single;
                single.appendDiagonal( Eigen::VectorXd::Ones( dimension ) );
                std::iota( indices.begin( ), indices.end( ), 0 );
                single.setBlock( indices, indices, settings.weightBlock_ );
            }
            for( std::size_t i = 0; i < count; ++i )
            {
                std::iota( indices.begin( ), indices.end( ), i * dimension );
                result.setBlock( indices, indices, settings.type_ == Type::constant_block ? settings.weightBlock_ : settings.weightBlocks_.at( i ) );
            }
        }
        return result;
    }

    std::size_t size( ) const
    {
        return diagonal_.size( );
    }

    bool hasOffDiagonalWeights( ) const
    {
        return !offDiagonal_.empty( );
    }

    static void validateDiagonal( const Eigen::VectorXd& diagonal )
    {
        if( !diagonal.allFinite( ) || ( diagonal.array( ) < 0.0 ).any( ) )
        {
            throw std::runtime_error( "Observation weight diagonals must be finite and nonnegative." );
        }
    }

    void appendDiagonal( const Eigen::VectorXd& diagonal )
    {
        validateDiagonal( diagonal );
        if( diagonal.size( ) > 0 )
        {
            diagonal_.insert( diagonal_.end( ), diagonal.data( ), diagonal.data( ) + diagonal.size( ) );
        }
    }

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
    void setBlock( const std::vector< Index >& rows,
                   const std::vector< Index >& columns,
                   const Eigen::MatrixXd& block )
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

    Eigen::VectorXd diagonalVector( ) const
    {
        Eigen::VectorXd result( diagonal_.size( ) );
        std::copy( diagonal_.begin( ), diagonal_.end( ), result.data( ) );
        return result;
    }

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

    bool operator==( const ObservationWeights& other ) const
    {
        return diagonal_ == other.diagonal_ && offDiagonal_ == other.offDiagonal_;
    }

    template< class Archive >
    void save( Archive& archive ) const
    {
        archive( diagonal_, offDiagonal_ );
    }

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

    std::vector< double > diagonal_;
    std::map< Entry, double > offDiagonal_;
};

}  // namespace observation_models
}  // namespace tudat

#endif  // TUDAT_OBSERVATION_WEIGHTS_H
