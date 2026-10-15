/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/astro/gravitation/harmonicGravityFieldVariations.h"
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <utility>

namespace tudat
{
namespace gravitation
{
namespace
{
//! Check canonical phase identities, finite matrix dimensions and the physical degree/order triangle.
void validateConstituents( const std::vector< HarmonicGravityConstituent >& constituents )
{
    if( constituents.empty( ) )
    {
        throw std::invalid_argument( "Empty harmonic gravity inventory" );
    }
    std::set< std::string > identifiers;

    for( const auto& wave : constituents )
    {
        if( wave.doodson.size( ) != 7 || input_output::parseDoodsonMultipliers( wave.doodson ) != wave.multipliers ||
            !identifiers.insert( wave.doodson ).second )
        {
            throw std::invalid_argument( "Invalid or duplicate canonical harmonic gravity identifier" );
        }
        if( wave.cCos.rows( ) == 0 || wave.cCos.cols( ) == 0 )
        {
            throw std::invalid_argument( "Empty harmonic gravity matrix" );
        }
        for( const auto* matrix : { &wave.cCos, &wave.cSin, &wave.sCos, &wave.sSin } )
        {
            if( matrix->rows( ) != wave.cCos.rows( ) || matrix->cols( ) != wave.cCos.cols( ) || !matrix->allFinite( ) )
            {
                throw std::invalid_argument( "Inconsistent or nonfinite harmonic gravity matrices" );
            }
        }
        // Canonical data must explicitly respect physical Sn0 and the degree/order triangle.
        if( !wave.sCos.col( 0 ).isZero( 0.0 ) || !wave.sSin.col( 0 ).isZero( 0.0 ) )
        {
            throw std::invalid_argument( "Canonical harmonic gravity Sn0 must be zero" );
        }
        for( int n = 0; n < wave.cCos.rows( ); ++n )
        {
            for( int m = n + 1; m < wave.cCos.cols( ); ++m )
            {
                for( const auto* matrix : { &wave.cCos, &wave.cSin, &wave.sCos, &wave.sSin } )
                {
                    if( ( *matrix )( n, m ) != 0.0 )
                    {
                        throw std::invalid_argument( "Nonzero harmonic gravity coefficient above degree/order triangle" );
                    }
                }
            }
        }
    }
}
//! Validate the selected matrix block before constructing the base correction dimensions.
int checkedMinimumDegree( const std::vector< HarmonicGravityConstituent >& waves,
                          int minimumDegree,
                          int minimumOrder,
                          int maximumDegree,
                          int maximumOrder )
{
    validateConstituents( waves );
    if( minimumDegree < 2 || maximumDegree < minimumDegree || minimumOrder < 0 || maximumOrder < minimumOrder ||
        maximumOrder > maximumDegree )
    {
        throw std::invalid_argument( "Invalid harmonic gravity degree/order selection" );
    }
    for( const auto& wave : waves )
    {
        if( maximumDegree >= wave.cCos.rows( ) || maximumOrder >= wave.cCos.cols( ) )
        {
            throw std::invalid_argument( "Harmonic gravity selection exceeds the explicit canonical matrices" );
        }
    }
    return minimumDegree;
}
}  // namespace
//! Convert the selected published Cp/Sp/Cm/Sm terms with IERS 6.15, retaining the requested constituent order.
std::vector< HarmonicGravityConstituent > convertHarmonicPotentialData( const input_output::HarmonicPotentialData& data,
                                                                        int maximumDegree,
                                                                        int maximumOrder,
                                                                        bool allowUnprofiledZeroFill,
                                                                        const std::vector< std::string >& selectedDoodsonIds )
{
    input_output::validateHarmonicPotentialData( data );
    if( maximumDegree < 2 || maximumOrder < 0 || maximumOrder > maximumDegree )
    {
        throw std::invalid_argument( "Invalid harmonic potential selection" );
    }
    std::vector< HarmonicGravityConstituent > result;

    const bool identified = data.profile != input_output::HarmonicPotentialProfile::iersPotential;

    std::vector< const input_output::HarmonicPotentialConstituent* > selected;

    if( selectedDoodsonIds.empty( ) )
    {
        for( const auto& source : data.constituents )
        {
            selected.push_back( &source );
        }
    }
    else
    {
        std::set< std::string > ids;

        for( const auto& id : selectedDoodsonIds )
        {
            const auto found = std::find_if(
                    data.constituents.begin( ), data.constituents.end( ), [ & ]( const auto& source ) { return source.doodson == id; } );
            if( !ids.insert( id ).second || found == data.constituents.end( ) )
            {
                throw std::invalid_argument( "Unknown or duplicate selected harmonic potential Doodson identifier" );
            }
            selected.push_back( &*found );
        }
    }
    for( const auto* sourcePointer : selected )
    {
        const auto& source = *sourcePointer;

        HarmonicGravityConstituent wave;

        wave.doodson = source.doodson;
        wave.multipliers = source.multipliers;
        wave.cCos = wave.cSin = wave.sCos = wave.sSin = Eigen::MatrixXd::Zero( maximumDegree + 1, maximumOrder + 1 );
        for( int n = 2; n <= maximumDegree; ++n )
        {
            for( int m = 0; m <= std::min( n, maximumOrder ); ++m )
            {
                const auto found = source.coefficients.find( { n, m } );

                if( found == source.coefficients.end( ) )
                {
                    if( identified && ( n < source.minimumDegree || n > source.maximumDegree || m > source.maximumOrder ) )
                    {
                        continue;
                    }
                    if( !identified && allowUnprofiledZeroFill )
                    {
                        continue;
                    }
                    throw std::invalid_argument( "Missing harmonic potential term inside selected wave extent: " + source.doodson + "/" +
                                                 std::to_string( n ) + "/" + std::to_string( m ) );
                }
                // IERS 6.15 fixes the sine signs; zonal sine coefficients remain identically zero.
                const auto& c = found->second.values;

                wave.cCos( n, m ) = data.stokesUnitScale * ( c[ 0 ] + c[ 2 ] );
                wave.cSin( n, m ) = data.stokesUnitScale * ( c[ 1 ] + c[ 3 ] );
                if( m != 0 )
                {
                    wave.sCos( n, m ) = data.stokesUnitScale * ( c[ 1 ] - c[ 3 ] );
                    wave.sSin( n, m ) = data.stokesUnitScale * ( c[ 2 ] - c[ 0 ] );
                }
            }
        }
        result.push_back( std::move( wave ) );
    }
    validateConstituents( result );
    return result;
}
//! Apply the explicit major-to-target map once to each canonical quadrature.
std::vector< HarmonicGravityConstituent > applyHarmonicGravityAdmittance( const std::vector< HarmonicGravityConstituent >& majors,
                                                                          const HarmonicGravityAdmittance& admittance )
{
    validateConstituents( majors );
    if( admittance.majorDoodsonIds.size( ) != majors.size( ) || admittance.targetDoodsonIds.empty( ) ||
        admittance.matrix.cols( ) != static_cast< Eigen::Index >( majors.size( ) ) ||
        admittance.matrix.rows( ) != static_cast< Eigen::Index >( admittance.targetDoodsonIds.size( ) ) || !admittance.matrix.allFinite( ) )
    {
        throw std::invalid_argument( "Invalid ordered harmonic gravity admittance dimensions" );
    }
    for( std::size_t i = 0; i < majors.size( ); ++i )
    {
        if( majors[ i ].doodson != admittance.majorDoodsonIds[ i ] || majors[ i ].cCos.rows( ) != majors.front( ).cCos.rows( ) ||
            majors[ i ].cCos.cols( ) != majors.front( ).cCos.cols( ) )
        {
            throw std::invalid_argument( "Admittance majors do not match their exact ordered canonical inventory" );
        }
    }
    std::set< std::string > identifiers;

    std::vector< HarmonicGravityConstituent > result;

    for( std::size_t i = 0; i < admittance.targetDoodsonIds.size( ); ++i )
    {
        HarmonicGravityConstituent target;

        target.doodson = admittance.targetDoodsonIds[ i ];
        if( target.doodson.size( ) != 7 || !identifiers.insert( target.doodson ).second )
        {
            throw std::invalid_argument( "Admittance targets must have unique canonical Doodson identifiers" );
        }
        target.multipliers = input_output::parseDoodsonMultipliers( target.doodson );
        target.cCos = target.cSin = target.sCos = target.sSin =
                Eigen::MatrixXd::Zero( majors.front( ).cCos.rows( ), majors.front( ).cCos.cols( ) );
        for( std::size_t j = 0; j < majors.size( ); ++j )
        {
            const double scale = admittance.matrix( i, j );

            // Published maps are sparse; avoid dense coefficient operations for absent input weights.
            if( scale == 0.0 )
            {
                continue;
            }
            target.cCos += scale * majors[ j ].cCos;
            target.cSin += scale * majors[ j ].cSin;
            target.sCos += scale * majors[ j ].sCos;
            target.sSin += scale * majors[ j ].sSin;
        }
        result.push_back( std::move( target ) );
    }
    validateConstituents( result );
    return result;
}
//! Retain the selected block and convert each degree from the source reference constants to the consuming field.
HarmonicGravityFieldVariations::HarmonicGravityFieldVariations( const std::vector< HarmonicGravityConstituent >& constituents,
                                                                const HarmonicGravityMetadata& metadata,
                                                                const std::function< Eigen::Vector6d( double ) >& doodsonArguments,
                                                                double outputGravitationalParameter,
                                                                double outputReferenceRadius,
                                                                int minimumDegree,
                                                                int minimumOrder,
                                                                int maximumDegree,
                                                                int maximumOrder ):
    GravityFieldVariations( checkedMinimumDegree( constituents, minimumDegree, minimumOrder, maximumDegree, maximumOrder ),
                            minimumOrder,
                            maximumDegree,
                            maximumOrder ),
    constituents_( constituents ), doodsonArguments_( doodsonArguments )
{
    if( !doodsonArguments_ ||
        ( metadata.component != HarmonicGravityComponent::oceanMass &&
          metadata.component != HarmonicGravityComponent::directAtmosphericMass ) ||
        !metadata.fullyNormalized || metadata.revision.empty( ) ||
        ( metadata.loading != HarmonicGravityLoading::included && metadata.loading != HarmonicGravityLoading::excluded ) )
    {
        throw std::invalid_argument( "Harmonic gravity requires explicit component, normalization, revision, loading and clock callback" );
    }
    for( double value : { metadata.gravitationalParameter, metadata.referenceRadius, outputGravitationalParameter, outputReferenceRadius } )
    {
        if( !std::isfinite( value ) || value <= 0.0 )
        {
            throw std::invalid_argument( "Harmonic gravity reference constants must be finite and positive" );
        }
    }
    for( auto& wave : constituents_ )
    {
        for( auto* matrix : { &wave.cCos, &wave.cSin, &wave.sCos, &wave.sSin } )
        {
            *matrix = matrix->block( minimumDegree, minimumOrder, numberOfDegrees_, numberOfOrders_ ).eval( );
            for( int n = minimumDegree; n <= maximumDegree; ++n )
            {
                const double scale = metadata.gravitationalParameter / outputGravitationalParameter *
                        std::pow( metadata.referenceRadius / outputReferenceRadius, n );
                if( !std::isfinite( scale ) )
                {
                    throw std::invalid_argument( "Nonfinite harmonic gravity GM/R scaling" );
                }
                matrix->row( n - minimumDegree ) *= scale;
            }
            if( !matrix->allFinite( ) )
            {
                throw std::invalid_argument( "Nonfinite scaled harmonic gravity coefficient" );
            }
        }
    }
}
//! Sum the canonical quadratures at the supplied physical phase, without extra nodal or loading factors.
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > HarmonicGravityFieldVariations::calculateSphericalHarmonicsCorrections( double time )
{
    if( !std::isfinite( time ) )
    {
        throw std::invalid_argument( "Nonfinite harmonic gravity epoch" );
    }
    const Eigen::Vector6d arguments = doodsonArguments_( time );

    if( !arguments.allFinite( ) )
    {
        throw std::invalid_argument( "Nonfinite harmonic gravity arguments" );
    }
    Eigen::MatrixXd cosine = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ ), sine = cosine;

    for( const auto& wave : constituents_ )
    {
        double phase = 0.0;

        for( int i = 0; i < 6; ++i )
        {
            phase += wave.multipliers[ i ] * arguments[ i ];
        }
        if( !std::isfinite( phase ) )
        {
            throw std::invalid_argument( "Nonfinite harmonic gravity phase" );
        }
        const double c = std::cos( phase ), s = std::sin( phase );

        cosine += c * wave.cCos + s * wave.cSin;
        sine += c * wave.sCos + s * wave.sSin;
    }
    return { cosine, sine };
}
}  // namespace gravitation
}  // namespace tudat
