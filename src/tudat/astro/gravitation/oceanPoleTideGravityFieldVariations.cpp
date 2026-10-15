/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/astro/gravitation/oceanPoleTideGravityFieldVariations.h"
#include <cmath>
#include <stdexcept>

namespace tudat
{
namespace gravitation
{
namespace
{
//! Validate source coverage and output bounds before the base allocates correction matrices.
int checkedMinimumDegree( const input_output::OceanPoleGravityData& data,
                          int minimumDegree,
                          int minimumOrder,
                          int maximumDegree,
                          int maximumOrder )
{
    input_output::validateOceanPoleGravityData( data );
    if( minimumDegree < 2 || maximumDegree < minimumDegree || maximumDegree > data.maximumDegree || minimumOrder < 0 ||
        maximumOrder < minimumOrder || maximumOrder > maximumDegree )
    {
        throw std::invalid_argument( "Invalid ocean-pole gravity degree/order selection" );
    }
    return minimumDegree;
}
}  // namespace

//! Retain the selected IERS 6.23 equilibrium integrals with explicit load Love numbers and physical reference constants.
OceanPoleTideGravityFieldVariations::OceanPoleTideGravityFieldVariations( const input_output::OceanPoleGravityData& data,
                                                                          const std::function< Eigen::Vector2d( double ) >& wobbleRadians,
                                                                          const std::vector< double >& loadLoveNumbers,
                                                                          const OceanPoleGravityConstants& constants,
                                                                          double outputGravitationalParameter,
                                                                          double outputReferenceRadius,
                                                                          int minimumDegree,
                                                                          int minimumOrder,
                                                                          int maximumDegree,
                                                                          int maximumOrder ):
    GravityFieldVariations( checkedMinimumDegree( data, minimumDegree, minimumOrder, maximumDegree, maximumOrder ),
                            minimumOrder,
                            maximumDegree,
                            maximumOrder ),
    wobbleRadians_( wobbleRadians ), gamma2_( constants.gamma2 )
{
    if( !wobbleRadians_ )
    {
        throw std::invalid_argument( "Missing ocean-pole gravity wobble callback" );
    }
    for( double value : { constants.angularVelocity,
                          constants.gravitationalConstant,
                          constants.oceanDensity,
                          constants.equatorialGravity,
                          constants.gravitationalParameter,
                          constants.referenceRadius,
                          outputGravitationalParameter,
                          outputReferenceRadius } )
    {
        if( !std::isfinite( value ) || value <= 0.0 )
        {
            throw std::invalid_argument( "Ocean-pole gravity constants must be finite and positive" );
        }
    }
    if( !std::isfinite( gamma2_.real( ) ) || !std::isfinite( gamma2_.imag( ) ) ||
        loadLoveNumbers.size( ) <= static_cast< std::size_t >( maximumDegree ) )
    {
        throw std::invalid_argument( "Missing ocean-pole gravity response constants" );
    }
    aReal_ = aImaginary_ = bReal_ = bImaginary_ = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ );
    // Eq. 6.23, followed by exactly one source-to-output GM/R conversion per degree.
    const double common = constants.angularVelocity * constants.angularVelocity * std::pow( constants.referenceRadius, 4 ) * 4.0 *
            mathematical_constants::PI * constants.gravitationalConstant * constants.oceanDensity /
            ( constants.gravitationalParameter * constants.equatorialGravity );

    for( int n = minimumDegree; n <= maximumDegree; ++n )
    {
        if( !std::isfinite( loadLoveNumbers[ n ] ) )
        {
            throw std::invalid_argument( "Nonfinite ocean-pole load Love number" );
        }
        const double factor = common * ( 1.0 + loadLoveNumbers[ n ] ) / ( 2.0 * n + 1.0 ) * constants.gravitationalParameter /
                outputGravitationalParameter * std::pow( constants.referenceRadius / outputReferenceRadius, n );

        if( !std::isfinite( factor ) )
        {
            throw std::invalid_argument( "Nonfinite ocean-pole gravity scaling" );
        }
        for( int m = minimumOrder; m <= std::min( n, maximumOrder ); ++m )
        {
            const auto found = data.coefficients.find( { n, m } );

            if( found == data.coefficients.end( ) )
            {
                throw std::invalid_argument( "Missing selected ocean-pole gravity coefficient" );
            }
            const int i = n - minimumDegree, j = m - minimumOrder;

            aReal_( i, j ) = factor * found->second.a.real( );
            aImaginary_( i, j ) = factor * found->second.a.imag( );
            bReal_( i, j ) = factor * found->second.b.real( );
            bImaginary_( i, j ) = factor * found->second.b.imag( );
        }
    }
}

//! Combine the scaled equilibrium integrals with the selected polar wobble and complex response.
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > OceanPoleTideGravityFieldVariations::calculateSphericalHarmonicsCorrections( double time )
{
    if( !std::isfinite( time ) )
    {
        throw std::invalid_argument( "Nonfinite ocean-pole gravity epoch" );
    }
    const Eigen::Vector2d wobble = wobbleRadians_( time );

    if( !wobble.allFinite( ) )
    {
        throw std::invalid_argument( "Nonfinite ocean-pole wobble" );
    }
    const double real = wobble.x( ) * gamma2_.real( ) + wobble.y( ) * gamma2_.imag( );

    const double imaginary = wobble.y( ) * gamma2_.real( ) - wobble.x( ) * gamma2_.imag( );

    return { aReal_ * real + aImaginary_ * imaginary, bReal_ * real + bImaginary_ * imaginary };
}
}  // namespace gravitation
}  // namespace tudat
