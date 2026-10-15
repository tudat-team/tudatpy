#include "tudat/astro/ground_stations/atmosphericTideLoading.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace tudat
{
namespace basic_astrodynamics
{
namespace
{
double phase( const Eigen::Matrix< double, 3, 4 >& coefficients, const double ut1 )
{
    if( !coefficients.allFinite( ) || !std::isfinite( ut1 ) )
    {
        throw std::runtime_error( "Nonfinite atmospheric tide synthesis input" );
    }

    // RP03 Eq.2 uses Universal Time; p1905 identifies its quadrature component with06:00 local time.
    // IERS2010 Eqs.19a/b prescribe the same convention in UT1. J2000's time origin is noon, not midnight.
    // Reduce before adding half a day, preserving phase resolution for large epochs.
    return 2 * mathematical_constants::PI * ( std::remainder( ut1, 86400.0 ) / 86400.0 + 0.5 );
}
Eigen::Matrix3d radialNorthEastBasis( double latitude, double longitude )
{
    const double a = latitude * mathematical_constants::PI / 180;
    const double b = longitude * mathematical_constants::PI / 180;
    Eigen::Matrix3d result;
    result.col( 0 ) << std::cos( a ) * std::cos( b ), std::cos( a ) * std::sin( b ), std::sin( a );
    result.col( 1 ) << -std::sin( a ) * std::cos( b ), -std::sin( a ) * std::sin( b ), std::cos( a );
    result.col( 2 ) << -std::sin( b ), std::cos( b ), 0;
    return result;
}
}  // namespace
Eigen::Vector3d synthesizeAtmosphericTideLoading( const Eigen::Matrix< double, 3, 4 >& coefficients, const double ut1 )
{
    const double angle = phase( coefficients, ut1 );
    return coefficients.col( 0 ) * std::cos( angle ) + coefficients.col( 1 ) * std::sin( angle ) +
            coefficients.col( 2 ) * std::cos( 2 * angle ) + coefficients.col( 3 ) * std::sin( 2 * angle );
}
Eigen::Vector3d differentiateAtmosphericTideLoadingWrtUt1( const Eigen::Matrix< double, 3, 4 >& coefficients, const double ut1 )
{
    const double angle = phase( coefficients, ut1 );
    return ( 2 * mathematical_constants::PI / 86400.0 ) *
            ( -coefficients.col( 0 ) * std::sin( angle ) + coefficients.col( 1 ) * std::cos( angle ) -
              2 * coefficients.col( 2 ) * std::sin( 2 * angle ) + 2 * coefficients.col( 3 ) * std::cos( 2 * angle ) );
}

namespace
{
Eigen::Matrix< double, 3, 4 > rp03GeocentreCoefficients( )
{
    // IERS2010 Chapter7 Eq.20/Table7.6: Cartesian X/Y/Z coefficients in metres.
    Eigen::Matrix< double, 3, 4 > coefficients;
    coefficients << 2.1188e-4, -7.6861e-4, 1.4472e-4, -1.7844e-4, -7.2766e-4, -2.3582e-4, -3.2691e-4, -1.5878e-4, -1.2176e-5, 3.2243e-5,
            -9.6271e-5, 1.6976e-5;
    return coefficients;
}
}  // namespace
Eigen::Vector3d calculateIers2010Rp03AtmosphericGeocentreTranslation( const double ut1SecondsSinceJ2000 )
{
    // Eq.18: X_crust=X_CM-dX; CE-to-CM adds this Cartesian vector after local rotation.
    return synthesizeAtmosphericTideLoading( rp03GeocentreCoefficients( ), ut1SecondsSinceJ2000 );
}
Eigen::Vector3d differentiateIers2010Rp03AtmosphericGeocentreTranslationWrtUt1( const double ut1SecondsSinceJ2000 )
{
    return differentiateAtmosphericTideLoadingWrtUt1( rp03GeocentreCoefficients( ), ut1SecondsSinceJ2000 );
}

AtmosphericTideGridInterpolator::AtmosphericTideGridInterpolator( const input_output::AtmosphericTideGrid& grid ):
    providerGrid_( grid ), modelCoefficients_( grid.coefficientsMetres )
{
    input_output::validateAtmosphericTideGrid( grid );
    const auto& lon = grid.longitudesDegrees;
    const auto& lat = grid.latitudesDegrees;
    const std::size_t nlon = lon.size( ), nlat = lat.size( );
    const double step = std::abs( lon[ 1 ] - lon[ 0 ] );
    duplicatedSeam_ = std::abs( std::abs( lon.back( ) - lon.front( ) ) - 360 ) < 1e-9;
    periodic_ = duplicatedSeam_ || std::abs( nlon * step - 360 ) < 1e-9;
    const std::size_t uniqueLongitudes = nlon - ( duplicatedSeam_ ? 1 : 0 );
    if( duplicatedSeam_ )
    {
        for( std::size_t j = 0; j < nlat; ++j )
        {
            auto& first = modelCoefficients_[ j ];
            auto& last = modelCoefficients_[ ( nlon - 1 ) * nlat + j ];
            seamDiscrepancy_ = std::max( seamDiscrepancy_, ( first - last ).cwiseAbs( ).maxCoeff( ) );
            const Eigen::Matrix< double, 3, 4 > mean = ( first + last ) / 2;
            first = mean;
            last = mean;
        }
    }
    for( std::size_t j = 0; j < nlat; ++j )
    {
        if( std::abs( lat[ j ] ) != 90 )
        {
            continue;
        }
        Eigen::Matrix< double, 3, 4 > mean = Eigen::Matrix< double, 3, 4 >::Zero( );
        Eigen::Matrix< double, 3, 4 > minimum, maximum;
        for( std::size_t i = 0; i < uniqueLongitudes; ++i )
        {
            const Eigen::Matrix< double, 3, 4 > cartesian =
                    radialNorthEastBasis( lat[ j ], lon[ i ] ) * grid.coefficientsMetres[ i * nlat + j ];
            if( i == 0 )
            {
                minimum = cartesian;
                maximum = cartesian;
            }
            else
            {
                minimum = minimum.cwiseMin( cartesian );
                maximum = maximum.cwiseMax( cartesian );
            }
            mean += radialNorthEastBasis( lat[ j ], lon[ i ] ) * modelCoefficients_[ i * nlat + j ];
        }
        polarDiscrepancy_ = std::max( polarDiscrepancy_, ( maximum - minimum ).cwiseAbs( ).maxCoeff( ) );
        mean /= uniqueLongitudes;
        for( std::size_t i = 0; i < nlon; ++i )
        {
            modelCoefficients_[ i * nlat + j ] = radialNorthEastBasis( lat[ j ], lon[ i ] ).transpose( ) * mean;
        }
    }
    if( seamDiscrepancy_ > maximumEndpointDiscrepancyMetres || polarDiscrepancy_ > maximumEndpointDiscrepancyMetres )
    {
        throw std::runtime_error( "Atmospheric grid seam/polar coefficients differ by more than provider rounding" );
    }
}

//! Interpolate provider metre coefficients, distinguishing invalid coordinates from absent grid coverage.
Eigen::Matrix< double, 3, 4 > AtmosphericTideGridInterpolator::interpolate( const double latitude, const double longitude ) const
{
    if( !std::isfinite( latitude ) || !std::isfinite( longitude ) || std::abs( latitude ) > 90 )
    {
        throw std::runtime_error( "Invalid atmospheric tide interpolation coordinate" );
    }
    const auto& lon = providerGrid_.longitudesDegrees;
    const auto& lat = providerGrid_.latitudesDegrees;
    const double latStep = lat[ 1 ] - lat[ 0 ], lonStep = lon[ 1 ] - lon[ 0 ];
    double y = ( latitude - lat.front( ) ) / latStep;
    if( y < -1e-10 || y > lat.size( ) - 1 + 1e-10 )
    {
        throw std::out_of_range( "Latitude outside atmospheric grid coverage" );
    }
    y = std::clamp( y, 0.0, static_cast< double >( lat.size( ) - 1 ) );
    double degrees = std::fmod( ( longitude - lon.front( ) ) * ( lonStep > 0 ? 1 : -1 ), 360.0 );
    if( degrees < 0 )
    {
        degrees += 360;
    }
    if( degrees < 1e-10 || 360 - degrees < 1e-10 )
    {
        degrees = 0;
    }
    double x = degrees / std::abs( lonStep );
    const std::size_t uniqueLongitudes = lon.size( ) - ( duplicatedSeam_ ? 1 : 0 );
    if( !periodic_ && x > lon.size( ) - 1 + 1e-10 )
    {
        throw std::out_of_range( "Longitude outside atmospheric grid coverage" );
    }
    if( !periodic_ )
    {
        x = std::min( x, static_cast< double >( lon.size( ) - 1 ) );
    }
    if( std::abs( latitude ) == 90 )
    {
        // At a pole the averaged/consistent Cartesian coefficient is unique. Interpolating local N/E components
        // between different longitude bases would reintroduce artificial longitude dependence.
        const std::size_t row = std::abs( latitude - lat.front( ) ) < 1e-10 ? 0 : lat.size( ) - 1;
        const Eigen::Matrix< double, 3, 4 > cartesian = radialNorthEastBasis( lat[ row ], lon.front( ) ) * modelCoefficients_[ row ];
        return radialNorthEastBasis( latitude, longitude ).transpose( ) * cartesian;
    }
    const std::size_t i = periodic_ ? std::min( static_cast< std::size_t >( std::floor( x ) ), uniqueLongitudes - 1 )
                                    : std::min( static_cast< std::size_t >( std::floor( x ) ), lon.size( ) - 2 );
    const std::size_t j = std::min( static_cast< std::size_t >( std::floor( y ) ), lat.size( ) - 2 );
    const std::size_t nextI = ( i + 1 ) % uniqueLongitudes;
    const double a = x - i, b = y - j;
    return ( 1 - a ) * ( ( 1 - b ) * modelCoefficients_[ i * lat.size( ) + j ] + b * modelCoefficients_[ i * lat.size( ) + j + 1 ] ) +
            a * ( ( 1 - b ) * modelCoefficients_[ nextI * lat.size( ) + j ] + b * modelCoefficients_[ nextI * lat.size( ) + j + 1 ] );
}
}  // namespace basic_astrodynamics
}  // namespace tudat
