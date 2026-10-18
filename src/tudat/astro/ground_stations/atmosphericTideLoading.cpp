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

input_output::AtmosphericTideCoefficients interpolateAtmosphericTideGrid( const input_output::AtmosphericTideGrid& grid,
                                                                          const double latitude,
                                                                          const double longitude )
{
    const auto& lon = grid.longitudesDegrees;
    const auto& lat = grid.latitudesDegrees;
    const double latStep = lat[ 1 ] - lat[ 0 ], lonStep = lon[ 1 ] - lon[ 0 ];
    const double y = ( latitude - lat.front( ) ) / latStep;
    if( y < -1e-10 || y > lat.size( ) - 1 + 1e-10 )
    {
        throw std::out_of_range( "Latitude outside atmospheric grid coverage" );
    }

    // Longitude offset from the first column, in [0, 360) degrees.
    double degrees = std::fmod( ( longitude - lon.front( ) ) * ( lonStep > 0 ? 1 : -1 ), 360.0 );
    if( degrees < 0 )
    {
        degrees += 360;
    }
    const double x = degrees / std::abs( lonStep );
    const bool duplicatedSeam = std::abs( std::abs( lon.back( ) - lon.front( ) ) - 360 ) < 1e-9;
    const bool periodic = duplicatedSeam || std::abs( lon.size( ) * std::abs( lonStep ) - 360 ) < 1e-9;
    const std::size_t columns = lon.size( ) - ( duplicatedSeam ? 1 : 0 );
    if( !periodic && x > lon.size( ) - 1 + 1e-10 )
    {
        throw std::out_of_range( "Longitude outside atmospheric grid coverage" );
    }
    const std::size_t i = std::min( static_cast< std::size_t >( std::floor( x ) ), periodic ? columns - 1 : lon.size( ) - 2 );
    const std::size_t j = std::min( static_cast< std::size_t >( std::floor( std::max( y, 0.0 ) ) ), lat.size( ) - 2 );
    const std::size_t nextI = ( i + 1 ) % columns;
    const double a = std::min( x - i, 1.0 ), b = std::clamp( y - j, 0.0, 1.0 );
    const auto& c = grid.coefficientsMetres;
    const std::size_t n = lat.size( );
    return ( 1 - a ) * ( ( 1 - b ) * c[ i * n + j ] + b * c[ i * n + j + 1 ] ) +
            a * ( ( 1 - b ) * c[ nextI * n + j ] + b * c[ nextI * n + j + 1 ] );
}
}  // namespace basic_astrodynamics
}  // namespace tudat
