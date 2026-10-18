#include "tudat/astro/ground_stations/oceanPoleTideDeformation.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace tudat
{
namespace basic_astrodynamics
{
Eigen::Vector3d calculateOceanPoleTideLoading( const Eigen::Matrix< double, 3, 2 >& coefficients,
                                               const Eigen::Vector2d& wobbleRadians,
                                               const double scaleMetresPerRadian )
{
    if( !coefficients.allFinite( ) || !wobbleRadians.allFinite( ) || !std::isfinite( scaleMetresPerRadian ) || scaleMetresPerRadian <= 0 )
    {
        throw std::runtime_error( "Nonfinite or invalid ocean pole-tide loading input" );
    }

    // IERS Eq. (29), gamma = 0.6870 + 0.0036i. Coefficients are dimensionless; m1/m2 are RADIANS.
    // Chapter 1: G=6.67428e-11, a=6378136.6, GM=3.986004418e14, ge=9.7803278, Omega=7.292115e-5;
    // rho=1025 (Section 7.1.5). Hp=sqrt(8*pi/15)*Omega^2*a^4/GM; K=4*pi*G*a*rho*Hp/(3*ge).
    return scaleMetresPerRadian *
            ( ( wobbleRadians.x( ) * 0.6870 + wobbleRadians.y( ) * 0.0036 ) * coefficients.col( 0 ) +
              ( wobbleRadians.y( ) * 0.6870 - wobbleRadians.x( ) * 0.0036 ) * coefficients.col( 1 ) );
}

//! Interpolate a checked grid, reporting absent latitude/longitude coverage separately from invalid coordinates.
static Eigen::Matrix< double, 3, 2 > interpolateValidatedOceanPoleTideCoefficients( const input_output::OceanPoleTideGrid& grid,
                                                                                    const double latitudeDegrees,
                                                                                    const double longitudeDegrees )
{
    if( !std::isfinite( latitudeDegrees ) || !std::isfinite( longitudeDegrees ) || latitudeDegrees < -90 || latitudeDegrees > 90 )
    {
        throw std::runtime_error( "Invalid ocean pole-tide interpolation coordinates" );
    }
    const double lastLatitude = grid.firstLatitudeDegrees + ( grid.latitudeCount - 1 ) * grid.latitudeStepDegrees;
    double latitude = latitudeDegrees;
    if( std::abs( latitude - grid.firstLatitudeDegrees ) < 1e-10 )
    {
        latitude = grid.firstLatitudeDegrees;
    }
    if( std::abs( latitude - lastLatitude ) < 1e-10 )
    {
        latitude = lastLatitude;
    }

    // The official grid stops half a step short of the poles; use the edge row there, but reject latitudes
    // further outside a (regional) grid.
    if( latitude < grid.firstLatitudeDegrees - grid.latitudeStepDegrees || latitude > lastLatitude + grid.latitudeStepDegrees )
    {
        throw std::out_of_range( "Ocean pole-tide latitude lies outside the supplied grid coverage" );
    }
    latitude = std::clamp( latitude, grid.firstLatitudeDegrees, lastLatitude );
    double longitude = longitudeDegrees - grid.firstLongitudeDegrees;

    // Even a regional extract accepts equivalent east longitudes, but it never interpolates across missing longitudes.
    longitude = std::fmod( longitude, 360.0 );
    if( longitude < 0 )
    {
        longitude += 360.0;
    }
    if( longitude < 1e-10 || 360.0 - longitude < 1e-10 )
    {
        longitude = 0.0;
    }
    const bool periodic = std::abs( grid.longitudeCount * grid.longitudeStepDegrees - 360.0 ) < 1e-9;
    const double lastLongitudeOffset = ( grid.longitudeCount - 1 ) * grid.longitudeStepDegrees;
    if( !periodic && longitude > lastLongitudeOffset + 1e-9 )
    {
        throw std::out_of_range( "Ocean pole-tide longitude lies outside the supplied grid coverage" );
    }
    if( !periodic )
    {
        longitude = std::min( longitude, lastLongitudeOffset );
    }
    const double y = ( latitude - grid.firstLatitudeDegrees ) / grid.latitudeStepDegrees;
    const std::size_t row = std::min( static_cast< std::size_t >( std::floor( y ) ), grid.latitudeCount - 2 );
    const double northWeight = y - row;
    const double x = longitude / grid.longitudeStepDegrees;
    const std::size_t column = periodic ? std::min( static_cast< std::size_t >( std::floor( x ) ), grid.longitudeCount - 1 )
                                        : std::min( static_cast< std::size_t >( std::floor( x ) ), grid.longitudeCount - 2 );
    const std::size_t nextColumn = ( column + 1 ) % grid.longitudeCount;
    const double eastWeight = x - column;
    const auto& sw = grid.coefficients.at( row * grid.longitudeCount + column );
    const auto& se = grid.coefficients.at( row * grid.longitudeCount + nextColumn );
    const auto& nw = grid.coefficients.at( ( row + 1 ) * grid.longitudeCount + column );
    const auto& ne = grid.coefficients.at( ( row + 1 ) * grid.longitudeCount + nextColumn );
    return ( 1 - northWeight ) * ( ( 1 - eastWeight ) * sw + eastWeight * se ) +
            northWeight * ( ( 1 - eastWeight ) * nw + eastWeight * ne );
}

Eigen::Matrix< double, 3, 2 > interpolateOceanPoleTideCoefficients( const input_output::OceanPoleTideGrid& grid,
                                                                    const double latitudeDegrees,
                                                                    const double longitudeDegrees )
{
    return interpolateValidatedOceanPoleTideCoefficients( grid, latitudeDegrees, longitudeDegrees );
}

OceanPoleTideDeformation::OceanPoleTideDeformation( const input_output::OceanPoleTideGrid& grid,
                                                    const std::shared_ptr< earth_orientation::PolarMotionCalculator > polarMotionCalculator,
                                                    const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter ):
    BodyDeformationModel( converter ), grid_( grid ), polarMotionCalculator_( polarMotionCalculator )
{
    if( !polarMotionCalculator_ )
    {
        throw std::runtime_error( "Ocean pole-tide loading requires a polar motion calculator" );
    }
}

//! Evaluate ocean pole-tide displacement at TT using valid nominal geocentric station geometry.
Eigen::Vector3d OceanPoleTideDeformation::calculateDisplacementAtTt(
        const Time& tt,
        const std::shared_ptr< ground_stations::GroundStationState > stationState )
{
    const double time = static_cast< double >( tt );
    if( !std::isfinite( time ) || !stationState )
    {
        throw std::runtime_error( "Invalid ocean pole-tide epoch or station" );
    }

    // A station at the geocentre has no local radial/north/east basis, regardless of grid coverage.
    const Eigen::Vector3d sphericalPosition = stationState->getNominalSphericalPosition( );

    if( !sphericalPosition.allFinite( ) || sphericalPosition.x( ) <= 0.0 )
    {
        throw std::runtime_error( "Invalid ocean pole-tide station geometry" );
    }
    const double radiansPerDegree = mathematical_constants::PI / 180.0;

    const auto coefficients = interpolateValidatedOceanPoleTideCoefficients(
            grid_, sphericalPosition.y( ) / radiansPerDegree, sphericalPosition.z( ) / radiansPerDegree );

    const Eigen::Vector2d pole = polarMotionCalculator_->getPositionOfCipInItrsAtTt( tt );
    const Eigen::Vector2d secularPole = earth_orientation::getSecularPolePositionInMas( time ) * ( radiansPerDegree / 3600000.0 );
    const Eigen::Vector2d wobble( pole.x( ) - secularPole.x( ), -( pole.y( ) - secularPole.y( ) ) );
    const Eigen::Vector3d radialNorthEast = calculateOceanPoleTideLoading( coefficients, wobble );
    const auto basis = stationState->getEnuGeocentricUnitVectors( );
    return radialNorthEast.z( ) * basis[ 0 ] + radialNorthEast.y( ) * basis[ 1 ] + radialNorthEast.x( ) * basis[ 2 ];
}

Eigen::Vector3d OceanPoleTideDeformation::evaluateDisplacement( const Time& tdb,
                                                                const std::shared_ptr< ground_stations::GroundStationState > station )
{
    return calculateDisplacementAtTt( timeConverter_->getCurrentTime< Time >( tdb_scale, tt_scale, tdb ), station );
}
Eigen::Vector6d OceanPoleTideDeformation::evaluateDisplacementState( const Time& tdb,
                                                                     const std::shared_ptr< ground_stations::GroundStationState > station )
{
    const Time tt = timeConverter_->getCurrentTime< Time >( tdb_scale, tt_scale, tdb );
    const auto position = [ this, station ]( const Time& time ) { return calculateDisplacementAtTt( time, station ); };
    Eigen::Vector6d state;
    state.head< 3 >( ) = position( tt );
    state.tail< 3 >( ) = differentiateDisplacement( tt, position ) * timeConverter_->getGeocentricTtRateWithRespectToTdb( tdb );
    return state;
}
}  // namespace basic_astrodynamics
}  // namespace tudat
