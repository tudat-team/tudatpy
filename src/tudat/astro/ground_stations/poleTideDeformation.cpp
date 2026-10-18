#include "tudat/astro/ground_stations/poleTideDeformation.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/earth_orientation/polarMotionCalculator.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"

namespace tudat
{

namespace basic_astrodynamics
{

namespace
{
Eigen::Vector3d poleTideGeocentricEnu( const Eigen::Vector2d& polarMas,
                                       const double tt,
                                       const std::shared_ptr< ground_stations::GroundStationState > nominalSiteState )
{
    const Eigen::Vector2d secularMas = earth_orientation::getSecularPolePositionInMas( tt );

    // IERS Eq. (25): wobble variables in arcseconds, as required by Eq. (26).
    const double m1 = ( polarMas( 0 ) - secularMas( 0 ) ) / 1000.0;
    const double m2 = -( polarMas( 1 ) - secularMas( 1 ) ) / 1000.0;

    // Retrieve geocentric latitude and east longitude, in radians.
    double latitude = nominalSiteState->getNominalLatitude( );
    double longitude = nominalSiteState->getNominalLongitude( );

    // IERS Eq. (26), with colatitude theta = pi/2 - latitude and north = -south.
    // The coefficients below convert the published millimetres to metres.
    // Up is radial here; it is not the ellipsoid normal of a geodetic ENU frame.
    Eigen::Vector3d enuDisplacement;
    enuDisplacement( 0 ) = 9.0E-3 * sin( latitude ) * ( m1 * sin( longitude ) - m2 * cos( longitude ) );
    enuDisplacement( 1 ) = -9.0E-3 * cos( 2.0 * latitude ) * ( m1 * cos( longitude ) + m2 * sin( longitude ) );
    enuDisplacement( 2 ) = -33.0E-3 * sin( 2.0 * latitude ) * ( m1 * cos( longitude ) + m2 * sin( longitude ) );

    return enuDisplacement;
}
Eigen::Vector3d geocentricEnuToCartesian( const Eigen::Vector3d& enu, const std::shared_ptr< ground_stations::GroundStationState > station )
{
    const auto basis = station->getEnuGeocentricUnitVectors( );
    return enu.x( ) * basis[ 0 ] + enu.y( ) * basis[ 1 ] + enu.z( ) * basis[ 2 ];
}
}  // namespace

Eigen::Vector3d PoleTideDeformation::displacementAtTt( const Time& tt,
                                                       const std::shared_ptr< ground_stations::GroundStationState > station )
{
    const Eigen::Vector2d polarMas =
            unit_conversions::convertRadiansToDegrees( polarMotionCalculator_->getPositionOfCipInItrsAtTt( tt ) ) * 3600000.0;
    return geocentricEnuToCartesian( poleTideGeocentricEnu( polarMas, static_cast< double >( tt ), station ), station );
}

Eigen::Vector3d PoleTideDeformation::evaluateDisplacement( const Time& tdb,
                                                           const std::shared_ptr< ground_stations::GroundStationState > station )
{
    return displacementAtTt( timeConverter_->getCurrentTime< Time >( tdb_scale, tt_scale, tdb ), station );
}

Eigen::Vector6d PoleTideDeformation::evaluateDisplacementState( const Time& tdb,
                                                                const std::shared_ptr< ground_stations::GroundStationState > station )
{
    const Time tt = timeConverter_->getCurrentTime< Time >( tdb_scale, tt_scale, tdb );
    const double rate = timeConverter_->getGeocentricTtRateWithRespectToTdb( tdb );
    const Eigen::Vector3d velocity =
            differentiateDisplacement( tt, [ & ]( const Time& sample ) { return displacementAtTt( sample, station ); } ) * rate;
    return ( Eigen::Vector6d( ) << displacementAtTt( tt, station ), velocity ).finished( );
}

}  // namespace basic_astrodynamics

}  // namespace tudat
