#include "tudat/astro/ground_stations/poleTideDeformation.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/earth_orientation/polarMotionCalculator.h"

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

Eigen::Vector3d PoleTideDeformation::calculateGeocentricEnuDisplacement(
        const double tt,
        const std::shared_ptr< ground_stations::GroundStationState > station )
{
    const Eigen::Vector2d polarMas =
            unit_conversions::convertRadiansToDegrees( polarMotionCalculator_->getPositionOfCipInItrsAtTt( tt ) ) * 3600000.0;
    return poleTideGeocentricEnu( polarMas, tt, station );
}

Eigen::Vector3d PoleTideDeformation::calculateDisplacement( const EarthDeformationEpoch& epoch, const Eigen::Vector3d& position )
{
    return calculateDisplacement(
            epoch, std::make_shared< ground_stations::GroundStationState >( position, coordinate_conversions::cartesian_position ) );
}
Eigen::Vector3d PoleTideDeformation::calculateDisplacement( const EarthDeformationEpoch& epoch,
                                                            const std::shared_ptr< ground_stations::GroundStationState > station )
{
    const auto native = contextClock_ ? epoch.withTimeConverter( contextClock_ ) : epoch;
    const Eigen::Vector2d polarMas =
            unit_conversions::convertRadiansToDegrees( polarMotionCalculator_->getPositionOfCipInItrsAtTt( native.getTt( ) ) ) * 3600000.0;
    return geocentricEnuToCartesian( poleTideGeocentricEnu( polarMas, static_cast< double >( native.getTt( ) ), station ), station );
}
Eigen::Vector3d PoleTideDeformation::calculateDisplacementVelocity( const EarthDeformationEpoch& epoch, const Eigen::Vector3d& position )
{
    return calculateTerrestrialDisplacementVelocity( contextClock_ ? epoch.withTimeConverter( contextClock_ ) : epoch, position );
}
Eigen::Vector3d PoleTideDeformation::calculateDisplacementVelocity( const EarthDeformationEpoch& epoch,
                                                                    const std::shared_ptr< ground_stations::GroundStationState > station )
{
    return calculateTerrestrialDisplacementVelocity( contextClock_ ? epoch.withTimeConverter( contextClock_ ) : epoch, station );
}
Eigen::Vector6d PoleTideDeformation::calculateDisplacementState( const EarthDeformationEpoch& epoch, const Eigen::Vector3d& position )
{
    return calculateNativeDisplacementState( epoch, position );
}
Eigen::Vector6d PoleTideDeformation::calculateDisplacementState( const EarthDeformationEpoch& epoch,
                                                                 const std::shared_ptr< ground_stations::GroundStationState > station )
{
    return calculateNativeDisplacementState( epoch, station );
}

Eigen::Vector3d PoleTideDeformation::calculateEnuDisplacement(
        const double ephemerisTime,
        const std::shared_ptr< ground_stations::GroundStationState > nominalSiteState )
{
    return nominalSiteState->getRotationFromBodyFixedToTopocentricFrame( ephemerisTime ) *
            calculateDisplacement( ephemerisTime, nominalSiteState );
}

//! Function to calculate site displacement due to pole tide deformation.
Eigen::Vector3d PoleTideDeformation::calculateDisplacement( const double ephemerisTime,
                                                            const std::shared_ptr< ground_stations::GroundStationState > nominalSiteState )
{
    // IERS Eqs. (27)-(28): rotate the spherical local displacement to ITRS.
    // Using the geocentric basis makes the Cartesian result independent of the
    // station's shape model (whose topocentric frame can be geodetic).
    const Eigen::Vector3d geocentricEnuDisplacement = calculateGeocentricEnuDisplacement( ephemerisTime, nominalSiteState );
    const std::vector< Eigen::Vector3d > unitVectors = nominalSiteState->getEnuGeocentricUnitVectors( );
    return geocentricEnuDisplacement( 0 ) * unitVectors[ 0 ] + geocentricEnuDisplacement( 1 ) * unitVectors[ 1 ] +
            geocentricEnuDisplacement( 2 ) * unitVectors[ 2 ];
}

}  // namespace basic_astrodynamics

}  // namespace tudat
