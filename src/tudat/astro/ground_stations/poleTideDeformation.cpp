#include "tudat/astro/ground_stations/poleTideDeformation.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/earth_orientation/polarMotionCalculator.h"

namespace tudat
{

namespace basic_astrodynamics
{

Eigen::Vector3d PoleTideDeformation::calculateGeocentricEnuDisplacement(
        const double ephemerisTime,
        const std::shared_ptr< ground_stations::GroundStationState > nominalSiteState )
{
    // Retrieve position of CIP in ITRS from polar motion calculator.
    double tt = ephemerisTime;
    Eigen::Vector2d cipPositionInItrs =
            unit_conversions::convertRadiansToDegrees( polarMotionCalculator_->getPositionOfCipInItrsAtTt( tt ) ) * 3600.0 * 1000.0;

    // Retrieve the conventional secular pole (IERS update of 1 February 2018), in mas.
    Eigen::Vector2d secularPolePositionInItrs = earth_orientation::getSecularPolePositionInMas( ephemerisTime );

    // IERS Eq. (25): wobble variables in arcseconds, as required by Eq. (26).
    const double m1 = ( cipPositionInItrs( 0 ) - secularPolePositionInItrs( 0 ) ) / 1000.0;
    const double m2 = -( cipPositionInItrs( 1 ) - secularPolePositionInItrs( 1 ) ) / 1000.0;

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
