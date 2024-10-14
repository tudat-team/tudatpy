#include "tudat/astro/ground_stations/poleTideDeformation.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/earth_orientation/polarMotionCalculator.h"

namespace tudat
{

namespace basic_astrodynamics
{

Eigen::Vector3d PoleTideDeformation::calculateEnuDisplacement(
    const double ephemerisTime,
    const std::shared_ptr< ground_stations::GroundStationState > nominalSiteState )
{
    // Retrieve position of CIP in ITRS from polar motion calculator.
    double tt = ephemerisTime;
    double utc = sofa_interface::convertTTtoUTC< double >( tt );
    Eigen::Vector2d cipPositionInItrs =
        unit_conversions::convertRadiansToDegrees( polarMotionCalculator_->getPositionOfCipInItrs( tt, utc ) ) * 3600.0 * 1000.0;

    // Retrieve secular pole
    Eigen::Vector2d secularPolePositionInItrs = earth_orientation::getSecularPolePositionInMas( ephemerisTime );

    // Calculate deviation of CIP position from mean.
    double m1 = cipPositionInItrs( 0 ) - secularPolePositionInItrs( 0 );
    double m2 = -( cipPositionInItrs( 1 ) - secularPolePositionInItrs( 1 ) );

    // Retrieve station latitude and longitude.
    double latitude = nominalSiteState->getNominalLatitude( );
    double longitude = nominalSiteState->getNominalLongitude( );

    double radiansToArcSeconds = unit_conversions::convertArcSecondsToRadians( 1.0 );

    // Calculate site displacement in local coordinates (east, north, up), IERS 2010 Conventions Eq. 7.26
    Eigen::Vector3d enuDisplacement;
    enuDisplacement( 0 ) = 1000.0 * radiansToArcSeconds * 9.0E-3 * cos( latitude ) * ( m1 * sin( longitude ) - m2 * cos( longitude ) );
    enuDisplacement( 1 ) = 1000.0 * -radiansToArcSeconds * 9.0E-3 * cos( 2.0 * latitude ) * ( m1 * cos( longitude ) + m2 * sin( longitude ) );
    enuDisplacement( 2 ) = 1000.0 * -radiansToArcSeconds * 33.0E-3 * sin( 2.0 * latitude ) * ( m1 * cos( longitude ) - m2 * sin( longitude ) );

    return enuDisplacement;
}
//! Function to calculate site displacement due to pole tide deformation.
Eigen::Vector3d PoleTideDeformation::calculateDisplacement(
    const double ephemerisTime,
    const std::shared_ptr< ground_stations::GroundStationState > nominalSiteState )
{
    // Transform to ITRS (geocentric) frame
    return  nominalSiteState->getRotationFromBodyFixedToTopocentricFrame( ephemerisTime ).inverse( ) *
        calculateEnuDisplacement( ephemerisTime, nominalSiteState );
}


}

}
