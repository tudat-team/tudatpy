#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Astrodynamics/SiteDisplacements/poleTideDeformation.h"

namespace tudat
{

namespace site_displacements
{

//! Function to calculate site displacement due to pole tide deformation.
Eigen::Vector3d PoleTideDeformation::calculateSiteDisplacement( const double terrestrialTime,
                                                                const boost::shared_ptr< NominalGroundStationState > nominalSiteState )
{
    using namespace tudat::earth_orientation;
    using mathematical_constants::PI;

    // Retrieve position of CIP in ITRS from polar motion calculator.
    Eigen::Vector2d cipPositionInItrs =
            unit_conversions::convertDegreesToRadians< Eigen::Vector2d >( polarMotionFunction_( terrestrialTime ) ) / 3600.0;

    // Determine years since J2000 and, subsequently, the mean CIP position at the given time.
    double yearsSinceJ2000 = terrestrialTime / physical_constants::JULIAN_DAY /
            physical_constants::JULIAN_YEAR_IN_DAYS;
    Eigen::Vector2d meanCipPositionInItrs =
            calculateMeanCipPositionInItrs( yearsSinceJ2000 );

    // Calculate deviation of CIP position from mean.
    Eigen::Vector2d cipPositionDeviationInItrs =
            cipPositionInItrs - meanCipPositionInItrs;

    // Retrieve station latitude and longitude.
    double latitude = PI / 2.0 - nominalSiteState->getNominalSphericalPosition( ).y( );
    double longitude = nominalSiteState->getNominalSphericalPosition( ).z( );

    double radiansToArcSeconds = 180.0 / mathematical_constants::PI / 3600.0;

    // Calculate site displacement in local coordinates (east, north, up), IERS 2010 Conventions Eq. 7.26
    Eigen::Vector3d siteDisplacement = Eigen::Vector3d::Zero( );
    siteDisplacement.x( ) = -radiansToArcSeconds * 9.0E-3 * cos( latitude ) * ( cipPositionDeviationInItrs.x( ) * sin( longitude ) -
                                                     cipPositionDeviationInItrs.y( ) * cos( longitude ) );
    siteDisplacement.y( ) = radiansToArcSeconds * 9.0E-3 * cos( 2.0 * latitude ) * ( cipPositionDeviationInItrs.x( ) * cos( longitude ) +
                                                          cipPositionDeviationInItrs.y( ) * sin( longitude ) );
    siteDisplacement.z( ) = radiansToArcSeconds * 33.0E-3 * sin( 2.0 * latitude ) * ( cipPositionDeviationInItrs.x( ) * cos( longitude ) +
                                                           cipPositionDeviationInItrs.y( ) * sin( longitude ) );

    // Transform to ITRS (geocentric) frame
    return  nominalSiteState->getEastNorthRadialGeocentricUnitVectors( )[ 0 ] * siteDisplacement[ 0 ] +
            nominalSiteState->getEastNorthRadialGeocentricUnitVectors( )[ 1 ] * siteDisplacement[ 1 ] +
            nominalSiteState->getEastNorthRadialGeocentricUnitVectors( )[ 2 ] * siteDisplacement[ 2 ];
}


}

}
