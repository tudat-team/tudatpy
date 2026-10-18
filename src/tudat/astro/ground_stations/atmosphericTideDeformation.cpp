#include "tudat/astro/ground_stations/atmosphericTideDeformation.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include <cmath>
#include <stdexcept>

namespace tudat
{
namespace basic_astrodynamics
{

Eigen::Vector3d AtmosphericTideDeformation::evaluateDisplacement( const Time& tdb,
                                                                  const std::shared_ptr< ground_stations::GroundStationState > station )
{
    return evaluate( tdb, station, false ).head< 3 >( );
}

Eigen::Vector6d AtmosphericTideDeformation::evaluateDisplacementState(
        const Time& tdb,
        const std::shared_ptr< ground_stations::GroundStationState > station )
{
    return evaluate( tdb, station, true );
}

Eigen::Vector6d AtmosphericTideDeformation::evaluate( const Time& tdb,
                                                      const std::shared_ptr< ground_stations::GroundStationState >& station,
                                                      bool withRate )
{
    // Geocentric latitude and longitude of the station, and its coefficients.
    const Eigen::Vector3d position = station->getNominalCartesianPosition( );
    const double longitude = std::atan2( position.y( ), position.x( ) );
    const double latitude = std::atan2( position.z( ), std::hypot( position.x( ), position.y( ) ) );
    input_output::AtmosphericTideCoefficients coefficients;
    if( grid_.coefficientsMetres.empty( ) )
    {
        const auto entry = stationCoefficients_.find( station->getSiteId( ) );
        if( entry == stationCoefficients_.end( ) )
        {
            throw std::out_of_range( "Atmospheric tide station not found: " + station->getSiteId( ) );
        }
        coefficients = entry->second;
    }
    else
    {
        coefficients = interpolateAtmosphericTideGrid(
                grid_, latitude * 180 / mathematical_constants::PI, longitude * 180 / mathematical_constants::PI );
    }

    // UT1 phase, and (for the rate) the derivative of UT1 with respect to TDB.
    const Time tai = timeConverter_->getTaiTime( tdb_scale, tdb );
    double ut1, ut1RateWrtTdb = 1;
    if( withRate )
    {
        const auto ut1WithRate = timeConverter_->getUt1FromTaiWithRate( tai );
        ut1 = static_cast< double >( ut1WithRate.ut1 );
        ut1RateWrtTdb = ut1WithRate.dUt1DTai * timeConverter_->getGeocentricTtRateWithRespectToTdb( tdb );
    }
    else
    {
        ut1 = static_cast< double >( timeConverter_->getUt1FromTai( tai ) );
    }

    // Rotate the radial, north and east components (as east, north, up) to body-fixed axes.
    const auto rotation = reference_frames::getEnuLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion( longitude, latitude );
    const Eigen::Vector3d local = synthesizeAtmosphericTideLoading( coefficients, ut1 );
    Eigen::Vector6d state = Eigen::Vector6d::Zero( );
    state.head< 3 >( ) = rotation * Eigen::Vector3d( local.z( ), local.y( ), local.x( ) );
    if( withRate )
    {
        const Eigen::Vector3d localRate = differentiateAtmosphericTideLoadingWrtUt1( coefficients, ut1 );
        state.tail< 3 >( ) = rotation * Eigen::Vector3d( localRate.z( ), localRate.y( ), localRate.x( ) );
    }
    if( convertToCentreOfMass_ )
    {
        state.head< 3 >( ) += calculateIers2010Rp03AtmosphericGeocentreTranslation( ut1 );
        if( withRate )
        {
            state.tail< 3 >( ) += differentiateIers2010Rp03AtmosphericGeocentreTranslationWrtUt1( ut1 );
        }
    }
    state.tail< 3 >( ) *= ut1RateWrtTdb;
    return state;
}

}  // namespace basic_astrodynamics
}  // namespace tudat
