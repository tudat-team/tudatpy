#include "tudat/astro/ground_stations/bodyDeformationModel.h"
#include "tudat/math/basic/numericalDerivative.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>

namespace tudat
{
namespace basic_astrodynamics
{
namespace
{
//! Displacement velocity from a fourth-order central difference.
template< typename Epoch >
Eigen::Vector3d differentiateDisplacementAtEpoch( const Epoch& time,
                                                  const double requestedStep,
                                                  const std::function< Eigen::Vector3d( const Epoch& ) >& displacement )
{
    // Use the actually representable step, which differs from the requested step only for double epochs.
    const double step = static_cast< double >( ( time + requestedStep ) - time );
    const Eigen::Vector3d velocity = numerical_derivatives::computeFourthOrderCentralDifference( displacement, time, step );
    if( !velocity.allFinite( ) )
    {
        throw std::runtime_error( "Non-finite displacement derivative." );
    }
    return velocity;
}
}  // namespace

BodyDeformationModel::BodyDeformationModel( const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter ):
    timeConverter_( converter ? converter : getGeocentricTimeConverter( ) )
{}

Eigen::Vector3d BodyDeformationModel::calculateDisplacement( const Time& tdb,
                                                             const std::shared_ptr< ground_stations::GroundStationState > station )
{
    if( !std::isfinite( static_cast< double >( tdb ) ) )
    {
        throw std::invalid_argument( "Body deformation requires a finite TDB epoch." );
    }
    if( !station )
    {
        throw std::invalid_argument( "Body deformation requires a station." );
    }
    const Eigen::Vector3d displacement = evaluateDisplacement( tdb, station );
    if( !displacement.allFinite( ) )
    {
        throw std::runtime_error( "Non-finite body deformation displacement." );
    }
    return displacement;
}

Eigen::Vector3d BodyDeformationModel::calculateDisplacement( const Time& tdb, const Eigen::Vector3d& bodyFixedPosition )
{
    return calculateDisplacement(
            tdb, std::make_shared< ground_stations::GroundStationState >( bodyFixedPosition, coordinate_conversions::cartesian_position ) );
}

Eigen::Vector6d BodyDeformationModel::calculateDisplacementState( const Time& tdb,
                                                                  const std::shared_ptr< ground_stations::GroundStationState > station )
{
    if( !std::isfinite( static_cast< double >( tdb ) ) )
    {
        throw std::invalid_argument( "Body deformation requires a finite TDB epoch." );
    }
    if( !station )
    {
        throw std::invalid_argument( "Body deformation requires a station." );
    }
    const Eigen::Vector6d state = evaluateDisplacementState( tdb, station );
    if( !state.allFinite( ) )
    {
        throw std::runtime_error( "Non-finite body deformation state." );
    }
    return state;
}

Eigen::Vector6d BodyDeformationModel::evaluateDisplacementState( const Time& tdb,
                                                                 const std::shared_ptr< ground_stations::GroundStationState > station )
{
    Eigen::Vector6d state;
    state.tail< 3 >( ) = differentiateDisplacement( tdb, [ & ]( const Time& time ) { return calculateDisplacement( time, station ); } );
    state.head< 3 >( ) = calculateDisplacement( tdb, station );
    return state;
}

Eigen::Vector3d BodyDeformationModel::differentiateDisplacement( const Time& time,
                                                                 const std::function< Eigen::Vector3d( const Time& ) >& displacement ) const
{
    return differentiateDisplacementAtEpoch< Time >( time, displacementDerivativeTimeStep_, displacement );
}

Eigen::Vector3d BodyDeformationModel::differentiateDisplacement(
        double time,
        const std::function< Eigen::Vector3d( const double& ) >& displacement ) const
{
    return differentiateDisplacementAtEpoch< double >( time, displacementDerivativeTimeStep_, displacement );
}

void BodyDeformationModel::setDisplacementDerivativeTimeStep( const double timeStep )
{
    if( !std::isfinite( timeStep ) || timeStep <= 0 )
    {
        throw std::invalid_argument( "Displacement derivative time step must be finite and positive." );
    }
    displacementDerivativeTimeStep_ = timeStep;
}

}  // namespace basic_astrodynamics
}  // namespace tudat
