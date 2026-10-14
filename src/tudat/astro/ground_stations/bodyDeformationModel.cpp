#include "tudat/astro/ground_stations/bodyDeformationModel.h"

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
template< typename Epoch >
Eigen::Vector3d differentiateDisplacementAtEpoch( const Epoch& time,
                                                  const double requestedStep,
                                                  const std::function< Eigen::Vector3d( const Epoch& ) >& displacement,
                                                  const std::function< void( const Epoch&, const Epoch& ) >& validateInterval )
{
    const auto derivative = [ & ]( double requested ) -> Eigen::Vector3d {
        // Use the actual representable first-point spacing in the denominator.
        const double step = static_cast< double >( ( time + requested ) - time );
        const std::array< Epoch, 4 > epochs = { time - 2.0 * step, time - step, time + step, time + 2.0 * step };
        if( !std::isfinite( static_cast< double >( time ) ) || !std::isfinite( step ) || step <= 0.0 ||
            !( epochs[ 0 ] < epochs[ 1 ] && epochs[ 1 ] < time && time < epochs[ 2 ] && epochs[ 2 ] < epochs[ 3 ] ) )
        {
            throw std::invalid_argument( "Displacement derivative step does not give finite distinct stencil epochs." );
        }
        // Floating epochs must support an approximately uniform stencil. Do not silently
        // differentiate across a sampling grid too coarse to represent the selected step.
        const double spacingTolerance = 16.0 * std::numeric_limits< double >::epsilon( ) * step;
        for( std::size_t i = 0; i < epochs.size( ); ++i )
        {
            const double expectedOffset = ( i < 2 ? static_cast< double >( i ) - 2.0 : static_cast< double >( i ) - 1.0 ) * step;
            if( !std::isfinite( static_cast< double >( epochs[ i ] ) ) ||
                std::fabs( static_cast< double >( epochs[ i ] - time ) - expectedOffset ) > spacingTolerance )
            {
                throw std::invalid_argument( "Displacement derivative epochs cannot represent a uniform stencil." );
            }
        }

        if( validateInterval ) validateInterval( epochs.front( ), epochs.back( ) );
        std::array< Eigen::Vector3d, 4 > values;
        for( std::size_t i = 0; i < epochs.size( ); ++i )
        {
            values[ i ] = displacement( epochs[ i ] );
            if( !values[ i ].allFinite( ) )
            {
                throw std::runtime_error( "Non-finite displacement in derivative stencil." );
            }
        }
        const Eigen::Vector3d result = ( 8.0 * ( values[ 2 ] - values[ 1 ] ) - ( values[ 3 ] - values[ 0 ] ) ) / ( 12.0 * step );
        if( !result.allFinite( ) )
        {
            throw std::runtime_error( "Non-finite displacement derivative." );
        }
        return result;
    };
    Eigen::Vector3d velocity;
    try
    {
        velocity = derivative( requestedStep );
        if constexpr( std::is_same< Epoch, Time >::value )
        {
            const double scalarCentre = static_cast< double >( time );
            const double resolutionStep = 4 *
                    std::max( std::nextafter( scalarCentre, INFINITY ) - scalarCentre,
                              scalarCentre - std::nextafter( scalarCentre, -INFINITY ) );
            if( requestedStep < resolutionStep )
            {
                const Eigen::Vector3d resolved = derivative( resolutionStep );
                if( ( resolved - velocity ).norm( ) > 1e-12 + 1e-6 * std::max( resolved.norm( ), velocity.norm( ) ) )
                    throw std::invalid_argument( "Displacement derivative step cannot resolve the supplied callback." );
            }
        }
    }
    catch( ... )
    {
        const std::exception_ptr original = std::current_exception( );
        try
        {
            displacement( time );
        }
        catch( ... )
        {
            // Restoration is best effort; preserve the original dynamic exception type.
        }
        std::rethrow_exception( original );
    }
    const Eigen::Vector3d central = displacement( time );
    if( !central.allFinite( ) )
    {
        throw std::runtime_error( "Non-finite central displacement after derivative evaluation." );
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
    if( !std::isfinite( static_cast< double >( tdb ) ) ) throw std::invalid_argument( "Body deformation requires a finite TDB epoch." );
    if( !station ) throw std::invalid_argument( "Body deformation requires a station." );
    const Eigen::Vector3d displacement = evaluateDisplacement( tdb, station );
    if( !displacement.allFinite( ) ) throw std::runtime_error( "Non-finite body deformation displacement." );
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
    if( !std::isfinite( static_cast< double >( tdb ) ) ) throw std::invalid_argument( "Body deformation requires a finite TDB epoch." );
    if( !station ) throw std::invalid_argument( "Body deformation requires a station." );
    const Eigen::Vector6d state = evaluateDisplacementState( tdb, station );
    if( !state.allFinite( ) ) throw std::runtime_error( "Non-finite body deformation state." );
    return state;
}

Eigen::Vector6d BodyDeformationModel::evaluateDisplacementState( const Time& tdb,
                                                                 const std::shared_ptr< ground_stations::GroundStationState > station )
{
    Eigen::Vector6d state;
    state.tail< 3 >( ) = differentiateDisplacement(
            tdb,
            [ & ]( const Time& time ) { return calculateDisplacement( time, station ); },
            [ & ]( const Time& lower, const Time& upper ) { validateDisplacementDerivativeInterval( lower, upper ); } );
    state.head< 3 >( ) = calculateDisplacement( tdb, station );
    return state;
}

Eigen::Vector3d BodyDeformationModel::differentiateDisplacement(
        const Time& time,
        const std::function< Eigen::Vector3d( const Time& ) >& displacement,
        const std::function< void( const Time&, const Time& ) >& validateInterval ) const
{
    return differentiateDisplacementAtEpoch< Time >( time, displacementDerivativeTimeStep_, displacement, validateInterval );
}

Eigen::Vector3d BodyDeformationModel::differentiateDisplacement(
        double time,
        const std::function< Eigen::Vector3d( const double& ) >& displacement,
        const std::function< void( const double&, const double& ) >& validateInterval ) const
{
    return differentiateDisplacementAtEpoch< double >( time, displacementDerivativeTimeStep_, displacement, validateInterval );
}

void BodyDeformationModel::setDisplacementDerivativeTimeStep( const double timeStep )
{
    if( !std::isfinite( timeStep ) || timeStep <= 0 )
        throw std::invalid_argument( "Displacement derivative time step must be finite and positive." );
    displacementDerivativeTimeStep_ = timeStep;
}

}  // namespace basic_astrodynamics
}  // namespace tudat
