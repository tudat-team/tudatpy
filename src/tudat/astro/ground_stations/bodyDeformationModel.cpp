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
Eigen::Vector3d differentiateDisplacement( const Epoch& time,
                                           const double requestedStep,
                                           const std::function< Eigen::Vector3d( const Epoch& ) >& displacement,
                                           const std::function< void( const Epoch&, const Epoch& ) >& validateInterval )
{
    // Use the actual representable first-point spacing in the denominator.
    const double step = static_cast< double >( ( time + requestedStep ) - time );
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

    validateInterval( epochs.front( ), epochs.back( ) );

    Eigen::Vector3d velocity;
    try
    {
        std::array< Eigen::Vector3d, 4 > values;
        for( std::size_t i = 0; i < epochs.size( ); ++i )
        {
            values[ i ] = displacement( epochs[ i ] );
            if( !values[ i ].allFinite( ) )
            {
                throw std::runtime_error( "Non-finite displacement in derivative stencil." );
            }
        }
        velocity = ( 8.0 * ( values[ 2 ] - values[ 1 ] ) - ( values[ 3 ] - values[ 0 ] ) ) / ( 12.0 * step );
        if( !velocity.allFinite( ) )
        {
            throw std::runtime_error( "Non-finite displacement derivative." );
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

Eigen::Vector3d BodyDeformationModel::calculateDisplacement( const EarthDeformationEpoch& epoch, const Eigen::Vector3d& position )
{
    return calculateDisplacement( epoch.getTdbSeconds( ), position );
}

Eigen::Vector3d BodyDeformationModel::calculateDisplacement( const EarthDeformationEpoch& epoch,
                                                             const std::shared_ptr< ground_stations::GroundStationState > station )
{
    return calculateDisplacement( epoch.getTdbSeconds( ), station );
}

Eigen::Vector3d BodyDeformationModel::calculateDisplacementVelocity( const EarthDeformationEpoch& epoch, const Eigen::Vector3d& position )
{
    return calculateDisplacementVelocity( epoch.getTdbSeconds( ), position );
}

Eigen::Vector3d BodyDeformationModel::calculateDisplacementVelocity( const EarthDeformationEpoch& epoch,
                                                                     const std::shared_ptr< ground_stations::GroundStationState > station )
{
    return calculateDisplacementVelocity( epoch.getTdbSeconds( ), station );
}

Eigen::Vector6d BodyDeformationModel::calculateDisplacementState( const EarthDeformationEpoch& epoch, const Eigen::Vector3d& position )
{
    return calculateDisplacementState( epoch.getTdbSeconds( ), position );
}

Eigen::Vector6d BodyDeformationModel::calculateDisplacementState( const EarthDeformationEpoch& epoch,
                                                                  const std::shared_ptr< ground_stations::GroundStationState > station )
{
    return calculateDisplacementState( epoch.getTdbSeconds( ), station );
}

namespace
{
template< typename Position >
Eigen::Vector6d nativeState( BodyDeformationModel& model, const EarthDeformationEpoch& epoch, const Position& position )
{
    Eigen::Vector6d state;
    try
    {
        state.tail< 3 >( ) = model.calculateDisplacementVelocity( epoch, position );
    }
    catch( ... )
    {
        const auto original = std::current_exception( );
        try
        {
            model.calculateDisplacement( epoch, position );
        }
        catch( ... )
        {}
        std::rethrow_exception( original );
    }
    state.head< 3 >( ) = model.calculateDisplacement( epoch, position );
    if( !state.allFinite( ) ) throw std::runtime_error( "Non-finite native deformation state" );
    return state;
}
}  // namespace

Eigen::Vector6d BodyDeformationModel::calculateNativeDisplacementState( const EarthDeformationEpoch& epoch,
                                                                        const Eigen::Vector3d& position )
{
    return nativeState( *this, epoch, position );
}
Eigen::Vector6d BodyDeformationModel::calculateNativeDisplacementState(
        const EarthDeformationEpoch& epoch,
        const std::shared_ptr< ground_stations::GroundStationState > station )
{
    return nativeState( *this, epoch, station );
}

Eigen::Vector3d BodyDeformationModel::calculateTerrestrialDisplacementVelocity( const EarthDeformationEpoch& epoch,
                                                                                const Eigen::Vector3d& position )
{
    const double rate = epoch.getTtRateWithRespectToTdb( );
    return calculateDisplacementVelocity( static_cast< double >( epoch.getTt( ) ), position ) * rate;
}
Eigen::Vector3d BodyDeformationModel::calculateTerrestrialDisplacementVelocity(
        const EarthDeformationEpoch& epoch,
        const std::shared_ptr< ground_stations::GroundStationState > station )
{
    const double rate = epoch.getTtRateWithRespectToTdb( );
    return calculateDisplacementVelocity( static_cast< double >( epoch.getTt( ) ), station ) * rate;
}
Eigen::Vector3d BodyDeformationModel::calculateRoutedDisplacementVelocity( const EarthDeformationEpoch& epoch,
                                                                           const Eigen::Vector3d& position )
{
    return differentiateDisplacement< Time >(
            epoch.getTdb( ),
            displacementDerivativeTimeStep_,
            [ & ]( const Time& time ) { return calculateDisplacement( epoch.withTdb( time ), position ); },
            [ & ]( const Time& lower, const Time& upper ) {
                validateContextDerivativeInterval( epoch, epoch.withTdb( lower ), epoch.withTdb( upper ) );
            } );
}
Eigen::Vector3d BodyDeformationModel::calculateRoutedDisplacementVelocity(
        const EarthDeformationEpoch& epoch,
        const std::shared_ptr< ground_stations::GroundStationState > station )
{
    return differentiateDisplacement< Time >(
            epoch.getTdb( ),
            displacementDerivativeTimeStep_,
            [ & ]( const Time& time ) { return calculateDisplacement( epoch.withTdb( time ), station ); },
            [ & ]( const Time& lower, const Time& upper ) {
                validateContextDerivativeInterval( epoch, epoch.withTdb( lower ), epoch.withTdb( upper ) );
            } );
}

void BodyDeformationModel::setDisplacementDerivativeTimeStep( const double timeStep )
{
    if( !std::isfinite( timeStep ) || timeStep <= 0.0 )
    {
        throw std::invalid_argument( "Displacement derivative time step must be finite and positive." );
    }
    displacementDerivativeTimeStep_ = timeStep;
}

Eigen::Vector3d BodyDeformationModel::calculateDisplacementVelocity( const double time, const Eigen::Vector3d& bodyFixedPosition )
{
    return differentiateDisplacement< double >(
            time,
            displacementDerivativeTimeStep_,
            [ & ]( const double epoch ) { return calculateDisplacement( epoch, bodyFixedPosition ); },
            [ & ]( const double lower, const double upper ) { validateDisplacementDerivativeInterval( time, lower, upper ); } );
}

Eigen::Vector3d BodyDeformationModel::calculateDisplacementVelocity(
        const double time,
        const std::shared_ptr< ground_stations::GroundStationState > stationState )
{
    return differentiateDisplacement< double >(
            time,
            displacementDerivativeTimeStep_,
            [ & ]( const double epoch ) { return calculateDisplacement( epoch, stationState ); },
            [ & ]( const double lower, const double upper ) { validateDisplacementDerivativeInterval( time, lower, upper ); } );
}

Eigen::Vector6d BodyDeformationModel::calculateDisplacementState( const double time, const Eigen::Vector3d& bodyFixedPosition )
{
    Eigen::Vector6d state;
    state.head< 3 >( ) = calculateDisplacement( time, bodyFixedPosition );
    state.tail< 3 >( ) = calculateDisplacementVelocity( time, bodyFixedPosition );
    if( !state.allFinite( ) )
    {
        throw std::runtime_error( "Non-finite displacement state." );
    }
    return state;
}

Eigen::Vector6d BodyDeformationModel::calculateDisplacementState(
        const double time,
        const std::shared_ptr< ground_stations::GroundStationState > stationState )
{
    Eigen::Vector6d state;
    state.head< 3 >( ) = calculateDisplacement( time, stationState );
    state.tail< 3 >( ) = calculateDisplacementVelocity( time, stationState );
    if( !state.allFinite( ) )
    {
        throw std::runtime_error( "Non-finite displacement state." );
    }
    return state;
}

}  // namespace basic_astrodynamics
}  // namespace tudat
