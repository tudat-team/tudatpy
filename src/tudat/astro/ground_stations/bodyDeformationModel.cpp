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
Eigen::Vector3d differentiateDisplacement( const double time,
                                           const double requestedStep,
                                           const std::function< Eigen::Vector3d( const double ) >& displacement,
                                           const std::function< void( const double, const double ) >& validateInterval )
{
    // Use the actual representable first-point spacing in the denominator.
    const double step = ( time + requestedStep ) - time;
    const std::array< double, 4 > epochs = { time - 2.0 * step, time - step, time + step, time + 2.0 * step };
    if( !std::isfinite( time ) || !std::isfinite( step ) || step <= 0.0 ||
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
        if( !std::isfinite( epochs[ i ] ) || std::fabs( ( epochs[ i ] - time ) - expectedOffset ) > spacingTolerance )
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
    return differentiateDisplacement(
            time,
            displacementDerivativeTimeStep_,
            [ & ]( const double epoch ) { return calculateDisplacement( epoch, bodyFixedPosition ); },
            [ & ]( const double lower, const double upper ) { validateDisplacementDerivativeInterval( time, lower, upper ); } );
}

Eigen::Vector3d BodyDeformationModel::calculateDisplacementVelocity(
        const double time,
        const std::shared_ptr< ground_stations::GroundStationState > stationState )
{
    return differentiateDisplacement(
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
