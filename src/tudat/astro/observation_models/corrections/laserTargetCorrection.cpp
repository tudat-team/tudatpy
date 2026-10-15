/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/astro/observation_models/corrections/laserTargetCorrection.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace tudat
{
namespace observation_models
{

//! Store the published additive observation correction with the opposite computed-surface sign.
LaserTargetSignature::LaserTargetSignature( const double centerOfMassDistance ):
    centerOfMassDistance_( centerOfMassDistance ), internalOpticalPathLength_( 0.0 )
{
    if( !std::isfinite( centerOfMassDistance_ ) || centerOfMassDistance_ < 0.0 )
    {
        throw std::invalid_argument( "Laser target CoM distance must be finite and nonnegative" );
    }
}

//! Store explicit optical points and validate their fixed relative return weights.
LaserTargetSignature::LaserTargetSignature( std::vector< Eigen::Vector3d > opticalPointOffsets,
                                            std::vector< double > returnWeights,
                                            const double internalOpticalPathLength,
                                            ReturnWeightFunction returnWeightFunction ):
    centerOfMassDistance_( 0.0 ), internalOpticalPathLength_( internalOpticalPathLength ),
    opticalPointOffsets_( std::move( opticalPointOffsets ) ), returnWeights_( std::move( returnWeights ) ),
    returnWeightFunction_( std::move( returnWeightFunction ) )
{
    if( !std::isfinite( internalOpticalPathLength_ ) || internalOpticalPathLength_ < 0.0 )
    {
        throw std::invalid_argument( "Laser internal optical path must be finite and nonnegative" );
    }
    if( opticalPointOffsets_.empty( ) || opticalPointOffsets_.size( ) != returnWeights_.size( ) )
    {
        throw std::invalid_argument( "Laser optical points require one return weight per point" );
    }
    for( const auto& point : opticalPointOffsets_ )
    {
        if( !point.allFinite( ) )
        {
            throw std::invalid_argument( "Laser optical point coordinates must be finite" );
        }
    }
    bool hasPositiveWeight = false;

    for( const double weight : returnWeights_ )
    {
        if( !std::isfinite( weight ) || weight < 0.0 )
        {
            throw std::invalid_argument( "Laser return weights must be finite and nonnegative" );
        }
        hasPositiveWeight = hasPositiveWeight || weight > 0.0;
    }
    if( !hasPositiveWeight )
    {
        throw std::invalid_argument( "Laser return weights must include a positive weight" );
    }
}

//! Evaluate a shared weighted return distribution without subtracting nearly equal geometric ranges.
LaserTargetResponse LaserTargetSignature::calculateResponse( const Eigen::Vector3d& targetToTransmitter,
                                                             const Eigen::Vector3d& targetToReceiver,
                                                             const double reflectionTime ) const
{
    const double uplinkRange = targetToTransmitter.norm( );

    const double downlinkRange = targetToReceiver.norm( );

    if( !targetToTransmitter.allFinite( ) || !targetToReceiver.allFinite( ) || !std::isfinite( uplinkRange ) ||
        !std::isfinite( downlinkRange ) || uplinkRange <= 0.0 || downlinkRange <= 0.0 || !std::isfinite( reflectionTime ) )
    {
        throw std::invalid_argument( "Laser target response requires finite separated link ends and reflection time" );
    }
    if( opticalPointOffsets_.empty( ) )
    {
        return { Eigen::Vector2d::Constant( -centerOfMassDistance_ ), 0.0 };
    }

    // Both legs see the same set of return weights, evaluated with their distinct retarded directions.
    std::vector< double > weights = returnWeightFunction_
            ? returnWeightFunction_( reflectionTime, targetToTransmitter / uplinkRange, targetToReceiver / downlinkRange )
            : returnWeights_;

    if( weights.size( ) != opticalPointOffsets_.size( ) )
    {
        throw std::runtime_error( "Laser directional response returned the wrong number of weights" );
    }
    double maximumWeight = 0.0;

    for( const double weight : weights )
    {
        if( !std::isfinite( weight ) || weight < 0.0 )
        {
            throw std::runtime_error( "Laser directional return weights must be finite and nonnegative" );
        }
        maximumWeight = std::max( maximumWeight, weight );
    }
    if( maximumWeight <= 0.0 )
    {
        throw std::runtime_error( "Laser directional response has no positive return weight" );
    }

    // Rationalizing |r-p|-|r| preserves small optical offsets even at lunar or larger ranges.
    std::vector< Eigen::Vector2d > pointPathCorrections( opticalPointOffsets_.size( ) );

    Eigen::Vector2d weightedMean = Eigen::Vector2d::Zero( );

    double weightSum = 0.0;

    for( unsigned int i = 0; i < opticalPointOffsets_.size( ); ++i )
    {
        const Eigen::Vector3d& point = opticalPointOffsets_[ i ];

        pointPathCorrections[ i ].x( ) = ( point.squaredNorm( ) - 2.0 * targetToTransmitter.dot( point ) ) /
                ( ( targetToTransmitter - point ).norm( ) + uplinkRange );
        pointPathCorrections[ i ].y( ) =
                ( point.squaredNorm( ) - 2.0 * targetToReceiver.dot( point ) ) / ( ( targetToReceiver - point ).norm( ) + downlinkRange );
        weights[ i ] /= maximumWeight;
        weightedMean += weights[ i ] * pointPathCorrections[ i ];
        weightSum += weights[ i ];
    }
    weightedMean /= weightSum;

    // Width describes individual half-round-trip returns about their centroid; it is never added to the mean.
    double variance = 0.0;

    for( unsigned int i = 0; i < pointPathCorrections.size( ); ++i )
    {
        const double centeredRange = 0.5 * ( pointPathCorrections[ i ] - weightedMean ).sum( );

        variance += weights[ i ] * centeredRange * centeredRange;
    }

    // A common internal optical path changes the centroid, while the geometric return spread stays the same.
    weightedMean.array( ) += internalOpticalPathLength_;
    if( !weightedMean.allFinite( ) || !std::isfinite( variance ) )
    {
        throw std::runtime_error( "Non-finite laser target response" );
    }
    return { weightedMean, std::sqrt( variance / weightSum ) };
}

//! Bind one shared response definition to the target's orientation.
LaserTargetLightTimeCorrection::LaserTargetLightTimeCorrection( std::shared_ptr< LaserTargetSignature > signature,
                                                                std::function< Eigen::Quaterniond( double ) > rotationToTargetFrame ):
    LightTimeCorrection( laser_target_signature ), signature_( std::move( signature ) ),
    rotationToTargetFrame_( std::move( rotationToTargetFrame ) )
{
    if( !signature_ || ( signature_->requiresTargetRotation( ) && !rotationToTargetFrame_ ) )
    {
        throw std::invalid_argument( "Laser target correction requires a signature and target rotation" );
    }
}

//! Correct only the current leg, keeping the optical point's uplink and downlink response consistent.
double LaserTargetLightTimeCorrection::calculateLightTimeCorrectionWithMultiLegLinkEndStates(
        const std::vector< Eigen::Vector6d >& linkEndsStates,
        const std::vector< double >& linkEndsTimes,
        const unsigned int currentMultiLegTransmitterIndex,
        const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings )
{
    if( linkEndsStates.size( ) != 4 || linkEndsTimes.size( ) != 4 ||
        ( currentMultiLegTransmitterIndex != 0 && currentMultiLegTransmitterIndex != 2 ) )
    {
        throw std::invalid_argument( "Laser target correction requires a station-target-station two-leg link" );
    }
    const Eigen::Quaterniond uplinkRotation =
            rotationToTargetFrame_ ? rotationToTargetFrame_( linkEndsTimes[ 1 ] ) : Eigen::Quaterniond::Identity( );

    const Eigen::Quaterniond downlinkRotation =
            rotationToTargetFrame_ ? rotationToTargetFrame_( linkEndsTimes[ 2 ] ) : Eigen::Quaterniond::Identity( );

    if( !uplinkRotation.coeffs( ).allFinite( ) || !downlinkRotation.coeffs( ).allFinite( ) ||
        std::abs( uplinkRotation.norm( ) - 1.0 ) > 1.0E-12 || std::abs( downlinkRotation.norm( ) - 1.0 ) > 1.0E-12 )
    {
        throw std::runtime_error( "Laser target rotation must be a finite unit quaternion" );
    }
    const LaserTargetResponse response =
            signature_->calculateResponse( uplinkRotation * ( linkEndsStates[ 0 ].head< 3 >( ) - linkEndsStates[ 1 ].head< 3 >( ) ),
                                           downlinkRotation * ( linkEndsStates[ 3 ].head< 3 >( ) - linkEndsStates[ 2 ].head< 3 >( ) ),
                                           linkEndsTimes[ 1 ] );

    return response.meanPathCorrections[ currentMultiLegTransmitterIndex / 2 ] / physical_constants::SPEED_OF_LIGHT;
}

}  // namespace observation_models
}  // namespace tudat
