/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/ground_stations/transmittingFrequencies.h"

namespace tudat
{

namespace ground_stations
{

template<>
double StationFrequencyInterpolator::getTemplatedCurrentFrequency( const double& lookupTime )
{
    return getCurrentFrequency( lookupTime );
}

template<>
double StationFrequencyInterpolator::getTemplatedCurrentFrequency( const Time& lookupTime )
{
    return getCurrentFrequency( lookupTime );
}

template<>
HighPrecisionStateScalar StationFrequencyInterpolator::getTemplatedCurrentFrequency( const double& lookupTime )
{
    return getCurrentLongFrequency( lookupTime );
}

template<>
HighPrecisionStateScalar StationFrequencyInterpolator::getTemplatedCurrentFrequency( const Time& lookupTime )
{
    return getCurrentLongFrequency( lookupTime );
}

template<>
double StationFrequencyInterpolator::getTemplatedFrequencyIntegral( const double& quadratureStartTime, const double& quadratureEndTime )
{
    return getFrequencyIntegral( quadratureStartTime, quadratureEndTime );
}

template<>
double StationFrequencyInterpolator::getTemplatedFrequencyIntegral( const Time& quadratureStartTime, const Time& quadratureEndTime )
{
    return getFrequencyIntegral( quadratureStartTime, quadratureEndTime );
}

template<>
HighPrecisionStateScalar StationFrequencyInterpolator::getTemplatedFrequencyIntegral( const double& quadratureStartTime,
                                                                                      const double& quadratureEndTime )
{
    return getLongFrequencyIntegral( quadratureStartTime, quadratureEndTime );
}

template<>
HighPrecisionStateScalar StationFrequencyInterpolator::getTemplatedFrequencyIntegral( const Time& quadratureStartTime,
                                                                                      const Time& quadratureEndTime )
{
    return getLongFrequencyIntegral( quadratureStartTime, quadratureEndTime );
}

void PiecewiseLinearFrequencyInterpolator::initialize( )
{
    invalidTimeBlocksStartTimes_.clear( );
    invalidTimeBlocksEndTimes_.clear( );
    invalidStartTimeLookupScheme_ = nullptr;

    // Check if dimensions of all vectors are consistent
    if( startTimes_.size( ) != endTimes_.size( ) || startTimes_.size( ) != rampRates_.size( ) ||
        startTimes_.size( ) != startFrequencies_.size( ) )
    {
        throw std::runtime_error(
                "Error when creating piecewise linear frequency interpolator: the dimensions of the specified vectors "
                "are not consistent: start times (" +
                std::to_string( startTimes_.size( ) ) + "), end times (" + std::to_string( endTimes_.size( ) ) + "), ramp rates (" +
                std::to_string( rampRates_.size( ) ) + "), start frequencies (" + std::to_string( startFrequencies_.size( ) ) + ")." );
    }

    std::map< Time, Time > endTimesMap;
    std::map< Time, double > rampRatesMap;
    std::map< Time, double > startFrequenciesMap;

    for( unsigned int i = 0; i < startTimes_.size( ); i++ )
    {
        endTimesMap[ startTimes_.at( i ) ] = endTimes_.at( i );
        rampRatesMap[ startTimes_.at( i ) ] = rampRates_.at( i );
        startFrequenciesMap[ startTimes_.at( i ) ] = startFrequencies_.at( i );
    }

    startTimes_ = utilities::createVectorFromMapKeys( endTimesMap );
    endTimes_ = utilities::createVectorFromMapValues( endTimesMap );
    rampRates_ = utilities::createVectorFromMapValues( rampRatesMap );
    startFrequencies_ = utilities::createVectorFromMapValues( startFrequenciesMap );

    for( unsigned int i = 0; i < startTimes_.size( ); i++ )
    {
        if( endTimes_.at( i ) <= startTimes_.at( i ) )
        {
            throw std::runtime_error( "Error when creating piecewise linear frequency interpolator: ramp end time (" +
                                      std::to_string( double( endTimes_.at( i ) ) ) + ") is not after ramp start time (" +
                                      std::to_string( double( startTimes_.at( i ) ) ) + ")." );
        }

        if( i > 0 )
        {
            if( startTimes_.at( i ) < endTimes_.at( i - 1 ) )
            {
                endTimes_.at( i - 1 ) = startTimes_.at( i );
            }
            else if( gapHandling_ != extrapolate_at_gaps && startTimes_.at( i ) > endTimes_.at( i - 1 ) )
            {
                invalidTimeBlocksStartTimes_.push_back( endTimes_.at( i - 1 ) );
                invalidTimeBlocksEndTimes_.push_back( startTimes_.at( i ) );
            }
        }
    }

    startTimeLookupScheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< Time > >( startTimes_ );

    if( !invalidTimeBlocksStartTimes_.empty( ) )
    {
        invalidStartTimeLookupScheme_ = std::make_shared< interpolators::BinarySearchLookupScheme< Time > >( invalidTimeBlocksStartTimes_ );
    }
}

void PiecewiseLinearFrequencyInterpolator::addFrequencyInterpolator(
        const std::shared_ptr< PiecewiseLinearFrequencyInterpolator > rampsToAdd )
{
    if( rampsToAdd->getGapHandling( ) != gapHandling_ )
    {
        std::cerr << "Warning: gap handling of ramp table to add (" << rampsToAdd->getGapHandling( )
                  << ") differs from base ramp table gap handling (" << gapHandling_ << "). Using base gap handling." << std::endl;
    }

    std::vector< Time > startTimesToAdd = rampsToAdd->getStartTimes( );
    std::vector< Time > endTimesToAdd = rampsToAdd->getEndTimes( );
    std::vector< double > rampRatesToAdd = rampsToAdd->getRampRates( );
    std::vector< double > startFrequenciesToAdd = rampsToAdd->getStartFrequencies( );

    std::map< Time, Time > endTimesMap;
    std::map< Time, double > rampRatesMap;
    std::map< Time, double > startFrequenciesMap;

    for( unsigned int i = 0; i < startTimes_.size( ); i++ )
    {
        endTimesMap[ startTimes_.at( i ) ] = endTimes_.at( i );
        rampRatesMap[ startTimes_.at( i ) ] = rampRates_.at( i );
        startFrequenciesMap[ startTimes_.at( i ) ] = startFrequencies_.at( i );
    }

    for( unsigned int i = 0; i < startTimesToAdd.size( ); i++ )
    {
        endTimesMap[ startTimesToAdd.at( i ) ] = endTimesToAdd.at( i );
        rampRatesMap[ startTimesToAdd.at( i ) ] = rampRatesToAdd.at( i );
        startFrequenciesMap[ startTimesToAdd.at( i ) ] = startFrequenciesToAdd.at( i );
    }

    startTimes_ = utilities::createVectorFromMapKeys( endTimesMap );
    endTimes_ = utilities::createVectorFromMapValues( endTimesMap );
    rampRates_ = utilities::createVectorFromMapValues( rampRatesMap );
    startFrequencies_ = utilities::createVectorFromMapValues( startFrequenciesMap );

    initialize( );
}

}  // namespace ground_stations

}  // namespace tudat
