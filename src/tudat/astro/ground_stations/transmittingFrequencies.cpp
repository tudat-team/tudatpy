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
long double StationFrequencyInterpolator::getTemplatedCurrentFrequency( const double& lookupTime )
{
    return getCurrentLongFrequency( lookupTime );
}

template<>
long double StationFrequencyInterpolator::getTemplatedCurrentFrequency( const Time& lookupTime )
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
long double StationFrequencyInterpolator::getTemplatedFrequencyIntegral( const double& quadratureStartTime,
                                                                         const double& quadratureEndTime )
{
    return getLongFrequencyIntegral( quadratureStartTime, quadratureEndTime );
}

template<>
long double StationFrequencyInterpolator::getTemplatedFrequencyIntegral( const Time& quadratureStartTime, const Time& quadratureEndTime )
{
    return getLongFrequencyIntegral( quadratureStartTime, quadratureEndTime );
}

void PiecewiseLinearFrequencyInterpolator::initialize( )
{
    invalidTimeBlocksStartTimes_.clear( );
    invalidTimeBlocksEndTimes_.clear( );

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

void MultiBandFrequencyInterpolator::addBandInterpolator( const basic_astrodynamics::FrequencyBands frequencyBand,
                                                          const std::shared_ptr< StationFrequencyInterpolator > bandInterpolator )
{
    if( bandInterpolator == nullptr )
    {
        throw std::runtime_error( "Error when adding frequency interpolator for " +
                                  basic_astrodynamics::getFrequencyBandString( frequencyBand ) +
                                  " to multi-band frequency interpolator: interpolator is null." );
    }

    if( std::dynamic_pointer_cast< MultiBandFrequencyInterpolator >( bandInterpolator ) != nullptr )
    {
        throw std::runtime_error( "Error when adding frequency interpolator for " +
                                  basic_astrodynamics::getFrequencyBandString( frequencyBand ) +
                                  " to multi-band frequency interpolator: interpolator to add is itself a multi-band interpolator." );
    }

    if( bandInterpolators_.count( frequencyBand ) == 0 )
    {
        bandInterpolators_[ frequencyBand ] = bandInterpolator;
    }
    else
    {
        std::shared_ptr< PiecewiseLinearFrequencyInterpolator > existingPiecewiseLinearInterpolator =
                std::dynamic_pointer_cast< PiecewiseLinearFrequencyInterpolator >( bandInterpolators_.at( frequencyBand ) );
        std::shared_ptr< PiecewiseLinearFrequencyInterpolator > piecewiseLinearInterpolatorToAdd =
                std::dynamic_pointer_cast< PiecewiseLinearFrequencyInterpolator >( bandInterpolator );
        if( existingPiecewiseLinearInterpolator == nullptr || piecewiseLinearInterpolatorToAdd == nullptr )
        {
            throw std::runtime_error( "Error when adding frequency interpolator for " +
                                      basic_astrodynamics::getFrequencyBandString( frequencyBand ) +
                                      " to multi-band frequency interpolator: an interpolator is already defined for this band, and "
                                      "only piecewise linear interpolators can be merged." );
        }
        existingPiecewiseLinearInterpolator->addFrequencyInterpolator( piecewiseLinearInterpolatorToAdd );
    }
}

std::shared_ptr< StationFrequencyInterpolator > MultiBandFrequencyInterpolator::getBandInterpolator(
        const basic_astrodynamics::FrequencyBands frequencyBand )
{
    auto bandIterator = bandInterpolators_.find( frequencyBand );
    if( bandIterator == bandInterpolators_.end( ) )
    {
        std::string availableBands;
        for( auto const& [ availableBand, bandInterpolator ] : bandInterpolators_ )
        {
            availableBands += ( availableBands.empty( ) ? "" : ", " ) + basic_astrodynamics::getFrequencyBandString( availableBand );
        }
        throw std::runtime_error( "Error when retrieving transmitted frequency: no frequency interpolator defined for " +
                                  basic_astrodynamics::getFrequencyBandString( frequencyBand ) +
                                  ", available bands are: " + availableBands + "." );
    }
    return bandIterator->second;
}

std::vector< basic_astrodynamics::FrequencyBands > MultiBandFrequencyInterpolator::getAvailableBands( ) const
{
    std::vector< basic_astrodynamics::FrequencyBands > availableBands;
    for( auto const& [ frequencyBand, bandInterpolator ] : bandInterpolators_ )
    {
        availableBands.push_back( frequencyBand );
    }
    return availableBands;
}

std::shared_ptr< StationFrequencyInterpolator > MultiBandFrequencyInterpolator::getSingleBandInterpolator( )
{
    if( bandInterpolators_.size( ) != 1 )
    {
        throw std::runtime_error(
                "Error when retrieving transmitted frequency: no frequency band was specified, but multi-band frequency "
                "interpolator contains " +
                std::to_string( bandInterpolators_.size( ) ) +
                " bands. Provide the frequency bands of the observation in the ancillary settings." );
    }
    return bandInterpolators_.begin( )->second;
}

double MultiBandFrequencyInterpolator::getCurrentFrequency( const double lookupTime )
{
    return getSingleBandInterpolator( )->getTemplatedCurrentFrequency< double, double >( lookupTime );
}

double MultiBandFrequencyInterpolator::getCurrentFrequency( const Time& lookupTime )
{
    return getSingleBandInterpolator( )->getTemplatedCurrentFrequency< double, Time >( lookupTime );
}

long double MultiBandFrequencyInterpolator::getCurrentLongFrequency( const double lookupTime )
{
    return getSingleBandInterpolator( )->getTemplatedCurrentFrequency< long double, double >( lookupTime );
}

long double MultiBandFrequencyInterpolator::getCurrentLongFrequency( const Time& lookupTime )
{
    return getSingleBandInterpolator( )->getTemplatedCurrentFrequency< long double, Time >( lookupTime );
}

double MultiBandFrequencyInterpolator::getFrequencyIntegral( const double quadratureStartTime, const double quadratureEndTime )
{
    return getSingleBandInterpolator( )->getTemplatedFrequencyIntegral< double, double >( quadratureStartTime, quadratureEndTime );
}

double MultiBandFrequencyInterpolator::getFrequencyIntegral( const Time& quadratureStartTime, const Time& quadratureEndTime )
{
    return getSingleBandInterpolator( )->getTemplatedFrequencyIntegral< double, Time >( quadratureStartTime, quadratureEndTime );
}

long double MultiBandFrequencyInterpolator::getLongFrequencyIntegral( const double quadratureStartTime, const double quadratureEndTime )
{
    return getSingleBandInterpolator( )->getTemplatedFrequencyIntegral< long double, double >( quadratureStartTime, quadratureEndTime );
}

long double MultiBandFrequencyInterpolator::getLongFrequencyIntegral( const Time& quadratureStartTime, const Time& quadratureEndTime )
{
    return getSingleBandInterpolator( )->getTemplatedFrequencyIntegral< long double, Time >( quadratureStartTime, quadratureEndTime );
}

std::shared_ptr< StationFrequencyInterpolator > combineTransmittingFrequencyCalculators(
        const std::shared_ptr< StationFrequencyInterpolator > existingCalculator,
        const std::shared_ptr< StationFrequencyInterpolator > calculatorToAdd,
        const std::optional< basic_astrodynamics::FrequencyBands >& frequencyBand,
        const std::string& errorMessageContext )
{
    if( calculatorToAdd == nullptr )
    {
        throw std::runtime_error( "Error when adding transmitted frequency calculator to " + errorMessageContext +
                                  ": calculator is null." );
    }

    std::shared_ptr< MultiBandFrequencyInterpolator > existingMultiBandCalculator =
            std::dynamic_pointer_cast< MultiBandFrequencyInterpolator >( existingCalculator );

    // Band-specific calculator: store in (new or existing) multi-band calculator
    if( frequencyBand.has_value( ) )
    {
        if( existingCalculator == nullptr )
        {
            existingMultiBandCalculator = std::make_shared< MultiBandFrequencyInterpolator >( );
        }
        else if( existingMultiBandCalculator == nullptr )
        {
            throw std::runtime_error( "Error when adding transmitted frequency calculator for " +
                                      basic_astrodynamics::getFrequencyBandString( frequencyBand.value( ) ) + " to " + errorMessageContext +
                                      ": existing calculator is not associated with a frequency band. Band-specific and band-independent "
                                      "frequency calculators cannot be combined." );
        }
        existingMultiBandCalculator->addBandInterpolator( frequencyBand.value( ), calculatorToAdd );
        return existingMultiBandCalculator;
    }

    // Band-independent calculator
    if( existingCalculator == nullptr )
    {
        return calculatorToAdd;
    }

    std::shared_ptr< MultiBandFrequencyInterpolator > multiBandCalculatorToAdd =
            std::dynamic_pointer_cast< MultiBandFrequencyInterpolator >( calculatorToAdd );
    if( existingMultiBandCalculator != nullptr && multiBandCalculatorToAdd != nullptr )
    {
        for( auto const& [ bandToAdd, bandInterpolatorToAdd ] : multiBandCalculatorToAdd->getBandInterpolators( ) )
        {
            existingMultiBandCalculator->addBandInterpolator( bandToAdd, bandInterpolatorToAdd );
        }
        return existingMultiBandCalculator;
    }
    else if( existingMultiBandCalculator != nullptr || multiBandCalculatorToAdd != nullptr )
    {
        throw std::runtime_error( "Error when adding transmitted frequency calculator to " + errorMessageContext +
                                  ": band-specific and band-independent frequency calculators cannot be combined." );
    }

    std::shared_ptr< PiecewiseLinearFrequencyInterpolator > existingPiecewiseLinearCalculator =
            std::dynamic_pointer_cast< PiecewiseLinearFrequencyInterpolator >( existingCalculator );
    std::shared_ptr< PiecewiseLinearFrequencyInterpolator > piecewiseLinearCalculatorToAdd =
            std::dynamic_pointer_cast< PiecewiseLinearFrequencyInterpolator >( calculatorToAdd );
    if( existingPiecewiseLinearCalculator == nullptr || piecewiseLinearCalculatorToAdd == nullptr )
    {
        throw std::runtime_error( "Error when adding transmitted frequency calculator to " + errorMessageContext +
                                  ": a calculator is already defined, and only piecewise linear calculators can be merged." );
    }
    existingPiecewiseLinearCalculator->addFrequencyInterpolator( piecewiseLinearCalculatorToAdd );
    return existingPiecewiseLinearCalculator;
}

}  // namespace ground_stations

}  // namespace tudat
