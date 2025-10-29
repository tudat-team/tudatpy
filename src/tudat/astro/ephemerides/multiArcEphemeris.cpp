/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#include "tudat/astro/ephemerides/multiArcEphemeris.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"

namespace tudat
{

namespace ephemerides
{


void MultiArcEphemeris::resetSingleArcEphemerides(
        const std::vector< std::shared_ptr< Ephemeris > >& singleArcEphemerides,
        const std::vector< double >& arcStartTimes )
{
    if( singleArcEphemerides.size( ) != arcStartTimes.size( ) )
    {
        throw std::runtime_error( "Error when updating multi-arc ephemeris, input is inconsistent" );
    }

    singleArcEphemerides_ = singleArcEphemerides;
    arcStartTimes_ = arcStartTimes;
    arcEndTimes_.resize( arcStartTimes_.size( ) );

    for( int i = 0; i < static_cast< int > ( singleArcEphemerides_.size( ) ) - 1 ; i++ )
    {
        std::pair< double, double > safeInterval = getSafeEphemerisEvaluationInterval( singleArcEphemerides_.at( i ) );
        if( safeInterval.second < arcStartTimes_.at( i + 1 ) )
        {
            arcEndTimes_[ i ] = safeInterval.second;
        }
        else
        {
            arcEndTimes_[ i ] = arcStartTimes_.at( i + 1 );
        }
    }
    if( singleArcEphemerides_.size( ) > 0 )
    {
        arcEndTimes_[ singleArcEphemerides_.size( ) - 1 ] = std::numeric_limits< double >::max( );
    }
    // Create times at which the look up changes from one arc to the other.
    std::vector< double > arcSplitTimes = arcStartTimes_;
    arcSplitTimes.push_back( std::numeric_limits< double >::max( ) );
    lookUpscheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >( arcSplitTimes );
}

Eigen::Vector6d MultiArcEphemeris::getCartesianState( const double currentTime )
{

    std::pair< bool, int > currentArc = getCurrentEphemerisArc( currentTime );
    if( currentArc.first )
    {
        return singleArcEphemerides_.at( currentArc.second )->getCartesianState( currentTime );
    }
    else
    {
        if( defaultEphemeris_ == nullptr )
        {
            throw exceptions::MultiArcEphemerisError< double >(
                    currentTime, std::make_pair( arcStartTimes_.at( currentArc.second ), arcEndTimes_.at( currentArc.second ) ), currentArc.second );
        }
        else
        {
            return defaultEphemeris_->getCartesianState( currentTime );
        }
    }
}

Eigen::Matrix< long double, 6, 1 > MultiArcEphemeris::getCartesianLongState( const double currentTime )
{
    std::pair< bool, int > currentArc = getCurrentEphemerisArc( currentTime );
    if( currentArc.first )
    {
        return singleArcEphemerides_.at( currentArc.second )->getCartesianLongState( double( currentTime ) );
    }
    else
    {
        if( defaultEphemeris_ == nullptr )
        {
            throw exceptions::MultiArcEphemerisError< double >(
                    currentTime, std::make_pair( arcStartTimes_.at( currentArc.second ), arcEndTimes_.at( currentArc.second ) ), currentArc.second );
        }
        else
        {
            return defaultEphemeris_->getCartesianLongState( currentTime );
        }
    }
}

Eigen::Vector6d MultiArcEphemeris::getCartesianStateFromExtendedTime( const Time& currentTime )
{
    std::pair< bool, int > currentArc = getCurrentEphemerisArc( currentTime );
    if( currentArc.first )
    {
        return singleArcEphemerides_.at( currentArc.second )->getCartesianStateFromExtendedTime( currentTime );
    }
    else
    {
        if( defaultEphemeris_ == nullptr )
        {
            throw exceptions::MultiArcEphemerisError< Time >(
                    currentTime, std::make_pair( arcStartTimes_.at( currentArc.second ), arcEndTimes_.at( currentArc.second ) ), currentArc.second );
        }
        else
        {
            return defaultEphemeris_->getCartesianStateFromExtendedTime( currentTime );
        }
    }
}


Eigen::Matrix< long double, 6, 1 > MultiArcEphemeris::getCartesianLongStateFromExtendedTime( const Time& currentTime )
{
    std::pair< bool, int > currentArc = getCurrentEphemerisArc( currentTime );
    if( currentArc.first )
    {
        return singleArcEphemerides_.at( currentArc.second )->getCartesianLongStateFromExtendedTime( double( currentTime ) );
    }
    else
    {
        if( defaultEphemeris_ == nullptr )
        {
            throw exceptions::MultiArcEphemerisError< Time >(
                    currentTime, std::make_pair( arcStartTimes_.at( currentArc.second ), arcEndTimes_.at( currentArc.second ) ), currentArc.second );
        }
        else
        {
            return defaultEphemeris_->getCartesianLongStateFromExtendedTime( currentTime );
        }
    }
}


std::pair< bool, int > MultiArcEphemeris::getCurrentEphemerisArc( const double currentTime )
{
    if( singleArcEphemerides_.size( ) == 0 )
    {
        throw std::runtime_error(
                "Error when retrieving state from multi-arc ephemeris; no constituent single-arc ephemerides are set" );
    }

    int arcIndex = lookUpscheme_->findNearestLowerNeighbour( currentTime );
    if( arcEndTimes_.at( arcIndex ) < currentTime )
    {
        return std::make_pair( false, arcIndex );
    }
    else
    {
        return std::make_pair( true, arcIndex );
    }
}

//! Function that retrieves the time interval at which an ephemeris can be safely interrogated
std::pair< double, double > getSafeEphemerisEvaluationInterval( const std::shared_ptr< ephemerides::Ephemeris > ephemerisModel )
{
    // Make default output pair
    std::pair< double, double > safeInterval =
            std::make_pair( std::numeric_limits< double >::lowest( ), std::numeric_limits< double >::max( ) );

    // Check if model is tabulated, and retrieve safe interval from model
    if( isTabulatedEphemeris( ephemerisModel ) )
    {
        safeInterval = getTabulatedEphemerisSafeInterval( ephemerisModel );
    }
    // Check if model is multi-arc, and retrieve safe intervals from first and last arc.
    else if( std::dynamic_pointer_cast< ephemerides::MultiArcEphemeris >( ephemerisModel ) != nullptr )
    {
        std::shared_ptr< ephemerides::MultiArcEphemeris > multiArcEphemerisModel =
                std::dynamic_pointer_cast< ephemerides::MultiArcEphemeris >( ephemerisModel );
        safeInterval.first = getSafeEphemerisEvaluationInterval( multiArcEphemerisModel->getSingleArcEphemerides( ).at( 0 ) ).first;
        safeInterval.second = getSafeEphemerisEvaluationInterval(
                                      multiArcEphemerisModel->getSingleArcEphemerides( ).at(
                                              multiArcEphemerisModel->getSingleArcEphemerides( ).size( ) - 1 ) ).second;
    }
    return safeInterval;
}


}  // namespace ephemerides

}  // namespace tudat
