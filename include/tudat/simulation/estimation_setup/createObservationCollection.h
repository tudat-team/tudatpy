/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATE_OBSERVATION_COLLECTION_H
#define TUDAT_CREATE_OBSERVATION_COLLECTION_H

#include <Eigen/Core>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "tudat/astro/ephemerides/tabulatedEphemeris.h"
#include "tudat/astro/ephemerides/tabulatedRotationalEphemeris.h"
#include "tudat/basics/timeType.h"
#include "tudat/basics/utilities.h"
#include "tudat/io/trackingSupplementaryData.h"
#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/simulation/environment_setup/body.h"

namespace tudat
{

namespace observation_models
{

template< typename EphemerisScalarType, typename EphemerisTimeType >
inline void resetTabulatedEphemerisFromTrackingSupplementaryStateHistory(
        const std::map< double, Eigen::Vector6d >& stateHistory,
        const std::shared_ptr< ephemerides::TabulatedCartesianEphemeris< EphemerisScalarType, EphemerisTimeType > > tabulatedEphemeris )
{
    std::map< EphemerisTimeType, Eigen::Matrix< EphemerisScalarType, 6, 1 > > castStateHistory;
    utilities::castMatrixMap< double, double, EphemerisTimeType, EphemerisScalarType, 6, 1 >( stateHistory, castStateHistory );

    tabulatedEphemeris->resetInterpolator(
            interpolators::createOneDimensionalInterpolator( castStateHistory, interpolators::linearInterpolation( ) ) );
}

inline void resetTabulatedEphemerisFromTrackingSupplementaryStateHistory(
        const std::map< double, Eigen::Vector6d >& stateHistory,
        const std::shared_ptr< ephemerides::Ephemeris > ephemeris,
        const std::string& bodyName )
{
    if( std::dynamic_pointer_cast< ephemerides::TabulatedCartesianEphemeris< double, double > >( ephemeris ) != nullptr )
    {
        resetTabulatedEphemerisFromTrackingSupplementaryStateHistory(
                stateHistory,
                std::dynamic_pointer_cast< ephemerides::TabulatedCartesianEphemeris< double, double > >( ephemeris ) );
    }
    else if( std::dynamic_pointer_cast< ephemerides::TabulatedCartesianEphemeris< long double, double > >( ephemeris ) != nullptr )
    {
        resetTabulatedEphemerisFromTrackingSupplementaryStateHistory(
                stateHistory,
                std::dynamic_pointer_cast< ephemerides::TabulatedCartesianEphemeris< long double, double > >( ephemeris ) );
    }
    else if( std::dynamic_pointer_cast< ephemerides::TabulatedCartesianEphemeris< double, Time > >( ephemeris ) != nullptr )
    {
        resetTabulatedEphemerisFromTrackingSupplementaryStateHistory(
                stateHistory,
                std::dynamic_pointer_cast< ephemerides::TabulatedCartesianEphemeris< double, Time > >( ephemeris ) );
    }
    else if( std::dynamic_pointer_cast< ephemerides::TabulatedCartesianEphemeris< long double, Time > >( ephemeris ) != nullptr )
    {
        resetTabulatedEphemerisFromTrackingSupplementaryStateHistory(
                stateHistory,
                std::dynamic_pointer_cast< ephemerides::TabulatedCartesianEphemeris< long double, Time > >( ephemeris ) );
    }
    else
    {
        throw std::runtime_error( "Error when setting tracking supplementary data in body " + bodyName +
                                  ", tabulated ephemeris type was not recognized." );
    }
}

template< typename EphemerisScalarType, typename EphemerisTimeType >
inline void resetTabulatedRotationalEphemerisFromTrackingSupplementaryStateHistory(
        const std::map< double, Eigen::Vector7d >& rotationalStateHistory,
        const std::shared_ptr< ephemerides::TabulatedRotationalEphemeris< EphemerisScalarType, EphemerisTimeType > >
                tabulatedRotationalEphemeris )
{
    std::map< EphemerisTimeType, Eigen::Matrix< EphemerisScalarType, 7, 1 > > castRotationalStateHistory;
    utilities::castMatrixMap< double, double, EphemerisTimeType, EphemerisScalarType, 7, 1 >( rotationalStateHistory,
                                                                                              castRotationalStateHistory );

    tabulatedRotationalEphemeris->reset(
            interpolators::createOneDimensionalInterpolator( castRotationalStateHistory, interpolators::linearInterpolation( ) ) );
}

inline void resetTabulatedRotationalEphemerisFromTrackingSupplementaryStateHistory(
        const std::map< double, Eigen::Vector7d >& rotationalStateHistory,
        const std::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris,
        const std::string& bodyName )
{
    if( std::dynamic_pointer_cast< ephemerides::TabulatedRotationalEphemeris< double, double > >( rotationalEphemeris ) != nullptr )
    {
        resetTabulatedRotationalEphemerisFromTrackingSupplementaryStateHistory(
                rotationalStateHistory,
                std::dynamic_pointer_cast< ephemerides::TabulatedRotationalEphemeris< double, double > >( rotationalEphemeris ) );
    }
    else if( std::dynamic_pointer_cast< ephemerides::TabulatedRotationalEphemeris< long double, double > >( rotationalEphemeris ) !=
             nullptr )
    {
        resetTabulatedRotationalEphemerisFromTrackingSupplementaryStateHistory(
                rotationalStateHistory,
                std::dynamic_pointer_cast< ephemerides::TabulatedRotationalEphemeris< long double, double > >( rotationalEphemeris ) );
    }
    else if( std::dynamic_pointer_cast< ephemerides::TabulatedRotationalEphemeris< double, Time > >( rotationalEphemeris ) != nullptr )
    {
        resetTabulatedRotationalEphemerisFromTrackingSupplementaryStateHistory(
                rotationalStateHistory,
                std::dynamic_pointer_cast< ephemerides::TabulatedRotationalEphemeris< double, Time > >( rotationalEphemeris ) );
    }
    else if( std::dynamic_pointer_cast< ephemerides::TabulatedRotationalEphemeris< long double, Time > >( rotationalEphemeris ) !=
             nullptr )
    {
        resetTabulatedRotationalEphemerisFromTrackingSupplementaryStateHistory(
                rotationalStateHistory,
                std::dynamic_pointer_cast< ephemerides::TabulatedRotationalEphemeris< long double, Time > >( rotationalEphemeris ) );
    }
    else
    {
        throw std::runtime_error( "Error when setting tracking supplementary data in body " + bodyName +
                                  ", tabulated rotational ephemeris type was not recognized." );
    }
}

inline void setTranslationalStateSupplementaryDataInBodies(
    simulation_setup::SystemOfBodies& bodies,
    const std::map< std::pair< std::string, std::string >, std::vector< data::TranslationalStateSupplementaryData > >&
            translationalStateSupplementaryData )
{
    for( auto it = translationalStateSupplementaryData.begin( ); it != translationalStateSupplementaryData.end( ); ++it )
    {
        std::string bodyName = it->first.first;
        std::string referencePointName = it->first.second;
        if( referencePointName != "" )
        {
            throw std::runtime_error( "Error, reference point ID for setting ephemeris from tracking supplementary data must be empty, but found " + referencePointName );
        }
        std::map< double, Eigen::Vector6d > stateHistory;
        std::vector< std::pair< double, std::pair< Eigen::Vector6d, Eigen::Vector6d > > > inconsistentDuplicateStateHistoryEntries;
        std::string frameOrigin;

        for( unsigned int i = 0; i < it->second.size( ); ++i )
        {
            if( i == 0 )
            {
                frameOrigin = it->second.at( i ).getFrameOrigin( );
            }
            else if( it->second.at( i ).getFrameOrigin( ) != frameOrigin )
            {
                throw std::runtime_error( "Error, inconsistent frame origins found when setting translational state from tracking "
                                          "supplementary data for body " +
                                          bodyName + ". Found " + it->second.at( i ).getFrameOrigin( ) + " but expected " +
                                          frameOrigin + "." );
            }

            const std::map< double, Eigen::Vector6d >& currentStateHistory = it->second.at( i ).getStateHistory( );
            for( auto stateIterator = currentStateHistory.begin( ); stateIterator != currentStateHistory.end( ); ++stateIterator )
            {
                if( stateHistory.count( stateIterator->first ) == 0 )
                {
                    stateHistory[ stateIterator->first ] = stateIterator->second;
                }
                else if( !stateHistory.at( stateIterator->first ).isApprox( stateIterator->second, 0.0 ) )
                {
                    inconsistentDuplicateStateHistoryEntries.push_back(
                            std::make_pair( stateIterator->first,
                                            std::make_pair( stateHistory.at( stateIterator->first ), stateIterator->second ) ) );
                }
            }
        }

        std::shared_ptr< ephemerides::Ephemeris > ephemeris = bodies.at( bodyName )->getEphemeris( );
        if( ephemeris != nullptr )
        {
            if( ephemeris->getReferenceFrameOrigin( ) != frameOrigin )
            {
                throw std::runtime_error( "Error when setting tracking supplementary data in body " + bodyName +
                                          ", existing ephemeris frame origin is " + ephemeris->getReferenceFrameOrigin( ) +
                                          " but supplementary data frame origin is " + frameOrigin + "." );
            }

            if( !ephemerides::isTabulatedEphemeris( ephemeris ) )
            {
                throw std::runtime_error( "Error when setting tracking supplementary data in body " + bodyName +
                                          ", existing ephemeris is not tabulated." );
            }

            resetTabulatedEphemerisFromTrackingSupplementaryStateHistory( stateHistory, ephemeris, bodyName );
        }
        else
        {
            std::map< Time, Eigen::Vector6d > timeStateHistory;
            utilities::castMatrixMap< double, double, Time, double, 6, 1 >( stateHistory, timeStateHistory );
            bodies.at( bodyName )
                    ->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< double, Time > >(
                            interpolators::createOneDimensionalInterpolator( timeStateHistory, interpolators::linearInterpolation( ) ),
                            frameOrigin ) );
        }
    }
}

inline void setRotationalStateSupplementaryDataInBodies(
    simulation_setup::SystemOfBodies& bodies,
    const std::map< std::pair< std::string, std::string >, std::vector< data::RotationalStateSupplementaryData > >&
            rotationalStateSupplementaryData )
{
    for( auto it = rotationalStateSupplementaryData.begin( ); it != rotationalStateSupplementaryData.end( ); ++it )
    {
        std::string bodyName = it->first.first;
        std::string referencePointName = it->first.second;
        if( referencePointName != "" )
        {
            throw std::runtime_error( "Error, reference point ID for setting rotational ephemeris from tracking supplementary data must be empty, but found " + referencePointName );
        }
        std::map< double, Eigen::Vector7d > rotationalStateHistory;
        std::vector< std::pair< double, std::pair< Eigen::Vector7d, Eigen::Vector7d > > > inconsistentDuplicateRotationalStateHistoryEntries;
        std::string baseFrameOrientation;

        for( unsigned int i = 0; i < it->second.size( ); ++i )
        {
            if( i == 0 )
            {
                baseFrameOrientation = it->second.at( i ).getBaseFrameOrientation( );
            }
            else if( it->second.at( i ).getBaseFrameOrientation( ) != baseFrameOrientation )
            {
                throw std::runtime_error( "Error, inconsistent base frame orientations found when setting rotational state from tracking "
                                          "supplementary data for body " +
                                          bodyName + ". Found " + it->second.at( i ).getBaseFrameOrientation( ) + " but expected " +
                                          baseFrameOrientation + "." );
            }

            const std::map< double, Eigen::Vector7d >& currentRotationalStateHistory =
                    it->second.at( i ).getRotationalStateHistory( );
            for( auto stateIterator = currentRotationalStateHistory.begin( );
                 stateIterator != currentRotationalStateHistory.end( );
                 ++stateIterator )
            {
                if( rotationalStateHistory.count( stateIterator->first ) == 0 )
                {
                    rotationalStateHistory[ stateIterator->first ] = stateIterator->second;
                }
                else if( !rotationalStateHistory.at( stateIterator->first ).isApprox( stateIterator->second, 0.0 ) )
                {
                    inconsistentDuplicateRotationalStateHistoryEntries.push_back(
                            std::make_pair( stateIterator->first,
                                            std::make_pair( rotationalStateHistory.at( stateIterator->first ),
                                                            stateIterator->second ) ) );
                }
            }
        }

        std::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris = bodies.at( bodyName )->getRotationalEphemeris( );
        if( rotationalEphemeris != nullptr )
        {
            if( rotationalEphemeris->getBaseFrameOrientation( ) != baseFrameOrientation )
            {
                throw std::runtime_error( "Error when setting tracking supplementary rotational data in body " + bodyName +
                                          ", existing rotational ephemeris base frame orientation is " +
                                          rotationalEphemeris->getBaseFrameOrientation( ) +
                                          " but supplementary data base frame orientation is " + baseFrameOrientation + "." );
            }

            if( !ephemerides::isTabulatedRotationalEphemeris( rotationalEphemeris ) )
            {
                throw std::runtime_error( "Error when setting tracking supplementary data in body " + bodyName +
                                          ", existing rotational ephemeris is not tabulated." );
            }

            resetTabulatedRotationalEphemerisFromTrackingSupplementaryStateHistory(
                    rotationalStateHistory, rotationalEphemeris, bodyName );
        }
        else
        {
            std::map< Time, Eigen::Vector7d > timeRotationalStateHistory;
            utilities::castMatrixMap< double, double, Time, double, 7, 1 >( rotationalStateHistory, timeRotationalStateHistory );
            bodies.at( bodyName )
                    ->setRotationalEphemeris(
                            std::make_shared< ephemerides::TabulatedRotationalEphemeris< double, Time > >(
                                    interpolators::createOneDimensionalInterpolator( timeRotationalStateHistory,
                                                                                     interpolators::linearInterpolation( ) ),
                                    baseFrameOrientation,
                                    bodyName + "_Fixed" ) );
        }
    }
}

inline void setFrequencySupplementaryDataInBodies(
    simulation_setup::SystemOfBodies& bodies,
    const std::map< std::pair< std::string, std::string >, std::vector< std::shared_ptr< data::FrequencySupplementaryData > > >&
            frequencySupplementaryData )
{
    for( auto it = frequencySupplementaryData.begin( ); it != frequencySupplementaryData.end( ); ++it )
    {
        std::string bodyName = it->first.first;
        std::string referencePointName = it->first.second;

        std::vector< data::RampedFrequencySupplementaryData::FrequencyRamp > frequencyRamps;
        std::vector< std::map< double, double > > piecewiseConstantFrequencyHistories;

        for( unsigned int i = 0; i < it->second.size( ); ++i )
        {
            if( it->second.at( i ) == nullptr )
            {
                throw std::runtime_error( "Error when setting frequency supplementary data in body " + bodyName +
                                          ", reference point " + referencePointName + ": frequency data entry is null." );
            }

            if( i > 0 &&
                it->second.at( i )->getFrequencySupplementaryDataType( ) != it->second.at( 0 )->getFrequencySupplementaryDataType( ) )
            {
                throw std::runtime_error( "Error when setting frequency supplementary data in body " + bodyName +
                                          ", reference point " + referencePointName +
                                          ": all frequency supplementary data entries must have the same type." );
            }

            if( it->second.at( i )->getFrequencySupplementaryDataType( ) == data::FrequencySupplementaryDataType::ramped_frequency )
            {
                std::shared_ptr< data::RampedFrequencySupplementaryData > rampedFrequencySupplementaryData =
                        std::dynamic_pointer_cast< data::RampedFrequencySupplementaryData >( it->second.at( i ) );
                if( rampedFrequencySupplementaryData == nullptr )
                {
                    throw std::runtime_error( "Error when setting frequency supplementary data in body " + bodyName +
                                              ", reference point " + referencePointName +
                                              ": frequency data type is ramped, but derived object type is inconsistent." );
                }
                const std::vector< data::RampedFrequencySupplementaryData::FrequencyRamp >& currentFrequencyRamps =
                        rampedFrequencySupplementaryData->getFrequencyRamps( );
                frequencyRamps.insert( frequencyRamps.end( ), currentFrequencyRamps.begin( ), currentFrequencyRamps.end( ) );
            }
            else if( it->second.at( i )->getFrequencySupplementaryDataType( ) ==
                     data::FrequencySupplementaryDataType::piecewise_constant_frequency )
            {
                std::shared_ptr< data::PiecewiseConstantFrequencySupplementaryData > piecewiseConstantFrequencySupplementaryData =
                        std::dynamic_pointer_cast< data::PiecewiseConstantFrequencySupplementaryData >( it->second.at( i ) );
                if( piecewiseConstantFrequencySupplementaryData == nullptr )
                {
                    throw std::runtime_error( "Error when setting frequency supplementary data in body " + bodyName +
                                              ", reference point " + referencePointName +
                                              ": frequency data type is piecewise constant, but derived object type is inconsistent." );
                }
                piecewiseConstantFrequencyHistories.push_back( piecewiseConstantFrequencySupplementaryData->getFrequencyHistory( ) );
            }
        }
    }
}

inline void setInstrumentSupplementaryDataInBodies(
    simulation_setup::SystemOfBodies& bodies,
    const std::map< std::pair< std::string, std::string >, std::vector< std::shared_ptr< data::InstrumentSupplementaryData > > >&
            instrumentSupplementaryData )
{
    for( auto it = instrumentSupplementaryData.begin( ); it != instrumentSupplementaryData.end( ); ++it )
    { }
}

inline void setTrackingSupplementaryDataInBodies(
    simulation_setup::SystemOfBodies& bodies,
    const std::vector< data::TrackingSupplementaryData >& supplementaryData )
{
    std::map< std::pair< std::string, std::string >, std::vector< data::TranslationalStateSupplementaryData > >
            translationalStateSupplementaryData;
    std::map< std::pair< std::string, std::string >, std::vector< data::RotationalStateSupplementaryData > >
            rotationalStateSupplementaryData;
    std::map< std::pair< std::string, std::string >, std::vector< std::shared_ptr< data::FrequencySupplementaryData > > >
            frequencySupplementaryData;
    std::map< std::pair< std::string, std::string >, std::vector< std::shared_ptr< data::InstrumentSupplementaryData > > >
            instrumentSupplementaryData;

    for( const data::TrackingSupplementaryData& currentSupplementaryData: supplementaryData )
    {
        const std::pair< std::string, std::string > bodyReferencePoint = std::make_pair(
                currentSupplementaryData.getBodyName( ), currentSupplementaryData.getReferencePointName( ) );

        translationalStateSupplementaryData[ bodyReferencePoint ].push_back(
                currentSupplementaryData.getTranslationalStateSupplementaryData( ) );
        rotationalStateSupplementaryData[ bodyReferencePoint ].push_back( currentSupplementaryData.getRotationalStateSupplementaryData( ) );

        const std::vector< std::shared_ptr< data::FrequencySupplementaryData > >& currentFrequencySupplementaryData =
                currentSupplementaryData.getFrequencySupplementaryData( );
        frequencySupplementaryData[ bodyReferencePoint ].insert( frequencySupplementaryData[ bodyReferencePoint ].end( ),
                                                                 currentFrequencySupplementaryData.begin( ),
                                                                 currentFrequencySupplementaryData.end( ) );

        const std::vector< std::shared_ptr< data::InstrumentSupplementaryData > >& currentInstrumentSupplementaryData =
                currentSupplementaryData.getInstrumentSupplementaryData( );
        instrumentSupplementaryData[ bodyReferencePoint ].insert( instrumentSupplementaryData[ bodyReferencePoint ].end( ),
                                                                  currentInstrumentSupplementaryData.begin( ),
                                                                  currentInstrumentSupplementaryData.end( ) );
    }

    setTranslationalStateSupplementaryDataInBodies( bodies, translationalStateSupplementaryData );
    setRotationalStateSupplementaryDataInBodies( bodies, rotationalStateSupplementaryData );
    setFrequencySupplementaryDataInBodies( bodies, frequencySupplementaryData );
    setInstrumentSupplementaryDataInBodies( bodies, instrumentSupplementaryData );
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_CREATE_OBSERVATION_COLLECTION_H
