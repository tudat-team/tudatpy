/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATION_LEGACY_PARSER_ADAPTER_H
#define TUDAT_OBSERVATION_LEGACY_PARSER_ADAPTER_H

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "tudat/simulation/estimation_setup/observationDataset.h"

namespace tudat
{

namespace observation_models
{

//! Evaluate a legacy ObservationCollectionParser against one dataset set.
/*!
 * This adapter preserves the old set-level parser behavior for compatibility
 * APIs. New row-level selection should use ObservationCondition instead.
 */
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type >
bool isObservationSetSelectedByLegacyParser( const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                                             const unsigned int setId,
                                             const std::shared_ptr< ObservationCollectionParser >& observationParser )
{
    const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = dataset.getObservationSetMetadata( setId );
    const LinkEnds linkEnds = dataset.getLinkDefinition( metadata.linkDefinitionId_ ).linkEnds_;
    const bool useOppositeCondition = observationParser->useOppositeCondition( );

    bool isSelected = false;
    switch( observationParser->getObservationParserType( ) )
    {
        case empty_parser: {
            isSelected = true;
            break;
        }
        case observable_type_parser: {
            const std::vector< ObservableType > observableTypes =
                    std::dynamic_pointer_cast< ObservationCollectionObservableTypeParser >( observationParser )->getObservableTypes( );
            isSelected = std::count( observableTypes.begin( ), observableTypes.end( ), metadata.observableType_ ) > 0;
            break;
        }
        case link_ends_parser: {
            const std::vector< LinkEnds > linkEndsVector =
                    std::dynamic_pointer_cast< ObservationCollectionLinkEndsParser >( observationParser )->getLinkEndsVector( );
            isSelected = std::count( linkEndsVector.begin( ), linkEndsVector.end( ), linkEnds ) > 0;
            break;
        }
        case link_end_string_parser: {
            const std::shared_ptr< ObservationCollectionLinkEndStringParser > linkEndStringObservationParser =
                    std::dynamic_pointer_cast< ObservationCollectionLinkEndStringParser >( observationParser );
            const std::vector< std::string > linkEndNames = linkEndStringObservationParser->getLinkEndNames( );
            const bool isReferencePoint = linkEndStringObservationParser->isReferencePoint( );

            for( const auto& linkEnd : linkEnds )
            {
                const std::string name = isReferencePoint ? linkEnd.second.getReferencePointName( ) : linkEnd.second.bodyName_;
                if( std::count( linkEndNames.begin( ), linkEndNames.end( ), name ) > 0 )
                {
                    isSelected = true;
                }
            }
            break;
        }
        case link_end_id_parser: {
            const std::vector< LinkEndId > linkEndIds =
                    std::dynamic_pointer_cast< ObservationCollectionLinkEndIdParser >( observationParser )->getLinkEndIds( );
            for( const auto& linkEnd : linkEnds )
            {
                if( std::count( linkEndIds.begin( ), linkEndIds.end( ), linkEnd.second ) > 0 )
                {
                    isSelected = true;
                }
            }
            break;
        }
        case link_end_type_parser: {
            const std::vector< LinkEndType > linkEndTypes =
                    std::dynamic_pointer_cast< ObservationCollectionLinkEndTypeParser >( observationParser )->getLinkEndTypes( );
            for( const auto& linkEnd : linkEnds )
            {
                if( std::count( linkEndTypes.begin( ), linkEndTypes.end( ), linkEnd.first ) > 0 )
                {
                    isSelected = true;
                }
            }
            break;
        }
        case single_link_end_parser: {
            const std::vector< std::pair< LinkEndType, LinkEndId > > singleLinkEnds =
                    std::dynamic_pointer_cast< ObservationCollectionSingleLinkEndParser >( observationParser )->getSingleLinkEnds( );
            for( const auto& singleLinkEnd : singleLinkEnds )
            {
                if( linkEnds.count( singleLinkEnd.first ) > 0 && linkEnds.at( singleLinkEnd.first ) == singleLinkEnd.second )
                {
                    isSelected = true;
                }
            }
            break;
        }
        case time_bounds_parser: {
            const std::vector< std::pair< double, double > > timeBoundsVector =
                    std::dynamic_pointer_cast< ObservationCollectionTimeBoundsParser >( observationParser )->getTimeBoundsVector( );
            const std::pair< TimeType, TimeType > setTimeBounds = dataset.getTimeBoundsForSet( setId );
            for( const auto& timeBounds : timeBoundsVector )
            {
                if( ( setTimeBounds.first >= timeBounds.first ) && ( setTimeBounds.second <= timeBounds.second ) )
                {
                    isSelected = true;
                }
            }
            break;
        }
        case ancillary_settings_parser: {
            const std::vector< std::shared_ptr< ObservationAncillarySimulationSettings > > ancillarySettings =
                    std::dynamic_pointer_cast< ObservationCollectionAncillarySettingsParser >( observationParser )->getAncillarySettings( );
            const std::shared_ptr< ObservationAncillarySimulationSettings >& currentSettings =
                    dataset.getAncillarySettings( metadata.ancillarySettingsId_ );
            for( const auto& settings : ancillarySettings )
            {
                if( currentSettings == settings )
                {
                    isSelected = true;
                }
                else if( currentSettings != nullptr && settings != nullptr &&
                         currentSettings->getDoubleData( ) == settings->getDoubleData( ) &&
                         currentSettings->getDoubleVectorData( ) == settings->getDoubleVectorData( ) )
                {
                    isSelected = true;
                }
            }
            break;
        }
        case multi_type_parser: {
            const std::shared_ptr< ObservationCollectionMultiTypeParser > multiTypeParser =
                    std::dynamic_pointer_cast< ObservationCollectionMultiTypeParser >( observationParser );
            const std::vector< std::shared_ptr< ObservationCollectionParser > > observationParsers =
                    multiTypeParser->getObservationParsers_( );

            if( multiTypeParser->areConditionsCombined( ) )
            {
                isSelected = true;
                for( const std::shared_ptr< ObservationCollectionParser >& parser : observationParsers )
                {
                    if( !isObservationSetSelectedByLegacyParser( dataset, setId, parser ) )
                    {
                        isSelected = false;
                    }
                }
            }
            else
            {
                for( const std::shared_ptr< ObservationCollectionParser >& parser : observationParsers )
                {
                    if( isObservationSetSelectedByLegacyParser( dataset, setId, parser ) )
                    {
                        isSelected = true;
                    }
                }
            }
            break;
        }
        default:
            throw std::runtime_error( "Observation parser type not recognised." );
    }

    return useOppositeCondition ? !isSelected : isSelected;
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATION_LEGACY_PARSER_ADAPTER_H
