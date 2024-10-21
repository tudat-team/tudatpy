/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/estimation_setup/observationsProcessing.h"

namespace tudat
{

namespace observation_models
{

//std::shared_ptr< ObservationCollectionParser > getObservationParserFromDependentVariableSettings(
//        const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > dependentVariableSettings )
//{
//    std::shared_ptr< ObservationCollectionParser > observationParser;
//
//    std::vector< std::shared_ptr< ObservationCollectionParser > > parserList;
//
//    // Check if relevant link end id and type are both specified
//    if ( ( dependentVariableSettings->linkEndId_ != LinkEndId( "", "" ) ) && ( dependentVariableSettings->linkEndType_ != unidentified_link_end ) )
//    {
//        parserList.push_back( std::make_shared< ObservationCollectionSingleLinkEndParser >( std::make_pair( dependentVariableSettings->linkEndType_, dependentVariableSettings->linkEndId_ ) ) );
//    }
//    // if only relevant link end id is specified
//    else if ( dependentVariableSettings->linkEndId_ != LinkEndId( "", "" ) )
//    {
//        parserList.push_back( std::make_shared< ObservationCollectionLinkEndIdParser >( dependentVariableSettings->linkEndId_ ) );
//    }
//    // if only relevant link end type is specified
//    else if ( dependentVariableSettings->linkEndType_ != unidentified_link_end )
//    {
//        parserList.push_back( std::make_shared< ObservationCollectionLinkEndTypeParser >( dependentVariableSettings->linkEndType_ ) );
//    }
//
//    // Check if originating link end id and type are both specified
//    if ( ( dependentVariableSettings->originatingLinkEndId_ != LinkEndId( "", "" ) ) && ( dependentVariableSettings->originatingLinkEndType_ != unidentified_link_end ) )
//    {
//        parserList.push_back( std::make_shared< ObservationCollectionSingleLinkEndParser >( std::make_pair( dependentVariableSettings->originatingLinkEndType_, dependentVariableSettings->originatingLinkEndId_ ) ) );
//    }
//        // if only originating link end id is specified
//    else if ( dependentVariableSettings->originatingLinkEndId_ != LinkEndId( "", "" ) )
//    {
//        parserList.push_back( std::make_shared< ObservationCollectionLinkEndIdParser >( dependentVariableSettings->originatingLinkEndId_ ) );
//    }
//        // if only originating link end type is specified
//    else if ( dependentVariableSettings->originatingLinkEndType_ != unidentified_link_end )
//    {
//        parserList.push_back( std::make_shared< ObservationCollectionLinkEndTypeParser >( dependentVariableSettings->originatingLinkEndType_ ) );
//    }
//
//    // Check if indirect condition on observable type via ancillary settings
//    if ( std::dynamic_pointer_cast< simulation_setup::AncillaryObservationDependentVariableSettings >( dependentVariableSettings ) != nullptr )
//    {
//        ObservableType observableType = std::dynamic_pointer_cast< simulation_setup::AncillaryObservationDependentVariableSettings >( dependentVariableSettings )->observableType_;
//        if ( observableType != undefined_observation_model )
//        {
//            parserList.push_back( std::make_shared< ObservationCollectionObservableTypeParser >( observableType ) );
//        }
//    }
//
//    // Create multi-type observation collection parser
//    if ( parserList.size( ) > 0 )
//    {
//        observationParser = std::make_shared< ObservationCollectionMultiTypeParser >( parserList, true );
//    }
//    else
//    {
//        observationParser = std::make_shared< ObservationCollectionParser >( );
//    }
//
//    return observationParser;
//}

}

}
