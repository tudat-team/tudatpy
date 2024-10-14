/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONSPROCESSING_H
#define TUDAT_OBSERVATIONSPROCESSING_H

#include <vector>

#include <memory>
#include <functional>

#include <Eigen/Core>

#include "tudat/basics/basicTypedefs.h"
#include "tudat/basics/timeType.h"
#include "tudat/basics/tudatTypeTraits.h"
#include "tudat/basics/utilities.h"

#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"

namespace tudat
{

namespace observation_models
{

enum ObservationFilterType
{
    residual_filtering,
    absolute_value_filtering,
    epochs_filtering,
    time_bounds_filtering,
    dependent_variable_filtering
};

struct ObservationFilterBase
{
public:

    ObservationFilterBase( const ObservationFilterType filterType, const bool filterOut = true, const bool useOppositeCondition = false ) :
            filterType_( filterType ), filterOut_( filterOut ), useOppositeCondition_( useOppositeCondition )
    { }

    virtual ~ObservationFilterBase( ){ }

    ObservationFilterType getFilterType( ) const
    {
        return filterType_;
    }

    bool filterOut( ) const
    {
        return filterOut_;
    }

    bool useOppositeCondition( ) const
    {
        return useOppositeCondition_;
    }

protected:

    ObservationFilterType filterType_;
    const bool filterOut_;
    const bool useOppositeCondition_;
};

template< typename FilterValueType >
bool checkFilterTypeConsistency( ObservationFilterType filterType )
{
    bool consistentType = true;
    switch( filterType )
    {
        case residual_filtering:
        case absolute_value_filtering:
        {
            if ( !( std::is_same< FilterValueType, double>::value ) ){ consistentType = false; }
            break;
        }
        case epochs_filtering:
        {
            if ( !( std::is_same< FilterValueType, std::vector< double > >::value ) ){ consistentType = false; }
            break;
        }
        case time_bounds_filtering:
        {
            if ( !( std::is_same< FilterValueType, std::pair< double, double > >::value ) ) { consistentType = false; }
            break;
        }
        default:
            break;
    }
    return consistentType;
}


template< typename FilterValueType >
struct ObservationFilter : public ObservationFilterBase
{
public:
    ObservationFilter( ObservationFilterType filterType, const FilterValueType& filterValue, const bool filterOut = true, const bool useOppositeCondition = false ) :
            ObservationFilterBase( filterType, filterOut, useOppositeCondition ), filterValue_( filterValue )
    {
        if ( !checkFilterTypeConsistency< FilterValueType >( filterType ) )
        {
            throw std::runtime_error( "Error when creating observation filter of type " + std::to_string( filterType ) + " filter value type is inconsistent." );
        }
    }

    virtual ~ObservationFilter( ){ }

    FilterValueType getFilterValue( ) const
    {
        return filterValue_;
    }

protected:

    FilterValueType filterValue_;
};


inline std::shared_ptr< ObservationFilterBase > observationFilter(
        const ObservationFilterType filterType, const double filterValue, const bool filterOut = true, const bool useOppositeCondition = false )
{
    return std::make_shared< ObservationFilter< double > >( filterType, filterValue, filterOut, useOppositeCondition );
}

inline std::shared_ptr< ObservationFilterBase > observationFilter(
        const ObservationFilterType filterType, const std::vector< double > filterValue, const bool filterOut = true, const bool useOppositeCondition = false )
{
    return std::make_shared< ObservationFilter< std::vector< double > > >( filterType, filterValue, filterOut, useOppositeCondition );
}

inline std::shared_ptr< ObservationFilterBase > observationFilter(
        const ObservationFilterType filterType, const std::pair< double, double > filterValue, const bool filterOut = true, const bool useOppositeCondition = false )
{
    return std::make_shared< ObservationFilter< std::pair< double, double > > >( filterType, filterValue, filterOut, useOppositeCondition );
}

enum ObservationSetSplitterType
{
    time_tags_splitter,
    time_interval_splitter,
    time_span_splitter,
    nb_observations_splitter,
    dependent_variables_splitter
};

struct ObservationSetSplitterBase
{
public:
    ObservationSetSplitterBase( const ObservationSetSplitterType splitterType, const int minNumberObservations = 0 ) :
    splitterType_( splitterType ), minNumberObservations_( minNumberObservations ){ }

    virtual ~ObservationSetSplitterBase( ){ }

    ObservationSetSplitterType getSplitterType( ) const
    {
        return splitterType_;
    }

    int getMinNumberObservations( ) const
    {
        return minNumberObservations_;
    }

protected:
    ObservationSetSplitterType splitterType_;
    int minNumberObservations_;
};

template< typename SetSplitType >
bool checkObsSetSplitTypeConsistency( ObservationSetSplitterType setSplitType )
{
    bool consistentType = true;
    switch( setSplitType )
    {
        case time_tags_splitter:
        {
            if ( !( std::is_same< SetSplitType, std::vector< double > >::value ) ){ consistentType = false; }
            break;
        }
        case time_interval_splitter:
        case time_span_splitter:
        {
            if ( !( std::is_same< SetSplitType, double >::value ) ){ consistentType = false; }
            break;
        }
        case nb_observations_splitter:
        {
            if ( !( std::is_same< SetSplitType, int >::value ) ) { consistentType = false; }
            break;
        }
        default:
            break;
    }
    return consistentType;
}

template< typename SplitterValueType >
struct ObservationSetSplitter : public ObservationSetSplitterBase
{
public:
    ObservationSetSplitter( ObservationSetSplitterType splitterType, const SplitterValueType& splitterValue, const int minNumberObservations ) :
            ObservationSetSplitterBase( splitterType, minNumberObservations ), splitterValue_( splitterValue )
    {
        if ( !checkObsSetSplitTypeConsistency< SplitterValueType >( splitterType ) )
        {
            throw std::runtime_error( "Error when creating observation set splitter of type " + std::to_string( splitterType ) + " split value type is inconsistent." );
        }
    }

    virtual ~ObservationSetSplitter( ){ }

    SplitterValueType getSplitterValue( ) const
    {
        return splitterValue_;
    }

protected:

    SplitterValueType splitterValue_;
};

inline std::shared_ptr< ObservationSetSplitterBase > observationSetSplitter(
        const ObservationSetSplitterType splitterType, const std::vector< double > splitterValue, const int minNumberObservations = 0 )
{
    return std::make_shared< ObservationSetSplitter< std::vector< double > > >( splitterType, splitterValue, minNumberObservations );
}

inline std::shared_ptr< ObservationSetSplitterBase > observationSetSplitter(
        const ObservationSetSplitterType splitterType, const double splitterValue, const int minNumberObservations = 0 )
{
    return std::make_shared< ObservationSetSplitter< double > >( splitterType, splitterValue, minNumberObservations );
}

inline std::shared_ptr< ObservationSetSplitterBase > observationSetSplitter(
        const ObservationSetSplitterType splitterType, const int splitterValue, const int minNumberObservations = 0 )
{
    return std::make_shared< ObservationSetSplitter< int > >( splitterType, splitterValue, minNumberObservations );
}

enum ObservationParserType
{
    empty_parser,
    observable_type_parser,
    link_ends_parser,
    link_end_string_parser,
    time_bounds_parser,
    ancillary_settings_parser,
    multi_type_parser
};

struct ObservationCollectionParser
{

public:

    ObservationCollectionParser( ) : parserType_( empty_parser ), useOppositeCondition_( false ){ }

    ObservationCollectionParser( const ObservationParserType parserType, const bool useOppositeCondition = false ) :
            parserType_( parserType ), useOppositeCondition_( useOppositeCondition ){ }

    virtual ~ObservationCollectionParser( ){ }

    ObservationParserType getObservationParserType( ) const
    {
        return parserType_;
    }

    bool useOppositeCondition( ) const
    {
        return useOppositeCondition_;
    }

protected:

    const ObservationParserType parserType_;

    const bool useOppositeCondition_;

};

struct ObservationCollectionObservableTypeParser : public ObservationCollectionParser
{
public:

    ObservationCollectionObservableTypeParser( const ObservableType observableType, const bool useOppositeCondition = false ) :
            ObservationCollectionParser( observable_type_parser, useOppositeCondition ), observableTypes_( std::vector< ObservableType >( { observableType } ) ){ }

    ObservationCollectionObservableTypeParser( const std::vector< ObservableType > observableTypes, const bool useOppositeCondition = false ) :
            ObservationCollectionParser( observable_type_parser, useOppositeCondition ), observableTypes_( observableTypes ){ }

    virtual ~ObservationCollectionObservableTypeParser( ){ }

    std::vector< ObservableType > getObservableTypes( ) const
    {
        return observableTypes_;
    }

protected:

    const std::vector< ObservableType > observableTypes_;

};

struct ObservationCollectionLinkEndsParser : public ObservationCollectionParser
{
public:

    ObservationCollectionLinkEndsParser( const LinkEnds linkEnds, const bool useOppositeCondition = false ) :
            ObservationCollectionParser( link_ends_parser, useOppositeCondition ), linkEndsVector_( std::vector< LinkEnds >( { linkEnds } ) ){ }

    ObservationCollectionLinkEndsParser( const std::vector< LinkEnds > linkEndsVector, const bool useOppositeCondition = false ) :
            ObservationCollectionParser( link_ends_parser, useOppositeCondition ), linkEndsVector_( linkEndsVector ){ }

    virtual ~ObservationCollectionLinkEndsParser( ){ }

    std::vector< LinkEnds > getLinkEndsVector( ) const
    {
        return linkEndsVector_;
    }

protected:

    const std::vector< LinkEnds > linkEndsVector_;

};

struct ObservationCollectionLinkEndStringParser : public ObservationCollectionParser
{
public:

    ObservationCollectionLinkEndStringParser( const std::string linkEndsNames,
                                          const bool isReferencePoint = false,
                                          const bool useOppositeCondition = false ) :
            ObservationCollectionParser( link_end_string_parser, useOppositeCondition ), linkEndsNames_( std::vector< std::string >( { linkEndsNames } ) ),
            isReferencePoint_( isReferencePoint ){ }

    ObservationCollectionLinkEndStringParser( const std::vector< std::string > linkEndsNames,
                                          const bool isReferencePoint = false,
                                          const bool useOppositeCondition = false ) :
            ObservationCollectionParser( link_end_string_parser, useOppositeCondition ), linkEndsNames_( linkEndsNames ), isReferencePoint_( isReferencePoint ){ }

    virtual ~ObservationCollectionLinkEndStringParser( ){ }

    std::vector< std::string > getLinkEndNames( ) const
    {
        return linkEndsNames_;
    }

    bool isReferencePoint( ) const
    {
        return isReferencePoint_;
    }

protected:

    const std::vector< std::string > linkEndsNames_;

    const bool isReferencePoint_;

};

struct ObservationCollectionTimeBoundsParser : public ObservationCollectionParser
{
public:

    ObservationCollectionTimeBoundsParser( const std::pair< double, double > timeBounds, const bool useOppositeCondition = false ) :
            ObservationCollectionParser( time_bounds_parser, useOppositeCondition ), timeBoundsVector_( std::vector< std::pair< double, double > >( { timeBounds } ) ){ }

    ObservationCollectionTimeBoundsParser( const std::vector< std::pair< double, double > > timeBoundsVector, const bool useOppositeCondition = false ) :
            ObservationCollectionParser( time_bounds_parser, useOppositeCondition ), timeBoundsVector_( timeBoundsVector ){ }

    virtual ~ObservationCollectionTimeBoundsParser( ){ }

    std::vector< std::pair< double, double > > getTimeBoundsVector( ) const
    {
        return timeBoundsVector_;
    }

protected:

    const std::vector< std::pair< double, double > > timeBoundsVector_;

};

struct ObservationCollectionAncillarySettingsParser : public ObservationCollectionParser
{
public:

    ObservationCollectionAncillarySettingsParser( const std::shared_ptr< ObservationAncilliarySimulationSettings > ancillarySettings, const bool useOppositeCondition = false ) :
            ObservationCollectionParser( ancillary_settings_parser, useOppositeCondition ),
            ancillarySettings_( std::vector< std::shared_ptr< ObservationAncilliarySimulationSettings > >( { ancillarySettings } ) ){ }

    ObservationCollectionAncillarySettingsParser( const std::vector< std::shared_ptr< ObservationAncilliarySimulationSettings > > ancillarySettings, const bool useOppositeCondition = false ) :
            ObservationCollectionParser( ancillary_settings_parser, useOppositeCondition ), ancillarySettings_( ancillarySettings ){ }

    virtual ~ObservationCollectionAncillarySettingsParser( ){ }

    std::vector< std::shared_ptr< ObservationAncilliarySimulationSettings > > getAncillarySettings( ) const
    {
        return ancillarySettings_;
    }

protected:

    const std::vector< std::shared_ptr< ObservationAncilliarySimulationSettings > > ancillarySettings_;

};

struct ObservationCollectionMultiTypeParser : public ObservationCollectionParser
{
public:

    ObservationCollectionMultiTypeParser( const std::vector< std::shared_ptr< ObservationCollectionParser > >& observationParsers,
                                          const bool combineConditions = false ) :
            ObservationCollectionParser( multi_type_parser, false ), observationParsers_( observationParsers ), combineConditions_( combineConditions ){ }

    virtual ~ObservationCollectionMultiTypeParser( ){ }

    std::vector< std::shared_ptr< ObservationCollectionParser > > getObservationParsers_(  ) const
    {
        return observationParsers_;
    }

    bool areConditionsCombined( ) const
    {
        return combineConditions_;
    }

protected:

    const std::vector< std::shared_ptr< ObservationCollectionParser > > observationParsers_;

    const bool combineConditions_;

};

inline std::shared_ptr< ObservationCollectionParser > observationParser( )
{
    return std::make_shared< ObservationCollectionParser >( );
}

inline std::shared_ptr< ObservationCollectionParser > observationParser( const ObservableType observableType, const bool useOppositeCondition = false )
{
    return std::make_shared< ObservationCollectionObservableTypeParser >( observableType, useOppositeCondition );
}

inline std::shared_ptr< ObservationCollectionParser > observationParser( const std::vector< ObservableType >& observableTypes, const bool useOppositeCondition = false )
{
    return std::make_shared< ObservationCollectionObservableTypeParser >( observableTypes, useOppositeCondition );
}

inline std::shared_ptr< ObservationCollectionParser > observationParser( const LinkEnds linkEnds, const bool useOppositeCondition = false )
{
    return std::make_shared< ObservationCollectionLinkEndsParser >( linkEnds, useOppositeCondition );
}

inline std::shared_ptr< ObservationCollectionParser > observationParser( const std::vector< LinkEnds >& linkEndsVector, const bool useOppositeCondition = false )
{
    return std::make_shared< ObservationCollectionLinkEndsParser >( linkEndsVector, useOppositeCondition );
}

inline std::shared_ptr< ObservationCollectionParser > observationParser( const std::string bodyName, const bool isReferencePoint = false, const bool useOppositeCondition = false )
{
    return std::make_shared< ObservationCollectionLinkEndStringParser >( bodyName, isReferencePoint, useOppositeCondition );
}

inline std::shared_ptr< ObservationCollectionParser > observationParser( const std::vector< std::string >& bodyNames,
                                                                         const bool isReferencePoint = false, const bool useOppositeCondition = false )
{
    return std::make_shared< ObservationCollectionLinkEndStringParser >( bodyNames, isReferencePoint, useOppositeCondition );
}

inline std::shared_ptr< ObservationCollectionParser > observationParser( const std::pair< double, double >& timeBounds, const bool useOppositeCondition = false )
{
    return std::make_shared< ObservationCollectionTimeBoundsParser >( timeBounds, useOppositeCondition );
}

inline std::shared_ptr< ObservationCollectionParser > observationParser( const std::vector< std::pair< double, double > >& timeBoundsVector, const bool useOppositeCondition = false )
{
    return std::make_shared< ObservationCollectionTimeBoundsParser >( timeBoundsVector, useOppositeCondition );
}

inline std::shared_ptr< ObservationCollectionParser > observationParser(
        const std::shared_ptr< ObservationAncilliarySimulationSettings > ancillarySettings, const bool useOppositeCondition = false )
{
    return std::make_shared< ObservationCollectionAncillarySettingsParser >( ancillarySettings, useOppositeCondition );
}

inline std::shared_ptr< ObservationCollectionParser > observationParser(
        const std::vector< std::shared_ptr< ObservationAncilliarySimulationSettings > >& ancillarySettings, const bool useOppositeCondition = false )
{
    return std::make_shared< ObservationCollectionAncillarySettingsParser >( ancillarySettings, useOppositeCondition );
}

inline std::shared_ptr< ObservationCollectionParser > observationParser( const std::vector< std::shared_ptr< ObservationCollectionParser > >& observationParsers,
                                                                         const bool combineConditions = false )
{
    return std::make_shared< ObservationCollectionMultiTypeParser >( observationParsers, combineConditions );
}


} // namespace observation_models

} // namespace tudat

#endif // TUDAT_OBSERVATIONSPROCESSING_H
