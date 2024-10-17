/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONOUTPUTSETTINGS
#define TUDAT_OBSERVATIONOUTPUTSETTINGS

#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/simulation/estimation_setup/observationsProcessing.h"

namespace tudat
{

namespace simulation_setup
{

using namespace observation_models;

enum ObservationDependentVariables
{
    station_elevation_angle,
    station_azimuth_angle,
    target_range,
    body_avoidance_angle_variable,
    link_body_center_distance,
    link_limb_distance,
    link_angle_with_orbital_plane,
    doppler_integration_time_dependent_variable,
    retransmission_delays_dependent_variable
};


class ObservationDependentVariableSettings
{
public:
    ObservationDependentVariableSettings(
            const ObservationDependentVariables variableType,
            const LinkEndId linkEndId = LinkEndId( "", "" ),
            const LinkEndType linkEndType = unidentified_link_end,
            const LinkEndId originatingLinkEndId = LinkEndId( "", "" ),
            const LinkEndType originatingLinkEndType = unidentified_link_end ):
        variableType_( variableType ), linkEndId_( linkEndId ), linkEndType_( linkEndType ), originatingLinkEndId_( originatingLinkEndId ), originatingLinkEndType_( originatingLinkEndType )
    {

    }

    virtual ~ObservationDependentVariableSettings( ){ }

    ObservationDependentVariables variableType_;

    virtual std::string getIdentifier( ) = 0;

    virtual bool areSettingsCompatible( const std::shared_ptr< ObservationDependentVariableSettings > otherSettings ) = 0;

    bool areBaseSettingsCompatible( const std::shared_ptr< ObservationDependentVariableSettings > otherSettings )
    {
        bool isCompatible = true;
        if ( variableType_ != otherSettings->variableType_ )
        {
            isCompatible = false;
        }
        else
        {
            if ( ( linkEndType_ != otherSettings->linkEndType_ && linkEndType_ != unidentified_link_end && otherSettings->linkEndType_ != unidentified_link_end )
            || ( ( linkEndId_ != otherSettings->linkEndId_ && linkEndId_ != LinkEndId( "", "" ) && otherSettings->linkEndId_ != LinkEndId( "", "" ) ) )
            || ( ( originatingLinkEndType_ != otherSettings->originatingLinkEndType_ &&
                originatingLinkEndType_ != unidentified_link_end && otherSettings->originatingLinkEndType_ != unidentified_link_end ) )
            || ( ( originatingLinkEndId_ != otherSettings->originatingLinkEndId_ &&
                originatingLinkEndId_ != LinkEndId( "", "" ) && otherSettings->originatingLinkEndId_ != LinkEndId( "", "" ) ) ) )
            {
                isCompatible = false;
            }
        }

        return isCompatible;
    }

    LinkEndId linkEndId_;
    LinkEndId originatingLinkEndId_;
    LinkEndType linkEndType_;
    LinkEndType originatingLinkEndType_;

};


enum IntegratedObservationPropertyHandling
{
    interval_start,
    interval_end,
    interval_undefined
};

std::string getIntegrationHandlingString( const IntegratedObservationPropertyHandling integratedObservableHandling );

class StationAngleObservationDependentVariableSettings: public ObservationDependentVariableSettings
{
public:
    StationAngleObservationDependentVariableSettings(
            const ObservationDependentVariables variableType,
            const LinkEndId relevantLinkEnd = LinkEndId( "", "" ),
            const LinkEndType linkEndRole = unidentified_link_end,
            const IntegratedObservationPropertyHandling integratedObservableHandling = interval_start, // interval_undefined,
            const LinkEndType originatingLinkEndRole = unidentified_link_end ):
        ObservationDependentVariableSettings( variableType, relevantLinkEnd, linkEndRole, LinkEndId( "", "" ), originatingLinkEndRole ),
//        relevantLinkEnd_( relevantLinkEnd ), linkEndRole_( linkEndRole ),
        integratedObservableHandling_( integratedObservableHandling ),
//        originatingLinkEndRole_( originatingLinkEndRole ),
        isLinkEndDefined_( true )
    {

    }

    StationAngleObservationDependentVariableSettings(
        const ObservationDependentVariables variableType,
        const LinkEndType linkEndRole,
        const IntegratedObservationPropertyHandling integratedObservableHandling = interval_start, // interval_undefined,
        const LinkEndType originatingLinkEndRole = unidentified_link_end ):
        ObservationDependentVariableSettings( variableType, LinkEndId( "", "" ), linkEndRole, LinkEndId( "", "" ), originatingLinkEndRole ),
//        linkEndRole_( linkEndRole ),
        integratedObservableHandling_( integratedObservableHandling ),
//        originatingLinkEndRole_( originatingLinkEndRole ),
        isLinkEndDefined_( false )
    {

    }

    std::string getIdentifier( )
    {
        std::string identifier;
        if( isLinkEndDefined_ )
        {
            identifier = ", station: (" + linkEndId_.bodyName_ + ", " + linkEndId_.stationName_ + ")";
            if ( linkEndType_ != unidentified_link_end )
            {
                identifier += " as " + getLinkEndTypeString( linkEndType_ );
            }
            if ( originatingLinkEndType_ != unidentified_link_end )
            {
                identifier += " link to " + getLinkEndTypeString( originatingLinkEndType_ );
            }
            identifier += getIntegrationHandlingString( integratedObservableHandling_ );
        }
        else
        {
            identifier += " link end " + getLinkEndTypeString( linkEndType_ );
            if ( originatingLinkEndType_ != unidentified_link_end )
            {
                identifier += " link to " + getLinkEndTypeString( originatingLinkEndType_ );
            }
            identifier += getIntegrationHandlingString( integratedObservableHandling_ );
        }
        return identifier;
    }

    bool areSettingsCompatible( const std::shared_ptr< ObservationDependentVariableSettings > otherSettings )
    {
        bool isCompatible = true;
        std::shared_ptr< StationAngleObservationDependentVariableSettings > stationAngleSettings =
                std::dynamic_pointer_cast< StationAngleObservationDependentVariableSettings >( otherSettings );
        if ( stationAngleSettings == nullptr )
        {
            isCompatible = false;
        }
        else
        {
            if ( !areBaseSettingsCompatible( otherSettings ) )
            {
                isCompatible = false;
            }
            else
            {
                if ( ( integratedObservableHandling_ != stationAngleSettings->integratedObservableHandling_ )
                && ( integratedObservableHandling_ != interval_undefined )
                && ( stationAngleSettings->integratedObservableHandling_ != interval_undefined ) )
                {
                    isCompatible = false;
                }
            }
        }

        return isCompatible;
    }

//    LinkEndId relevantLinkEnd_;

//    LinkEndType linkEndRole_;

    IntegratedObservationPropertyHandling integratedObservableHandling_;

//    LinkEndType originatingLinkEndRole_;

    bool isLinkEndDefined_;

};

class InterlinkObservationDependentVariableSettings: public ObservationDependentVariableSettings
{
public:
    InterlinkObservationDependentVariableSettings(
            const ObservationDependentVariables variableType,
            const LinkEndType startLinkEnd = unidentified_link_end,
            const LinkEndType endLinkEnd = unidentified_link_end,
            const IntegratedObservationPropertyHandling integratedObservableHandling = interval_start, // interval_undefined,
            const std::string relativeBody = "" ):
    ObservationDependentVariableSettings( variableType, LinkEndId( "", "" ), endLinkEnd, LinkEndId( "", "" ), startLinkEnd ),
//        startLinkEnd_( startLinkEnd ), endLinkEnd_( endLinkEnd ),
        integratedObservableHandling_( integratedObservableHandling ),
        relativeBody_( relativeBody ){ }

    ~InterlinkObservationDependentVariableSettings( ){ }

    std::string getIdentifier( )
    {
        std::string identifier = ", link from " + getLinkEndTypeString( originatingLinkEndType_ ) + " to " +
            getLinkEndTypeString( linkEndType_ );
        if( relativeBody_ != "" )
        {
            identifier += " with " + relativeBody_ + " as relative body";
        }
        identifier += getIntegrationHandlingString( integratedObservableHandling_ );

        return identifier;
    }

    bool areSettingsCompatible( const std::shared_ptr< ObservationDependentVariableSettings > otherSettings )
    {
        bool isCompatible = true;
        std::shared_ptr< InterlinkObservationDependentVariableSettings > interlinkSettings =
                std::dynamic_pointer_cast< InterlinkObservationDependentVariableSettings >( otherSettings );
        if ( interlinkSettings == nullptr )
        {
            isCompatible = false;
        }
        else
        {
            if ( !areBaseSettingsCompatible( otherSettings ) )
            {
                isCompatible = false;
            }
            else
            {
                if ( ( integratedObservableHandling_ != interlinkSettings->integratedObservableHandling_ )
                     && ( integratedObservableHandling_ != interval_undefined )
                     && ( interlinkSettings->integratedObservableHandling_ != interval_undefined ) )
                {
                    isCompatible = false;
                }
                if ( ( relativeBody_ != interlinkSettings->relativeBody_ ) && ( relativeBody_ != "" ) && ( interlinkSettings->relativeBody_ != "" ) )
                {
                    isCompatible = false;
                }
            }
        }

        return isCompatible;
    }

//    LinkEndType startLinkEnd_;

//    LinkEndType endLinkEnd_;

    IntegratedObservationPropertyHandling integratedObservableHandling_;

    std::string relativeBody_;
};


std::string getObservationDependentVariableName(
        const ObservationDependentVariables variableType );

std::string getObservationDependentVariableId(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );

bool isObservationDependentVariableVectorial(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );


bool isObservationDependentVariableAncilliarySetting(
    const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );


bool isObservationDependentVariableGroundStationProperty(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );

bool isObservationDependentVariableSimpleLinkProperty(
    const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );

int getObservationDependentVariableSize(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );

bool doesStationAngleVariableExistForGivenLink(
        const ObservableType observableType,
        const LinkEnds& linkEnds,
        const std::shared_ptr< StationAngleObservationDependentVariableSettings > variableSettings );

bool doesInterlinkVariableExistForGivenLink(
    const ObservableType observableType,
    const LinkEnds& linkEnds,
    const std::shared_ptr< InterlinkObservationDependentVariableSettings > variableSettings );

bool doesObservationDependentVariableExistForGivenLink(
        const ObservableType observableType,
        const LinkEnds& linkEnds,
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );


inline std::shared_ptr< ObservationDependentVariableSettings > elevationAngleAtLinkEndTypeDependentVariable(
        const LinkEndType linkEndRole,
        const IntegratedObservationPropertyHandling integratedObservableHandling = interval_undefined,
        const LinkEndType originatingLinkEndRole = unidentified_link_end )
{
    return std::make_shared< StationAngleObservationDependentVariableSettings >( station_elevation_angle, linkEndRole, integratedObservableHandling, originatingLinkEndRole );
}

inline std::shared_ptr< ObservationDependentVariableSettings > azimuthAngleAtLinkEndTypeDependentVariable(
    const LinkEndType linkEndRole,
    const IntegratedObservationPropertyHandling integratedObservableHandling = interval_undefined,
    const LinkEndType originatingLinkEndRole = unidentified_link_end )
{
    return std::make_shared< StationAngleObservationDependentVariableSettings >( station_azimuth_angle, linkEndRole, integratedObservableHandling, originatingLinkEndRole );
}

inline std::shared_ptr< ObservationDependentVariableSettings > targetRangeBetweenLinkEndsDependentVariable(
        const LinkEndType startLinkEnd,
        const LinkEndType endLinkEnd,
        const IntegratedObservationPropertyHandling integratedObservableHandling = interval_undefined )
{
    return std::make_shared< InterlinkObservationDependentVariableSettings >( target_range, startLinkEnd, endLinkEnd, integratedObservableHandling );
}

inline std::shared_ptr< ObservationDependentVariableSettings > bodyAvoidanceAngleDependentVariable(
        const LinkEndType startLinkEnd,
        const LinkEndType endLinkEnd,
        const std::string relativeBody,
        const IntegratedObservationPropertyHandling integratedObservableHandling = interval_undefined )
{
    return std::make_shared< InterlinkObservationDependentVariableSettings >( body_avoidance_angle_variable, startLinkEnd, endLinkEnd, integratedObservableHandling, relativeBody );
}

inline std::shared_ptr< ObservationDependentVariableSettings > linkBodyCenterDistanceDependentVariable(
        const LinkEndType startLinkEnd,
        const LinkEndType endLinkEnd,
        const std::string relativeBody,
        const IntegratedObservationPropertyHandling integratedObservableHandling = interval_undefined )
{
    return std::make_shared< InterlinkObservationDependentVariableSettings >( link_body_center_distance, startLinkEnd, endLinkEnd, integratedObservableHandling, relativeBody );
}

inline std::shared_ptr< ObservationDependentVariableSettings > linkLimbDistanceDependentVariable(
        const LinkEndType startLinkEnd,
        const LinkEndType endLinkEnd,
        const std::string relativeBody,
        const IntegratedObservationPropertyHandling integratedObservableHandling = interval_undefined )
{
    return std::make_shared< InterlinkObservationDependentVariableSettings >( link_limb_distance, startLinkEnd, endLinkEnd, integratedObservableHandling, relativeBody );
}

inline std::shared_ptr< ObservationDependentVariableSettings > linkAngleWrtOrbitalPlaneDependentVariable(
        const LinkEndType startLinkEnd,
        const LinkEndType endLinkEnd,
        const std::string relativeBody,
        const IntegratedObservationPropertyHandling integratedObservableHandling = interval_undefined )
{
    return std::make_shared< InterlinkObservationDependentVariableSettings >( link_angle_with_orbital_plane, startLinkEnd, endLinkEnd, integratedObservableHandling, relativeBody );
}

//inline std::shared_ptr< ObservationDependentVariableSettings > dopplerIntegrationTimeDependentVariable( )
//{
//    return std::make_shared< ObservationDependentVariableSettings >( doppler_integration_time_dependent_variable );
//}
//
//inline std::shared_ptr< ObservationDependentVariableSettings > retransmissionDelaysDependentVariable( )
//{
//    return std::make_shared< ObservationDependentVariableSettings >( retransmission_delays_dependent_variable );
//}

//// for interlink
//inline std::shared_ptr< ObservationDependentVariableSettings > dependentVariable(
//        const ObservationDependentVariables variableType,
//        const std::pair< LinkEndId, LinkEndId > linkEnds =
//                std::make_pair( LinkEndId( "", "" ), LinkEndId( "", "" ) ),/*
//        const std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( )*/ )
//{
//
//}
//
//// for interlink
//inline std::shared_ptr< ObservationDependentVariableSettings > dependentVariable(
//        const ObservationDependentVariables variableType,
//        const std::pair< LinkEndType, LinkEndType > linkEnds,
//        const std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) )
//{
//
//}
//
//
//// for station angle
//inline std::shared_ptr< ObservationDependentVariableSettings > dependentVariable(
//        const ObservationDependentVariables variableType,
//        const LinkEndId linkEnds,
//        const std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) )
//{
//
//}
//
//inline std::shared_ptr< ObservationDependentVariableSettings > dependentVariable(
//        const ObservationDependentVariables variableType,
//        const LinkEndType linkEnds,
//        const std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) )
//{
//
//}

//std::shared_ptr< ObservationCollectionParser > getObservationParserFromDependentVariableSettings(
//        const std::shared_ptr< ObservationDependentVariableSettings > dependentVariableSettings )
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
//    // if only originating link end id is specified
//    else if ( dependentVariableSettings->originatingLinkEndId_ != LinkEndId( "", "" ) )
//    {
//        parserList.push_back( std::make_shared< ObservationCollectionLinkEndIdParser >( dependentVariableSettings->originatingLinkEndId_ ) );
//    }
//    // if only originating link end type is specified
//    else if ( dependentVariableSettings->originatingLinkEndType_ != unidentified_link_end )
//    {
//        parserList.push_back( std::make_shared< ObservationCollectionLinkEndTypeParser >( dependentVariableSettings->originatingLinkEndType_ ) );
//    }
//
//    // Create multi-type observation collection parser
//    observationParser = std::make_shared< ObservationCollectionMultiTypeParser >( parserList, true );
//
//    return observationParser;
//}

//obsCollection->addDependentVariable( elevation_angle( receiver ) );
//Eigen::MatrixXd = obsCollection->getDependentVariable( elevation_angle( "sTATION1" ) );
//std::vector< Eigen::MatrixXd > test = obsCollection->getDependentVariable( elevation_angle( "sTATION1" ) );
//
//dependentVariable( target_range );
//dependentVariable( target_range, std::make_pair( receiver, reflector1 ) );
//dependentVariable( target_range, std::make_pair( "station1", "spacecraft" ) );
//
//dependentVariable( elevation_angle );
//dependentVariable( elevation_angle, receiver );
//dependentVariable( elevation_angle, "station1" );
//
//observationCollection->addDependentVariable( elevationAngle( receiver ), observationParser( one_way_doppler ) );




// type -> always needs to be provided
// linkEndsId

}

}
#endif // TUDAT_OBSERVATIONOUTPUTSETTINGS
