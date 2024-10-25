/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/estimation_setup/observationOutputSettings.h"
#include "tudat/astro/observation_models/observationViabilityCalculator.h"

namespace tudat
{

namespace simulation_setup
{


std::string getIntegrationHandlingString( const IntegratedObservationPropertyHandling integratedObservableHandling )
{
    std::string identifier = "";
    if( integratedObservableHandling != interval_undefined )
    {
        if( integratedObservableHandling == interval_start )
        {
            identifier += ", start of integration interval";
        }
        else if( integratedObservableHandling == interval_end )
        {
            identifier += ", end of integration interval";
        }
    }
    return identifier;

}

//! Function checking whether the interlinks between two link ends are compatible (i.e., for both the originating and receiving ends of the interlink,
//! the link end's types and IDs should be either identical, or undefined).
bool areInterlinksCompatible(
        const std::pair< LinkEndType, LinkEndId >& firstReceivingLinkEnd,
        const std::pair< LinkEndType, LinkEndId >& firstOriginatingLinkEnd,
        const std::pair< LinkEndType, LinkEndId >& secondReceivingLinkEnd,
        const std::pair< LinkEndType, LinkEndId >& secondOriginatingLinkEnd )
{
    // Check if the receiving link ends (both types and link IDs) are either identical or undefined (in which case the two settings would be considered compatible).
    bool compatibleReceivingEndTypes = ( firstReceivingLinkEnd.first == secondReceivingLinkEnd.first
            || firstReceivingLinkEnd.first == unidentified_link_end || secondReceivingLinkEnd.first == unidentified_link_end );
    bool compatibleReceivingEndIds = ( firstReceivingLinkEnd.second == secondReceivingLinkEnd.second
            || firstReceivingLinkEnd.second == LinkEndId( "", "" ) || secondReceivingLinkEnd.second == LinkEndId( "", "" ) );

    // Check if the originating link ends (both types and link IDs) are either identical or undefined (in which case the two settings would be considered compatible).
    bool compatibleOriginatingEndTypes = ( firstOriginatingLinkEnd.first == secondOriginatingLinkEnd.first
            || firstOriginatingLinkEnd.first == unidentified_link_end || secondOriginatingLinkEnd.first == unidentified_link_end );
    bool compatibleOriginatingEndIds = ( firstOriginatingLinkEnd.second == secondOriginatingLinkEnd.second
            || firstOriginatingLinkEnd.second == LinkEndId( "", "" ) || secondOriginatingLinkEnd.second == LinkEndId( "", "" ) );

    // Check whether both link ends' types and IDs match
    bool areLinkEndsCompatible = ( compatibleReceivingEndTypes && compatibleReceivingEndIds && compatibleOriginatingEndTypes && compatibleOriginatingEndIds );

    return areLinkEndsCompatible;
}

//! Function returning the dependent variable name
std::string getObservationDependentVariableName(
        const ObservationDependentVariables variableType )
{
    std::string dependentVariableName;
    switch( variableType )
    {
    case station_elevation_angle:
    {
        dependentVariableName = "Station elevation angle ";
        break;
    }
    case station_azimuth_angle:
    {
        dependentVariableName = "Station azimuth angle ";
        break;
    }
    case target_range:
    {
        dependentVariableName = "Range between link ends ";
        break;
    }
    case body_avoidance_angle_variable:
    {
        dependentVariableName = "Body avoidance angle ";
        break;
    }
    case link_body_center_distance:
    {
        dependentVariableName = "Link to body center distance ";
        break;
    }
    case link_limb_distance:
    {
        dependentVariableName = "Link to body limb distance ";
        break;
    }
    case link_angle_with_orbital_plane:
    {
        dependentVariableName = "Angle between link vector and orbital plane ";
        break;
    }
    case integration_time_dependent_variable:
    {
        dependentVariableName = "Integration count time ";
        break;
    }
    case retransmission_delays_dependent_variable:
    {
        dependentVariableName = "Retransmission delays ";
        break;
    }
    default:
        throw std::runtime_error( "Error when checking observation dependent variable. Type " +
                                  std::to_string( variableType ) +
                                  " not found when retrieving variable name." );
    }
    return dependentVariableName;
}

//! Function that returns a string uniquely describing a dependent variable settings object
std::string getObservationDependentVariableId(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings )
{
    return getObservationDependentVariableName( variableSettings->variableType_ ) +
            variableSettings->getIdentifier( );
}

//! Function checking whether a given dependent variable should be given as a vector of size > 1
bool isObservationDependentVariableVectorial( const ObservationDependentVariables variableType )
{
    bool isVariableVectorial = false;
    switch( variableType )
    {
    case station_elevation_angle:
        isVariableVectorial = false;
        break;
    case station_azimuth_angle:
        isVariableVectorial = false;
        break;
    case target_range:
        isVariableVectorial = false;
        break;
    case body_avoidance_angle_variable:
        isVariableVectorial = false;
        break;
    case integration_time_dependent_variable:
        isVariableVectorial = false;
        break;
    case link_body_center_distance:
        isVariableVectorial = false;
        break;
    case link_limb_distance:
        isVariableVectorial = false;
        break;
    case link_angle_with_orbital_plane:
        isVariableVectorial = false;
        break;
    case retransmission_delays_dependent_variable:
        isVariableVectorial = true;
        break;
    default:
        throw std::runtime_error( "Error when checking observation dependent variable. Type " + getObservationDependentVariableName( variableType ) +
                                  " not found when checking if variable is vectorial." );
    }
    return isVariableVectorial;
}

//! Function checking whether a given dependent variable type is related to some ancillary settings
bool isObservationDependentVariableAncilliarySetting(
    const ObservationDependentVariables variableType )
{
    bool isAncilliarySetting = false;
    switch( variableType )
    {
    case station_elevation_angle:
        break;
    case station_azimuth_angle:
        break;
    case target_range:
        break;
    case body_avoidance_angle_variable:
        break;
    case link_body_center_distance:
        break;
    case link_limb_distance:
        break;
    case link_angle_with_orbital_plane:
        break;
    case integration_time_dependent_variable:
        isAncilliarySetting = true;
        break;
    case retransmission_delays_dependent_variable:
        isAncilliarySetting = true;
        break;
    default:
        throw std::runtime_error( "Error when checking observation dependent variable. Type " + getObservationDependentVariableName( variableType ) +
                                  " not found when checking for ancilliary setting." );

    }
    return isAncilliarySetting;
}

//! Function checking whether a given dependent variable is related to a ground station property
bool isObservationDependentVariableGroundStationProperty( const ObservationDependentVariables variableType )
{
    bool isGroundStationProperty = false;
    switch( variableType )
    {
    case station_elevation_angle:
        isGroundStationProperty = true;
        break;
    case station_azimuth_angle:
        isGroundStationProperty = true;
        break;
    case target_range:
        break;
    case body_avoidance_angle_variable:
        break;
    case integration_time_dependent_variable:
        break;
    case link_body_center_distance:
        break;
    case link_limb_distance:
        break;
    case link_angle_with_orbital_plane:
        break;
    case retransmission_delays_dependent_variable:
        break;
    default:
        throw std::runtime_error( "Error when checking observation dependent variable. Type " + getObservationDependentVariableName( variableType ) +
                                  " not found when checking for ground station dependency." );
    }
    return isGroundStationProperty;
}

//! Function checking whether a given dependent variable type is an interlink property
bool isObservationDependentVariableInterlinkProperty( const ObservationDependentVariables variableType )
{
    bool isInterlinkProperty = false;
    switch( variableType )
    {
        case station_elevation_angle:
            break;
        case station_azimuth_angle:
            break;
        case target_range:
            isInterlinkProperty = true;
            break;
        case body_avoidance_angle_variable:
            isInterlinkProperty = true;
            break;
        case integration_time_dependent_variable:
            break;
        case link_body_center_distance:
            isInterlinkProperty = true;
            break;
        case link_limb_distance:
            isInterlinkProperty = true;
            break;
        case link_angle_with_orbital_plane:
            isInterlinkProperty = true;
            break;
        case retransmission_delays_dependent_variable:
            break;
        default:
            throw std::runtime_error( "Error when checking observation dependent variable. Type " + getObservationDependentVariableName( variableType ) +
                                      " not found when checking if interlink property." );
    }
    return isInterlinkProperty;
}

//! Function that returns the size of a given dependent variable (can be link end-dependent for ancillary settings dependent variables)
int getObservationDependentVariableSize(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings,
        const LinkEnds linkEnds )
{
    int variableSize = 0;
    if( !isObservationDependentVariableVectorial( variableSettings->variableType_ ) )
    {
        variableSize = 1;
    }
    else
    {
        switch( variableSettings->variableType_ )
        {
            case retransmission_delays_dependent_variable:
                variableSize = linkEnds.size( ) - 2;
                break;
            default:
                throw std::runtime_error( "Error when checking observation dependent variable. Type " + getObservationDependentVariableId( variableSettings ) +
                " not found when determining parameter size." );

        }
    }
    return variableSize;
}

//! Function that checks whether a given station angle dependent variable can be computed for a given observable type and link ends.
bool doesStationAngleVariableExistForGivenLink(
        const observation_models::ObservableType observableType,
        const observation_models::LinkEnds& linkEnds,
        const std::shared_ptr< StationAngleObservationDependentVariableSettings > variableSettings )
{
    bool doesLinkHaveDependency = false;

    if( linkEnds.size( ) > 1 )
    {
        if( variableSettings->isLinkEndDefined_ )
        {
            std::vector<observation_models::LinkEndType> linkEndTypeList = getLinkEndTypesForGivenLinkEndId(
                linkEnds, variableSettings->linkEndId_ );
            if ( linkEndTypeList.size( ) > 0 )
            {
                doesLinkHaveDependency = true;
            }
        }
        else
        {
            if( linkEnds.count( variableSettings->linkEndType_ ) > 0 )
            {
                doesLinkHaveDependency = true;
            }
        }
    }
    return doesLinkHaveDependency;
}

//! Function that checks whether a given interlink dependent variable can be computed for a given observable type and link ends.
bool doesInterlinkVariableExistForGivenLink(
    const observation_models::ObservableType observableType,
    const observation_models::LinkEnds& linkEnds,
    const std::shared_ptr< InterlinkObservationDependentVariableSettings > variableSettings )
{
    bool doesLinkHaveDependency = true;
    if( variableSettings->originatingLinkEndType_ != observation_models::unidentified_link_end )
    {
        if( linkEnds.count( variableSettings->originatingLinkEndType_ ) == 0 )
        {
            doesLinkHaveDependency = false;
        }
    }

    if( variableSettings->linkEndType_ != observation_models::unidentified_link_end )
    {
        if( linkEnds.count( variableSettings->linkEndType_ ) == 0 )
        {
            doesLinkHaveDependency = false;
        }
    }

    return doesLinkHaveDependency;
}

//! Function that checks whether a given dependent variable can be computed for a given observable type and link ends.
bool doesObservationDependentVariableExistForGivenLink(
        const observation_models::ObservableType observableType,
        const observation_models::LinkEnds& linkEnds,
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings )
{
    bool doesLinkHaveDependency = false;
    switch( variableSettings->variableType_ )
    {
    case station_elevation_angle:
        doesLinkHaveDependency = doesStationAngleVariableExistForGivenLink(
                   observableType, linkEnds, std::dynamic_pointer_cast< StationAngleObservationDependentVariableSettings >(
                        variableSettings ) );
        break;
    case station_azimuth_angle:
        doesLinkHaveDependency = doesStationAngleVariableExistForGivenLink(
                   observableType, linkEnds, std::dynamic_pointer_cast< StationAngleObservationDependentVariableSettings >(
                        variableSettings ) );
        break;
    case target_range:
        doesLinkHaveDependency = doesInterlinkVariableExistForGivenLink(
            observableType, linkEnds, std::dynamic_pointer_cast< InterlinkObservationDependentVariableSettings >(
                variableSettings ) );
        break;
    case body_avoidance_angle_variable:
        doesLinkHaveDependency = doesInterlinkVariableExistForGivenLink(
            observableType, linkEnds, std::dynamic_pointer_cast< InterlinkObservationDependentVariableSettings >(
                variableSettings ) );
        break;
    case link_body_center_distance:
        doesLinkHaveDependency = doesInterlinkVariableExistForGivenLink(
            observableType, linkEnds, std::dynamic_pointer_cast< InterlinkObservationDependentVariableSettings >(
                variableSettings ) );
        break;
    case link_limb_distance:
        doesLinkHaveDependency = doesInterlinkVariableExistForGivenLink(
            observableType, linkEnds, std::dynamic_pointer_cast< InterlinkObservationDependentVariableSettings >(
                variableSettings ) );
        break;
    case link_angle_with_orbital_plane:
        doesLinkHaveDependency = doesInterlinkVariableExistForGivenLink(
            observableType, linkEnds, std::dynamic_pointer_cast< InterlinkObservationDependentVariableSettings >(
                variableSettings ) );
        break;
    case integration_time_dependent_variable:
        doesLinkHaveDependency = true;
        break;
    case retransmission_delays_dependent_variable:
        doesLinkHaveDependency = true;
        break;
    default:
        throw std::runtime_error( "Error when checking observation dependent variable. Type " +
                                  getObservationDependentVariableId( variableSettings ) +
                                  " not found when checking if variable exists for given link." );

    }
    return doesLinkHaveDependency;
}

//! Function checking whether a given dependent variable is link end-dependent (false for ancillary settings dependent variables)
bool isObservationDependentVariableLinkEndDependent( const ObservationDependentVariables variableType )
{
    bool linkEndDependent = true;
    switch( variableType )
    {
        case station_elevation_angle:
        case station_azimuth_angle:
        case target_range:
        case body_avoidance_angle_variable:
        case link_body_center_distance:
        case link_limb_distance:
        case link_angle_with_orbital_plane:
            break;
        case integration_time_dependent_variable:
            linkEndDependent = false;
            break;
        case retransmission_delays_dependent_variable:
            linkEndDependent = false;
            break;
        default:
            throw std::runtime_error( "Error when checking observation dependent variable. Type " + getObservationDependentVariableName( variableType ) +
                                      " not found when checking if variable is link end dependent." );
    }
    return linkEndDependent;
}

//! Function checking whether a given interlink dependent variable depends on the link direction (i.e., depends on which link end is the receiving/originating one)
bool isInterlinkPropertyDirectionAgnostic( const ObservationDependentVariables variableType )
{
    if ( !isObservationDependentVariableInterlinkProperty( variableType ) )
    {
        throw std::runtime_error( "Error when checking if interlink dependent variable is link direction-agnostic, type " +
                                          getObservationDependentVariableName( variableType ) + " is not interlink property." );
    }
    bool isDirectionAgnostic = false;
    switch( variableType )
    {
        case target_range:
            isDirectionAgnostic = true;
            break;
        case body_avoidance_angle_variable:
            break;
        case link_body_center_distance:
            isDirectionAgnostic = true;
            break;
        case link_limb_distance:
            isDirectionAgnostic = true;
            break;
        case link_angle_with_orbital_plane:
            break;
        default:
            throw std::runtime_error( "Error when checking observation dependent variable. Type " + getObservationDependentVariableName( variableType ) +
                                      " not found when checking if interlink variable is link direction-agnostic." );
    }
    return isDirectionAgnostic;
}

//! Returns a function which checks whether an ancillary settings dependent variable exists for a given observable type (dependent on whether
//! such ancillary settings are defined for a specific observable type)
std::function< bool( const ObservableType observableType ) > getIsObservableTypeCompatibleFunction(
        const ObservationDependentVariables variableType )
{
    std::function< bool( const ObservableType observableType ) > isObservableTypeCompatibleFunction;
    if ( !isObservationDependentVariableAncilliarySetting( variableType ) )
    {
        throw std::runtime_error( "Error when retrieving function defining whether observable type is compatible with ancillary settings, the input dependent variable "
                                  "is not of ancillary settings type." );
    }
    switch( variableType )
    {
        case integration_time_dependent_variable:
        {
            isObservableTypeCompatibleFunction = std::bind( &observation_models::isObservableOfIntegratedType, std::placeholders::_1 );
            break;
        }
        case retransmission_delays_dependent_variable:
        {
            isObservableTypeCompatibleFunction = std::bind( &observation_models::observableCanHaveRetransmissionDelay, std::placeholders::_1 );
            break;
        }
        default:
            throw std::runtime_error( "Error when checking observation dependent variable. Type " + std::to_string( variableType ) +
                                      " not found when retriving function defining whether observable type is compatible with requested ancillary settings." );
    }

    return isObservableTypeCompatibleFunction;
}

//! Function that returns a fully defined ObservationDependentVariableSettings object, using the original settings (possibly incomplete, i.e.
//! with missing link ends information) and specific link ends information for both the receiving and originating ends of the link
std::shared_ptr< ObservationDependentVariableSettings > createCompleteObservationDependentVariableSettings(
        const std::shared_ptr< ObservationDependentVariableSettings > originalSettings,
        const LinkEndType& linkEndType, const LinkEndId& linkEndId, const LinkEndType& originatingLinkEndType, const LinkEndId& originatingLinkEndId )
{
    std::shared_ptr< ObservationDependentVariableSettings > completeSettings;

    if ( std::dynamic_pointer_cast< StationAngleObservationDependentVariableSettings >( originalSettings ) != nullptr )
    {
        // Create fully defined station angle dependent variable settings
        std::shared_ptr< StationAngleObservationDependentVariableSettings > stationAngleSettings =
                std::dynamic_pointer_cast< StationAngleObservationDependentVariableSettings >( originalSettings );
        completeSettings = std::make_shared< StationAngleObservationDependentVariableSettings >(
                originalSettings->variableType_, linkEndId, linkEndType, originatingLinkEndId, originatingLinkEndType,
                stationAngleSettings->integratedObservableHandling_ );
    }
    else if ( std::dynamic_pointer_cast< InterlinkObservationDependentVariableSettings >( originalSettings ) != nullptr )
    {
        // Create fully defined interlink property dependent variable settings
        std::shared_ptr< InterlinkObservationDependentVariableSettings > interlinkSettings =
                std::dynamic_pointer_cast< InterlinkObservationDependentVariableSettings >( originalSettings );
        completeSettings = std::make_shared< InterlinkObservationDependentVariableSettings >(
                originalSettings->variableType_, originatingLinkEndType, linkEndType, originatingLinkEndId, linkEndId,
                interlinkSettings->integratedObservableHandling_, interlinkSettings->relativeBody_ );
    }
    else
    {
        completeSettings = std::make_shared< ObservationDependentVariableSettings >( originalSettings->variableType_, linkEndId, linkEndType,
                                                                                     originatingLinkEndId, originatingLinkEndType );
    }

    return completeSettings;
}

//! Function that returns a list of all compatible dependent variable settings that can be created for a given observable type and link ends, from a base
//! dependent variable settings that might not be entirely defined (i.e., some link ends not specified, etc.)
std::vector< std::shared_ptr< ObservationDependentVariableSettings > > createAllCompatibleDependentVariableSettings(
        const ObservableType observableType,
        const LinkEnds& linkEnds,
        std::shared_ptr< ObservationDependentVariableSettings > dependentVariableSettings )
{
    // Retrieve relevant information about all interlinks for current observable type and link ends
    std::vector< std::pair< std::pair< LinkEndType, LinkEndId >, std::pair< LinkEndType, LinkEndId > > > interlinksInSet = getInterlinks( observableType, linkEnds );

    // Get interlink requirements from dependent variable settings
    std::pair< std::pair< LinkEndType, LinkEndId >, std::pair< LinkEndType, LinkEndId > > interlinksSettings =
            std::make_pair( std::make_pair( dependentVariableSettings->linkEndType_, dependentVariableSettings->linkEndId_ ),
                            std::make_pair( dependentVariableSettings->originatingLinkEndType_, dependentVariableSettings->originatingLinkEndId_ ) );

    std::vector< std::shared_ptr< ObservationDependentVariableSettings > > allDependentVariablesSettings;
    // If the dependent variable is independent of the link ends
    if ( !isObservationDependentVariableLinkEndDependent( dependentVariableSettings->variableType_ ) )
    {
        // If the dependent variable refers to an observation ancillary data, check whether it is compatible with the current observable type
        if ( std::dynamic_pointer_cast< AncillaryObservationDependentVariableSettings >( dependentVariableSettings ) != nullptr )
        {
            std::shared_ptr< AncillaryObservationDependentVariableSettings > ancillaryDependentVariableSettings =
                    std::dynamic_pointer_cast< AncillaryObservationDependentVariableSettings >( dependentVariableSettings );
            if ( ancillaryDependentVariableSettings->isObservableTypeCompatible_( observableType ) )
            {
                std::shared_ptr< AncillaryObservationDependentVariableSettings > completeAncillaryDependentVariableSettings =
                        std::make_shared< AncillaryObservationDependentVariableSettings >( ancillaryDependentVariableSettings->variableType_, observableType );
                allDependentVariablesSettings.push_back( completeAncillaryDependentVariableSettings );
            }
        }
        else
        {
            allDependentVariablesSettings.push_back( dependentVariableSettings );
        }
    }
    // For link end-dependent variables, all possible links must be considered and a separate dependent variable settings object is created for each of them.
    else
    {
        std::vector< std::pair< std::pair< LinkEndType, LinkEndId >, std::pair< LinkEndType, LinkEndId > > > interlinksToCreateList;
        for ( auto interlink : interlinksInSet )
        {
            std::pair< LinkEndType, LinkEndId > receivingLinkEnd = interlink.first;
            std::pair< LinkEndType, LinkEndId > originatingLinkEnd = interlink.second;
            std::pair< LinkEndType, LinkEndId > receivingLinkEndInSettings = interlinksSettings.first;
            std::pair< LinkEndType, LinkEndId > originatingLinkEndInSettings = interlinksSettings.second;

            // Check if the link ends (both receiving and originating ends) are either identical or undefined (in which case the two settings would be considered compatible).
            bool directLinksMatch = areInterlinksCompatible( receivingLinkEnd, originatingLinkEnd, receivingLinkEndInSettings, originatingLinkEndInSettings );

            // Check if inverting the receiving/originating ends of the link would lead to compatible link definitions (for dependent variables that are independent of the
            // link "direction").
            bool revertedLinksMatch = areInterlinksCompatible( receivingLinkEnd, originatingLinkEnd, originatingLinkEndInSettings, receivingLinkEndInSettings );
            if ( revertedLinksMatch )
            {
                interlink = std::make_pair( interlink.second, interlink.first );
            }


            if ( directLinksMatch || revertedLinksMatch )
            {
                // If dependent variable is not a station property, include all interlinks
                if ( !isObservationDependentVariableGroundStationProperty( dependentVariableSettings->variableType_ ) )
                {
                    interlinksToCreateList.push_back( interlink );
                }
                // If dependent variable is a station property, check that the first link end has a ground station defined
                else
                {
                    // if station defined for start link end
                    if ( interlink.first.second.stationName_ != "" )
                    {
                        interlinksToCreateList.push_back( interlink );
                    }
                    // if station only defined for end link end
                    else if ( interlink.second.second.stationName_ != "" )
                    {
                        // Reverse link order
                        interlinksToCreateList.push_back( std::make_pair( interlink.second, interlink.first ) );
                    }
                }
            }
        }

        // Create fully defined dependent variable settings for all relevant links.
        for ( auto interlink : interlinksToCreateList )
        {
            std::shared_ptr< ObservationDependentVariableSettings > completeSettings = createCompleteObservationDependentVariableSettings(
                    dependentVariableSettings, interlink.first.first, interlink.first.second, interlink.second.first, interlink.second.second );

            allDependentVariablesSettings.push_back( completeSettings );
        }
    }

    return allDependentVariablesSettings;

}

}

}
