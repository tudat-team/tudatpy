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

std::string getObservationDependentVariableName( const ObservationDependentVariables variableType );

bool isObservationDependentVariableAncilliarySetting( const ObservationDependentVariables variableType );

bool isObservationDependentVariableInterlinkProperty( const ObservationDependentVariables variableType );

bool isInterlinkPropertyDirectionAgnostic( const ObservationDependentVariables variableType );

bool isObservationDependentVariableVectorial( const ObservationDependentVariables variableType );

bool isObservationDependentVariableGroundStationProperty( const ObservationDependentVariables variableType );

bool isObservationDependentVariableLinkEndDependent( const ObservationDependentVariables variableType );


class ObservationDependentVariableSettings
{
public:
    ObservationDependentVariableSettings(
            const ObservationDependentVariables variableType,
            const LinkEndId linkEndId = LinkEndId( "", "" ),
            const LinkEndType linkEndType = unidentified_link_end,
            const LinkEndId originatingLinkEndId = LinkEndId( "", "" ),
            const LinkEndType originatingLinkEndType = unidentified_link_end ):
        variableType_( variableType ), linkEndId_( linkEndId ), linkEndType_( linkEndType ),
        originatingLinkEndId_( originatingLinkEndId ), originatingLinkEndType_( originatingLinkEndType ){ }

    virtual ~ObservationDependentVariableSettings( ){ }

    ObservationDependentVariables variableType_;

    std::string getBaseIdentifier( )
    {
        std::string identifier = ", link end";
        if( linkEndId_ != LinkEndId( "", "" ) && linkEndType_ != unidentified_link_end )
        {
            identifier += ": (" + linkEndId_.bodyName_ + ", " + linkEndId_.stationName_ + ") as " + getLinkEndTypeString( linkEndType_ );
        }
        else if ( linkEndId_ != LinkEndId( "", "" ) )
        {
            identifier = ": (" + linkEndId_.bodyName_ + ", " + linkEndId_.stationName_ + ")";
        }
        else if ( linkEndType_ != unidentified_link_end )
        {
            identifier = " of type " + getLinkEndTypeString( linkEndType_ );
        }

        if( originatingLinkEndId_ != LinkEndId( "", "" ) && originatingLinkEndType_ != unidentified_link_end )
        {
            identifier += " to link end: (" + originatingLinkEndId_.bodyName_ + ", " + originatingLinkEndId_.stationName_ + ") as " + getLinkEndTypeString( originatingLinkEndType_ );
        }
        else if ( originatingLinkEndId_ != LinkEndId( "", "" ) )
        {
            identifier += " to link end: (" + originatingLinkEndId_.bodyName_ + ", " + originatingLinkEndId_.stationName_ + ")";
        }
        else if ( originatingLinkEndType_ != unidentified_link_end )
        {
            identifier += " to link end of type " + getLinkEndTypeString( originatingLinkEndType_ );
        }
        return identifier;
    }

    virtual std::string getIdentifier( )
    {
        return getBaseIdentifier( );
    }

    virtual bool areSettingsCompatible( const std::shared_ptr< ObservationDependentVariableSettings > otherSettings )
    {
        return areBaseSettingsCompatible( otherSettings, true );
    }

    bool areBaseSettingsCompatible( const std::shared_ptr< ObservationDependentVariableSettings > otherSettings,
                                    const bool revertedLinksAllowed = true )
    {
        bool isCompatible = true;
        if ( variableType_ != otherSettings->variableType_ )
        {
            isCompatible = false;
        }
        else
        {
            bool directLinksMatch =
                    ( ( linkEndType_ == otherSettings->linkEndType_ || linkEndType_ == unidentified_link_end || otherSettings->linkEndType_ == unidentified_link_end )
                    && ( ( linkEndId_ == otherSettings->linkEndId_ || linkEndId_ == LinkEndId( "", "" ) || otherSettings->linkEndId_ == LinkEndId( "", "" ) ) )
                    && ( ( originatingLinkEndType_ == otherSettings->originatingLinkEndType_ ||
                    originatingLinkEndType_ == unidentified_link_end || otherSettings->originatingLinkEndType_ == unidentified_link_end ) )
                    && ( ( originatingLinkEndId_ == otherSettings->originatingLinkEndId_ ||
                    originatingLinkEndId_ == LinkEndId( "", "" ) || otherSettings->originatingLinkEndId_ == LinkEndId( "", "" ) ) ) );

            bool revertedLinksMatch =
                    ( ( linkEndType_ == otherSettings->originatingLinkEndType_ || linkEndType_ == unidentified_link_end || otherSettings->originatingLinkEndType_ == unidentified_link_end )
                    && ( ( linkEndId_ == otherSettings->originatingLinkEndId_ || linkEndId_ == LinkEndId( "", "" ) || otherSettings->originatingLinkEndId_ == LinkEndId( "", "" ) ) )
                    && ( ( originatingLinkEndType_ == otherSettings->linkEndType_ ||
                    originatingLinkEndType_ == unidentified_link_end || otherSettings->linkEndType_ == unidentified_link_end ) )
                    && ( ( originatingLinkEndId_ == otherSettings->linkEndId_ ||
                    originatingLinkEndId_ == LinkEndId( "", "" ) || otherSettings->linkEndId_ == LinkEndId( "", "" ) ) ) );

            if ( !directLinksMatch && ( !revertedLinksAllowed || !revertedLinksMatch ) )
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
            const LinkEndId originatingLinkEndId = LinkEndId( "", "" ),
            const LinkEndType originatingLinkEndRole = unidentified_link_end,
            const IntegratedObservationPropertyHandling integratedObservableHandling = interval_start ):
        ObservationDependentVariableSettings( variableType, relevantLinkEnd, linkEndRole, originatingLinkEndId, originatingLinkEndRole ),
        integratedObservableHandling_( integratedObservableHandling ),
        isLinkEndDefined_( ( relevantLinkEnd != LinkEndId( "", "" ) ? true : false ) ){ }

    std::string getIdentifier( )
    {
        return getBaseIdentifier( ) + getIntegrationHandlingString( integratedObservableHandling_ );
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
            if ( !areBaseSettingsCompatible( otherSettings, false ) )
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

    IntegratedObservationPropertyHandling integratedObservableHandling_;

    bool isLinkEndDefined_;

};

class InterlinkObservationDependentVariableSettings: public ObservationDependentVariableSettings
{
public:
    InterlinkObservationDependentVariableSettings(
            const ObservationDependentVariables variableType,
            const LinkEndType startLinkEndType = unidentified_link_end,
            const LinkEndType endLinkEndType = unidentified_link_end,
            const LinkEndId startLinkEndId = LinkEndId( "", "" ),
            const LinkEndId endLinkEndId = LinkEndId( "", "" ),
            const IntegratedObservationPropertyHandling integratedObservableHandling = interval_start,
            const std::string relativeBody = "" ):
    ObservationDependentVariableSettings( variableType, endLinkEndId, endLinkEndType, startLinkEndId, startLinkEndType ),
        integratedObservableHandling_( integratedObservableHandling ),
        relativeBody_( relativeBody ){ }

    ~InterlinkObservationDependentVariableSettings( ){ }

    std::string getIdentifier( )
    {
        std::string identifier = getBaseIdentifier( );
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
            if ( !areBaseSettingsCompatible( otherSettings, isInterlinkPropertyDirectionAgnostic( variableType_ ) ) )
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

    IntegratedObservationPropertyHandling integratedObservableHandling_;

    std::string relativeBody_;
};

std::function< bool( const ObservableType observableType ) > getIsObservableTypeCompatibleFunction(
        const ObservationDependentVariables variableType );

class AncillaryObservationDependentVariableSettings: public ObservationDependentVariableSettings
{
public:
    AncillaryObservationDependentVariableSettings(
            const ObservationDependentVariables variableType,
            const ObservableType observableType = undefined_observation_model ):
            ObservationDependentVariableSettings( variableType ), observableType_( observableType )
    {
        isObservableTypeCompatible_ = getIsObservableTypeCompatibleFunction( variableType );
    }

    ~AncillaryObservationDependentVariableSettings( ){ }

    std::string getIdentifier( )
    {
        std::string identifier = getBaseIdentifier( );
        if ( observableType_ != undefined_observation_model )
        {
            identifier += ", for observable of type " + std::to_string( observableType_ );
        }
        return identifier;
    }

    bool areSettingsCompatible( const std::shared_ptr< ObservationDependentVariableSettings > otherSettings )
    {
        bool isCompatible = true;
        std::shared_ptr< AncillaryObservationDependentVariableSettings > ancillarySettings =
                std::dynamic_pointer_cast< AncillaryObservationDependentVariableSettings >( otherSettings );
        if ( ancillarySettings == nullptr )
        {
            isCompatible = false;
        }
        else
        {
            if ( variableType_ != otherSettings->variableType_ )
            {
                isCompatible = false;
            }
            if ( ( observableType_ != ancillarySettings->observableType_ ) && ( observableType_ != undefined_observation_model )
            && ( ancillarySettings->observableType_ != undefined_observation_model ) )
            {
                isCompatible = false;
            }
        }

        return isCompatible;
    }

    ObservableType observableType_;

    std::function< bool( const ObservableType observableType ) > isObservableTypeCompatible_;

};


std::string getObservationDependentVariableId(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );

int getObservationDependentVariableSize(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings,
        const LinkEnds linkEnds );

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

std::shared_ptr< ObservationDependentVariableSettings > createCompleteObservationDependentVariableSettings(
        const std::shared_ptr< ObservationDependentVariableSettings > originalSettings,
        const LinkEndType& linkEndType, const LinkEndId& linkEndId, const LinkEndType& originatingLinkEndType, const LinkEndId& originatingLinkEndId );

std::vector< std::shared_ptr< ObservationDependentVariableSettings > > createAllCompatibleDependentVariableSettings(
        const ObservableType observableType,
        const LinkEnds& linkEnds,
        std::shared_ptr< ObservationDependentVariableSettings > dependentVariableSettings );

inline std::shared_ptr< ObservationDependentVariableSettings > elevationAngleDependentVariable(
        const LinkEndType linkEndRole = unidentified_link_end,
        const LinkEndId linkEndId = LinkEndId( "", "" ),
        const LinkEndType originatingLinkEndRole = unidentified_link_end,
        const LinkEndId originatingLinkEndId = LinkEndId( "", "" ),
        const IntegratedObservationPropertyHandling integratedObservableHandling = interval_start )
{
    return std::make_shared< StationAngleObservationDependentVariableSettings >(
            station_elevation_angle, linkEndId, linkEndRole, originatingLinkEndId, originatingLinkEndRole, integratedObservableHandling );
}

inline std::shared_ptr< ObservationDependentVariableSettings > azimuthAngleDependentVariable(
    const LinkEndType linkEndRole = unidentified_link_end,
    const LinkEndId linkEndId = LinkEndId( "", "" ),
    const LinkEndType originatingLinkEndRole = unidentified_link_end,
    const LinkEndId originatingLinkEndId = LinkEndId( "", "" ),
    const IntegratedObservationPropertyHandling integratedObservableHandling = interval_start )
{
    return std::make_shared< StationAngleObservationDependentVariableSettings >(
            station_azimuth_angle, linkEndId, linkEndRole, originatingLinkEndId, originatingLinkEndRole, integratedObservableHandling );
}

inline std::shared_ptr< ObservationDependentVariableSettings > targetRangeBetweenLinkEndsDependentVariable(
        const LinkEndType startLinkEndType = unidentified_link_end,
        const LinkEndType endLinkEndType = unidentified_link_end,
        const LinkEndId startLinkEndId = LinkEndId( "", "" ),
        const LinkEndId endLinkEndId = LinkEndId( "", "" ),
        const IntegratedObservationPropertyHandling integratedObservableHandling = interval_start )
{
    return std::make_shared< InterlinkObservationDependentVariableSettings >(
            target_range, startLinkEndType, endLinkEndType, startLinkEndId, endLinkEndId, integratedObservableHandling );
}

inline std::shared_ptr< ObservationDependentVariableSettings > bodyAvoidanceAngleDependentVariable(
        const std::string relativeBody,
        const LinkEndType startLinkEndType = unidentified_link_end,
        const LinkEndType endLinkEndType = unidentified_link_end,
        const LinkEndId startLinkEndId = LinkEndId( "", "" ),
        const LinkEndId endLinkEndId = LinkEndId( "", "" ),
        const IntegratedObservationPropertyHandling integratedObservableHandling = interval_start )
{
    return std::make_shared< InterlinkObservationDependentVariableSettings >(
            body_avoidance_angle_variable, startLinkEndType, endLinkEndType, startLinkEndId, endLinkEndId, integratedObservableHandling, relativeBody );
}

inline std::shared_ptr< ObservationDependentVariableSettings > linkBodyCenterDistanceDependentVariable(
        const std::string relativeBody,
        const LinkEndType startLinkEndType = unidentified_link_end,
        const LinkEndType endLinkEndType = unidentified_link_end,
        const LinkEndId startLinkEndId = LinkEndId( "", "" ),
        const LinkEndId endLinkEndId = LinkEndId( "", "" ),
        const IntegratedObservationPropertyHandling integratedObservableHandling = interval_start )
{
    return std::make_shared< InterlinkObservationDependentVariableSettings >(
            link_body_center_distance, startLinkEndType, endLinkEndType, startLinkEndId, endLinkEndId, integratedObservableHandling, relativeBody );
}

inline std::shared_ptr< ObservationDependentVariableSettings > linkLimbDistanceDependentVariable(
        const std::string relativeBody,
        const LinkEndType startLinkEndType = unidentified_link_end,
        const LinkEndType endLinkEndType = unidentified_link_end,
        const LinkEndId startLinkEndId = LinkEndId( "", "" ),
        const LinkEndId endLinkEndId = LinkEndId( "", "" ),
        const IntegratedObservationPropertyHandling integratedObservableHandling = interval_start )
{
    return std::make_shared< InterlinkObservationDependentVariableSettings >(
            link_limb_distance, startLinkEndType, endLinkEndType, startLinkEndId, endLinkEndId, integratedObservableHandling, relativeBody );
}

inline std::shared_ptr< ObservationDependentVariableSettings > linkAngleWrtOrbitalPlaneDependentVariable(
        const std::string relativeBody,
        const LinkEndType startLinkEndType = unidentified_link_end,
        const LinkEndType endLinkEndType = unidentified_link_end,
        const LinkEndId startLinkEndId = LinkEndId( "", "" ),
        const LinkEndId endLinkEndId = LinkEndId( "", "" ),
        const IntegratedObservationPropertyHandling integratedObservableHandling = interval_start )
{
    return std::make_shared< InterlinkObservationDependentVariableSettings >(
            link_angle_with_orbital_plane, startLinkEndType, endLinkEndType, startLinkEndId, endLinkEndId, integratedObservableHandling, relativeBody );
}

inline std::shared_ptr< ObservationDependentVariableSettings > dopplerIntegrationTimeDependentVariable(
        const ObservableType observableType = undefined_observation_model )
{
    return std::make_shared< AncillaryObservationDependentVariableSettings >( doppler_integration_time_dependent_variable, observableType );
}

inline std::shared_ptr< ObservationDependentVariableSettings > retransmissionDelaysDependentVariable(
        const ObservableType observableType = undefined_observation_model )
{
    return std::make_shared< AncillaryObservationDependentVariableSettings >( retransmission_delays_dependent_variable, observableType );
}

}

}
#endif // TUDAT_OBSERVATIONOUTPUTSETTINGS
