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

#include <cereal/cereal.hpp>
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>

#include <memory>
#include <functional>
#include <string>
#include <vector>

#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/corrections/lightTimeCorrection.h"
#include "tudat/io/serialization/base.h"

namespace tudat
{

namespace simulation_setup
{

using namespace observation_models;

enum ObservationDependentVariables {
    station_elevation_angle,
    station_azimuth_angle,
    target_range,
    body_avoidance_angle_variable,
    link_body_center_distance,
    link_limb_distance,
    link_angle_with_orbital_plane,
    integration_time_dependent_variable,
    retransmission_delays_dependent_variable,
    link_end_epochs_dependent_variable,
    light_time_correction_components
};

//! Function checking whether the interlinks between two link ends are compatible (i.e., for both the originating and receiving ends of the interlink,
//! the link end's types and IDs should be either identical, or undefined).
bool areInterlinksCompatible( const std::pair< LinkEndType, LinkEndId >& firstReceivingLinkEnd,
                              const std::pair< LinkEndType, LinkEndId >& firstOriginatingLinkEnd,
                              const std::pair< LinkEndType, LinkEndId >& secondReceivingLinkEnd,
                              const std::pair< LinkEndType, LinkEndId >& secondOriginatingLinkEnd );

//! Function returning the dependent variable name
std::string getObservationDependentVariableName( const ObservationDependentVariables variableType );

//! Function checking whether a given dependent variable type is related to some ancillary settings
bool isObservationDependentVariableAncillarySetting( const ObservationDependentVariables variableType );

//! Function checking whether a given dependent variable type is an interlink property
bool isObservationDependentVariableInterlinkProperty( const ObservationDependentVariables variableType );

//! Function checking whether a given interlink dependent variable depends on the link direction (i.e., depends on which link end is the receiving/originating one)
bool isInterlinkPropertyDirectionAgnostic( const ObservationDependentVariables variableType );

//! Function checking whether a given dependent variable should be given as a vector of size > 1
bool isObservationDependentVariableVectorial( const ObservationDependentVariables variableType );

//! Function checking whether a given dependent variable is related to a ground station property
bool isObservationDependentVariableGroundStationProperty( const ObservationDependentVariables variableType );

//! Function checking whether a given dependent variable is link end-dependent (false for ancillary settings dependent variables)
bool isObservationDependentVariableLinkEndDependent( const ObservationDependentVariables variableType );

//! Base class for observation dependent variable settings
class ObservationDependentVariableSettings
{
public:
    ObservationDependentVariableSettings( const ObservationDependentVariables variableType,
                                          const LinkEndId linkEndId = LinkEndId( "", "" ),
                                          const LinkEndType linkEndType = unidentified_link_end,
                                          const LinkEndId originatingLinkEndId = LinkEndId( "", "" ),
                                          const LinkEndType originatingLinkEndType = unidentified_link_end ):
        variableType_( variableType ), linkEndId_( linkEndId ), linkEndType_( linkEndType ), originatingLinkEndId_( originatingLinkEndId ),
        originatingLinkEndType_( originatingLinkEndType )
    {}

    virtual ~ObservationDependentVariableSettings( ) {}

    ObservationDependentVariables variableType_;

    // Used for serialization testing
    bool operator==( const ObservationDependentVariableSettings& rhs ) const
    {
        return equals( rhs );
    }

    bool operator!=( const ObservationDependentVariableSettings& rhs ) const
    {
        return !equals( rhs );
    }

    //! Get identifier for base dependent variable settings
    std::string getBaseIdentifier( )
    {
        std::string identifier = ", link end";
        if( linkEndId_ != LinkEndId( "", "" ) && linkEndType_ != unidentified_link_end )
        {
            identifier += ": (" + linkEndId_.bodyName_ + ", " + linkEndId_.getReferencePointName( ) + ") as " +
                    getLinkEndTypeString( linkEndType_ );
        }
        else if( linkEndId_ != LinkEndId( "", "" ) )
        {
            identifier = ": (" + linkEndId_.bodyName_ + ", " + linkEndId_.getReferencePointName( ) + ")";
        }
        else if( linkEndType_ != unidentified_link_end )
        {
            identifier = " of type " + getLinkEndTypeString( linkEndType_ );
        }

        if( originatingLinkEndId_ != LinkEndId( "", "" ) && originatingLinkEndType_ != unidentified_link_end )
        {
            identifier += " to link end: (" + originatingLinkEndId_.bodyName_ + ", " + originatingLinkEndId_.getReferencePointName( ) +
                    ") as " + getLinkEndTypeString( originatingLinkEndType_ );
        }
        else if( originatingLinkEndId_ != LinkEndId( "", "" ) )
        {
            identifier += " to link end: (" + originatingLinkEndId_.bodyName_ + ", " + originatingLinkEndId_.getReferencePointName( ) + ")";
        }
        else if( originatingLinkEndType_ != unidentified_link_end )
        {
            identifier += " to link end of type " + getLinkEndTypeString( originatingLinkEndType_ );
        }
        return identifier;
    }

    virtual std::string getIdentifier( )
    {
        return getBaseIdentifier( );
    }

    //! Function that checks whether two dependent variable settings are compatible (i.e., they might refer to the same dependent variable).
    virtual bool areSettingsCompatible( const std::shared_ptr< ObservationDependentVariableSettings > otherSettings )
    {
        return areBaseSettingsCompatible( otherSettings, true );
    }

    //! Function that checks whether the base elements of two dependent variable settings objects are compatible. This implies that the dependent variable type should be identical,
    //! and that the originating and receiving parts of the link for which the dependent variable is computed should be either identical or not specified. The function checks whether
    //! the dependent variable under consideration depends on the link "direction" (i.e., which link end is the originating/receiving one) or is agnostic to it. In the latter case,
    //! dependent variable settings with "reverted" links (i.e., originating/receiving ends are inverted) are considered compatible.
    bool areBaseSettingsCompatible( const std::shared_ptr< ObservationDependentVariableSettings > otherSettings,
                                    const bool revertedLinksAllowed = true )
    {
        bool isCompatible = true;
        // Check if the dependent variables are identical
        if( variableType_ != otherSettings->variableType_ )
        {
            isCompatible = false;
        }
        else
        {
            std::pair< LinkEndType, LinkEndId > receivingLinkEnd = std::make_pair( linkEndType_, linkEndId_ );
            std::pair< LinkEndType, LinkEndId > originatingLinkEnd = std::make_pair( originatingLinkEndType_, originatingLinkEndId_ );
            std::pair< LinkEndType, LinkEndId > otherSettingsReceivingLinkEnd =
                    std::make_pair( otherSettings->linkEndType_, otherSettings->linkEndId_ );
            std::pair< LinkEndType, LinkEndId > otherSettingsOriginatingLinkEnd =
                    std::make_pair( otherSettings->originatingLinkEndType_, otherSettings->originatingLinkEndId_ );

            // Check if the link ends (both receiving and originating ends) are either identical or undefined (in which case the two
            // settings would be considered compatible).
            bool directLinksMatch = areInterlinksCompatible(
                    receivingLinkEnd, originatingLinkEnd, otherSettingsReceivingLinkEnd, otherSettingsOriginatingLinkEnd );

            // Check if inverting the receiving/originating ends of the link would lead to compatible link definitions (for dependent
            // variables that are independent of the link "direction").
            bool revertedLinksMatch = areInterlinksCompatible(
                    receivingLinkEnd, originatingLinkEnd, otherSettingsOriginatingLinkEnd, otherSettingsReceivingLinkEnd );

            // Check if the links are either directly compatible, or compatible once reverted if allowed for the specific dependent variable
            // under consideration.
            if( !directLinksMatch && ( !revertedLinksAllowed || !revertedLinksMatch ) )
            {
                isCompatible = false;
            }
        }

        return isCompatible;
    }

    //! Link end ID (receiving end of the link)
    LinkEndId linkEndId_;

    //! Link end type (receiving end of the link)
    LinkEndType linkEndType_;

    //! Link end ID (originating end of the link)
    LinkEndId originatingLinkEndId_;

    //! Link end type (originating end of the link)
    LinkEndType originatingLinkEndType_;

    //! Save dependent variable settings to a JSON file
    TUDAT_DEFINE_FILE_IO_POLYMORPHIC( ObservationDependentVariableSettings )

protected:
    // Default constructor for serialization
    ObservationDependentVariableSettings( ): variableType_( station_elevation_angle ) {}

    // Each derived class should implement this function such that it returns true if a deserialized object is
    // equal to the original object.
    virtual bool equals( const ObservationDependentVariableSettings& rhs ) const
    {
        return variableType_ == rhs.variableType_ && linkEndId_ == rhs.linkEndId_ && linkEndType_ == rhs.linkEndType_ &&
                originatingLinkEndId_ == rhs.originatingLinkEndId_ && originatingLinkEndType_ == rhs.originatingLinkEndType_;
    }

private:
    friend class cereal::access;

    template< class Archive >
    void save( Archive& ar ) const
    {
        ar( CEREAL_NVP( variableType_ ) );
        ar( CEREAL_NVP( linkEndId_ ) );
        ar( CEREAL_NVP( linkEndType_ ) );
        ar( CEREAL_NVP( originatingLinkEndId_ ) );
        ar( CEREAL_NVP( originatingLinkEndType_ ) );
    }

    template< class Archive >
    void load( Archive& ar )
    {
        ar( CEREAL_NVP( variableType_ ) );
        ar( CEREAL_NVP( linkEndId_ ) );
        ar( CEREAL_NVP( linkEndType_ ) );
        ar( CEREAL_NVP( originatingLinkEndId_ ) );
        ar( CEREAL_NVP( originatingLinkEndType_ ) );
    }
};

enum IntegratedObservationPropertyHandling { interval_start, interval_end, interval_undefined };

std::string getIntegrationHandlingString( const IntegratedObservationPropertyHandling integratedObservableHandling );

class StationAngleObservationDependentVariableSettings : public ObservationDependentVariableSettings
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
        isLinkEndDefined_( ( relevantLinkEnd != LinkEndId( "", "" ) ? true : false ) )
    {}

    std::string getIdentifier( )
    {
        return getBaseIdentifier( ) + getIntegrationHandlingString( integratedObservableHandling_ );
    }

    //! Function that checks whether two dependent variable settings are compatible (i.e., they might refer to the same dependent variable).
    bool areSettingsCompatible( const std::shared_ptr< ObservationDependentVariableSettings > otherSettings )
    {
        bool isCompatible = true;
        std::shared_ptr< StationAngleObservationDependentVariableSettings > stationAngleSettings =
                std::dynamic_pointer_cast< StationAngleObservationDependentVariableSettings >( otherSettings );

        //! Check that both settings are station angle properties
        if( stationAngleSettings == nullptr )
        {
            isCompatible = false;
        }
        else
        {
            // Check whether base settings  are compatible (i.e., same variable type and compatible originating/receiving link ends)
            if( !areBaseSettingsCompatible( otherSettings, false ) )
            {
                isCompatible = false;
            }
            else
            {
                // Check that the time to consider for an integrated observable, if defined, is consistent
                if( ( integratedObservableHandling_ != stationAngleSettings->integratedObservableHandling_ ) &&
                    ( integratedObservableHandling_ != interval_undefined ) &&
                    ( stationAngleSettings->integratedObservableHandling_ != interval_undefined ) )
                {
                    isCompatible = false;
                }
            }
        }

        return isCompatible;
    }

    IntegratedObservationPropertyHandling integratedObservableHandling_;

    bool isLinkEndDefined_;

protected:
    // Default constructor for serialization
    StationAngleObservationDependentVariableSettings( ): integratedObservableHandling_( interval_undefined ), isLinkEndDefined_( false ) {}

    // Used for serialization testing
    bool equals( const ObservationDependentVariableSettings& other ) const override
    {
        const auto* rhs = dynamic_cast< const StationAngleObservationDependentVariableSettings* >( &other );
        if( !rhs )
        {
            return false;
        }
        return ObservationDependentVariableSettings::equals( other ) &&
                integratedObservableHandling_ == rhs->integratedObservableHandling_ && isLinkEndDefined_ == rhs->isLinkEndDefined_;
    }

private:
    friend class cereal::access;

    template< class Archive >
    void save( Archive& ar ) const
    {
        ar( cereal::base_class< ObservationDependentVariableSettings >( this ) );
        ar( CEREAL_NVP( integratedObservableHandling_ ) );
        ar( CEREAL_NVP( isLinkEndDefined_ ) );
    }

    template< class Archive >
    void load( Archive& ar )
    {
        ar( cereal::base_class< ObservationDependentVariableSettings >( this ) );
        ar( CEREAL_NVP( integratedObservableHandling_ ) );
        ar( CEREAL_NVP( isLinkEndDefined_ ) );
    }
};

class InterlinkObservationDependentVariableSettings : public ObservationDependentVariableSettings
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
        integratedObservableHandling_( integratedObservableHandling ), relativeBody_( relativeBody )
    {}

    ~InterlinkObservationDependentVariableSettings( ) {}

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

    //! Function that checks whether two dependent variable settings are compatible (i.e., they might refer to the same dependent variable).
    bool areSettingsCompatible( const std::shared_ptr< ObservationDependentVariableSettings > otherSettings )
    {
        bool isCompatible = true;
        std::shared_ptr< InterlinkObservationDependentVariableSettings > interlinkSettings =
                std::dynamic_pointer_cast< InterlinkObservationDependentVariableSettings >( otherSettings );

        //! Check that both settings are interlink properties
        if( interlinkSettings == nullptr )
        {
            isCompatible = false;
        }
        else
        {
            // Check whether base settings  are compatible (i.e., same variable type and compatible originating/receiving link ends)
            if( !areBaseSettingsCompatible( otherSettings, isInterlinkPropertyDirectionAgnostic( variableType_ ) ) )
            {
                isCompatible = false;
            }
            else
            {
                // Check that the time to consider for an integrated observable, if defined, is consistent
                if( ( integratedObservableHandling_ != interlinkSettings->integratedObservableHandling_ ) &&
                    ( integratedObservableHandling_ != interval_undefined ) &&
                    ( interlinkSettings->integratedObservableHandling_ != interval_undefined ) )
                {
                    isCompatible = false;
                }
                //! Check that relative body, if defined, is identical.
                if( ( relativeBody_ != interlinkSettings->relativeBody_ ) && ( relativeBody_ != "" ) &&
                    ( interlinkSettings->relativeBody_ != "" ) )
                {
                    isCompatible = false;
                }
            }
        }

        return isCompatible;
    }

    IntegratedObservationPropertyHandling integratedObservableHandling_;

    std::string relativeBody_;

protected:
    // Default constructor for serialization
    InterlinkObservationDependentVariableSettings( ): integratedObservableHandling_( interval_undefined ) {}

    // Used for serialization testing
    bool equals( const ObservationDependentVariableSettings& other ) const override
    {
        const auto* rhs = dynamic_cast< const InterlinkObservationDependentVariableSettings* >( &other );
        if( !rhs )
        {
            return false;
        }
        return ObservationDependentVariableSettings::equals( other ) &&
                integratedObservableHandling_ == rhs->integratedObservableHandling_ && relativeBody_ == rhs->relativeBody_;
    }

private:
    friend class cereal::access;

    template< class Archive >
    void save( Archive& ar ) const
    {
        ar( cereal::base_class< ObservationDependentVariableSettings >( this ) );
        ar( CEREAL_NVP( integratedObservableHandling_ ) );
        ar( CEREAL_NVP( relativeBody_ ) );
    }

    template< class Archive >
    void load( Archive& ar )
    {
        ar( cereal::base_class< ObservationDependentVariableSettings >( this ) );
        ar( CEREAL_NVP( integratedObservableHandling_ ) );
        ar( CEREAL_NVP( relativeBody_ ) );
    }
};

//! Returns a function which checks whether an ancillary settings dependent variable exists for a given observable type (dependent on whether
//! such ancillary settings are defined for a specific observable type)
std::function< bool( const ObservableType observableType ) > getIsObservableTypeCompatibleFunction(
        const ObservationDependentVariables variableType );

class AncillaryObservationDependentVariableSettings : public ObservationDependentVariableSettings
{
public:
    AncillaryObservationDependentVariableSettings( const ObservationDependentVariables variableType,
                                                   const ObservableType observableType = undefined_observation_model ):
        ObservationDependentVariableSettings( variableType ), observableType_( observableType )
    {
        // Check whether the ancillary settings exists for the required observable type
        isObservableTypeCompatible_ = getIsObservableTypeCompatibleFunction( variableType );
    }

    ~AncillaryObservationDependentVariableSettings( ) {}

    std::string getIdentifier( )
    {
        std::string identifier = getBaseIdentifier( );
        if( observableType_ != undefined_observation_model )
        {
            identifier += ", for observable of type " + std::to_string( observableType_ );
        }
        return identifier;
    }

    //! Function that checks whether two dependent variable settings are compatible (i.e., they might refer to the same dependent variable).
    bool areSettingsCompatible( const std::shared_ptr< ObservationDependentVariableSettings > otherSettings )
    {
        bool isCompatible = true;
        std::shared_ptr< AncillaryObservationDependentVariableSettings > ancillarySettings =
                std::dynamic_pointer_cast< AncillaryObservationDependentVariableSettings >( otherSettings );

        //! Check that both settings refer to observation ancillary settings
        if( ancillarySettings == nullptr )
        {
            isCompatible = false;
        }
        else
        {
            // Check that the variable type is consistent
            if( variableType_ != otherSettings->variableType_ )
            {
                isCompatible = false;
            }
            // Check that the observable type for which the ancillary data should be retrieved, if defined, is consistent
            if( ( observableType_ != ancillarySettings->observableType_ ) && ( observableType_ != undefined_observation_model ) &&
                ( ancillarySettings->observableType_ != undefined_observation_model ) )
            {
                isCompatible = false;
            }
        }

        return isCompatible;
    }

    ObservableType observableType_;

    std::function< bool( const ObservableType observableType ) > isObservableTypeCompatible_;

protected:
    // Default constructor for serialization
    AncillaryObservationDependentVariableSettings( ): observableType_( undefined_observation_model ) {}

    // Used for serialization testing
    bool equals( const ObservationDependentVariableSettings& other ) const override
    {
        const auto* rhs = dynamic_cast< const AncillaryObservationDependentVariableSettings* >( &other );
        if( !rhs )
        {
            return false;
        }
        // isObservableTypeCompatible_ is a std::function (not comparable) — compare only scalar members
        return ObservationDependentVariableSettings::equals( other ) && observableType_ == rhs->observableType_;
    }

private:
    friend class cereal::access;

    template< class Archive >
    void save( Archive& ar ) const
    {
        ar( cereal::base_class< ObservationDependentVariableSettings >( this ) );
        ar( CEREAL_NVP( observableType_ ) );
    }

    template< class Archive >
    void load( Archive& ar )
    {
        ar( cereal::base_class< ObservationDependentVariableSettings >( this ) );
        ar( CEREAL_NVP( observableType_ ) );
        // Reconstruct the function after loading
        isObservableTypeCompatible_ = getIsObservableTypeCompatibleFunction( variableType_ );
    }
};

//! Function that returns a string uniquely describing a dependent variable settings object
std::string getObservationDependentVariableId( const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );

//! Function that returns the size of a given dependent variable (can be link end-dependent for ancillary settings dependent variables)
int getObservationDependentVariableSize( const std::shared_ptr< ObservationDependentVariableSettings > variableSettings,
                                         const LinkEnds linkEnds );

//! Function that checks whether a given station angle dependent variable can be computed for a given observable type and link ends.
bool doesStationAngleVariableExistForGivenLink(
        const ObservableType observableType,
        const LinkEnds& linkEnds,
        const std::shared_ptr< StationAngleObservationDependentVariableSettings > variableSettings );

//! Function that checks whether a given interlink dependent variable can be computed for a given observable type and link ends.
bool doesInterlinkVariableExistForGivenLink( const ObservableType observableType,
                                             const LinkEnds& linkEnds,
                                             const std::shared_ptr< InterlinkObservationDependentVariableSettings > variableSettings );

//! Function that checks whether a given dependent variable can be computed for a given observable type and link ends.
bool doesObservationDependentVariableExistForGivenLink( const ObservableType observableType,
                                                        const LinkEnds& linkEnds,
                                                        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings );

//! Function that returns a fully defined ObservationDependentVariableSettings object, using the original settings (possibly incomplete, i.e.
//! with missing link ends information) and specific link ends information for both the receiving and originating ends of the link
std::shared_ptr< ObservationDependentVariableSettings > createCompleteObservationDependentVariableSettings(
        const std::shared_ptr< ObservationDependentVariableSettings > originalSettings,
        const LinkEndType& linkEndType,
        const LinkEndId& linkEndId,
        const LinkEndType& originatingLinkEndType,
        const LinkEndId& originatingLinkEndId );

//! Function that returns a list of all compatible dependent variable settings that can be created for a given observable type and link ends, from a base
//! dependent variable settings that might not be entirely defined (i.e., some link ends not specified, etc.)
std::vector< std::shared_ptr< ObservationDependentVariableSettings > > createAllCompatibleDependentVariableSettings(
        const ObservableType observableType,
        const LinkEnds& linkEnds,
        std::shared_ptr< ObservationDependentVariableSettings > dependentVariableSettings );

//! Function to create a dependent variable computing the elevation angle from a given station
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

//! Function to create a dependent variable computing the azimuth angle from a given station
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

//! Function to create a dependent variable computing the range between the two ends of the link
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

//! Function to create a dependent variable computing the link avoidance angle w.r.t. a given body
inline std::shared_ptr< ObservationDependentVariableSettings > bodyAvoidanceAngleDependentVariable(
        const std::string relativeBody,
        const LinkEndType startLinkEndType = unidentified_link_end,
        const LinkEndType endLinkEndType = unidentified_link_end,
        const LinkEndId startLinkEndId = LinkEndId( "", "" ),
        const LinkEndId endLinkEndId = LinkEndId( "", "" ),
        const IntegratedObservationPropertyHandling integratedObservableHandling = interval_start )
{
    return std::make_shared< InterlinkObservationDependentVariableSettings >( body_avoidance_angle_variable,
                                                                              startLinkEndType,
                                                                              endLinkEndType,
                                                                              startLinkEndId,
                                                                              endLinkEndId,
                                                                              integratedObservableHandling,
                                                                              relativeBody );
}

//! Function to create a dependent variable computing the minimum distance between the center of a body and the link
inline std::shared_ptr< ObservationDependentVariableSettings > linkBodyCenterDistanceDependentVariable(
        const std::string relativeBody,
        const LinkEndType startLinkEndType = unidentified_link_end,
        const LinkEndType endLinkEndType = unidentified_link_end,
        const LinkEndId startLinkEndId = LinkEndId( "", "" ),
        const LinkEndId endLinkEndId = LinkEndId( "", "" ),
        const IntegratedObservationPropertyHandling integratedObservableHandling = interval_start )
{
    return std::make_shared< InterlinkObservationDependentVariableSettings >( link_body_center_distance,
                                                                              startLinkEndType,
                                                                              endLinkEndType,
                                                                              startLinkEndId,
                                                                              endLinkEndId,
                                                                              integratedObservableHandling,
                                                                              relativeBody );
}

//! Function to create a dependent variable computing the minimum distance between the limb of a body and the link
inline std::shared_ptr< ObservationDependentVariableSettings > linkLimbDistanceDependentVariable(
        const std::string relativeBody,
        const LinkEndType startLinkEndType = unidentified_link_end,
        const LinkEndType endLinkEndType = unidentified_link_end,
        const LinkEndId startLinkEndId = LinkEndId( "", "" ),
        const LinkEndId endLinkEndId = LinkEndId( "", "" ),
        const IntegratedObservationPropertyHandling integratedObservableHandling = interval_start )
{
    return std::make_shared< InterlinkObservationDependentVariableSettings >( link_limb_distance,
                                                                              startLinkEndType,
                                                                              endLinkEndType,
                                                                              startLinkEndId,
                                                                              endLinkEndId,
                                                                              integratedObservableHandling,
                                                                              relativeBody );
}

//! Function to create a dependent variable computing the link - orbital plane angle
inline std::shared_ptr< ObservationDependentVariableSettings > linkAngleWrtOrbitalPlaneDependentVariable(
        const std::string relativeBody,
        const LinkEndType startLinkEndType = unidentified_link_end,
        const LinkEndType endLinkEndType = unidentified_link_end,
        const LinkEndId startLinkEndId = LinkEndId( "", "" ),
        const LinkEndId endLinkEndId = LinkEndId( "", "" ),
        const IntegratedObservationPropertyHandling integratedObservableHandling = interval_start )
{
    return std::make_shared< InterlinkObservationDependentVariableSettings >( link_angle_with_orbital_plane,
                                                                              startLinkEndType,
                                                                              endLinkEndType,
                                                                              startLinkEndId,
                                                                              endLinkEndId,
                                                                              integratedObservableHandling,
                                                                              relativeBody );
}

//! Function to create integration time dpeendent variable
inline std::shared_ptr< ObservationDependentVariableSettings > integrationTimeDependentVariable(
        const ObservableType observableType = undefined_observation_model )
{
    return std::make_shared< AncillaryObservationDependentVariableSettings >( integration_time_dependent_variable, observableType );
}

//! Function to create retransmission delays dependent variable
inline std::shared_ptr< ObservationDependentVariableSettings > retransmissionDelaysDependentVariable(
        const ObservableType observableType = undefined_observation_model )
{
    return std::make_shared< AncillaryObservationDependentVariableSettings >( retransmission_delays_dependent_variable, observableType );
}

//! Function to create link-end epochs dependent variable
inline std::shared_ptr< ObservationDependentVariableSettings > linkEndEpochsDependentVariable(
        const ObservableType observableType = undefined_observation_model )
{
    return std::make_shared< AncillaryObservationDependentVariableSettings >( link_end_epochs_dependent_variable, observableType );
}

//! Settings for saving the per-type light-time correction contributions on a single leg
//! (transmitter → receiver) of an observable. The receiving link end is stored in the base
//! `linkEndId_` / `linkEndType_` members; the originating (transmitting) end is stored in the
//! base `originatingLinkEndId_` / `originatingLinkEndType_` members. An optional filter selects
//! a subset of correction types by `LightTimeCorrectionType` (empty = keep all registered
//! corrections). Useful e.g. for separating code vs. carrier corrections, which can have opposite
//! signs (ionospheric correction flips sign between group delay and phase advance).
class LightTimeCorrectionComponentsDependentVariableSettings : public ObservationDependentVariableSettings
{
public:
    LightTimeCorrectionComponentsDependentVariableSettings(
            const LinkEndType transmitterLinkEndType,
            const LinkEndType receiverLinkEndType,
            const LinkEndId transmitterLinkEndId,
            const LinkEndId receiverLinkEndId,
            const std::vector< observation_models::LightTimeCorrectionType > correctionTypeFilter =
                    std::vector< observation_models::LightTimeCorrectionType >( ) ):
        ObservationDependentVariableSettings( light_time_correction_components,
                                              receiverLinkEndId,
                                              receiverLinkEndType,
                                              transmitterLinkEndId,
                                              transmitterLinkEndType ),
        correctionTypeFilter_( correctionTypeFilter ), resolvedSize_( -1 )
    {}

    ~LightTimeCorrectionComponentsDependentVariableSettings( ) {}

    bool areSettingsCompatible( const std::shared_ptr< ObservationDependentVariableSettings > otherSettings ) override
    {
        std::shared_ptr< LightTimeCorrectionComponentsDependentVariableSettings > castSettings =
                std::dynamic_pointer_cast< LightTimeCorrectionComponentsDependentVariableSettings >( otherSettings );
        if( castSettings == nullptr )
        {
            return false;
        }
        if( !areBaseSettingsCompatible( otherSettings, false ) )
        {
            return false;
        }
        return correctionTypeFilter_ == castSettings->correctionTypeFilter_;
    }

    std::vector< observation_models::LightTimeCorrectionType > correctionTypeFilter_;

    //! Size resolved from the actual LightTimeCalculator on the selected leg. A negative value
    //! means the settings have not yet been registered with an observation model.
    int resolvedSize_;

protected:
    // Default constructor for serialization
    LightTimeCorrectionComponentsDependentVariableSettings( ): resolvedSize_( -1 ) {}

    bool equals( const ObservationDependentVariableSettings& rhs ) const override
    {
        if( !ObservationDependentVariableSettings::equals( rhs ) )
        {
            return false;
        }
        const auto* derived = dynamic_cast< const LightTimeCorrectionComponentsDependentVariableSettings* >( &rhs );
        if( !derived )
        {
            return false;
        }
        return correctionTypeFilter_ == derived->correctionTypeFilter_ && resolvedSize_ == derived->resolvedSize_;
    }

private:
    friend class cereal::access;

    template< class Archive >
    void save( Archive& ar ) const
    {
        ar( cereal::base_class< ObservationDependentVariableSettings >( this ) );
        ar( CEREAL_NVP( correctionTypeFilter_ ) );
        ar( CEREAL_NVP( resolvedSize_ ) );
    }

    template< class Archive >
    void load( Archive& ar )
    {
        ar( cereal::base_class< ObservationDependentVariableSettings >( this ) );
        ar( CEREAL_NVP( correctionTypeFilter_ ) );
        ar( CEREAL_NVP( resolvedSize_ ) );
    }
};

//! Function to create a dependent variable saving the individual light-time correction
//! contributions on a single observation leg (transmitter → receiver). The output is a vector
//! with one entry per registered `LightTimeCorrection` on that leg, in registration order. If
//! `correctionTypeFilter` is supplied, all registered corrections whose type appears in the
//! filter are returned. The output is grouped by the first occurrence of each distinct filter
//! type, while preserving registration order within each type group, so multiple registered
//! corrections of the same type produce multiple output entries. Repeated types in
//! `correctionTypeFilter` do not create duplicate groups or additional entries.
inline std::shared_ptr< ObservationDependentVariableSettings > lightTimeCorrectionComponentsDependentVariable(
        const LinkEndType transmitterLinkEndType = unidentified_link_end,
        const LinkEndType receiverLinkEndType = unidentified_link_end,
        const LinkEndId transmitterLinkEndId = LinkEndId( "", "" ),
        const LinkEndId receiverLinkEndId = LinkEndId( "", "" ),
        const std::vector< observation_models::LightTimeCorrectionType > correctionTypeFilter =
                std::vector< observation_models::LightTimeCorrectionType >( ) )
{
    return std::make_shared< LightTimeCorrectionComponentsDependentVariableSettings >(
            transmitterLinkEndType, receiverLinkEndType, transmitterLinkEndId, receiverLinkEndId, correctionTypeFilter );
}

}  // namespace simulation_setup

}  // namespace tudat

// Register derived classes for polymorphic serialization
CEREAL_REGISTER_TYPE( tudat::simulation_setup::ObservationDependentVariableSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::StationAngleObservationDependentVariableSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::InterlinkObservationDependentVariableSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::AncillaryObservationDependentVariableSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::LightTimeCorrectionComponentsDependentVariableSettings )

CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::ObservationDependentVariableSettings,
                                      tudat::simulation_setup::StationAngleObservationDependentVariableSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::ObservationDependentVariableSettings,
                                      tudat::simulation_setup::InterlinkObservationDependentVariableSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::ObservationDependentVariableSettings,
                                      tudat::simulation_setup::AncillaryObservationDependentVariableSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::ObservationDependentVariableSettings,
                                      tudat::simulation_setup::LightTimeCorrectionComponentsDependentVariableSettings )

#endif  // TUDAT_OBSERVATIONOUTPUTSETTINGS
