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
    case doppler_integration_time_dependent_variable:
    {
        dependentVariableName = "Doppler integration count time ";
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

std::string getObservationDependentVariableId(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings )
{
    return getObservationDependentVariableName( variableSettings->variableType_ ) +
            variableSettings->getIdentifier( );
}


bool isObservationDependentVariableVectorial(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings )
{
    bool isVariableVectorial = false;
    switch( variableSettings->variableType_ )
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
    case doppler_integration_time_dependent_variable:
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
    default:
        throw std::runtime_error( "Error when checking observation dependent variable. Type " +
                                  getObservationDependentVariableId( variableSettings ) +
                                  " not found when checking if variable is vectorial." );

    }
    return isVariableVectorial;
}


bool isObservationDependentVariableAncilliarySetting(
    const std::shared_ptr< ObservationDependentVariableSettings > variableSettings )
{
    bool isAncilliarySetting = false;
    switch( variableSettings->variableType_ )
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
    case link_angle_with_orbital_plane:
        break;
    case doppler_integration_time_dependent_variable:
        isAncilliarySetting = true;
        break;
    case retransmission_delays_dependent_variable:
        isAncilliarySetting = true;
        break;
    default:
        throw std::runtime_error( "Error when checking observation dependent variable. Type " +
                                  getObservationDependentVariableId( variableSettings ) +
                                  " not found when checking for ancilliary setting." );

    }
    return isAncilliarySetting;
}


bool isObservationDependentVariableGroundStationProperty(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings )
{
    bool isGroundStationProperty = false;
    switch( variableSettings->variableType_ )
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
    case doppler_integration_time_dependent_variable:
        break;
    case link_body_center_distance:
        break;
    case link_limb_distance:
        break;
    case link_angle_with_orbital_plane:
        break;
    default:
        throw std::runtime_error( "Error when checking observation dependent variable. Type " +
                                  getObservationDependentVariableId( variableSettings ) +
                                  " not found when checking for ground station dependency." );

    }
    return isGroundStationProperty;
}

int getObservationDependentVariableSize(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings )
{
    int variableSize = 0;
    if( !isObservationDependentVariableVectorial( variableSettings ) )
    {
        variableSize = 1;
    }
    else
    {
        switch( variableSettings->variableType_ )
        {
        default:
            throw std::runtime_error( "Error when checking observation dependent variable. Type " +
                                      getObservationDependentVariableId( variableSettings ) +
                                      " not found when determining parameter size." );

        }
    }
    return variableSize;
}

bool doesStationAngleVariableExistForGivenLink(
        const observation_models::ObservableType observableType,
        const observation_models::LinkEnds& linkEnds,
        const std::shared_ptr< StationAngleObservationDependentVariableSettings > variableSettings )
{
    bool doesLinkHaveDependency = false;

    if( linkEnds.size( ) > 1 )
    {
        std::vector< observation_models::LinkEndType > linkEndTypeList = getLinkEndTypesForGivenLinkEndId(
                linkEnds, variableSettings->relevantLinkEnd_ );
        if( linkEndTypeList.size( ) > 0 )
        {
            doesLinkHaveDependency = true;
        }
    }
    return doesLinkHaveDependency;
}

bool doesInterlinkVariableExistForGivenLink(
    const observation_models::ObservableType observableType,
    const observation_models::LinkEnds& linkEnds,
    const std::shared_ptr< InterlinkObservationDependentVariableSettings > variableSettings )
{
    bool doesLinkHaveDependency = true;
    if( variableSettings->startLinkEnd_ != observation_models::unidentified_link_end )
    {
        if( linkEnds.count( variableSettings->startLinkEnd_ ) == 0 )
        {
            doesLinkHaveDependency = false;
        }
    }

    if( variableSettings->endLinkEnd_ != observation_models::unidentified_link_end )
    {
        if( linkEnds.count( variableSettings->endLinkEnd_ ) == 0 )
        {
            doesLinkHaveDependency = false;
        }
    }

    return doesLinkHaveDependency;
}

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
    default:
        throw std::runtime_error( "Error when checking observation dependent variable. Type " +
                                  getObservationDependentVariableId( variableSettings ) +
                                  " not found when checking if variable exists for given link." );

    }
    return doesLinkHaveDependency;
}

}

}
