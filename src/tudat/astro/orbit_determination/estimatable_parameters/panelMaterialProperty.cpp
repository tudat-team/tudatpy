/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/estimatable_parameters/panelMaterialProperty.h"

#include <algorithm>
#include <iostream>
#include <numeric>
#include <stdexcept>

namespace tudat
{

namespace estimatable_parameters
{

PanelMaterialPropertyParameter::PanelMaterialPropertyParameter(
        const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > vehiclePanels,
        const std::string& associatedBody,
        const std::string& panelTypeId,
        const EstimatebleParametersEnum parameterType ):
    EstimatableParameter< double >( parameterType, associatedBody, panelTypeId ), vehiclePanels_( vehiclePanels ),
    panelTypeId_( panelTypeId )
{
    if( !( ( parameterType == energy_accomodation_coefficient ) || ( parameterType == normal_accomodation_coefficient ) ||
           ( parameterType == tangential_accomodation_coefficient ) || ( parameterType == normal_velocity_at_wall_ratio ) ) )
    {
        throw std::runtime_error( "Error when creating estimated panel material property for " + panelTypeId + " of " + associatedBody +
                                  ", input type is inconsistent" );
    }

    if( vehiclePanels_.size( ) < 1 )
    {
        throw std::runtime_error( "Error when creating estimated panel material property for " + panelTypeId + " of " + associatedBody +
                                  ", no corresponding panels defined" );
    }

    normalizeValue( );
}

double PanelMaterialPropertyParameter::normalizeValue( )
{
    std::vector< double > values = getPanelMaterialPropertyValues( );
    double averageValue = std::accumulate( values.begin( ), values.end( ), 0.0 ) / values.size( );

    bool allValuesSame = std::all_of( values.begin( ), values.end( ), [ & ]( double value ) { return value == values[ 0 ]; } );

    if( !allValuesSame )
    {
        std::cerr << "Warning: material property values for panel group " << panelTypeId_
                  << " are not consistent. Resetting all to the average value." << std::endl;

        for( unsigned int i = 0; i < vehiclePanels_.size( ); i++ )
        {
            setPanelMaterialPropertyValue( vehiclePanels_.at( i ), averageValue );
        }
    }

    return averageValue;
}

double PanelMaterialPropertyParameter::getParameterValue( )
{
    return normalizeValue( );
}

void PanelMaterialPropertyParameter::setParameterValue( double parameterValue )
{
    for( unsigned int i = 0; i < vehiclePanels_.size( ); i++ )
    {
        setPanelMaterialPropertyValue( vehiclePanels_.at( i ), parameterValue );
    }
}

std::vector< double > PanelMaterialPropertyParameter::getPanelMaterialPropertyValues( )
{
    std::vector< double > values;
    for( unsigned int i = 0; i < vehiclePanels_.size( ); i++ )
    {
        values.push_back( getPanelMaterialPropertyValue( vehiclePanels_.at( i ) ) );
    }
    return values;
}

double PanelMaterialPropertyParameter::getPanelMaterialPropertyValue(
        const std::shared_ptr< system_models::VehicleExteriorPanel > vehiclePanel )
{
    switch( parameterName_.first )
    {
        case energy_accomodation_coefficient:
            return vehiclePanel->getEnergyAccomodationCoefficient( );
        case normal_accomodation_coefficient:
            return vehiclePanel->getNormalAccomodationCoefficient( );
        case tangential_accomodation_coefficient:
            return vehiclePanel->getTangentialAccomodationCoefficient( );
        case normal_velocity_at_wall_ratio:
            return vehiclePanel->getNormalVelocityAtWallRatio( );
        default:
            throw std::runtime_error( "Error when retrieving panel material property value, parameter type is inconsistent" );
    }
}

void PanelMaterialPropertyParameter::setPanelMaterialPropertyValue(
        const std::shared_ptr< system_models::VehicleExteriorPanel > vehiclePanel,
        const double parameterValue )
{
    switch( parameterName_.first )
    {
        case energy_accomodation_coefficient:
            vehiclePanel->setEnergyAccomodationCoefficient( parameterValue );
            break;
        case normal_accomodation_coefficient:
            vehiclePanel->setNormalAccomodationCoefficient( parameterValue );
            break;
        case tangential_accomodation_coefficient:
            vehiclePanel->setTangentialAccomodationCoefficient( parameterValue );
            break;
        case normal_velocity_at_wall_ratio:
            vehiclePanel->setNormalVelocityAtWallRatio( parameterValue );
            break;
        default:
            throw std::runtime_error( "Error when setting panel material property value, parameter type is inconsistent" );
    }
}

}  // namespace estimatable_parameters

}  // namespace tudat
