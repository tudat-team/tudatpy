/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */



#include "tudat/astro/orbit_determination/estimatable_parameters/specularReflectivity.h"

namespace tudat
{

namespace estimatable_parameters
{

SpecularReflectivity::SpecularReflectivity(
    const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > vehiclePanels,
    const std::string& associatedBody,
    const std::string panelTypeId):
    EstimatableParameter< double >( specular_reflectivity, associatedBody, panelTypeId ),
    panelTypeId_( panelTypeId )
{
    // check if the panelTypeID exists among the panels
    if(vehiclePanels.size() < 1)
    {
        throw std::runtime_error( "Error when creating estimated specular reflectivity for " +
                                  panelTypeId + " of " + associatedBody + ", no corresponding panels defined" );
    }

    for( unsigned int i = 0; i < vehiclePanels.size( ); i++ )
    {
        bool isNull = true;
        if(vehiclePanels.at ( i )->getReflectionLaw() != nullptr )
        {
            auto reflectionLaw = std::dynamic_pointer_cast<electromagnetism::SpecularDiffuseMixReflectionLaw>(vehiclePanels.at(i)->getReflectionLaw());
            if( reflectionLaw != nullptr )
            {
                reflectionLaws_.push_back( reflectionLaw );
                isNull = false;
            }
        }

        if( isNull )
        {
            throw std::runtime_error( "Error when creating estimated specular reflectivity for " +
                                      panelTypeId + " of " + associatedBody + ", detected incompatible panel reflection law" );
        }
    }
    normalizeValue( );
}

double SpecularReflectivity::normalizeValue( )
{

    // Retrieve all specular reflectivity values for the panels corresponding to the given panelTypeId
    std::vector<double> specularReflectivities = getPanelSpecularReflectivities(  );

    // Calculate the average diffuse reflectivity
    double averageSpecularReflectivity = std::accumulate(specularReflectivities.begin(), specularReflectivities.end(), 0.0) / specularReflectivities.size();

    // Check if all values are the same
    bool allValuesSame = std::all_of(specularReflectivities.begin(), specularReflectivities.end(),
                                     [&](double value) { return value == specularReflectivities[0]; });

    // If not all values are the same, print a warning and reset the values to the average
    if (!allValuesSame)
    {
        std::cerr << "Warning: Specular reflectivity values for panel group "
                  << panelTypeId_ << " are not consistent. Resetting all to the average value." << std::endl;

        // Set all panels' specular reflectivity to the average value
        for( unsigned int i  = 0; i < reflectionLaws_.size( ); i++ )
        {
            reflectionLaws_.at( i )->setSpecularReflectivity( averageSpecularReflectivity, true );
        }
    }

    // Return the average diffuse reflectivity in all cases
    return averageSpecularReflectivity;
}

double SpecularReflectivity::getParameterValue( )
{
    return normalizeValue( );
}

void SpecularReflectivity::setParameterValue( double parameterValue )
{
    for( unsigned int i  = 0; i < reflectionLaws_.size( ); i++ )
    {
        reflectionLaws_.at( i )->setSpecularReflectivity( parameterValue, true );
        double absorptivity = reflectionLaws_.at( i )->getAbsorptivity( );

        // Check if the absorptivity is negative, which indicates non-physical behavior
        if ( absorptivity < 0.0 )
        {
            // Print a warning with both the non-physical specular reflectivity and the resulting negative absorptivity
            std::cerr << "Warning: Non-physical behavior detected for panel group "
                      << panelTypeId_ << ". Specular reflectivity = " << parameterValue
                      << ", resulting absorptivity = " << absorptivity << "."
                      << std::endl;
        }
    }
}

std::vector< double > SpecularReflectivity::getPanelSpecularReflectivities( )
{
    std::vector< double > values;
    for( unsigned int i = 0; i < reflectionLaws_.size( ); i++ )
    {
        values.push_back( reflectionLaws_.at( i )->getSpecularReflectivity( ) );
    }
    return values;
}


} // namespace estimatable_parameters

} // namespace tudat

