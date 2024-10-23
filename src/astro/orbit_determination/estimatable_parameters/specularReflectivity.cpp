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


double SpecularReflectivity::normalizeValue( )
{

    // Retrieve all specular reflectivity values for the panels corresponding to the given panelTypeId
    std::vector<double> specularReflectivities = radiationPressureInterface_->getSpecularReflectivityForPanelTypeId(panelTypeId_);

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
        radiationPressureInterface_->setGroupSpecularReflectivity(panelTypeId_, averageSpecularReflectivity);

        // Adjust the absorptivity
        double averageDiffuseReflectivity = radiationPressureInterface_->getAverageDiffuseReflectivity(panelTypeId_);
        double apsorptivity = 1 - averageDiffuseReflectivity - averageSpecularReflectivity;
        radiationPressureInterface_->setGroupAbsorptivity(panelTypeId_, apsorptivity);
    }

    std::cout<<"Getting "<<averageSpecularReflectivity<<std::endl;

    // Return the average diffuse reflectivity in all cases
    return averageSpecularReflectivity;
}

double SpecularReflectivity::getParameterValue( )
{
    return normalizeValue( );
}

void SpecularReflectivity::setParameterValue( double parameterValue )
{
    std::cout<<"Setting "<<parameterValue<<std::endl;


    // Set the specular reflectivity, even if it's greater than 1 (non-physical case)
    radiationPressureInterface_->setGroupSpecularReflectivity(panelTypeId_, parameterValue);
    std::vector<double> specularReflectivities = radiationPressureInterface_->getSpecularReflectivityForPanelTypeId(panelTypeId_);
    for( int i = 0; i < specularReflectivities.size( ); i++ )
    {
        std::cout<<"Val "<<specularReflectivities.at( i )<<std::endl;
    }
    // Get the current diffuse reflectivity
    double diffuseReflectivity = radiationPressureInterface_->getAverageDiffuseReflectivity( panelTypeId_ );

    // Compute the absorptivity such that the sum of specular + diffuse + absorptivity = 1
    double absorptivity = 1.0 - parameterValue - diffuseReflectivity;

    // Check if the absorptivity is negative, which indicates non-physical behavior
    if ( absorptivity < 0.0 )
    {
        // Print a warning with both the non-physical specular reflectivity and the resulting negative absorptivity
        std::cerr << "Warning: Non-physical behavior detected for panel group "
                  << panelTypeId_ << ". Specular reflectivity = " << parameterValue
                  << ", resulting absorptivity = " << absorptivity << "."
                  << std::endl;
    }

    // Set the absorptivity for the panel group (even if it's negative)
    radiationPressureInterface_->setGroupAbsorptivity(panelTypeId_, absorptivity);
}

} // namespace estimatable_parameters

} // namespace tudat

