/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DIFFUSEREFLECTIVITY_H
#define TUDAT_DIFFUSEREFLECTIVITY_H




#include <cmath>

#include <Eigen/Core>

#include "tudat/basics/utilities.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/electromagnetism/radiationPressureInterface.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for the estimation of a diffuse reflectivity coefficient per group of panels
class DiffuseReflectivity: public EstimatableParameter< double >
{
    public:
        //! Constructor.
        /*!
        * Constructor DiffuseReflectivity
        * \param radiationPressureInterface Object containing the radiation pressure coefficient to be estimated.
        * \param associatedBody Name of body containing the radiationPressureInterface object
        * \param panelTypeID Name of panel group for which to estimate the coefficient
        */
        DiffuseReflectivity(
                const std::shared_ptr< electromagnetism::PaneledRadiationPressureTargetModel > radiationPressureInterface,
                const std::string& associatedBody,
                const std::string panelTypeId):
            EstimatableParameter< double >( diffuse_reflectivity, associatedBody, panelTypeId ),
            radiationPressureInterface_( radiationPressureInterface ),
            panelTypeId_( panelTypeId )
        {
            // check if the panelTypeID exists among the panels
            std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > panelsFromId =
            radiationPressureInterface_->getPanelsFromId( panelTypeId_);
            if(panelsFromId.size() < 1){
                throw std::runtime_error( "Error when creating estimated diffuse reflectivity for " +
                panelTypeId + " of " + associatedBody + ", no corresponding panels defined" );
            }

        }
        //! Destructor.
        ~DiffuseReflectivity( ) { }

        double getParameterValue( )
        {
            // Retrieve all diffuse reflectivity values for the panels corresponding to the given panelTypeId
            std::vector<double> diffuseReflectivities = radiationPressureInterface_->getDiffuseReflectivityForPanelTypeId(panelTypeId_);

            // Check if all values are the same
            bool allValuesSame = std::all_of(diffuseReflectivities.begin(), diffuseReflectivities.end(),
                                             [&](double value) { return value == diffuseReflectivities[0]; });

            // If not all values are the same, print a warning and reset the values to the average
            if (!allValuesSame)
            {
                std::cerr << "Warning: Diffuse reflectivity values for panel group "
                          << panelTypeId_ << " are not consistent. Resetting all to the average value." << std::endl;

                // Calculate the average specular reflectivity
                double averageDiffuseReflectivity = std::accumulate(diffuseReflectivities.begin(), diffuseReflectivities.end(), 0.0) / diffuseReflectivities.size();

                // Set all panels' diffuse reflectivity to the average value
                radiationPressureInterface_->setGroupDiffuseReflectivity(panelTypeId_, averageDiffuseReflectivity);
            }

            // Return the average value (or the original value if all values were the same)
            return diffuseReflectivities[0];
        }


        void setParameterValue( double parameterValue )
        {
            // Set the diffuse reflectivity, even if it's greater than 1 (non-physical case)
            radiationPressureInterface_->setGroupDiffuseReflectivity(panelTypeId_, parameterValue);

            // Get the current specular reflectivity
            double specularReflectivity = radiationPressureInterface_->getAverageSpecularReflectivity( panelTypeId_ );

            // Compute the absorptivity such that the sum of specular + diffuse + absorptivity = 1
            double absorptivity = 1.0 - parameterValue - specularReflectivity;

            // Check if the absorptivity is negative, which indicates non-physical behavior
            if ( absorptivity < 0.0 )
            {
                // Print a warning with both the non-physical specular reflectivity and the resulting negative absorptivity
                std::cerr << "Warning: Non-physical behavior detected for panel group "
                          << panelTypeId_ << ". Diffuse reflectivity = " << parameterValue
                          << ", resulting absorptivity = " << absorptivity << "."
                          << std::endl;
            }

            // Set the absorptivity for the panel group (even if it's negative)
            radiationPressureInterface_->setGroupAbsorptivity(panelTypeId_, absorptivity);

        }

        int getParameterSize( ){ return 1; }

    protected:

    private:

        std::shared_ptr< electromagnetism::PaneledRadiationPressureTargetModel > radiationPressureInterface_;

        std::string panelTypeId_;

};


} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_DIFFUSEREFLECTIVITY_H