/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SPECULARREFLECTIVITY_H
#define TUDAT_SPECULARREFLECTIVITY_H




#include <cmath>

#include <Eigen/Core>

#include "tudat/basics/utilities.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/electromagnetism/radiationPressureInterface.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for the estimation of a specular reflectivity coefficient per group of panels
class SpecularReflectivity: public EstimatableParameter< double >
{
    public:
        //! Constructor.
        /*!
        * ConstructorSpecularReflectivity
        * \param radiationPressureInterface Object containing the radiation pressure coefficient to be estimated.
        * \param associatedBody Name of body containing the radiationPressureInterface object
        * \param panelTypeID Name of panel group for which to estimate the coefficient
        */
        SpecularReflectivity(
                const std::shared_ptr< electromagnetism::PaneledRadiationPressureTargetModel > radiationPressureInterface,
                const std::string& associatedBody,
                const std::string panelTypeId):
            EstimatableParameter< double >( specular_reflectivity, associatedBody, panelTypeId ),
            radiationPressureInterface_( radiationPressureInterface ),
            panelTypeId_( panelTypeId )
        {
            // check if the panelTypeID exists among the panels
            std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > panelsFromId =
            radiationPressureInterface_->getPanelsFromId( panelTypeId_);
            if(panelsFromId.size() < 1){
                throw std::runtime_error( "Error when creating estimated specular reflectivity for " +
                panelTypeId + " of " + associatedBody + ", no corresponding panels defined" );
            }

        }
        //! Destructor.
        ~SpecularReflectivity( ) { }

        double getParameterValue( )
        {
            return radiationPressureInterface_->getAverageSpecularReflectivity( panelTypeId_ );
        }

        void setParameterValue( double parameterValue )
        {
            // set the diffuse reflectivity
            if ( parameterValue > 1.0 ){
                parameterValue = 1.0;
            }
            radiationPressureInterface_->setGroupSpecularReflectivity(panelTypeId_, parameterValue );

            // update the absorptivity
            double diffuseReflectivity = radiationPressureInterface_->getAverageDiffuseReflectivity( panelTypeId_ );
            if(parameterValue + diffuseReflectivity >= 1.0)
            {
                std::cerr << "Warning: When updating reflectivity coefficients for panel group "
                << panelTypeId_ << ", sum of specular and diffuse is larger than 1. Setting absorptivity to zero" << std::endl;
                double absorptivity = 0;
                radiationPressureInterface_->setGroupAbsorptivity(panelTypeId_, absorptivity);

            }
            else
            {
                double absorptivity = 1.0 - parameterValue - diffuseReflectivity;
                radiationPressureInterface_->setGroupAbsorptivity(panelTypeId_, absorptivity);
            }

        }

        int getParameterSize( ){ return 1; }

    protected:

    private:

        std::shared_ptr< electromagnetism::PaneledRadiationPressureTargetModel > radiationPressureInterface_;

        std::string panelTypeId_;

};


} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_SPECULARREFLECTIVITY_H