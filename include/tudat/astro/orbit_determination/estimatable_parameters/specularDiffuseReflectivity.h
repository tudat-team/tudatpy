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
class SpecularDiffuseReflectivityParameter: public EstimatableParameter< double >
{
    public:
        //! Constructor.
        /*!
        * ConstructorSpecularReflectivity
        * \param radiationPressureInterface Object containing the radiation pressure coefficient to be estimated.
        * \param associatedBody Name of body containing the radiationPressureInterface object
        * \param panelTypeID Name of panel group for which to estimate the coefficient
        */
        SpecularDiffuseReflectivityParameter(
                const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > vehiclePanels,
                const std::string& associatedBody,
                const std::string& panelTypeId,
                const EstimatebleParametersEnum parameterType );

        //! Destructor.
        ~SpecularDiffuseReflectivityParameter( ) { }

        double normalizeValue( );

        double getParameterValue( );

        void setParameterValue( double parameterValue );

        int getParameterSize( ){ return 1; }

    protected:

    private:

        std::vector< double > getPanelReflectivities( );

        std::vector< std::shared_ptr< electromagnetism::SpecularDiffuseMixReflectionLaw > > reflectionLaws_;

        std::string panelTypeId_;

};


} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_SPECULARREFLECTIVITY_H