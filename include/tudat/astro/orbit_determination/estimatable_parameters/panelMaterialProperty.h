/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PANELMATERIALPROPERTY_H
#define TUDAT_PANELMATERIALPROPERTY_H

#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/system_models/vehicleExteriorPanels.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for the estimation of a scalar material property per group of panels.
class PanelMaterialPropertyParameter : public EstimatableParameter< double >
{
public:
    PanelMaterialPropertyParameter( const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > vehiclePanels,
                                    const std::string& associatedBody,
                                    const std::string& panelTypeId,
                                    const EstimatebleParametersEnum parameterType );

    ~PanelMaterialPropertyParameter( ) {}

    double normalizeValue( );

    double getParameterValue( );

    void setParameterValue( double parameterValue );

    int getParameterSize( )
    {
        return 1;
    }

private:
    std::vector< double > getPanelMaterialPropertyValues( );

    double getPanelMaterialPropertyValue( const std::shared_ptr< system_models::VehicleExteriorPanel > vehiclePanel );

    void setPanelMaterialPropertyValue( const std::shared_ptr< system_models::VehicleExteriorPanel > vehiclePanel,
                                        const double parameterValue );

    std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > vehiclePanels_;

    std::string panelTypeId_;
};

}  // namespace estimatable_parameters

}  // namespace tudat

#endif  // TUDAT_PANELMATERIALPROPERTY_H
