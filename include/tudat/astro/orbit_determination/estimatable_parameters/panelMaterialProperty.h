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

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/system_models/vehicleExteriorPanels.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Estimatable scalar material property shared by a group of vehicle exterior panels.
/*!
 * The parameter is backed by all supplied panels, which are expected to have the same panel type identifier. Setting the
 * parameter writes the value to every panel. On construction and retrieval, inconsistent panel values are replaced by their
 * arithmetic mean so that the group remains represented by a single scalar.
 */
class PanelMaterialPropertyParameter : public EstimatableParameter< double >
{
public:
    //! Create an estimatable material property for one panel group.
    /*!
     * \param vehiclePanels Non-empty collection of panels sharing the estimated property. If the initial property values differ,
     * they are replaced by their arithmetic mean.
     * \param associatedBody Name of the body carrying the panels.
     * \param panelTypeId Panel type identifier of the group represented by this parameter.
     * \param parameterType Material property to estimate; must be an accommodation coefficient or the normal velocity at wall
     * ratio.
     */
    PanelMaterialPropertyParameter( const std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > >& vehiclePanels,
                                    const std::string& associatedBody,
                                    const std::string& panelTypeId,
                                    const EstimatebleParametersEnum parameterType );

    ~PanelMaterialPropertyParameter( ) override = default;

    //! Return the shared scalar value, first averaging and synchronizing inconsistent panel values.
    double getParameterValue( ) override;

    //! Set the material property to the supplied value on every panel in the group.
    void setParameterValue( const double parameterValue ) override;

    //! Return the scalar parameter size.
    int getParameterSize( ) override
    {
        return 1;
    }

private:
    //! Return the common value, replacing inconsistent panel values by their arithmetic mean.
    double normalizeValue( );

    //! Retrieve the selected material property from every panel.
    std::vector< double > getPanelMaterialPropertyValues( );

    //! Retrieve the selected material property from one panel.
    double getPanelMaterialPropertyValue( const std::shared_ptr< system_models::VehicleExteriorPanel >& vehiclePanel );

    //! Set the selected material property on one panel.
    void setPanelMaterialPropertyValue( const std::shared_ptr< system_models::VehicleExteriorPanel >& vehiclePanel,
                                        const double parameterValue );

    std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > vehiclePanels_;

    std::string panelTypeId_;
};

}  // namespace estimatable_parameters

}  // namespace tudat

#endif  // TUDAT_PANELMATERIALPROPERTY_H
