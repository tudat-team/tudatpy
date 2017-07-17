/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "aerodynamics.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `AerodynamicCoefficientSettings` object.
//! Called automatically by `nlohmann::json` when using a constructor such as
//! `json( aerodynamicCoefficientSettings )`.
void to_json( json& jsonObject,
              const boost::shared_ptr< AerodynamicCoefficientSettings >& aerodynamicCoefficientSettings )
{
    using namespace json_interface;

    // Initialise
    jsonObject = json( );

    // Get type
    jsonObject[ "type" ] = stringFromEnum( aerodynamicCoefficientSettings->getAerodynamicCoefficientType( ),
                                           aerodynamicCoefficientTypes );

    // Reference area
    jsonObject[ "referenceArea" ] = aerodynamicCoefficientSettings->getReferenceArea( );

    /// ConstantAerodynamicCoefficientSettings
    boost::shared_ptr< ConstantAerodynamicCoefficientSettings > constantAerodynamicCoefficientSettings =
            boost::dynamic_pointer_cast< ConstantAerodynamicCoefficientSettings >( aerodynamicCoefficientSettings );
    if ( constantAerodynamicCoefficientSettings )
    {
        // Force coefficients
        jsonObject[ "forceCoefficients" ] = stdFromEigen< double >(
                    constantAerodynamicCoefficientSettings->getConstantForceCoefficient( ) );
    }
}

} // namespace simulation_setup


namespace json_interface
{

//! Create a shared pointer to an `AerodynamicCoefficientSettings` object from a `json` object.
/*!
 * Create a shared pointer to an `AerodynamicCoefficientSettings` object from a `json` object.
 * \param settings `json` object containing only the settings for one aerodynamic coefficients settings.
 * \param fallbackArea Fallback reference area to be used when no reference area is speciefied in `settings`.
 * \return Shared pointer to an `AerodynamicCoefficientSettings` object.
 */
boost::shared_ptr< simulation_setup::AerodynamicCoefficientSettings > createAerodynamicCoefficientSettings(
        const json& settings, const double& fallbackArea )
{
    using namespace simulation_setup;

    // Get aerodynamic coefficient type (constant by default)
    const AerodynamicCoefficientTypes aerodynamicCoefficientType = enumFromString(
                getValue< std::string >( settings, "type", "constantAerodynamicCoefficients" ),
                aerodynamicCoefficientTypes );

    // Get reference area (use fallback value if not NaN when referenceArea is not provided)
    const double referenceArea = isnan( fallbackArea ) ? getValue< double >( settings, "referenceArea" )
                                                       : getValue( settings, "referenceArea", fallbackArea );

    if ( aerodynamicCoefficientType == constant_aerodynamic_coefficients )
    {
        // Get force coefficients
        // Read from key "forceCoefficients".
        // If not defined, use vector { CD, 0, 0 }, where CD is read from key "dragCoefficient".
        const Eigen::Vector3d forceCoefficients = eigenFromStd< double, 3 >(
                    getValue< std::vector< double > >( settings, "forceCoefficients",
        { getValue< double >( settings, "dragCoefficient" ), 0.0, 0.0 } ) );

        // Create and return settings
        return boost::make_shared< ConstantAerodynamicCoefficientSettings >( referenceArea, forceCoefficients );
    }
    else
    {
        throw std::runtime_error( stringFromEnum( aerodynamicCoefficientType, aerodynamicCoefficientTypes )
                                  + " not supported by json_interface." );
    }
}

} // namespace json_interfaces

} // namespace tudat
