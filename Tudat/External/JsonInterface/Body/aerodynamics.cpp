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

//! Convert `AerodynamicCoefficientTypes` to `json`.
void to_json( json& jsonObject, const AerodynamicCoefficientTypes& aerodynamicCoefficientType )
{
    jsonObject = json( json_interface::stringFromEnum( aerodynamicCoefficientType, aerodynamicCoefficientTypes ) );
}

//! Convert `json` to `AerodynamicCoefficientTypes`.
void from_json( const json& jsonObject, AerodynamicCoefficientTypes& aerodynamicCoefficientType )
{
    aerodynamicCoefficientType = json_interface::enumFromString(
                jsonObject.get< std::string >( ), aerodynamicCoefficientTypes );
}

//! Create a `json` object from a shared pointer to a `AerodynamicCoefficientSettings` object.
//! Called automatically by `nlohmann::json` when using a constructor such as
//! `json( aerodynamicCoefficientSettings )`.
void to_json( json& jsonObject,
              const boost::shared_ptr< AerodynamicCoefficientSettings >& aerodynamicCoefficientSettings )
{
    using namespace json_interface;
    using Keys = Keys::Body::Aerodynamics;

    // Initialise
    jsonObject = json( );

    // Get type
    jsonObject[ Keys::type ] = aerodynamicCoefficientSettings->getAerodynamicCoefficientType( );

    // Reference area
    jsonObject[ Keys::referenceArea ] = aerodynamicCoefficientSettings->getReferenceArea( );

    /// ConstantAerodynamicCoefficientSettings
    boost::shared_ptr< ConstantAerodynamicCoefficientSettings > constantAerodynamicCoefficientSettings =
            boost::dynamic_pointer_cast< ConstantAerodynamicCoefficientSettings >( aerodynamicCoefficientSettings );
    if ( constantAerodynamicCoefficientSettings )
    {
        jsonObject[ Keys::forceCoefficients ] = constantAerodynamicCoefficientSettings->constantForceCoefficient_;
        jsonObject[ Keys::momentCoefficients ] = constantAerodynamicCoefficientSettings->constantMomentCoefficient_;
        jsonObject[ Keys::areCoefficientsInAerodynamicFrame ] =
                constantAerodynamicCoefficientSettings->areCoefficientsInAerodynamicFrame_;
        jsonObject[ Keys::areCoefficientsInNegativeAxisDirection ] =
                constantAerodynamicCoefficientSettings->areCoefficientsInNegativeAxisDirection_;
        return;
    }

    // FIXME: derivered classes missing
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
        const json& settings, const KeyTree& keyTree, const double& fallbackArea )
{
    using namespace simulation_setup;
    using Keys = Keys::Body::Aerodynamics;

    // Aerodynamic coefficient type (constant by default)
    const AerodynamicCoefficientTypes aerodynamicCoefficientType =
            getValue( settings, keyTree + Keys::type, constant_aerodynamic_coefficients );

    // Reference area (use fallback area if reference area not provided, final value cannont be NaN)
    const double referenceArea = getNumeric( settings, keyTree + Keys::referenceArea, fallbackArea );

    switch ( aerodynamicCoefficientType )
    {
    case constant_aerodynamic_coefficients:
    {
        // Read forceCoefficients. If not defined, use [ dragCoefficient, 0, 0 ].
        Eigen::Vector3d forceCoefficients = Eigen::Vector3d::Zero( );
        const auto forceCoefficientsPointer =
                getValuePointer< Eigen::Vector3d >( settings, keyTree + Keys::forceCoefficients );
        if ( forceCoefficientsPointer )
        {
            forceCoefficients = *forceCoefficientsPointer;
        }
        else
        {
            forceCoefficients( 0 ) = getNumeric< double >( settings, keyTree + Keys::dragCoefficient );
        }

        // Create settings (with default arguments)
        ConstantAerodynamicCoefficientSettings constantAerodynamicCoefficientSettings(
                    referenceArea, forceCoefficients );

        // Moment coefficients
        const auto momentCoefficients =
                getValuePointer< Eigen::Vector3d >( settings, keyTree + Keys::momentCoefficients );
        if ( momentCoefficients )
        {
            constantAerodynamicCoefficientSettings.constantMomentCoefficient_ = *momentCoefficients ;
        }

        // areCoefficientsInAerodynamicFrame
        const auto areCoefficientsInAerodynamicFrame =
                getValuePointer< bool >( settings, keyTree + Keys::areCoefficientsInAerodynamicFrame );
        if ( areCoefficientsInAerodynamicFrame )
        {
            constantAerodynamicCoefficientSettings.areCoefficientsInAerodynamicFrame_ =
                    *areCoefficientsInAerodynamicFrame;
        }

        // areCoefficientsInNegativeAxisDirection
        const auto areCoefficientsInNegativeAxisDirection =
                getValuePointer< bool >( settings, keyTree + Keys::areCoefficientsInNegativeAxisDirection );
        if ( areCoefficientsInNegativeAxisDirection )
        {
            constantAerodynamicCoefficientSettings.areCoefficientsInNegativeAxisDirection_ =
                    *areCoefficientsInNegativeAxisDirection;
        }

        // Return shared pointer
        return boost::make_shared< ConstantAerodynamicCoefficientSettings >( constantAerodynamicCoefficientSettings );
    }
    default:
        throw std::runtime_error( stringFromEnum( aerodynamicCoefficientType, aerodynamicCoefficientTypes )
                                  + " not supported by json_interface." );
    }
}

} // namespace json_interface

} // namespace tudat
