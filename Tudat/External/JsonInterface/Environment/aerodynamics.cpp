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
    jsonObject = json_interface::stringFromEnum( aerodynamicCoefficientType, aerodynamicCoefficientTypes );
}

//! Convert `json` to `AerodynamicCoefficientTypes`.
void from_json( const json& jsonObject, AerodynamicCoefficientTypes& aerodynamicCoefficientType )
{
    aerodynamicCoefficientType = json_interface::enumFromString(
                jsonObject.get< std::string >( ), aerodynamicCoefficientTypes );
}

//! Create a `json` object from a shared pointer to a `AerodynamicCoefficientSettings` object.
void to_json( json& jsonObject,
              const boost::shared_ptr< AerodynamicCoefficientSettings >& aerodynamicCoefficientSettings )
{
    if ( aerodynamicCoefficientSettings )
    {
        using namespace json_interface;
        using K = Keys::Body::Aerodynamics;

        // Get type
        jsonObject[ K::type ] = aerodynamicCoefficientSettings->getAerodynamicCoefficientType( );

        // Reference area
        jsonObject[ K::referenceArea ] = aerodynamicCoefficientSettings->getReferenceArea( );

        /// ConstantAerodynamicCoefficientSettings
        boost::shared_ptr< ConstantAerodynamicCoefficientSettings > constantAerodynamicCoefficientSettings =
                boost::dynamic_pointer_cast< ConstantAerodynamicCoefficientSettings >( aerodynamicCoefficientSettings );
        if ( constantAerodynamicCoefficientSettings )
        {
            jsonObject[ K::forceCoefficients ] =
                    constantAerodynamicCoefficientSettings->getConstantForceCoefficient( );
            // FIXME: jsonObject[ K::momentCoefficients ] =
            //         constantAerodynamicCoefficientSettings->getConstantMomentCoefficient( );
            jsonObject[ K::areCoefficientsInAerodynamicFrame ] =
                    constantAerodynamicCoefficientSettings->getAreCoefficientsInAerodynamicFrame( );
            jsonObject[ K::areCoefficientsInNegativeAxisDirection ] =
                    constantAerodynamicCoefficientSettings->getAreCoefficientsInNegativeAxisDirection( );
            return;
        }

        // FIXME: derivered classes missing
    }
}

//! Create a `json` object from a shared pointer to a `AerodynamicCoefficientSettings` object.
void from_json( const json& jsonObject,
                boost::shared_ptr< AerodynamicCoefficientSettings >& aerodynamicCoefficientSettings )
{
    using namespace json_interface;
    using K = Keys::Body::Aerodynamics;

    // Aerodynamic coefficient type (constant by default)
    const AerodynamicCoefficientTypes aerodynamicCoefficientType =
            getValue( jsonObject, K::type, constant_aerodynamic_coefficients );

    // Reference area (use fallback area if reference area not provided, final value cannont be NaN)
    double fallbackReferenceArea = getNumeric< double >(
                jsonObject, SpecialKeys::up / K::referenceArea, TUDAT_NAN, true );

    switch ( aerodynamicCoefficientType )
    {
    case constant_aerodynamic_coefficients:
    {
        // Create settings (with default arguments)
        ConstantAerodynamicCoefficientSettings defaults( TUDAT_NAN, Eigen::Vector3d( ) );

        // Read forceCoefficients. If not defined, use [ dragCoefficient, 0, 0 ].
        Eigen::Vector3d forceCoefficients = Eigen::Vector3d::Zero( );
        const auto forceCoefficientsPointer =
                getValuePointer< Eigen::Vector3d >( jsonObject, K::forceCoefficients );
        if ( forceCoefficientsPointer )
        {
            forceCoefficients = *forceCoefficientsPointer;
        }
        else
        {
            forceCoefficients( 0 ) = getNumeric< double >( jsonObject, K::dragCoefficient );
        }

        // Return shared pointer
        aerodynamicCoefficientSettings = boost::make_shared< ConstantAerodynamicCoefficientSettings >(
                    getNumeric( jsonObject, K::referenceArea, fallbackReferenceArea ),
                    forceCoefficients,
                    getValue( jsonObject, K::areCoefficientsInAerodynamicFrame,
                              defaults.getAreCoefficientsInAerodynamicFrame( ) ),
                    getValue( jsonObject, K::areCoefficientsInNegativeAxisDirection,
                              defaults.getAreCoefficientsInNegativeAxisDirection( ) ) );
        return;
    }
    default:
        throw std::runtime_error( stringFromEnum( aerodynamicCoefficientType, aerodynamicCoefficientTypes )
                                  + " not supported by json_interface." );
    }
}

} // namespace simulation_setup

} // namespace tudat
