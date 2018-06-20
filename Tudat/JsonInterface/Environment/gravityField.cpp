/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "Tudat/JsonInterface/Environment/gravityField.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `GravityFieldSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< GravityFieldSettings >& gravityFieldSettings )
{
    if ( ! gravityFieldSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Body::GravityField;

    // Type
    const GravityFieldType gravityFieldType = gravityFieldSettings->getGravityFieldType( );
    jsonObject[ K::type ] = gravityFieldType;

    switch ( gravityFieldType )
    {
    case central:
    {
        std::shared_ptr< CentralGravityFieldSettings > centralGravityFieldSettings =
                std::dynamic_pointer_cast< CentralGravityFieldSettings >( gravityFieldSettings );
        assertNonnullptrPointer( centralGravityFieldSettings );
        jsonObject[ K::gravitationalParameter ] = centralGravityFieldSettings->getGravitationalParameter( );
        return;
    }
    case central_spice:
        return;
    case spherical_harmonic:
    {
        std::shared_ptr< SphericalHarmonicsGravityFieldSettings > shGravityFieldSettings =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >( gravityFieldSettings );
        assertNonnullptrPointer( shGravityFieldSettings );

        // FromFileSphericalHarmonicsGravityFieldSettings
        std::shared_ptr< FromFileSphericalHarmonicsGravityFieldSettings > shModelGravityFieldSettings =
                std::dynamic_pointer_cast< FromFileSphericalHarmonicsGravityFieldSettings >( gravityFieldSettings );
        if ( shModelGravityFieldSettings )
        {
            const SphericalHarmonicsModel model = shModelGravityFieldSettings->getSphericalHarmonicsModel( );
            if ( model == customModel )
            {
                jsonObject[ K::file ] = boost::filesystem::path( shModelGravityFieldSettings->getFilePath( ) );
                jsonObject[ K::associatedReferenceFrame ] = shModelGravityFieldSettings->getAssociatedReferenceFrame( );
                jsonObject[ K::maximumDegree ] = shModelGravityFieldSettings->getMaximumDegree( );
                jsonObject[ K::maximumOrder ] = shModelGravityFieldSettings->getMaximumOrder( );

                // Gravitational parameter (index)
                const int gmIndex = shModelGravityFieldSettings->getGravitationalParameterIndex( );
                jsonObject[ K::gravitationalParameterIndex ] = gmIndex;
                if ( gmIndex < 0 )
                {
                    jsonObject[ K::gravitationalParameter ] =
                            shModelGravityFieldSettings->getGravitationalParameter( );
                }

                // Reference radius (index)
                const int rIndex = shModelGravityFieldSettings->getReferenceRadiusIndex( );
                jsonObject[ K::referenceRadiusIndex ] = rIndex;
                if ( rIndex < 0 )
                {
                    jsonObject[ K::referenceRadius ] =
                            shModelGravityFieldSettings->getReferenceRadius( );
                }
            }
            else
            {
                jsonObject[ K::model ] = model;
            }
        }
        else
        {
            jsonObject[ K::gravitationalParameter ] = shGravityFieldSettings->getGravitationalParameter( );
            jsonObject[ K::referenceRadius ] = shGravityFieldSettings->getReferenceRadius( );
            jsonObject[ K::cosineCoefficients ] = shGravityFieldSettings->getCosineCoefficients( );
            jsonObject[ K::sineCoefficients ] = shGravityFieldSettings->getSineCoefficients( );
            jsonObject[ K::associatedReferenceFrame ] = shGravityFieldSettings->getAssociatedReferenceFrame( );
        }
        return;
    }
    default:
        handleUnimplementedEnumValue( gravityFieldType, gravityFieldTypes, unsupportedGravityFieldTypes );
    }
}

//! Create a shared pointer to a `GravityFieldSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< GravityFieldSettings >& gravityFieldSettings )
{
    using namespace json_interface;
    using K = Keys::Body::GravityField;

    // Get atmosphere model type
    const GravityFieldType gravityFieldType = getValue< GravityFieldType >( jsonObject, K::type );

    switch ( gravityFieldType ) {
    case central:
    {
        gravityFieldSettings = std::make_shared< CentralGravityFieldSettings >(
                    getValue< double >( jsonObject, K::gravitationalParameter ) );
        return;
    }
    case central_spice:
    {
        gravityFieldSettings = std::make_shared< GravityFieldSettings >( central_spice );
        return;
    }
    case spherical_harmonic:
    {
        // load coefficients from custom file
        if ( isDefined( jsonObject, K::file ) )
        {
            const int gmIndex = getValue( jsonObject, K::gravitationalParameterIndex, 0 );
            const int radiusIndex = getValue( jsonObject, K::referenceRadiusIndex, 1 );
            gravityFieldSettings = std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >(
                        getValue< boost::filesystem::path >( jsonObject, K::file ).string( ),
                        getValue< std::string >( jsonObject, K::associatedReferenceFrame ),
                        getValue< int >( jsonObject, K::maximumDegree ),
                        getValue< int >( jsonObject, K::maximumOrder ),
                        gmIndex,
                        radiusIndex,
                        getValue< double >( jsonObject, K::gravitationalParameter, TUDAT_NAN ),
                        getValue< double >( jsonObject, K::referenceRadius, TUDAT_NAN ) );
            return;
        }

        // load coefficients from model included in Tudat
        if ( isDefined( jsonObject, K::model ) )
        {
            gravityFieldSettings = std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >(
                        getValue< SphericalHarmonicsModel >( jsonObject, K::model ) );
            return;
        }

        // user-provided coefficients in JSON object
        gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                    getValue< double >( jsonObject, K::gravitationalParameter ),
                    getValue< double >( jsonObject, K::referenceRadius ),
                    getValue< Eigen::MatrixXd >( jsonObject, K::cosineCoefficients ),
                    getValue< Eigen::MatrixXd >( jsonObject, K::sineCoefficients ),
                    getValue< std::string >( jsonObject, K::associatedReferenceFrame ) );
        return;
    }
    default:
        handleUnimplementedEnumValue( gravityFieldType, gravityFieldTypes, unsupportedGravityFieldTypes );
    }
}

} // namespace simulation_setup

} // namespace tudat
