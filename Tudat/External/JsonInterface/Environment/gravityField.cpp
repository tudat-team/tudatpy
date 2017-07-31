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

#include "gravityField.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `GravityFieldSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< GravityFieldSettings >& gravityFieldSettings )
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
        boost::shared_ptr< CentralGravityFieldSettings > centralGravityFieldSettings =
                boost::dynamic_pointer_cast< CentralGravityFieldSettings >( gravityFieldSettings );
        enforceNonNullPointer( centralGravityFieldSettings );
        jsonObject[ K::gravitationalParameter ] = centralGravityFieldSettings->getGravitationalParameter( );
        return;
    }
    case central_spice:
        return;
    case spherical_harmonic:
    {
        boost::shared_ptr< SphericalHarmonicsGravityFieldSettings > sphericalHarmonicsGravityFieldSettings =
                boost::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >( gravityFieldSettings );
        enforceNonNullPointer( sphericalHarmonicsGravityFieldSettings );
        jsonObject[ K::associatedReferenceFrame ] =
                sphericalHarmonicsGravityFieldSettings->getAssociatedReferenceFrame( );

        /// SphericalHarmonicsFileGravityFieldSettings
        boost::shared_ptr< SphericalHarmonicsFileGravityFieldSettings > sphericalHarmonicsFileGravityFieldSettings =
                boost::dynamic_pointer_cast< SphericalHarmonicsFileGravityFieldSettings >( gravityFieldSettings );
        if ( sphericalHarmonicsFileGravityFieldSettings )
        {
            jsonObject[ K::file ] = path( sphericalHarmonicsFileGravityFieldSettings->fileName );
            jsonObject[ K::maximumDegree ] = sphericalHarmonicsFileGravityFieldSettings->maximumDegree;
            jsonObject[ K::maximumOrder ] = sphericalHarmonicsFileGravityFieldSettings->maximumOrder;

            // Gravitational parameter (index)
            const int gmIndex = sphericalHarmonicsFileGravityFieldSettings->gravitationalParameterIndex;
            jsonObject[ K::gravitationalParameterIndex ] = gmIndex;
            if ( ! ( gmIndex >= 0 ) )
            {
                jsonObject[ K::gravitationalParameter ] =
                        sphericalHarmonicsFileGravityFieldSettings->getGravitationalParameter( );
            }

            // Reference radius (index)
            const int rIndex = sphericalHarmonicsFileGravityFieldSettings->referenceRadiusIndex;
            jsonObject[ K::referenceRadiusIndex ] = rIndex;
            if ( ! ( rIndex >= 0 ) )
            {
                jsonObject[ K::referenceRadius ] =
                        sphericalHarmonicsFileGravityFieldSettings->getReferenceRadius( );
            }

            return;
        }

        jsonObject[ K::gravitationalParameter ] =
                sphericalHarmonicsGravityFieldSettings->getGravitationalParameter( );
        jsonObject[ K::referenceRadius ] =
                sphericalHarmonicsGravityFieldSettings->getReferenceRadius( );
        jsonObject[ K::cosineCoefficients ] =
                sphericalHarmonicsGravityFieldSettings->getCosineCoefficients( );
        jsonObject[ K::sineCoefficients ] =
                sphericalHarmonicsGravityFieldSettings->getSineCoefficients( );

        return;
    }
    default:
        jsonObject = handleUnimplementedEnumValueToJson( gravityFieldType, gravityFieldTypes,
                                                         unsupportedGravityFieldTypes );
    }
}

//! Create a shared pointer to a `GravityFieldSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< GravityFieldSettings >& gravityFieldSettings )
{
    using namespace json_interface;
    using K = Keys::Body::GravityField;

    // Get atmosphere model type
    const GravityFieldType gravityFieldType = getValue< GravityFieldType >( jsonObject, K::type );

    switch ( gravityFieldType ) {
    case central:
    {
        gravityFieldSettings = boost::make_shared< CentralGravityFieldSettings >(
                    getValue< double >( jsonObject, K::gravitationalParameter ) );
        return;
    }
    case central_spice:
    {
        gravityFieldSettings = boost::make_shared< GravityFieldSettings >( central_spice );
        return;
    }
    case spherical_harmonic:
    {
        const boost::shared_ptr< path > file = getOptional< path >( jsonObject, K::file );
        if ( file )  /// SphericalHarmonicsFileGravityFieldSettings
        {
            const int gmIndex = getNumeric( jsonObject, K::gravitationalParameterIndex, 0 );
            const int radiusIndex = getNumeric( jsonObject, K::referenceRadiusIndex, 1 );
            gravityFieldSettings = boost::make_shared< SphericalHarmonicsFileGravityFieldSettings >(
                        file->string( ),
                        getValue< std::string >( jsonObject, K::associatedReferenceFrame ),
                        getValue< int >( jsonObject, K::maximumDegree ),
                        getValue< int >( jsonObject, K::maximumOrder ),
                        gmIndex,
                        radiusIndex,
                        getNumeric< double >( jsonObject, K::gravitationalParameter, TUDAT_NAN, gmIndex >= 0 ),
                        getNumeric< double >( jsonObject, K::referenceRadius, TUDAT_NAN, radiusIndex >= 0 ) );
            return;
        }
        else   /// SphericalHarmonicsGravityFieldSettings
        {
            gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >(
                        getNumeric< double >( jsonObject, K::gravitationalParameter ),
                        getNumeric< double >( jsonObject, K::referenceRadius ),
                        getValue< Eigen::MatrixXd >( jsonObject, K::cosineCoefficients ),
                        getValue< Eigen::MatrixXd >( jsonObject, K::sineCoefficients ),
                        getValue< std::string >( jsonObject, K::associatedReferenceFrame ) );
            return;
        }
    }
    default:
        handleUnimplementedEnumValueFromJson( gravityFieldType, gravityFieldTypes, unsupportedGravityFieldTypes );
    }
}

} // namespace simulation_setup

} // namespace tudat
