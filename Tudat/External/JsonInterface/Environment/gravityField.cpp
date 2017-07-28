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

//! Convert `GravityFieldType`s to `json`.
void to_json( json& jsonObject, const GravityFieldType& gravityFieldType )
{
    jsonObject = json_interface::stringFromEnum( gravityFieldType, gravityFieldTypes );
}

//! Convert `json` to `GravityFieldType`.
void from_json( const json& jsonObject, GravityFieldType& gravityFieldType )
{
    gravityFieldType = json_interface::enumFromString( jsonObject.get< std::string >( ), gravityFieldTypes );
}

//! Create a `json` object from a shared pointer to a `GravityFieldSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< GravityFieldSettings >& gravityFieldSettings )
{
    if ( gravityFieldSettings )
    {
        using namespace json_interface;
        using K = Keys::Body::GravityField;

        // Type
        jsonObject[ K::type ] = gravityFieldSettings->getGravityFieldType( );

        /// central_spice
        if ( jsonObject[ K::type ] == central_spice )
        {
            return;
        }

        /// CentralGravityFieldSettings
        boost::shared_ptr< CentralGravityFieldSettings > centralGravityFieldSettings =
                boost::dynamic_pointer_cast< CentralGravityFieldSettings >( gravityFieldSettings );
        if ( centralGravityFieldSettings )
        {
            jsonObject[ K::gravitationalParameter ] = centralGravityFieldSettings->getGravitationalParameter( );
            return;
        }

        /// SphericalHarmonicsGravityFieldSettings
        boost::shared_ptr< SphericalHarmonicsGravityFieldSettings > sphericalHarmonicsGravityFieldSettings =
                boost::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >( gravityFieldSettings );
        if ( sphericalHarmonicsGravityFieldSettings )
        {
            jsonObject[ K::associatedReferenceFrame ] =
                    sphericalHarmonicsGravityFieldSettings->getAssociatedReferenceFrame( );

            /// SphericalHarmonicsFileGravityFieldSettings
            boost::shared_ptr< SphericalHarmonicsFileGravityFieldSettings > sphericalHarmonicsFileGravityFieldSettings =
                    boost::dynamic_pointer_cast< SphericalHarmonicsFileGravityFieldSettings >(
                        sphericalHarmonicsGravityFieldSettings );
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
        const auto file = getValuePointer< path >( jsonObject, K::file );
        if ( file )  /// SphericalHarmonicsFileGravityFieldSettings
        {
            const int gmIndex = getNumeric( jsonObject, K::gravitationalParameterIndex, 0 );
            const int radiusIndex = getNumeric( jsonObject, K::referenceRadiusIndex, 1 );
            gravityFieldSettings = boost::make_shared< SphericalHarmonicsFileGravityFieldSettings >(
                        file->string( ),
                        getValue< std::string >( jsonObject, K::associatedReferenceFrame ),
                        getNumeric< int >( jsonObject, K::maximumDegree ),
                        getNumeric< int >( jsonObject, K::maximumOrder ),
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
        throw std::runtime_error( stringFromEnum( gravityFieldType, gravityFieldTypes )
                                  + " not supported by json_interface." );
    }
}

} // namespace simulation_setup

} // namespace tudat
