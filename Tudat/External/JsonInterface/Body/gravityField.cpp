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
    jsonObject = json( json_interface::stringFromEnum( gravityFieldType, gravityFieldTypes ) );
}

//! Convert `json` to `GravityFieldType`.
void from_json( const json& jsonObject, GravityFieldType& gravityFieldType )
{
    gravityFieldType = json_interface::enumFromString( jsonObject.get< std::string >( ), gravityFieldTypes );
}

//! Create a `json` object from a shared pointer to a `GravityFieldSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< GravityFieldSettings >& gravityFieldSettings )
{
    using namespace json_interface;
    using Keys = Keys::Body::GravityField;

    // Initialise
    jsonObject = json( );

    // Type
    jsonObject[ Keys::type ] = gravityFieldSettings->getGravityFieldType( );

    /// central_spice
    if ( jsonObject[ Keys::type ] == central_spice )
    {
        return;
    }

    /// CentralGravityFieldSettings
    boost::shared_ptr< CentralGravityFieldSettings > centralGravityFieldSettings =
            boost::dynamic_pointer_cast< CentralGravityFieldSettings >( gravityFieldSettings );
    if ( centralGravityFieldSettings )
    {
        jsonObject[ Keys::gravitationalParameter ] = centralGravityFieldSettings->getGravitationalParameter( );
        return;
    }

    /// SphericalHarmonicsGravityFieldSettings
    boost::shared_ptr< SphericalHarmonicsGravityFieldSettings > sphericalHarmonicsGravityFieldSettings =
            boost::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >( gravityFieldSettings );
    if ( sphericalHarmonicsGravityFieldSettings )
    {
        jsonObject[ Keys::associatedReferenceFrame ] =
                sphericalHarmonicsGravityFieldSettings->getAssociatedReferenceFrame( );

        /// SphericalHarmonicsFileGravityFieldSettings
        boost::shared_ptr< SphericalHarmonicsFileGravityFieldSettings > sphericalHarmonicsFileGravityFieldSettings =
                boost::dynamic_pointer_cast< SphericalHarmonicsFileGravityFieldSettings >(
                    sphericalHarmonicsGravityFieldSettings );
        if ( sphericalHarmonicsFileGravityFieldSettings )
        {
            jsonObject[ Keys::file ] = path( sphericalHarmonicsFileGravityFieldSettings->fileName );
            jsonObject[ Keys::maximumDegree ] = sphericalHarmonicsFileGravityFieldSettings->maximumDegree;
            jsonObject[ Keys::maximumOrder ] = sphericalHarmonicsFileGravityFieldSettings->maximumOrder;

            // Gravitational parameter (index)
            const int gmIndex = sphericalHarmonicsFileGravityFieldSettings->gravitationalParameterIndex;
            jsonObject[ Keys::gravitationalParameterIndex ] = gmIndex;
            if ( ! ( gmIndex >= 0 ) )
            {
                jsonObject[ Keys::gravitationalParameter ] =
                        sphericalHarmonicsFileGravityFieldSettings->getGravitationalParameter( );
            }

            // Reference radius parameter (index)
            const int referenceRadiusIndex = sphericalHarmonicsFileGravityFieldSettings->referenceRadiusIndex;
            jsonObject[ Keys::referenceRadiusIndex ] = referenceRadiusIndex;
            if ( ! ( referenceRadiusIndex >= 0 ) )
            {
                jsonObject[ Keys::referenceRadius ] =
                        sphericalHarmonicsFileGravityFieldSettings->getReferenceRadius( );
            }

            return;
        }

        jsonObject[ Keys::gravitationalParameter ] =
                sphericalHarmonicsGravityFieldSettings->getGravitationalParameter( );
        jsonObject[ Keys::referenceRadius ] =
                sphericalHarmonicsGravityFieldSettings->getReferenceRadius( );
        jsonObject[ Keys::cosineCoefficients ] =
                sphericalHarmonicsGravityFieldSettings->getCosineCoefficients( );
        jsonObject[ Keys::sineCoefficients ] =
                sphericalHarmonicsGravityFieldSettings->getSineCoefficients( );
        return;
    }
}

} // namespace simulation_setup


namespace json_interface
{

//! Create a shared pointer to a `GravityFieldSettings` object from a `json` object.
boost::shared_ptr< simulation_setup::GravityFieldSettings > createGravityFieldSettings(
        const json& settings, const KeyTree& keyTree )
{
    using namespace simulation_setup;
    using Keys = Keys::Body::GravityField;

    // Get atmosphere model type
    const GravityFieldType gravityFieldType = getValue< GravityFieldType >( settings, keyTree + Keys::type );

    switch ( gravityFieldType ) {
    case central:
        return boost::make_shared< CentralGravityFieldSettings >(
                    getValue< double >( settings, keyTree + Keys::gravitationalParameter ) );
    case central_spice:
        return boost::make_shared< GravityFieldSettings >( central_spice );
    case spherical_harmonic:
    {
        const auto file = getValuePointer< path >( settings, keyTree + Keys::file );
        if ( file )  /// SphericalHarmonicsFileGravityFieldSettings
        {
            const int gmIndex = getNumeric( settings, keyTree + Keys::gravitationalParameterIndex, 0 );
            const int radiusIndex = getNumeric( settings, keyTree + Keys::referenceRadiusIndex, 1 );
            return boost::make_shared< SphericalHarmonicsFileGravityFieldSettings >(
                        file->string( ),
                        getValue< std::string >( settings, keyTree + Keys::associatedReferenceFrame ),
                        getNumeric< int >( settings, keyTree + Keys::maximumDegree ),
                        getNumeric< int >( settings, keyTree + Keys::maximumOrder ),
                        gmIndex,
                        radiusIndex,
                        getNumeric< double >( settings, keyTree + Keys::gravitationalParameter,
                                              TUDAT_NAN, gmIndex >= 0 ),
                        getNumeric< double >( settings, keyTree + Keys::referenceRadius,
                                              TUDAT_NAN, radiusIndex >= 0 ) );
        }
        else   /// SphericalHarmonicsGravityFieldSettings
        {
            return boost::make_shared< SphericalHarmonicsGravityFieldSettings >(
                        getNumeric< double >( settings, keyTree + Keys::gravitationalParameter ),
                        getNumeric< double >( settings, keyTree + Keys::referenceRadius ),
                        getValue< Eigen::MatrixXd >( settings, keyTree + Keys::cosineCoefficients ),
                        getValue< Eigen::MatrixXd >( settings, keyTree + Keys::sineCoefficients ),
                        getValue< std::string >( settings, keyTree + Keys::associatedReferenceFrame ) );
        }
    }
    default:
        throw std::runtime_error( stringFromEnum( gravityFieldType, gravityFieldTypes )
                                  + " not supported by json_interface." );
    }
}

} // namespace json_interface

} // namespace tudat
