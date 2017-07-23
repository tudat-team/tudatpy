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

#include "shapeModel.h"

namespace tudat
{

namespace simulation_setup
{

//! Convert `BodyShapeTypes`s to `json`.
void to_json( json& jsonObject, const BodyShapeTypes& bodyShapeType )
{
    jsonObject = json_interface::stringFromEnum( bodyShapeType, bodyShapeTypes );
}

//! Convert `json` to `BodyShapeTypes`.
void from_json( const json& jsonObject, BodyShapeTypes& bodyShapeType )
{
    bodyShapeType = json_interface::enumFromString( jsonObject.get< std::string >( ), bodyShapeTypes );
}

//! Create a `json` object from a shared pointer to a `BodyShapeSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< BodyShapeSettings >& bodyShapeSettings )
{
    if ( bodyShapeSettings )
    {
        using namespace json_interface;
        using Keys = Keys::Body::ShapeModel;

        // Base class settings
        jsonObject[ Keys::type ] = bodyShapeSettings->getBodyShapeType( );

        /// spherical_spice
        if ( jsonObject[ Keys::type ] == spherical_spice )
        {
            return;
        }

        /// SphericalBodyShapeSettings
        boost::shared_ptr< SphericalBodyShapeSettings > sphericalBodyShapeSettings =
                boost::dynamic_pointer_cast< SphericalBodyShapeSettings >( bodyShapeSettings );
        if ( sphericalBodyShapeSettings )
        {
            jsonObject[ Keys::radius ] = sphericalBodyShapeSettings->getRadius( );
            return;
        }

        /// OblateSphericalBodyShapeSettings
        boost::shared_ptr< OblateSphericalBodyShapeSettings > oblateSphericalBodyShapeSettings =
                boost::dynamic_pointer_cast< OblateSphericalBodyShapeSettings >( bodyShapeSettings );
        if ( oblateSphericalBodyShapeSettings )
        {
            jsonObject[ Keys::equatorialRadius ] = oblateSphericalBodyShapeSettings->getEquatorialRadius( );
            jsonObject[ Keys::flattening ] = oblateSphericalBodyShapeSettings->getFlattening( );
            return;
        }
    }
}

} // namespace simulation_setup


namespace json_interface
{

//! Create a shared pointer to a `BodyShapeSettings` object from a `json` object.
boost::shared_ptr< simulation_setup::BodyShapeSettings > createShapeModelSettings(
        const json& settings, const KeyTree& keyTree )
{
    using namespace simulation_setup;
    using Keys = Keys::Body::ShapeModel;

    // Base class settings
    const BodyShapeTypes bodyShapeType = getValue< BodyShapeTypes >( settings, keyTree + Keys::type );

    switch ( bodyShapeType ) {
    case spherical:
        return boost::make_shared< SphericalBodyShapeSettings >(
                    getNumeric< double >( settings, keyTree + Keys::radius ) );
    case spherical_spice:
        return boost::make_shared< BodyShapeSettings >( bodyShapeType );
    case oblate_spheroid:
        return boost::make_shared< OblateSphericalBodyShapeSettings >(
                    getNumeric< double >( settings, keyTree + Keys::equatorialRadius ),
                    getNumeric< double >( settings, keyTree + Keys::flattening ) );
    default:
        throw std::runtime_error( stringFromEnum( bodyShapeType, bodyShapeTypes )
                                  + " not supported by json_interface." );
    }
}

} // namespace json_interface

} // namespace tudat
