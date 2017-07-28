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
        using K = Keys::Body::ShapeModel;

        // Base class settings
        jsonObject[ K::type ] = bodyShapeSettings->getBodyShapeType( );

        /// spherical_spice
        if ( jsonObject[ K::type ] == spherical_spice )
        {
            return;
        }

        /// SphericalBodyShapeSettings
        boost::shared_ptr< SphericalBodyShapeSettings > sphericalBodyShapeSettings =
                boost::dynamic_pointer_cast< SphericalBodyShapeSettings >( bodyShapeSettings );
        if ( sphericalBodyShapeSettings )
        {
            jsonObject[ K::radius ] = sphericalBodyShapeSettings->getRadius( );
            return;
        }

        /// OblateSphericalBodyShapeSettings
        boost::shared_ptr< OblateSphericalBodyShapeSettings > oblateSphericalBodyShapeSettings =
                boost::dynamic_pointer_cast< OblateSphericalBodyShapeSettings >( bodyShapeSettings );
        if ( oblateSphericalBodyShapeSettings )
        {
            jsonObject[ K::equatorialRadius ] = oblateSphericalBodyShapeSettings->getEquatorialRadius( );
            jsonObject[ K::flattening ] = oblateSphericalBodyShapeSettings->getFlattening( );
            return;
        }
    }
}

//! Create a shared pointer to a `BodyShapeSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< BodyShapeSettings >& bodyShapeSettings )
{
    using namespace json_interface;
    using K = Keys::Body::ShapeModel;

    // Base class settings
    const BodyShapeTypes bodyShapeType = getValue< BodyShapeTypes >( jsonObject, K::type );

    switch ( bodyShapeType ) {
    case spherical:
    {
        bodyShapeSettings = boost::make_shared< SphericalBodyShapeSettings >(
                    getNumeric< double >( jsonObject, K::radius ) );
        return;
    }
    case spherical_spice:
    {
        bodyShapeSettings = boost::make_shared< BodyShapeSettings >( bodyShapeType );
        return;
    }
    case oblate_spheroid:
    {
        bodyShapeSettings = boost::make_shared< OblateSphericalBodyShapeSettings >(
                    getNumeric< double >( jsonObject, K::equatorialRadius ),
                    getNumeric< double >( jsonObject, K::flattening ) );
        return;
    }
    default:
        throw std::runtime_error( stringFromEnum( bodyShapeType, bodyShapeTypes )
                                  + " not supported by json_interface." );
    }
}

} // namespace simulation_setup

} // namespace tudat
