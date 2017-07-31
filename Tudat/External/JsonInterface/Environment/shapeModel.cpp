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

//! Create a `json` object from a shared pointer to a `BodyShapeSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< BodyShapeSettings >& bodyShapeSettings )
{
    if ( ! bodyShapeSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Body::ShapeModel;

    const BodyShapeTypes bodyShapeType = bodyShapeSettings->getBodyShapeType( );
    jsonObject[ K::type ] = bodyShapeType;

    switch ( bodyShapeType )
    {
    case spherical:
    {
        boost::shared_ptr< SphericalBodyShapeSettings > sphericalBodyShapeSettings =
                boost::dynamic_pointer_cast< SphericalBodyShapeSettings >( bodyShapeSettings );
        enforceNonNullPointer( sphericalBodyShapeSettings );
        jsonObject[ K::radius ] = sphericalBodyShapeSettings->getRadius( );
        return;
    }
    case spherical_spice:
        return;
    case oblate_spheroid:
    {
        boost::shared_ptr< OblateSphericalBodyShapeSettings > oblateSphericalBodyShapeSettings =
                boost::dynamic_pointer_cast< OblateSphericalBodyShapeSettings >( bodyShapeSettings );
        enforceNonNullPointer( oblateSphericalBodyShapeSettings );
        jsonObject[ K::equatorialRadius ] = oblateSphericalBodyShapeSettings->getEquatorialRadius( );
        jsonObject[ K::flattening ] = oblateSphericalBodyShapeSettings->getFlattening( );
        return;
    }
    default:
        jsonObject = handleUnimplementedEnumValueToJson( bodyShapeType, bodyShapeTypes, unsupportedBodyShapeTypes );
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
                    getValue< double >( jsonObject, K::flattening ) );
        return;
    }
    default:
        handleUnimplementedEnumValueFromJson( bodyShapeType, bodyShapeTypes, unsupportedBodyShapeTypes );
    }
}

} // namespace simulation_setup

} // namespace tudat
