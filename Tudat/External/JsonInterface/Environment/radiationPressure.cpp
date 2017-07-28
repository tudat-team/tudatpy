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

#include "radiationPressure.h"

namespace tudat
{

namespace simulation_setup
{

//! Convert `RadiationPressureType` to `json`.
void to_json( json& jsonObject, const RadiationPressureType& radiationPressureType )
{
    jsonObject = json_interface::stringFromEnum( radiationPressureType, radiationPressureTypes );
}

//! Convert `json` to `RadiationPressureType`.
void from_json( const json& jsonObject, RadiationPressureType& radiationPressureType )
{
    radiationPressureType = json_interface::enumFromString( jsonObject.get< std::string >( ), radiationPressureTypes );
}

//! Create a `json` object from a shared pointer to a `RadiationPressureInterfaceSettings` object.
void to_json( json& jsonObject,
              const boost::shared_ptr< RadiationPressureInterfaceSettings >& radiationPressureInterfaceSettings )
{
    if ( ! radiationPressureInterfaceSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Body::RadiationPressure;

    const RadiationPressureType radiationPressureType =
            radiationPressureInterfaceSettings->getRadiationPressureType( );
    jsonObject[ K::type ] = radiationPressureType;

    switch ( radiationPressureType )
    {
    case cannon_ball:
    {
        boost::shared_ptr< CannonBallRadiationPressureInterfaceSettings > cannonBallRadiationPressureInterfaceSettings =
                boost::dynamic_pointer_cast< CannonBallRadiationPressureInterfaceSettings >(
                    radiationPressureInterfaceSettings );
        enforceNonNullPointer( cannonBallRadiationPressureInterfaceSettings );
        jsonObject[ K::referenceArea ] = cannonBallRadiationPressureInterfaceSettings->getArea( );
        jsonObject[ K::radiationPressureCoefficient ] =
                cannonBallRadiationPressureInterfaceSettings->getRadiationPressureCoefficient( );
        jsonObject[ K::ocultingBodies ] =
                cannonBallRadiationPressureInterfaceSettings->getOccultingBodies( );
        return;
    }
    default:
        jsonObject = handleUnimplementedEnumValueToJson( radiationPressureType, radiationPressureTypes,
                                                         unsupportedRadiationPressureTypes );
    }
}

//! Create a `json` object from a shared pointer to a `RadiationPressureInterfaceSettings` object.
void from_json( const json& jsonObject,
                boost::shared_ptr< RadiationPressureInterfaceSettings >& radiationPressureInterfaceSettings )
{
    using namespace json_interface;
    using K = Keys::Body::RadiationPressure;

    // Get radiation pressure coefficient type (cannonBall by default)
    const RadiationPressureType radiationPressureType = getValue( jsonObject, K::type, cannon_ball );

    const std::string sourceBodyName = getParentKey(
                jsonObject, "Creating RadiationPressureInterfaceSettings out of context is not allowed "
                            "because the name of the radiating body cannot be inferred." );

    // Reference area (use fallback area if reference area not provided, final value cannont be NaN)
    double fallbackReferenceArea = getNumeric< double >(
                jsonObject, SpecialKeys::up / SpecialKeys::up / Keys::Body::referenceArea, TUDAT_NAN, true );

    switch ( radiationPressureType )
    {
    case cannon_ball:
    {
        CannonBallRadiationPressureInterfaceSettings defaults( "", TUDAT_NAN, TUDAT_NAN );
        radiationPressureInterfaceSettings = boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                    sourceBodyName,
                    getNumeric( jsonObject, K::referenceArea, fallbackReferenceArea ),
                    getValue< double >( jsonObject, K::radiationPressureCoefficient ),
                    getValue( jsonObject, K::ocultingBodies, defaults.getOccultingBodies( ) ) );
        return;
    }
    default:
        handleUnimplementedEnumValueFromJson( radiationPressureType, radiationPressureTypes,
                                              unsupportedRadiationPressureTypes );
    }
}

} // namespace simulation_setup

} // namespace tudat
