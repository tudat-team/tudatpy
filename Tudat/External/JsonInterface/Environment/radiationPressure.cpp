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
        jsonObject[ K::occultingBodies ] =
                cannonBallRadiationPressureInterfaceSettings->getOccultingBodies( );
        return;
    }
    default:
        handleUnimplementedEnumValue( radiationPressureType, radiationPressureTypes,
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

    // Reference area (either from the current object or from the current object's parent's parent, i.e. the body)
    const double referenceArea = getValue< double >(
                jsonObject, { K::referenceArea, SpecialKeys::up / SpecialKeys::up / Keys::Body::referenceArea } );

    switch ( radiationPressureType )
    {
    case cannon_ball:
    {
        std::string sourceBody;
        if ( ! defined( jsonObject, K::sourceBody ) )
        {
            try
            {
                sourceBody = getParentKey( jsonObject );
            }
            catch ( ... ) { }
        }
        if ( sourceBody.empty( ) )
        {
            sourceBody = getValue< std::string >( jsonObject, K::sourceBody );
        }

        CannonBallRadiationPressureInterfaceSettings defaults( "", TUDAT_NAN, TUDAT_NAN );
        radiationPressureInterfaceSettings = boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                    sourceBody,
                    referenceArea,
                    getValue< double >( jsonObject, K::radiationPressureCoefficient ),
                    getValue( jsonObject, K::occultingBodies, defaults.getOccultingBodies( ) ) );
        return;
    }
    default:
        handleUnimplementedEnumValue( radiationPressureType, radiationPressureTypes,
                                      unsupportedRadiationPressureTypes );
    }
}

} // namespace simulation_setup

} // namespace tudat
