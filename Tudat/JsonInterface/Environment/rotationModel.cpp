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

#include "Tudat/JsonInterface/Environment/rotationModel.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `RotationModelSettings` object.
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< RotationModelSettings >& rotationModelSettings )
{
    if ( ! rotationModelSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Body::RotationModel;

    const RotationModelType rotationModelType = rotationModelSettings->getRotationType( );
    jsonObject[ K::type ] = rotationModelType;
    jsonObject[ K::originalFrame ] = rotationModelSettings->getOriginalFrame( );
    jsonObject[ K::targetFrame ] = rotationModelSettings->getTargetFrame( );

    switch ( rotationModelType )
    {
    case simple_rotation_model:
    {
        boost::shared_ptr< SimpleRotationModelSettings > simpleRotationModelSettings =
                boost::dynamic_pointer_cast< SimpleRotationModelSettings >( rotationModelSettings );
        assertNonNullPointer( simpleRotationModelSettings );
        jsonObject[ K::initialOrientation ] = simpleRotationModelSettings->getInitialOrientation( );
        jsonObject[ K::initialTime ] = simpleRotationModelSettings->getInitialTime( );
        jsonObject[ K::rotationRate ] = simpleRotationModelSettings->getRotationRate( );
        return;
    }
    case spice_rotation_model:
        return;
    default:
        handleUnimplementedEnumValue( rotationModelType, rotationModelTypes, unsupportedRotationModelTypes );
    }
}

//! Create a shared pointer to a `RotationModelSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, boost::shared_ptr< RotationModelSettings >& rotationModelSettings )
{
    using namespace json_interface;
    using K = Keys::Body::RotationModel;

    // Base class settings
    const RotationModelType rotationModelType = getValue< RotationModelType >( jsonObject, K::type );
    const std::string originalFrame = getValue< std::string >( jsonObject, K::originalFrame );
    const std::string targetFrame = getValue< std::string >( jsonObject, K::targetFrame );

    switch ( rotationModelType ) {
    case simple_rotation_model:
    {
        const double initialTime = getValue< double >( jsonObject, K::initialTime );

        // Get JSON object for initialOrientation (or create it if not defined)
        nlohmann::json jsonInitialOrientation;
        if ( isDefined( jsonObject, K::initialOrientation ) )
        {
            jsonInitialOrientation = getValue< nlohmann::json >( jsonObject, K::initialOrientation );
        }
        else
        {
            jsonInitialOrientation[ K::originalFrame ] = originalFrame;
            jsonInitialOrientation[ K::targetFrame ] = targetFrame;
            jsonInitialOrientation[ K::initialTime ] = initialTime;
        }

        rotationModelSettings = boost::make_shared< SimpleRotationModelSettings >(
                    originalFrame,
                    targetFrame,
                    getAs< Eigen::Quaterniond >( jsonInitialOrientation ),
                    initialTime,
                    getValue< double >( jsonObject, K::rotationRate ) );
        return;
    }
    case spice_rotation_model:
    {
        rotationModelSettings = boost::make_shared< RotationModelSettings >(
                    rotationModelType, originalFrame, targetFrame );
        return;
    }
    default:
        handleUnimplementedEnumValue( rotationModelType, rotationModelTypes, unsupportedRotationModelTypes );
    }
}

} // namespace simulation_setup

} // namespace tudat
