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

#include "rotationModel.h"

namespace tudat
{

namespace simulation_setup
{

//! Convert `RotationModelType`s to `json`.
void to_json( json& jsonObject, const RotationModelType& rotationModelType )
{
    jsonObject = json_interface::stringFromEnum( rotationModelType, rotationModelTypes );
}

//! Convert `json` to `RotationModelType`.
void from_json( const json& jsonObject, RotationModelType& rotationModelType )
{
    rotationModelType = json_interface::enumFromString( jsonObject.get< std::string >( ), rotationModelTypes );
}

//! Create a `json` object from a shared pointer to a `RotationModelSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< RotationModelSettings >& rotationModelSettings )
{
    if ( rotationModelSettings )
    {
        using namespace json_interface;
        using K = Keys::Body::RotationModel;

        // Base class settings
        jsonObject[ K::type ] = rotationModelSettings->getRotationType( );
        jsonObject[ K::originalFrame ] = rotationModelSettings->getOriginalFrame( );
        jsonObject[ K::targetFrame ] = rotationModelSettings->getTargetFrame( );

        /// spice_rotation_model
        if ( jsonObject[ K::type ] == spice_rotation_model )
        {
            return;
        }

        /// SimpleRotationModelSettings
        boost::shared_ptr< SimpleRotationModelSettings > simpleRotationModelSettings =
                boost::dynamic_pointer_cast< SimpleRotationModelSettings >( rotationModelSettings );
        if ( simpleRotationModelSettings )
        {
            jsonObject[ K::initialOrientation ] = simpleRotationModelSettings->getInitialOrientation( );
            jsonObject[ K::initialTime ] = simpleRotationModelSettings->getInitialTime( );
            jsonObject[ K::rotationRate ] = simpleRotationModelSettings->getRotationRate( );
            return;
        }
    }
}

//! Create a shared pointer to a `RotationModelSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< RotationModelSettings >& rotationModelSettings )
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
        rotationModelSettings = boost::make_shared< SimpleRotationModelSettings >(
                    originalFrame,
                    targetFrame,
                    getValue< Eigen::Quaterniond >( jsonObject, K::initialOrientation ),
                    getEpoch< double >( jsonObject, K::initialTime ),
                    getNumeric< double >( jsonObject, K::rotationRate ) );
        return;
    }
    case spice_rotation_model:
    {
        rotationModelSettings = boost::make_shared< RotationModelSettings >(
                    rotationModelType, originalFrame, targetFrame );
        return;
    }
    default:
        throw std::runtime_error( stringFromEnum( rotationModelType, rotationModelTypes )
                                  + " not supported by json_interface." );
    }
}

} // namespace simulation_setup

} // namespace tudat
