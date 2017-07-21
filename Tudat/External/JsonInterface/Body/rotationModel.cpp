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
    jsonObject = json( json_interface::stringFromEnum( rotationModelType, rotationModelTypes ) );
}

//! Convert `json` to `RotationModelType`.
void from_json( const json& jsonObject, RotationModelType& rotationModelType )
{
    rotationModelType = json_interface::enumFromString( jsonObject.get< std::string >( ), rotationModelTypes );
}

//! Create a `json` object from a shared pointer to a `RotationModelSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< RotationModelSettings >& rotationModelSettings )
{
    using namespace json_interface;
    using Keys = Keys::Body::RotationModel;

    // Initialise
    jsonObject = json( );

    // Base class settings
    jsonObject[ Keys::type ] = rotationModelSettings->getRotationType( );
    jsonObject[ Keys::originalFrame ] = rotationModelSettings->getOriginalFrame( );
    jsonObject[ Keys::targetFrame ] = rotationModelSettings->getTargetFrame( );

    /// spice_rotation_model
    if ( jsonObject[ Keys::type ] == spice_rotation_model )
    {
        return;
    }

    /// SimpleRotationModelSettings
    boost::shared_ptr< SimpleRotationModelSettings > simpleRotationModelSettings =
        boost::dynamic_pointer_cast< SimpleRotationModelSettings >( rotationModelSettings );
    if ( simpleRotationModelSettings )
    {
    jsonObject[ Keys::initialOrientation ] = simpleRotationModelSettings->getInitialOrientation( );
    jsonObject[ Keys::initialTime ] = simpleRotationModelSettings->getInitialTime( );
    jsonObject[ Keys::rotationRate ] = simpleRotationModelSettings->getRotationRate( );
    return;
    }
}

} // namespace simulation_setup


namespace json_interface
{

//! Create a shared pointer to a `RotationModelSettings` object from a `json` object.
boost::shared_ptr< simulation_setup::RotationModelSettings > createRotationModelSettings(
        const json& settings, const KeyTree& keyTree )
{
    using namespace simulation_setup;
    using Keys = Keys::Body::RotationModel;

    // Base class settings
    const RotationModelType rotationModelType = getValue< RotationModelType >( settings, keyTree + Keys::type );
    const std::string originalFrame = getValue< std::string >( settings, keyTree + Keys::originalFrame );
    const std::string targetFrame = getValue< std::string >( settings, keyTree + Keys::targetFrame );

    switch ( rotationModelType ) {
    case simple_rotation_model:
    return boost::make_shared< SimpleRotationModelSettings >(
            originalFrame, targetFrame,
            getValue< Eigen::Quaterniond >( settings, keyTree + Keys::initialOrientation ),
            getEpoch< double >( settings, keyTree + Keys::initialTime ),
            getNumeric< double >( settings, keyTree + Keys::rotationRate ) );
    case spice_rotation_model:
    return boost::make_shared< RotationModelSettings >( rotationModelType, originalFrame, targetFrame );
    default:
    throw std::runtime_error( stringFromEnum( rotationModelType, rotationModelTypes )
                                  + " not supported by json_interface." );
    }
}

} // namespace json_interface

} // namespace tudat
