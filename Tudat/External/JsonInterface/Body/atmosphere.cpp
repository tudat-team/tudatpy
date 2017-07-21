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

#include "atmosphere.h"

namespace tudat
{

namespace simulation_setup
{

//! Convert `AtmosphereTypes` to `json`.
void to_json( json& jsonObject, const AtmosphereTypes& atmosphereType )
{
    jsonObject = json( json_interface::stringFromEnum( atmosphereType, atmosphereTypes ) );
}

//! Convert `json` to `AtmosphereTypes`.
void from_json( const json& jsonObject, AtmosphereTypes& atmosphereType )
{
    atmosphereType = json_interface::enumFromString( jsonObject.get< std::string >( ), atmosphereTypes );
}

//! Create a `json` object from a shared pointer to a `AtmosphereSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< AtmosphereSettings >& atmosphereSettings )
{
    using namespace json_interface;
    using Keys = Keys::Body::Atmosphere;

    // Initialise
    jsonObject = json( );

    // Type
    jsonObject[ Keys::type ] = atmosphereSettings->getAtmosphereType( );

    /// ExponentialAtmosphereSettings
    boost::shared_ptr< ExponentialAtmosphereSettings > exponentialAtmosphereSettings =
            boost::dynamic_pointer_cast< ExponentialAtmosphereSettings >( atmosphereSettings );
    if ( exponentialAtmosphereSettings )
    {
        jsonObject[ Keys::densityScaleHeight ] = exponentialAtmosphereSettings->getDensityScaleHeight( );
        jsonObject[ Keys::constantTemperature ] = exponentialAtmosphereSettings->getConstantTemperature( );
        jsonObject[ Keys::densityAtZeroAltitude ] = exponentialAtmosphereSettings->getDensityAtZeroAltitude( );
        jsonObject[ Keys::specificGasConstant ] = exponentialAtmosphereSettings->getSpecificGasConstant( );
        return;
    }

    /// TabulatedAtmosphereSettings
    boost::shared_ptr< TabulatedAtmosphereSettings > tabulatedAtmosphereSettings =
            boost::dynamic_pointer_cast< TabulatedAtmosphereSettings >( atmosphereSettings );
    if ( tabulatedAtmosphereSettings )
    {
        jsonObject[ Keys::file ] = path( tabulatedAtmosphereSettings->getAtmosphereFile( ) );
        return;
    }

    /// NRLMSISE00AtmosphereSettings
    boost::shared_ptr< NRLMSISE00AtmosphereSettings > nrlmsise00AtmosphereSettings =
            boost::dynamic_pointer_cast< NRLMSISE00AtmosphereSettings >( atmosphereSettings );
    if ( nrlmsise00AtmosphereSettings )
    {
        jsonObject[ Keys::spaceWeatherFile ] = path( nrlmsise00AtmosphereSettings->getSpaceWeatherFile( ) );
        return;
    }
}

} // namespace simulation_setup


namespace json_interface
{

//! Create a shared pointer to a `AtmosphereSettings` object from a `json` object.
boost::shared_ptr< simulation_setup::AtmosphereSettings > createAtmosphereSettings(
        const json& settings, const KeyTree& keyTree )
{
    using namespace simulation_setup;
    using Keys = Keys::Body::Atmosphere;

    // Get atmosphere model type
    const AtmosphereTypes atmosphereType = getValue< AtmosphereTypes >( settings, keyTree + Keys::type );

    switch ( atmosphereType ) {
    case exponential_atmosphere:
        return boost::make_shared< ExponentialAtmosphereSettings >(
                    getValue< double >( settings, keyTree + Keys::densityScaleHeight ),
                    getValue< double >( settings, keyTree + Keys::constantTemperature ),
                    getValue< double >( settings, keyTree + Keys::densityAtZeroAltitude ),
                    getValue< double >( settings, keyTree + Keys::specificGasConstant ) );
    case tabulated_atmosphere:
        return boost::make_shared< TabulatedAtmosphereSettings >(
                    getValue< path >( settings, keyTree + Keys::file ).string( ) );
    case nrlmsise00:
    {
        const boost::shared_ptr< path > spaceWeatherFile =
                getValuePointer< path >( settings, keyTree + Keys::spaceWeatherFile );
        if ( spaceWeatherFile )
        {
            return boost::make_shared< NRLMSISE00AtmosphereSettings >( spaceWeatherFile->string( ) );
        }
        else
        {
            return boost::make_shared< AtmosphereSettings >( nrlmsise00 );
        }
    }
    default:
        throw std::runtime_error( stringFromEnum( atmosphereType, atmosphereTypes )
                                  + " not supported by json_interface." );
    }
}

} // namespace json_interface

} // namespace tudat
