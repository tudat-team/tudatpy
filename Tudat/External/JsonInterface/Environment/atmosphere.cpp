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
    jsonObject = json_interface::stringFromEnum( atmosphereType, atmosphereTypes );
}

//! Convert `json` to `AtmosphereTypes`.
void from_json( const json& jsonObject, AtmosphereTypes& atmosphereType )
{
    atmosphereType = json_interface::enumFromString( jsonObject.get< std::string >( ), atmosphereTypes );
}

//! Create a `json` object from a shared pointer to a `AtmosphereSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< AtmosphereSettings >& atmosphereSettings )
{
    if ( atmosphereSettings )
    {
        using namespace json_interface;
        using K = Keys::Body::Atmosphere;

        // Type
        jsonObject[ K::type ] = atmosphereSettings->getAtmosphereType( );

        /// ExponentialAtmosphereSettings
        boost::shared_ptr< ExponentialAtmosphereSettings > exponentialAtmosphereSettings =
                boost::dynamic_pointer_cast< ExponentialAtmosphereSettings >( atmosphereSettings );
        if ( exponentialAtmosphereSettings )
        {
            jsonObject[ K::densityScaleHeight ] = exponentialAtmosphereSettings->getDensityScaleHeight( );
            jsonObject[ K::constantTemperature ] = exponentialAtmosphereSettings->getConstantTemperature( );
            jsonObject[ K::densityAtZeroAltitude ] = exponentialAtmosphereSettings->getDensityAtZeroAltitude( );
            jsonObject[ K::specificGasConstant ] = exponentialAtmosphereSettings->getSpecificGasConstant( );
            return;
        }

        /// TabulatedAtmosphereSettings
        boost::shared_ptr< TabulatedAtmosphereSettings > tabulatedAtmosphereSettings =
                boost::dynamic_pointer_cast< TabulatedAtmosphereSettings >( atmosphereSettings );
        if ( tabulatedAtmosphereSettings )
        {
            jsonObject[ K::file ] = path( tabulatedAtmosphereSettings->getAtmosphereFile( ) );
            return;
        }

        /// NRLMSISE00AtmosphereSettings
        boost::shared_ptr< NRLMSISE00AtmosphereSettings > nrlmsise00AtmosphereSettings =
                boost::dynamic_pointer_cast< NRLMSISE00AtmosphereSettings >( atmosphereSettings );
        if ( nrlmsise00AtmosphereSettings )
        {
            jsonObject[ K::spaceWeatherFile ] = path( nrlmsise00AtmosphereSettings->getSpaceWeatherFile( ) );
            return;
        }
    }
}

//! Create a shared pointer to a `AtmosphereSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< AtmosphereSettings >& atmosphereSettings )
{
    using namespace json_interface;
    using K = Keys::Body::Atmosphere;

    // Get atmosphere model type
    const AtmosphereTypes atmosphereType = getValue< AtmosphereTypes >( jsonObject, K::type );

    switch ( atmosphereType ) {
    case exponential_atmosphere:
    {
        atmosphereSettings = boost::make_shared< ExponentialAtmosphereSettings >(
                    getValue< double >( jsonObject, K::densityScaleHeight ),
                    getValue< double >( jsonObject, K::constantTemperature ),
                    getValue< double >( jsonObject, K::densityAtZeroAltitude ),
                    getValue< double >( jsonObject, K::specificGasConstant ) );
        return;
    }
    case tabulated_atmosphere:
    {
        atmosphereSettings = boost::make_shared< TabulatedAtmosphereSettings >(
                    getValue< path >( jsonObject, K::file ).string( ) );
        return;
    }
    case nrlmsise00:
    {
        const boost::shared_ptr< path > spaceWeatherFile =
                getValuePointer< path >( jsonObject, K::spaceWeatherFile );
        if ( spaceWeatherFile )
        {
            atmosphereSettings = boost::make_shared< NRLMSISE00AtmosphereSettings >( spaceWeatherFile->string( ) );
        }
        else
        {
            atmosphereSettings = boost::make_shared< AtmosphereSettings >( nrlmsise00 );
        }
        return;
    }
    default:
        throw std::runtime_error( stringFromEnum( atmosphereType, atmosphereTypes )
                                  + " not supported by json_interface." );
    }
}

} // namespace simulation_setup

} // namespace tudat
