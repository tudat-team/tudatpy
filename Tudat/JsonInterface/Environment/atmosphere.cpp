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

#include "Tudat/JsonInterface/Environment/atmosphere.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `AtmosphereSettings` object.
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< AtmosphereSettings >& atmosphereSettings )
{
    if ( ! atmosphereSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Body::Atmosphere;

    const AtmosphereTypes atmosphereType = atmosphereSettings->getAtmosphereType( );
    jsonObject[ K::type ] = atmosphereType;

    // ExponentialAtmosphereSettings
    switch ( atmosphereType )
    {
    case exponential_atmosphere:
    {
        boost::shared_ptr< ExponentialAtmosphereSettings > exponentialAtmosphereSettings =
                boost::dynamic_pointer_cast< ExponentialAtmosphereSettings >( atmosphereSettings );
        assertNonNullPointer( exponentialAtmosphereSettings );
        jsonObject[ K::densityScaleHeight ] = exponentialAtmosphereSettings->getDensityScaleHeight( );
        jsonObject[ K::constantTemperature ] = exponentialAtmosphereSettings->getConstantTemperature( );
        jsonObject[ K::densityAtZeroAltitude ] = exponentialAtmosphereSettings->getDensityAtZeroAltitude( );
        jsonObject[ K::specificGasConstant ] = exponentialAtmosphereSettings->getSpecificGasConstant( );
        return;
    }
    case tabulated_atmosphere:
    {
        boost::shared_ptr< TabulatedAtmosphereSettings > tabulatedAtmosphereSettings =
                boost::dynamic_pointer_cast< TabulatedAtmosphereSettings >( atmosphereSettings );
        assertNonNullPointer( tabulatedAtmosphereSettings );
        jsonObject[ K::file ] = boost::filesystem::path( tabulatedAtmosphereSettings->getAtmosphereFile( ) );
        return;
    }
    case nrlmsise00:
    {
        boost::shared_ptr< NRLMSISE00AtmosphereSettings > nrlmsise00AtmosphereSettings =
                boost::dynamic_pointer_cast< NRLMSISE00AtmosphereSettings >( atmosphereSettings );
        if ( nrlmsise00AtmosphereSettings )
        {
            jsonObject[ K::spaceWeatherFile ] = boost::filesystem::path( nrlmsise00AtmosphereSettings->getSpaceWeatherFile( ) );
            return;
        }
        // If not a NRLMSISE00AtmosphereSettings, it is a AtmosphereSettings with default space weather file.
        return;
    }
    default:
        handleUnimplementedEnumValue( atmosphereType, atmosphereTypes, unsupportedAtmosphereTypes );
    }
}

//! Create a shared pointer to a `AtmosphereSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, boost::shared_ptr< AtmosphereSettings >& atmosphereSettings )
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
                    getValue< boost::filesystem::path >( jsonObject, K::file ).string( ) );
        return;
    }
    case nrlmsise00:
    {
        if ( isDefined( jsonObject, K::spaceWeatherFile ) )
        {
            atmosphereSettings = boost::make_shared< NRLMSISE00AtmosphereSettings >(
                        getValue< boost::filesystem::path >( jsonObject, K::spaceWeatherFile ).string( ) );
        }
        else
        {
            atmosphereSettings = boost::make_shared< AtmosphereSettings >( nrlmsise00 );
        }
        return;
    }
    default:
        handleUnimplementedEnumValue( atmosphereType, atmosphereTypes, unsupportedAtmosphereTypes );
    }
}

} // namespace simulation_setup

} // namespace tudat
