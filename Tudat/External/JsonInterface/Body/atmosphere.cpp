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

#include <regex>

#include <boost/filesystem.hpp>

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
    }

    /// TabulatedAtmosphereSettings
    boost::shared_ptr< TabulatedAtmosphereSettings > tabulatedAtmosphereSettings =
            boost::dynamic_pointer_cast< TabulatedAtmosphereSettings >( atmosphereSettings );
    if ( tabulatedAtmosphereSettings )
    {
        jsonObject[ Keys::atmosphereFile ] =
                boost::filesystem::canonical( tabulatedAtmosphereSettings->getAtmosphereFile( ) ).string( );
    }

    /// NRLMSISE00AtmosphereSettings
    boost::shared_ptr< NRLMSISE00AtmosphereSettings > nrlmsise00AtmosphereSettings =
            boost::dynamic_pointer_cast< NRLMSISE00AtmosphereSettings >( atmosphereSettings );
    if ( nrlmsise00AtmosphereSettings )
    {
        jsonObject[ Keys::spaceWeatherFile ] =
                boost::filesystem::canonical( nrlmsise00AtmosphereSettings->getSpaceWeatherFile( ) ).string( );
    }
}

} // namespace simulation_setup


namespace json_interface
{

//! Create a shared pointer to a `AtmosphereSettings` object from a `json` object.
boost::shared_ptr< simulation_setup::AtmosphereSettings > createAtmosphereSettings( const json &settings )
{
    using namespace simulation_setup;
    using Keys = Keys::Body::Atmosphere;

    // Get atmosphere model type
    const AtmosphereTypes atmosphereType = getValue< AtmosphereTypes >( settings, Keys::type );

    if ( atmosphereType == exponential_atmosphere )
    {
        return boost::make_shared< ExponentialAtmosphereSettings >(
                    getValue< double >( settings, Keys::densityScaleHeight ),
                    getValue< double >( settings, Keys::constantTemperature ),
                    getValue< double >( settings, Keys::densityAtZeroAltitude ),
                    getValue< double >( settings, Keys::specificGasConstant ) );
    }
    else if ( atmosphereType == tabulated_atmosphere )
    {
        return boost::make_shared< TabulatedAtmosphereSettings >(
                    getValue< std::string >( settings, Keys::atmosphereFile ) );
    }
    else if ( atmosphereType == nrlmsise00 )
    {
        const auto spaceWeatherFile = getValuePointer< std::string >( settings, Keys::spaceWeatherFile );
        if ( spaceWeatherFile )
        {
            return boost::make_shared< NRLMSISE00AtmosphereSettings >( *spaceWeatherFile );
        }
        else
        {
            return boost::make_shared< AtmosphereSettings >( nrlmsise00 );
        }
    }
    else
    {
        throw std::runtime_error( stringFromEnum( atmosphereType, atmosphereTypes )
                                  + " not supported by json_interface." );
    }
}

} // namespace json_interface

} // namespace tudat
