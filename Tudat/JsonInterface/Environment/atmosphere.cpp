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

#include "Tudat/JsonInterface/Mathematics/interpolation.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `AtmosphereSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< AtmosphereSettings >& atmosphereSettings )
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
        std::shared_ptr< ExponentialAtmosphereSettings > exponentialAtmosphereSettings =
                std::dynamic_pointer_cast< ExponentialAtmosphereSettings >( atmosphereSettings );
        assertNonnullptrPointer( exponentialAtmosphereSettings );
        jsonObject[ K::densityScaleHeight ] = exponentialAtmosphereSettings->getDensityScaleHeight( );
        jsonObject[ K::constantTemperature ] = exponentialAtmosphereSettings->getConstantTemperature( );
        jsonObject[ K::densityAtZeroAltitude ] = exponentialAtmosphereSettings->getDensityAtZeroAltitude( );
        jsonObject[ K::specificGasConstant ] = exponentialAtmosphereSettings->getSpecificGasConstant( );
        jsonObject[ K::ratioOfSpecificHeats ] = exponentialAtmosphereSettings->getRatioOfSpecificHeats( );
        return;
    }
    case tabulated_atmosphere:
    {
        std::shared_ptr< TabulatedAtmosphereSettings > tabulatedAtmosphereSettings =
                std::dynamic_pointer_cast< TabulatedAtmosphereSettings >( atmosphereSettings );
        assertNonnullptrPointer( tabulatedAtmosphereSettings );
        jsonObject[ K::file ] = tabulatedAtmosphereSettings->getAtmosphereFile( );
        jsonObject[ K::independentVariablesNames ] = tabulatedAtmosphereSettings->getIndependentVariables( );
        jsonObject[ K::dependentVariablesNames ] = tabulatedAtmosphereSettings->getDependentVariables( );
        jsonObject[ K::specificGasConstant ] = tabulatedAtmosphereSettings->getSpecificGasConstant( );
        jsonObject[ K::ratioOfSpecificHeats ] = tabulatedAtmosphereSettings->getRatioOfSpecificHeats( );
        jsonObject[ K::boundaryHandling ] = tabulatedAtmosphereSettings->getBoundaryHandling( );
        return;
    }
    case nrlmsise00:
    {
        std::shared_ptr< NRLMSISE00AtmosphereSettings > nrlmsise00AtmosphereSettings =
                std::dynamic_pointer_cast< NRLMSISE00AtmosphereSettings >( atmosphereSettings );
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
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< AtmosphereSettings >& atmosphereSettings )
{
    using namespace json_interface;
    using K = Keys::Body::Atmosphere;

    // Get atmosphere model type
    const AtmosphereTypes atmosphereType = getValue< AtmosphereTypes >( jsonObject, K::type );

    switch ( atmosphereType ) {
    case exponential_atmosphere:
    {
        ExponentialAtmosphereSettings defaults( 0.0, 0.0, 0.0 );
        atmosphereSettings = std::make_shared< ExponentialAtmosphereSettings >(
                    getValue< double >( jsonObject, K::densityScaleHeight ),
                    getValue< double >( jsonObject, K::constantTemperature ),
                    getValue< double >( jsonObject, K::densityAtZeroAltitude ),
                    getValue< double >( jsonObject, K::specificGasConstant, defaults.getSpecificGasConstant( ) ),
                    getValue< double >( jsonObject, K::ratioOfSpecificHeats, defaults.getRatioOfSpecificHeats( ) ) );
        return;
    }
    case tabulated_atmosphere:
    {
        TabulatedAtmosphereSettings defaults( getValue< std::map< int, std::string > >( jsonObject, K::file ) );
        atmosphereSettings = std::make_shared< TabulatedAtmosphereSettings >(
                    getValue< std::map< int, std::string > >( jsonObject, K::file ),
                    getValue< std::vector< AtmosphereIndependentVariables > >( jsonObject, K::independentVariablesNames,
                                                                               defaults.getIndependentVariables( ) ),
                    getValue< std::vector< AtmosphereDependentVariables > >( jsonObject, K::dependentVariablesNames,
                                                                             defaults.getDependentVariables( ) ),
                    getValue< double >( jsonObject, K::specificGasConstant, defaults.getSpecificGasConstant( ) ),
                    getValue< double >( jsonObject, K::ratioOfSpecificHeats, defaults.getRatioOfSpecificHeats( ) ),
                    getValue< std::vector< interpolators::BoundaryInterpolationType > >( jsonObject, K::boundaryHandling,
                                                                                         defaults.getBoundaryHandling( ) ) );
        return;
    }
    case nrlmsise00:
    {
        if ( isDefined( jsonObject, K::spaceWeatherFile ) )
        {
            atmosphereSettings = std::make_shared< NRLMSISE00AtmosphereSettings >(
                        getValue< boost::filesystem::path >( jsonObject, K::spaceWeatherFile ).string( ) );
        }
        else
        {
            atmosphereSettings = std::make_shared< AtmosphereSettings >( nrlmsise00 );
        }
        return;
    }
    default:
        handleUnimplementedEnumValue( atmosphereType, atmosphereTypes, unsupportedAtmosphereTypes );
    }
}

} // namespace simulation_setup

} // namespace tudat
