/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>
#include <boost/bind/bind.hpp>

#include "tudat/astro/aerodynamics/exponentialAtmosphere.h"
#include "tudat/astro/aerodynamics/tabulatedAtmosphere.h"
#if TUDAT_BUILD_WITH_NRLMSISE
#include "tudat/astro/aerodynamics/nrlmsise00Atmosphere.h"
#include "tudat/astro/aerodynamics/nrlmsise00InputFunctions.h"
#endif
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/solarActivityData.h"
#include "tudat/simulation/environment_setup/createAtmosphereModel.h"

using namespace boost::placeholders;

namespace tudat
{

namespace simulation_setup
{

//! Function to create a wind model.
std::shared_ptr< aerodynamics::WindModel > createWindModel(
        const std::shared_ptr< WindModelSettings > windSettings,
        const std::string& body )
{
    std::shared_ptr< aerodynamics::WindModel > windModel;

    // Check wind model type and create requested model
    switch( windSettings->getWindModelType( ) )
    {
    case constant_wind_model:
    {
        // Check input consistency
        std::shared_ptr< ConstantWindModelSettings > customWindModelSettings =
                std::dynamic_pointer_cast< ConstantWindModelSettings >( windSettings );

        if( customWindModelSettings == nullptr )
        {
            throw std::runtime_error( "Error when making constant wind model for body " + body + ", input is incompatible" );
        }
        windModel = std::make_shared< aerodynamics::ConstantWindModel >(
                    customWindModelSettings->getConstantWindVelocity( ),
                    customWindModelSettings->getAssociatedFrame( ) );
        break;
    }
    case custom_wind_model:
    {
        // Check input consistency
        std::shared_ptr< CustomWindModelSettings > customWindModelSettings =
                std::dynamic_pointer_cast< CustomWindModelSettings >( windSettings );

        if( customWindModelSettings == nullptr )
        {
            throw std::runtime_error( "Error when making custom wind model for body " + body + ", input is incompatible" );
        }
        windModel = std::make_shared< aerodynamics::CustomWindModel >(
                    customWindModelSettings->getWindFunction( ),
                    customWindModelSettings->getAssociatedFrame( ) );

        break;
    }
    default:
        throw std::runtime_error( "Error when making wind model for body " + body + ", input type not recognized" );
    }

    return windModel;

}

//! Function to create an atmosphere model.
std::shared_ptr< aerodynamics::AtmosphereModel > createAtmosphereModel(
        const std::shared_ptr< AtmosphereSettings > atmosphereSettings,
        const std::string& body )
{
    using namespace tudat::aerodynamics;

    // Declare return object.
    std::shared_ptr< AtmosphereModel > atmosphereModel;

    // Check which type of atmosphere model is to be created.
    switch( atmosphereSettings->getAtmosphereType( ) )
    {
    case exponential_atmosphere:
    {
        // Check whether settings for atmosphere are consistent with its type.
        std::shared_ptr< ExponentialAtmosphereSettings > exponentialAtmosphereSettings =
                std::dynamic_pointer_cast< ExponentialAtmosphereSettings >( atmosphereSettings );

        if( exponentialAtmosphereSettings == nullptr )
        {
            throw std::runtime_error( "Error, expected exponential atmosphere settings for body " + body );
        }
        else
        {
            // Create and initialize exponential atmosphere model.
            std::shared_ptr< ExponentialAtmosphere > exponentialAtmosphereModel;
            if ( exponentialAtmosphereSettings->getBodyName( ) == undefined_body )
            {
                exponentialAtmosphereModel = std::make_shared< ExponentialAtmosphere >(
                            exponentialAtmosphereSettings->getDensityScaleHeight( ) ,
                            exponentialAtmosphereSettings->getConstantTemperature( ),
                            exponentialAtmosphereSettings->getDensityAtZeroAltitude( ),
                            exponentialAtmosphereSettings->getSpecificGasConstant( ),
                            exponentialAtmosphereSettings->getRatioOfSpecificHeats( ) );
            }
            else
            {
                exponentialAtmosphereModel = std::make_shared< ExponentialAtmosphere >(
                            exponentialAtmosphereSettings->getBodyName( ) );
            }
            atmosphereModel = exponentialAtmosphereModel;
        }
        break;
    }
    case custom_constant_temperature_atmosphere:
    {
        // Check whether settings for atmosphere are consistent with its type.
        std::shared_ptr< CustomConstantTemperatureAtmosphereSettings > customConstantTemperatureAtmosphereSettings =
                std::dynamic_pointer_cast< CustomConstantTemperatureAtmosphereSettings >( atmosphereSettings );
        if( customConstantTemperatureAtmosphereSettings == nullptr )
        {
            throw std::runtime_error( "Error, expected exponential atmosphere settings for body " + body );
        }
        else
        {
            // Create and initialize exponential atmosphere model.
            std::shared_ptr< CustomConstantTemperatureAtmosphere > customConstantTemperatureAtmosphereModel;
            if ( customConstantTemperatureAtmosphereSettings->getModelSpecificParameters( ).empty( ) )
            {
                customConstantTemperatureAtmosphereModel = std::make_shared< CustomConstantTemperatureAtmosphere >(
                            customConstantTemperatureAtmosphereSettings->getDensityFunction( ) ,
                            customConstantTemperatureAtmosphereSettings->getConstantTemperature( ),
                            customConstantTemperatureAtmosphereSettings->getSpecificGasConstant( ),
                            customConstantTemperatureAtmosphereSettings->getRatioOfSpecificHeats( ) );
            }
            else
            {
                customConstantTemperatureAtmosphereModel = std::make_shared< CustomConstantTemperatureAtmosphere >(
                            customConstantTemperatureAtmosphereSettings->getDensityFunctionType( ),
                            customConstantTemperatureAtmosphereSettings->getConstantTemperature( ),
                            customConstantTemperatureAtmosphereSettings->getSpecificGasConstant( ),
                            customConstantTemperatureAtmosphereSettings->getRatioOfSpecificHeats( ),
                            customConstantTemperatureAtmosphereSettings->getModelSpecificParameters( ) );
            }
            atmosphereModel = customConstantTemperatureAtmosphereModel;
        }
        break;
    }
    case tabulated_atmosphere:
    {
        // Check whether settings for atmosphere are consistent with its type
        std::shared_ptr< TabulatedAtmosphereSettings > tabulatedAtmosphereSettings =
                std::dynamic_pointer_cast< TabulatedAtmosphereSettings >( atmosphereSettings );
        if( tabulatedAtmosphereSettings == nullptr )
        {
            throw std::runtime_error( "Error, expected tabulated atmosphere settings for body " + body );
        }
        else
        {
            // Create and initialize tabulated atmosphere model.
            atmosphereModel = std::make_shared< TabulatedAtmosphere >(
                        tabulatedAtmosphereSettings->getAtmosphereFile( ),
                        tabulatedAtmosphereSettings->getIndependentVariables( ),
                        tabulatedAtmosphereSettings->getDependentVariables( ),
                        tabulatedAtmosphereSettings->getSpecificGasConstant( ),
                        tabulatedAtmosphereSettings->getRatioOfSpecificHeats( ),
                        tabulatedAtmosphereSettings->getBoundaryHandling( ),
                        tabulatedAtmosphereSettings->getDefaultExtrapolationValue( ) );
        }
        break;
    }
#if TUDAT_BUILD_WITH_NRLMSISE
    case nrlmsise00:
    {
        std::string spaceWeatherFilePath;
        std::shared_ptr< NRLMSISE00AtmosphereSettings > nrlmsise00AtmosphereSettings =
                std::dynamic_pointer_cast< NRLMSISE00AtmosphereSettings >( atmosphereSettings );

        if( nrlmsise00AtmosphereSettings == nullptr )
        {
            // Use default space weather file stored in tudatBundle.
            spaceWeatherFilePath = paths::getSpaceWeatherDataPath( ) + "/sw19571001.txt";
        }
        else
        {
            // Use space weather file specified by user.
            spaceWeatherFilePath = nrlmsise00AtmosphereSettings->getSpaceWeatherFile( );
        }

        std::cout<<"Reading space weather file "<<spaceWeatherFilePath<<std::endl;
        tudat::input_output::solar_activity::SolarActivityDataMap solarActivityData =
                tudat::input_output::solar_activity::readSolarActivityData( spaceWeatherFilePath ) ;

        // Create atmosphere model using NRLMISE00 input function
        std::function< tudat::aerodynamics::NRLMSISE00Input( double, double, double, double ) > inputFunction =
                std::bind( &tudat::aerodynamics::nrlmsiseInputFunction,
                           std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4,
                           solarActivityData, false, TUDAT_NAN );
        atmosphereModel = std::make_shared< aerodynamics::NRLMSISE00Atmosphere >( inputFunction );
        break;
    }
#endif
    case scaled_atmosphere:
    {
        // Check consistency of type and class.
        std::shared_ptr< ScaledAtmosphereSettings > scaledAtmosphereSettings =
                std::dynamic_pointer_cast< ScaledAtmosphereSettings >(
                    atmosphereSettings );
        if( scaledAtmosphereSettings == nullptr )
        {
            throw std::runtime_error(
                        "Error, expected scaled atmosphere settings for body " + body );
        }
        else
        {
            std::shared_ptr< AtmosphereModel > baseAtmosphere = createAtmosphereModel(
                        scaledAtmosphereSettings->getBaseSettings( ), body );
            atmosphereModel = std::make_shared< ScaledAtmosphereModel >(
                        baseAtmosphere, scaledAtmosphereSettings->getScaling( ), scaledAtmosphereSettings->getIsScalingAbsolute( ) );
        }
        break;
    }
    default:
        throw std::runtime_error( "Error, did not recognize atmosphere model settings type " +
                                  std::to_string( atmosphereSettings->getAtmosphereType( ) ) );
    }

    if( atmosphereSettings->getWindSettings( ) != nullptr )
    {
        atmosphereModel->setWindModel( createWindModel( atmosphereSettings->getWindSettings( ), body ) );
    }

    return atmosphereModel;
}

} // namespace simulation_setup

} // namespace tudat
