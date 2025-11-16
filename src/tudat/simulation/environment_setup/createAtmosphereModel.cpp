/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/aerodynamics/exponentialAtmosphere.h"
#include "tudat/astro/aerodynamics/tabulatedAtmosphere.h"
#include "tudat/astro/aerodynamics/marsDtmAtmosphereModel.h"
#if TUDAT_BUILD_WITH_MCD
#include "tudat/astro/aerodynamics/mcdAtmosphereModel.h"
#endif
#if TUDAT_BUILD_WITH_NRLMSISE
#include "tudat/astro/aerodynamics/nrlmsise00Atmosphere.h"
#include "tudat/astro/aerodynamics/nrlmsise00InputFunctions.h"
#endif
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/solarActivityData.h"
#include "tudat/simulation/environment_setup/createAtmosphereModel.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create a wind model.
std::shared_ptr< aerodynamics::WindModel > createWindModel( const std::shared_ptr< WindModelSettings > windSettings,
                                                            const std::string& body )
{
    std::shared_ptr< aerodynamics::WindModel > windModel;

    // Check wind model type and create requested model
    switch( windSettings->getWindModelType( ) )
    {
        case constant_wind_model: {
            // Check input consistency
            std::shared_ptr< ConstantWindModelSettings > customWindModelSettings =
                    std::dynamic_pointer_cast< ConstantWindModelSettings >( windSettings );

            if( customWindModelSettings == nullptr )
            {
                throw std::runtime_error( "Error when making constant wind model for body " + body + ", input is incompatible" );
            }
            windModel = std::make_shared< aerodynamics::ConstantWindModel >( customWindModelSettings->getConstantWindVelocity( ),
                                                                             customWindModelSettings->getAssociatedFrame( ) );
            break;
        }
        case custom_wind_model: {
            // Check input consistency
            std::shared_ptr< CustomWindModelSettings > customWindModelSettings =
                    std::dynamic_pointer_cast< CustomWindModelSettings >( windSettings );

            if( customWindModelSettings == nullptr )
            {
                throw std::runtime_error( "Error when making custom wind model for body " + body + ", input is incompatible" );
            }
            windModel = std::make_shared< aerodynamics::CustomWindModel >( customWindModelSettings->getWindFunction( ),
                                                                           customWindModelSettings->getAssociatedFrame( ) );

            break;
        }
        default:
            throw std::runtime_error( "Error when making wind model for body " + body + ", input type not recognized" );
    }

    return windModel;
}

//! Function to create an atmosphere model.
std::shared_ptr< aerodynamics::AtmosphereModel > createAtmosphereModel( const std::shared_ptr< AtmosphereSettings > atmosphereSettings,
                                                                        const std::string& body )
{
    using namespace tudat::aerodynamics;

    // Declare return object.
    std::shared_ptr< AtmosphereModel > atmosphereModel;

    // Check which type of atmosphere model is to be created.
    switch( atmosphereSettings->getAtmosphereType( ) )
    {
        case exponential_atmosphere: {
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
                if( exponentialAtmosphereSettings->getBodyName( ) == undefined_body )
                {
                    exponentialAtmosphereModel =
                            std::make_shared< ExponentialAtmosphere >( exponentialAtmosphereSettings->getDensityScaleHeight( ),
                                                                       exponentialAtmosphereSettings->getConstantTemperature( ),
                                                                       exponentialAtmosphereSettings->getDensityAtZeroAltitude( ),
                                                                       exponentialAtmosphereSettings->getSpecificGasConstant( ),
                                                                       exponentialAtmosphereSettings->getRatioOfSpecificHeats( ) );
                }
                else
                {
                    exponentialAtmosphereModel = std::make_shared< ExponentialAtmosphere >( exponentialAtmosphereSettings->getBodyName( ) );
                }
                atmosphereModel = exponentialAtmosphereModel;
            }
            break;
        }
        case custom_constant_temperature_atmosphere: {
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
                if( customConstantTemperatureAtmosphereSettings->getModelSpecificParameters( ).empty( ) )
                {
                    customConstantTemperatureAtmosphereModel = std::make_shared< CustomConstantTemperatureAtmosphere >(
                            customConstantTemperatureAtmosphereSettings->getDensityFunction( ),
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
        case tabulated_atmosphere: {
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
                atmosphereModel = std::make_shared< TabulatedAtmosphere >( tabulatedAtmosphereSettings->getAtmosphereFile( ),
                                                                           tabulatedAtmosphereSettings->getIndependentVariables( ),
                                                                           tabulatedAtmosphereSettings->getDependentVariables( ),
                                                                           tabulatedAtmosphereSettings->getSpecificGasConstant( ),
                                                                           tabulatedAtmosphereSettings->getRatioOfSpecificHeats( ),
                                                                           tabulatedAtmosphereSettings->getBoundaryHandling( ),
                                                                           tabulatedAtmosphereSettings->getDefaultExtrapolationValue( ) );
            }
            break;
        }
        case mars_dtm_atmosphere: {
            std::shared_ptr< MarsDtmAtmosphereSettings > marsDtmAtmosphereSettings =
                    std::dynamic_pointer_cast< MarsDtmAtmosphereSettings >( atmosphereSettings );

            if( marsDtmAtmosphereSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating Mars DTM atmosphere model, model settings are incompatible." );
            }

            std::string spaceWeatherFilePath;

            if( marsDtmAtmosphereSettings->getSpaceWeatherFile( ) == "" )
            {
                // Use default space weather file stored in tudatBundle.
                spaceWeatherFilePath = paths::getSpaceWeatherDataPath( ) + "/sw19571001.txt";
            }
            else
            {
                // Use space weather file specified by user.
                spaceWeatherFilePath = marsDtmAtmosphereSettings->getSpaceWeatherFile( );
            }

            tudat::input_output::solar_activity::SolarActivityDataMap solarActivityData =
                    tudat::input_output::solar_activity::readSolarActivityData( spaceWeatherFilePath );
            std::shared_ptr< input_output::solar_activity::SolarActivityContainer > solarActivityContainer =
                    std::make_shared< input_output::solar_activity::SolarActivityContainer >( solarActivityData );
            std::function< double( const double ) > f107Function = [ = ]( const double time ) {
                return solarActivityContainer->getSolarActivityData( time )->solarRadioFlux107Observed / 2.25;
            };

            // Create atmosphere model using NRLMISE00 input function
            atmosphereModel = std::make_shared< aerodynamics::MarsDtmAtmosphereModel >( f107Function );
            break;
        }
#if TUDAT_BUILD_WITH_MCD
        case mcd_atmosphere: {
            std::shared_ptr< McdAtmosphereSettings > mcdAtmosphereSettings =
                    std::dynamic_pointer_cast< McdAtmosphereSettings >( atmosphereSettings );

            if( mcdAtmosphereSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating MCD atmosphere model, model settings are incompatible." );
            }

            // Create atmosphere model using MCD
            atmosphereModel = std::make_shared< aerodynamics::McdAtmosphereModel >( mcdAtmosphereSettings->getMcdDataPath( ),
                                                                                    mcdAtmosphereSettings->getDustScenario( ),
                                                                                    mcdAtmosphereSettings->getPerturbationKey( ),
                                                                                    mcdAtmosphereSettings->getPerturbationSeed( ),
                                                                                    mcdAtmosphereSettings->getGravityWaveLength( ),
                                                                                    mcdAtmosphereSettings->getHighResolutionMode( ) );
            break;
        }
#endif

#if TUDAT_BUILD_WITH_NRLMSISE
        case nrlmsise00: {
            std::string spaceWeatherFilePath;
            bool useStormConditions;
            bool useAnomalousOxygen;

            // Attempt to cast the atmosphereSettings to NRLMSISE00AtmosphereSettings
            std::shared_ptr< NRLMSISE00AtmosphereSettings > nrlmsise00AtmosphereSettings =
                    std::dynamic_pointer_cast< NRLMSISE00AtmosphereSettings >( atmosphereSettings );

            if( nrlmsise00AtmosphereSettings == nullptr )
            {
                // Use default space weather file stored in tudatBundle and geomagnetic storm conditions when no settings provided
                spaceWeatherFilePath = paths::getSpaceWeatherDataPath( ) + "/sw19571001.txt";
                useStormConditions = false;  // Default geomagnetic activity when no settings provided
                useAnomalousOxygen = true;   // Default geomagnetic activity when no settings provided
            }
            else
            {
                useStormConditions = nrlmsise00AtmosphereSettings->getUseStormConditions( );
                useAnomalousOxygen = nrlmsise00AtmosphereSettings->getUseAnomalousOxygen( );
                if( nrlmsise00AtmosphereSettings->getSpaceWeatherFile( ).empty( ) )
                {
                    // Use default space weather file stored in tudatBundle when the file is not provided
                    spaceWeatherFilePath = paths::getSpaceWeatherDataPath( ) + "/sw19571001.txt";
                }
                else
                {
                    spaceWeatherFilePath = nrlmsise00AtmosphereSettings->getSpaceWeatherFile( );
                }
            }

            // Read the solar activity data from the specified file
            tudat::input_output::solar_activity::SolarActivityDataMap solarActivityData =
                    tudat::input_output::solar_activity::readSolarActivityData( spaceWeatherFilePath );

            // Create the atmosphere model using the NRLMSISE00 input function
            atmosphereModel = std::make_shared< aerodynamics::NRLMSISE00Atmosphere >(
                    solarActivityData, true, useStormConditions, useAnomalousOxygen );
            break;
        }
#endif
        case scaled_atmosphere: {
            // Check consistency of type and class.
            std::shared_ptr< ScaledAtmosphereSettings > scaledAtmosphereSettings =
                    std::dynamic_pointer_cast< ScaledAtmosphereSettings >( atmosphereSettings );
            if( scaledAtmosphereSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected scaled atmosphere settings for body " + body );
            }
            else
            {
                std::shared_ptr< AtmosphereModel > baseAtmosphere =
                        createAtmosphereModel( scaledAtmosphereSettings->getBaseSettings( ), body );
                atmosphereModel = std::make_shared< ScaledAtmosphereModel >(
                        baseAtmosphere, scaledAtmosphereSettings->getScaling( ), scaledAtmosphereSettings->getIsScalingAbsolute( ) );
            }
            break;
        }
        default:
            throw std::runtime_error( "Error, did not recognize atmosphere model settings type " +
                                      std::to_string( atmosphereSettings->getAtmosphereType( ) )
#if TUDAT_BUILD_WITH_MCD
                                      + " (MCD support: enabled)"
#else
                                      + " (MCD support: disabled)"
#endif
            );
    }

    if( atmosphereSettings->getWindSettings( ) != nullptr )
    {
        atmosphereModel->setWindModel( createWindModel( atmosphereSettings->getWindSettings( ), body ) );
    }

    return atmosphereModel;
}

}  // namespace simulation_setup

}  // namespace tudat
