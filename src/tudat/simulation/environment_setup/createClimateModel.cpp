/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <stdexcept>

#include "tudat/simulation/environment_setup/createClimateModel.h"
#if TUDAT_BUILD_WITH_MCD_INTERFACE
#include "tudat/interface/mcd/marsClimateDatabaseClimateModel.h"
#endif

namespace tudat
{

namespace simulation_setup
{

//! Function to create an atmosphere model.
std::shared_ptr< environment::ClimateModel > createClimateModel( std::shared_ptr< ClimateModelSettings > climateModelSettings,
                                                                 std::shared_ptr< simulation_setup::Body > body )
{
    std::shared_ptr< environment::ClimateModel > climateModel;

    switch( climateModelSettings->climateModelType_ )
    {
#if TUDAT_BUILD_WITH_MCD_INTERFACE
        case mars_climate_database: {
            std::shared_ptr< MarsClimateDatabaseClimateModelSettings > mcdClimateModelSettings =
                    std::dynamic_pointer_cast< MarsClimateDatabaseClimateModelSettings >( climateModelSettings );
            if( mcdClimateModelSettings == nullptr )
            {
                throw std::runtime_error( "Error in creating MCD climate model" );
            }

            if( body->getBodyName( ) != "Mars" )
            {
                throw std::runtime_error( "Error, trying to create MCD for a body different than Mars" );
            }

            climateModel = std::make_shared< mcd_interface::MarsClimateDatabaseClimateModel >( mcdClimateModelSettings->mcdDataPath_,
                                                                                               mcdClimateModelSettings->dustScenario_,
                                                                                               mcdClimateModelSettings->perturbationKey_,
                                                                                               mcdClimateModelSettings->perturbationSeed_,
                                                                                               mcdClimateModelSettings->gravityWaveLength_,
                                                                                               mcdClimateModelSettings->highResolutionMode_,
                                                                                               mcdClimateModelSettings->maximumCacheSize_ );
            break;
        }
#endif
        default:
            throw std::runtime_error( "Error when making climate model, input type not recognized" );
    }

    return climateModel;
};

}  // namespace simulation_setup

}  // namespace tudat
