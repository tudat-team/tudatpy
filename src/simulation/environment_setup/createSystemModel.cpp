/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/environment_setup/createSystemModel.h"
#include "tudat/simulation/environment_setup/createThrustModelGuidance.h"

namespace tudat
{

namespace simulation_setup
{



void addEngineModel(
        const std::string& bodyName,
        const std::string& engineName,
        const std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > thrustSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const Eigen::Vector3d bodyFixedThrustDirection )
{
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > magnitudeUpdateSettings;
    std::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustMagnitudeWrapper = createThrustMagnitudeWrapper(
                thrustSettings,
            bodies, bodyName, magnitudeUpdateSettings );


    std::shared_ptr< system_models::EngineModel > vehicleEngineModel =
            std::make_shared< system_models::EngineModel >( thrustMagnitudeWrapper, engineName, [=](const double){return bodyFixedThrustDirection; } );

    if( bodies.at( bodyName )->getVehicleSystems( ) == nullptr )
    {
        std::shared_ptr< system_models::VehicleSystems > vehicleSystems = std::make_shared<
                system_models::VehicleSystems >(  );
        bodies.at( bodyName )->setVehicleSystems( vehicleSystems );
    }
    bodies.at( bodyName )->getVehicleSystems( )->setEngineModel( vehicleEngineModel );

}

} // namespace simulation_setup

} // namespace tudat

