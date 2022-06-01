/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATESSYTEMMODEL_H
#define TUDAT_CREATESSYTEMMODEL_H

#include <memory>

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/thrustSettings.h"

namespace tudat
{

namespace simulation_setup
{


void addEngineModel(
        const std::string& bodyName,
        const std::string& engineName,
        const std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > thrustSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const Eigen::Vector3d bodyFixedThrustDirection = Eigen::Vector3d::UnitX( ) );


} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATESSYTEMMODEL_H
