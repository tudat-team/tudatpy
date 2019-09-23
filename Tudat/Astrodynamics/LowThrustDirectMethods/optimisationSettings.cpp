/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/LowThrustDirectMethods/lowThrustLegSettings.h"

namespace tudat
{

namespace transfer_trajectories
{

////! Function to create the object determining the direction of the thrust acceleration.
//std::shared_ptr< transfer_trajectories::LowThrustLeg  > createLowThrustLeg(
//        const std::shared_ptr< LowThrustLegSettings >& lowThrustLegSettings,
//        const Eigen::Vector6d& stateAtDeparture,
//        const Eigen::Vector6d& stateAtArrival,
//        const double& timeOfFlight,
//        const simulation_setup::NamedBodyMap& bodyMap,
//        const std::string& bodyToPropagate )
//{
//    // Declare LowThrustLeg object.
//    std::shared_ptr< LowThrustLeg > lowThrustLeg;

//    switch( lowThrustLegSettings->lowThrustLegType_ )
//    {
//    case hodographic_shaping_leg:
//        break;
//    case spherical_shaping_leg:
//        break;
//    case sims_flanagan_leg:
//        break;
//    case hybrid_method_leg:
//        break;
//    }

//    // Return low-thrust leg.
//    return lowThrustLeg;
//}

} // namespace transfer_trajectories

} // namespace tudat
