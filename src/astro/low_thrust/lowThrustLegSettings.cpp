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

#include "tudat/astro/low_thrust/lowThrustLegSettings.h"

namespace tudat
{

namespace low_thrust_trajectories
{

//std::shared_ptr< low_thrust_trajectories::LowThrustLeg  > createLowThrustLeg(
//        const std::shared_ptr< LowThrustLegSettings >& lowThrustLegSettings,
//        const Eigen::Vector6d& stateAtDeparture,
//        const Eigen::Vector6d& stateAtArrival,
//        const double& timeOfFlight )
//{
//    // Declare LowThrustLeg object.
//    std::shared_ptr< LowThrustLeg > lowThrustLeg;

//    switch( lowThrustLegSettings->lowThrustLegType_ )
//    {
//
//#if( TUDAT_BUILD_WITH_PAGMO )
//    case sims_flanagan_leg:
//    {
//        std::shared_ptr< SimsFlanaganLegSettings > simsFlanaganSettings =
//                std::dynamic_pointer_cast< SimsFlanaganLegSettings >( lowThrustLegSettings );

//        lowThrustLeg = std::make_shared< low_thrust_trajectories::SimsFlanagan >(
//                    stateAtDeparture, stateAtArrival, simsFlanaganSettings->centralBodyGravitationalParameter_,
//                    simsFlanaganSettings->vehicleInitialMass_, simsFlanaganSettings->maximumThrust_,
//                    simsFlanaganSettings->specificImpulseFunction_,
//                    simsFlanaganSettings->numberSegments_, timeOfFlight,
//                    simsFlanaganSettings->optimisationSettings_ );

//        break;
//    }
//    case hybrid_method_leg:
//    {
//        throw std::runtime_error( "Hybrid low-thrust leg disabled" );
////        std::shared_ptr< HybridMethodLegSettings > hybridMethodSettings =
////                std::dynamic_pointer_cast< HybridMethodLegSettings >( lowThrustLegSettings );

////        lowThrustLeg = std::make_shared< low_thrust_trajectories::HybridMethod >(
////                    stateAtDeparture, stateAtArrival, hybridMethodSettings->maximumThrust_, hybridMethodSettings->specificImpulse_,
////                    timeOfFlight,
////                    hybridMethodSettings->optimisationSettings_, hybridMethodSettings->initialAndFinalMEEcostatesBounds_ );

//        break;
//    }
//#endif
//    }
//
//    // Return low-thrust leg.
//    return lowThrustLeg;
//}

} // namespace transfer_trajectories

} // namespace tudat
