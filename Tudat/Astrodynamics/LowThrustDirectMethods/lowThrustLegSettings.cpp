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
//        simulation_setup::NamedBodyMap& bodyMap,
//        const std::string& bodyToPropagate,
//        const std::string& centralBody )
//{
//    // Declare LowThrustLeg object.
//    std::shared_ptr< LowThrustLeg > lowThrustLeg;

//    switch( lowThrustLegSettings->lowThrustLegType_ )
//    {
//    case hodographic_shaping_leg:
//    {
//        std::shared_ptr< HodographicShapingLegSettings > hodographicShapingSettings =
//                std::dynamic_pointer_cast< HodographicShapingLegSettings >( lowThrustLegSettings );

//        lowThrustLeg = std::make_shared< shape_based_methods::HodographicShaping >(
//                    stateAtDeparture, stateAtArrival, timeOfFlight, hodographicShapingSettings->numberOfRevolutions_,
//                    bodyMap, bodyToPropagate, centralBody, hodographicShapingSettings->radialVelocityFunctionComponents_,
//                    hodographicShapingSettings->normalVelocityFunctionComponents_, hodographicShapingSettings->axialVelocityFunctionComponents_,
//                    hodographicShapingSettings->freeCoefficientsNormalVelocityFunction_, hodographicShapingSettings->freeCoefficientsNormalVelocityFunction_,
//                    hodographicShapingSettings->freeCoefficientsAxialVelocityFunction_ );
//        break;
//    }
//    case spherical_shaping_leg:
//    {
//        std::shared_ptr< SphericalShapingLegSettings > sphericalShapingSettings =
//                std::dynamic_pointer_cast< SphericalShapingLegSettings >( lowThrustLegSettings );

//        lowThrustLeg = std::make_shared< shape_based_methods::SphericalShaping >(
//                    stateAtDeparture, stateAtArrival, timeOfFlight, sphericalShapingSettings->numberOfRevolutions_,
//                    bodyMap, bodyToPropagate, centralBody, sphericalShapingSettings->initialValueFreeCoefficient_,
//                    sphericalShapingSettings->rootFinderSettings_, sphericalShapingSettings->boundsFreeCoefficient_.first,
//                    sphericalShapingSettings->boundsFreeCoefficient_.second );
//        break;
//    }
//    case sims_flanagan_leg:
//    {
//        std::shared_ptr< SimsFlanaganLegSettings > SimsFlanaganSettings =
//                std::dynamic_pointer_cast< SimsFlanaganLegSettings >( lowThrustLegSettings );

//        lowThrustLeg = std::make_shared< low_thrust_direct_methods::SimsFlanagan >(
//                    stateAtDeparture, stateAtArrival, SimsFlanaganSettings->maximumThrust_, SimsFlanaganSettings->specificImpulseFunction_,
//                    SimsFlanaganSettings->numberSegments_, timeOfFlight, bodyMap, bodyToPropagate, centralBody,
//                    /*SimsFlanaganSettings->optimisationAlgorithm_, SimsFlanaganSettings->numberOfGenerations_,
//                    SimsFlanaganSettings->numberOfIndividualsPerPopulation_,*/ SimsFlanaganSettings->optimisationSettings_/*,
//                    SimsFlanaganSettings->relativeToleranceConstraints_,
//                    SimsFlanaganSettings->initialGuessThrustModel_*/ );

//        break;
//    }
//    case hybrid_method_leg:
//    {
//        std::shared_ptr< HybridMethodLegSettings > hybridMethodSettings =
//                std::dynamic_pointer_cast< HybridMethodLegSettings >( lowThrustLegSettings );

//        lowThrustLeg = std::make_shared< low_thrust_direct_methods::HybridMethod >(
//                    stateAtDeparture, stateAtArrival, hybridMethodSettings->maximumThrust_, hybridMethodSettings->specificImpulse_,
//                    timeOfFlight, bodyMap, bodyToPropagate, centralBody, hybridMethodSettings->integratorSettings_,
//                    /*hybridMethodSettings->optimisationAlgorithm_, hybridMethodSettings->numberOfGenerations_,
//                    hybridMethodSettings->numberOfIndividualsPerPopulation_,*/ hybridMethodSettings->optimisationSettings_/*, hybridMethodSettings->relativeToleranceConstraints_,
//                    hybridMethodSettings->initialGuessThrustModel_*/ );

//        break;
//    }
//    }

//    // Return low-thrust leg.
//    return lowThrustLeg;
//}

} // namespace transfer_trajectories

} // namespace tudat
