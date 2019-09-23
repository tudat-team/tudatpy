/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/LowThrustDirectMethods/lowThrustLeg.h"

namespace tudat
{

namespace transfer_trajectories
{

//! Retrieve acceleration model (thrust).
std::shared_ptr< propulsion::ThrustAcceleration > LowThrustLeg::getLowThrustAccelerationModel(
        std::function< double( const double ) > specificImpulseFunction )
{

    std::shared_ptr< simulation_setup::Body > vehicle = bodyMap_[ bodyToPropagate_ ];

    // Define thrust magnitude function from the shaped trajectory.
    std::function< double( const double ) > thrustMagnitudeFunction = [ = ]( const double currentTime )
    {

        // Compute current independent variable.
        double currentIndependentVariable = convertTimeToIndependentVariable( currentTime );

        // Compute current acceleration.
        double currentAcceleration = computeCurrentThrustAccelerationMagnitude( currentIndependentVariable );

        // Compute current mass of the vehicle.
        double currentMass = vehicle->getBodyMass();

        // Compute and return magnitude of the low-thrust force.
        return currentAcceleration * currentMass;
    };

    // Define thrust magnitude settings from thrust magnitude function.
    std::shared_ptr< simulation_setup::FromFunctionThrustMagnitudeSettings > thrustMagnitudeSettings
            = std::make_shared< simulation_setup::FromFunctionThrustMagnitudeSettings >(
                thrustMagnitudeFunction, specificImpulseFunction );


    // Define thrust direction function from the shaped trajectory.
    std::function< Eigen::Vector3d( const double ) > thrustDirectionFunction = [ = ]( const double currentTime )
    {
        // Compute current independent variable.
        double currentIndependentVariable = convertTimeToIndependentVariable( currentTime );

        // Compute current direction of the acceleration vector.
        Eigen::Vector3d currentAccelerationDirection = computeCurrentThrustAccelerationDirection( currentIndependentVariable );

        // Return direction of the low-thrust acceleration.
        return currentAccelerationDirection;
    };

    // Define thrust direction settings from the direction of thrust acceleration retrieved from the shaping method.
    std::shared_ptr< simulation_setup::CustomThrustDirectionSettings > thrustDirectionSettings =
            std::make_shared< simulation_setup::CustomThrustDirectionSettings >( thrustDirectionFunction );

    // Define thrust acceleration settings.
    std::shared_ptr< simulation_setup::ThrustAccelerationSettings > thrustAccelerationSettings =
            std::make_shared< simulation_setup::ThrustAccelerationSettings >(
                thrustDirectionSettings, thrustMagnitudeSettings );

    // Create low thrust acceleration model.
    std::shared_ptr< propulsion::ThrustAcceleration > lowThrustAccelerationModel =
            createThrustAcceleratioModel( thrustAccelerationSettings, bodyMap_, bodyToPropagate_ );

    return lowThrustAccelerationModel;

}


////! Full propagation.
//void LowThrustLeg::computeSemiAnalyticalAndFullPropagation(
//        std::function< double ( const double ) > specificImpulseFunction,
//        const std::shared_ptr<numerical_integrators::IntegratorSettings<double> > integratorSettings,
//        std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >,
//                std::shared_ptr< propagators::PropagatorSettings< double > > >& propagatorSettings,
//        std::map< double, Eigen::VectorXd >& fullPropagationResults,
//        std::map< double, Eigen::VectorXd >& semiAnalyticalResults,
//        std::map< double, Eigen::VectorXd>& dependentVariablesHistory )





} // namespace transfer_trajectories

} // namespace tudat
