/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/shapeBasedMethod.h"

namespace tudat
{

namespace shape_based_methods
{


basic_astrodynamics::AccelerationMap ShapeBasedMethod::retrieveLowThrustAccelerationMap(
        std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings )
{

    // Create low thrust acceleration model.
    std::shared_ptr< propulsion::ThrustAcceleration > lowThrustAccelerationModel =
            getLowThrustAccelerationModel( specificImpulseFunction, integratorSettings );

    // Acceleration from the central body.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > accelerationSettingsMap;
    accelerationSettingsMap[ centralBody_ ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );

    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ bodyToPropagate_ ] = accelerationSettingsMap;

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap_, accelerationMap, std::vector< std::string >{ bodyToPropagate_ }, std::vector< std::string >{ centralBody_ } );

    accelerationModelMap[ bodyToPropagate_ ][ bodyToPropagate_ ].push_back( lowThrustAccelerationModel );

    return accelerationModelMap;

}


//! Returns state history.
void ShapeBasedMethod::getTrajectory(
        std::vector< double >& epochsVector,
        std::map< double, Eigen::Vector6d >& propagatedTrajectory )
{
    propagatedTrajectory.clear( );

    for ( int i = 0 ; i < epochsVector.size() ; i++ )
    {
        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
        {
            throw std::runtime_error( "Error when retrieving the mass profile of a shape-based trajectory, "
                                      "epochs are not provided in increasing order." );
        }

        double currentIndependentVariable = convertTimeToIndependentVariable( epochsVector[ i ] );
        propagatedTrajectory[ epochsVector[ i ] ] = computeCurrentStateVector( currentIndependentVariable );
    }
}


//! Compute current thrust vector.
Eigen::Vector3d ShapeBasedMethod::computeCurrentThrust( double time,
                                                    std::function< double ( const double ) > specificImpulseFunction,
                                                    std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    double independentVariable = convertTimeToIndependentVariable( time );

    Eigen::Vector3d currentThrustVector = computeCurrentMass( 0.0, time, initialMass_, specificImpulseFunction, integratorSettings )
            * computeCurrentThrustAccelerationMagnitude( independentVariable, specificImpulseFunction, integratorSettings )
            * computeCurrentThrustAccelerationDirection( independentVariable, specificImpulseFunction, integratorSettings );

    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialMass_ );
    return currentThrustVector;
}


//! Return thrust profile.
void ShapeBasedMethod::getThrustProfile(
        std::vector< double >& epochsVector,
        std::map< double, Eigen::VectorXd >& thrustProfile,
        std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    thrustProfile.clear( );
    std::map< double, Eigen::VectorXd > massProfile;

    getMassProfile( epochsVector, massProfile, specificImpulseFunction, integratorSettings );

    for ( unsigned int i = 0 ; i < epochsVector.size() ; i++ )
    {
        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
        {
            throw std::runtime_error( "Error when retrieving the thrust profile of a shape-based trajectory, "
                                      "epochs are not provided in increasing order." );
        }

        double independentVariable = convertTimeToIndependentVariable( epochsVector[ i ] );

        double currentMass = massProfile[ epochsVector[ i ] ][ 0 ];
        thrustProfile[ epochsVector[ i ] ] = currentMass * computeCurrentThrustAccelerationMagnitude( independentVariable, specificImpulseFunction, integratorSettings )
                * computeCurrentThrustAccelerationDirection( independentVariable, specificImpulseFunction, integratorSettings );
    }
}


//! Define appropriate translational state propagator settings for the full propagation.
std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > ShapeBasedMethod::createLowThrustTranslationalStatePropagatorSettings(
        basic_astrodynamics::AccelerationMap accelerationModelMap,
        std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave )
{

    // Create termination conditions settings.
    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
            std::shared_ptr< propagators::PropagationTerminationSettings > > terminationConditions;

    terminationConditions.first = std::make_shared< propagators::PropagationTimeTerminationSettings >( 0.0, true );
    terminationConditions.second = std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight_, true );

    // Compute state vector at half of the time of flight.
    double independentVariableAtHalfTimeOfFlight = convertTimeToIndependentVariable( timeOfFlight_ / 2.0 );
    Eigen::Vector6d stateAtHalfOfTimeOfFlight = computeCurrentStateVector( independentVariableAtHalfTimeOfFlight );

    // Define translational state propagator settings.
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > translationalStatePropagatorSettings;

    // Define backward translational state propagation settings.
    translationalStatePropagatorSettings.first = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
            ( std::vector< std::string >{ centralBody_ }, accelerationModelMap, std::vector< std::string >{ bodyToPropagate_ },
              stateAtHalfOfTimeOfFlight, terminationConditions.first, propagators::cowell, dependentVariablesToSave );

    // Define forward translational state propagation settings.
    translationalStatePropagatorSettings.second = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
            ( std::vector< std::string >{ centralBody_ }, accelerationModelMap, std::vector< std::string >{ bodyToPropagate_ },
              stateAtHalfOfTimeOfFlight, terminationConditions.second, propagators::cowell, dependentVariablesToSave );

    return translationalStatePropagatorSettings;
}


} // namespace transfer_trajectories

} // namespace tudat
