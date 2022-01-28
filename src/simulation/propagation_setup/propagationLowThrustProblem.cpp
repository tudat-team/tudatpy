#include "tudat/simulation/propagation_setup/propagationLowThrustProblem.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"

namespace tudat
{

namespace simulation_setup
{

void computeLowThrustLegSemiAnalyticalAndFullPropagation(
        const std::shared_ptr< low_thrust_trajectories::LowThrustLeg > lowThrustLeg,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >,
        std::shared_ptr< propagators::PropagatorSettings< double > > >& propagatorSettings,
        std::map< double, Eigen::VectorXd >& fullPropagationResults,
        std::map< double, Eigen::Vector6d >& semiAnalyticalResults,
        std::map< double, Eigen::VectorXd >& dependentVariablesHistory )
{

    fullPropagationResults.clear();
    semiAnalyticalResults.clear();
    dependentVariablesHistory.clear();

    // Define backward propagator settings variables.
    integratorSettings->initialTimeStep_ = - std::fabs( integratorSettings->initialTimeStep_ );
    integratorSettings->initialTime_ = lowThrustLeg->getTimeOfFlight( ) / 2.0;

    // Perform the backward propagation.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationBackwards(
                bodies, integratorSettings, propagatorSettings.first );
    std::map< double, Eigen::VectorXd > stateHistoryFullProblemBackwardPropagation =
            dynamicsSimulatorIntegrationBackwards.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableHistoryBackwardPropagation =
            dynamicsSimulatorIntegrationBackwards.getDependentVariableHistory( );

    // Declare vector of epochs at which the trajectory has to be calculated.
    std::vector< double > epochsVector;

    // Compute and save full propagation and shaping method results along the backward propagation direction
    for( std::map< double, Eigen::VectorXd >::iterator itr = stateHistoryFullProblemBackwardPropagation.begin( );
         itr != stateHistoryFullProblemBackwardPropagation.end( ); itr++ )
    {
        epochsVector.push_back( itr->first );
        fullPropagationResults[ itr->first ] = itr->second;
        dependentVariablesHistory[ itr->first ] = dependentVariableHistoryBackwardPropagation[ itr->first ];
    }

    // Define forward propagator settings variables.
    integratorSettings->initialTimeStep_ = std::fabs( integratorSettings->initialTimeStep_ );
    integratorSettings->initialTime_ = lowThrustLeg->getTimeOfFlight( ) / 2.0;

    // Perform forward propagation.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationForwards(
                bodies, integratorSettings, propagatorSettings.second );
    std::map< double, Eigen::VectorXd > stateHistoryFullProblemForwardPropagation =
            dynamicsSimulatorIntegrationForwards.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableHistoryForwardPropagation =
            dynamicsSimulatorIntegrationForwards.getDependentVariableHistory( );

    // Compute and save full propagation and shaping method results along the forward propagation direction.
    for( std::map< double, Eigen::VectorXd >::iterator itr = stateHistoryFullProblemForwardPropagation.begin( );
         itr != stateHistoryFullProblemForwardPropagation.end( ); itr++ )
    {
        epochsVector.push_back( itr->first );
        fullPropagationResults[ itr->first ] = itr->second;
        dependentVariablesHistory[ itr->first ] = dependentVariableHistoryForwardPropagation[ itr->first ];
    }

    lowThrustLeg->getTrajectory( epochsVector, semiAnalyticalResults );

    // Reset initial integrator settings.
    integratorSettings->initialTimeStep_ = std::fabs( integratorSettings->initialTimeStep_ );

}


basic_astrodynamics::AccelerationMap retrieveLowThrustAccelerationMap(
        const std::shared_ptr< low_thrust_trajectories::LowThrustLeg > lowThrustLeg,
        const simulation_setup::SystemOfBodies& bodies,
        const std::string& bodyToPropagate,
        const std::string& centralBody,
        const std::function< double( const double ) > specificImpulseFunction,
        const double lowThrustLegInitialTime )
{

    // Acceleration from the central body.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > accelerationSettingsMap;
    accelerationSettingsMap[ centralBody ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                          basic_astrodynamics::point_mass_gravity ) );
    accelerationSettingsMap[ bodyToPropagate ].push_back(
                getLowThrustLegAccelerationSettings(
                    lowThrustLeg, bodies, bodyToPropagate, specificImpulseFunction, lowThrustLegInitialTime  ) );


    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ bodyToPropagate ] = accelerationSettingsMap;

    return createAccelerationModelsMap(
                bodies, accelerationMap, std::vector< std::string >{ bodyToPropagate }, std::vector< std::string >{ centralBody } );

}


//! Define appropriate translational state propagator settings for the full propagation.
std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > createLowThrustTranslationalStatePropagatorSettings(
        const std::shared_ptr< low_thrust_trajectories::LowThrustLeg > lowThrustLeg,
        const std::string& bodyToPropagate,
        const std::string& centralBody,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave )
{

    // Create termination conditions settings.
    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
            std::shared_ptr< propagators::PropagationTerminationSettings > > terminationConditions;

    terminationConditions.first = std::make_shared< propagators::PropagationTimeTerminationSettings >(
                0.0, true );
    terminationConditions.second = std::make_shared< propagators::PropagationTimeTerminationSettings >(
                lowThrustLeg->getTimeOfFlight( ), true );

    // Compute state vector at half of the time of flight.
    double independentVariableAtHalfTimeOfFlight = lowThrustLeg->convertTimeToIndependentVariable(
                lowThrustLeg->getTimeOfFlight( ) / 2.0 );
    Eigen::Vector6d stateAtHalfOfTimeOfFlight = lowThrustLeg->computeCurrentStateVector(
                independentVariableAtHalfTimeOfFlight );

    // Define translational state propagator settings.
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > translationalStatePropagatorSettings;

    // Define backward translational state propagation settings.
    translationalStatePropagatorSettings.first = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
            ( std::vector< std::string >{ centralBody }, accelerationModelMap, std::vector< std::string >{ bodyToPropagate },
              stateAtHalfOfTimeOfFlight, terminationConditions.first, propagators::cowell, dependentVariablesToSave );

    // Define forward translational state propagation settings.
    translationalStatePropagatorSettings.second = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
            ( std::vector< std::string >{ centralBody }, accelerationModelMap, std::vector< std::string >{ bodyToPropagate },
              stateAtHalfOfTimeOfFlight, terminationConditions.second, propagators::cowell, dependentVariablesToSave );

    return translationalStatePropagatorSettings;
}

//! Define appropriate propagator settings for the full propagation.
std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >,
std::shared_ptr< propagators::PropagatorSettings< double > > > createLowThrustPropagatorSettings(
        const std::shared_ptr< low_thrust_trajectories::LowThrustLeg > lowThrustLeg,
        const double bodyMassAtMidPoint,
        const simulation_setup::SystemOfBodies& bodies,
        const std::string& bodyToPropagate,
        const std::string& centralBody,
        const std::function< double( const double ) > specificImpulseFunction,
        const basic_astrodynamics::AccelerationMap perturbingAccelerationsMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings,
        const std::shared_ptr< propagators::DependentVariableSaveSettings >& dependentVariablesToSave )
{

    // Define propagator settings.
    std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >, std::shared_ptr< propagators::PropagatorSettings< double > > >
            propagatorSettings;

    // Create termination conditions settings.
    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
            std::shared_ptr< propagators::PropagationTerminationSettings > > terminationConditions;

    terminationConditions.first = std::make_shared< propagators::PropagationTimeTerminationSettings >(
                0.0, true );
    terminationConditions.second = std::make_shared< propagators::PropagationTimeTerminationSettings >(
                lowThrustLeg->getTimeOfFlight( ), true );

    // Re-initialise spacecraft mass in system of bodies.
    bodies.at( bodyToPropagate )->setConstantBodyMass( bodyMassAtMidPoint );

    // Retrieve low-thrust trajectory nominal accelerations (thrust + central gravity accelerations).
    basic_astrodynamics::AccelerationMap lowThrustTrajectoryAccelerations = retrieveLowThrustAccelerationMap(
               lowThrustLeg,  bodies, bodyToPropagate, centralBody, specificImpulseFunction, 0.0 );

    // Add perturbing accelerations given as input.
    basic_astrodynamics::AccelerationMap accelerationModelMap = perturbingAccelerationsMap;
    accelerationModelMap[ bodyToPropagate ][ bodyToPropagate ].push_back(
                lowThrustTrajectoryAccelerations.at( bodyToPropagate ).at( bodyToPropagate ).at( 0 ) );
    accelerationModelMap[ bodyToPropagate ][ centralBody ].push_back(
                lowThrustTrajectoryAccelerations.at( bodyToPropagate ).at( centralBody ).at( 0 ) );


    // Define translational state propagation settings
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > translationalStatePropagatorSettings =
            createLowThrustTranslationalStatePropagatorSettings(
                lowThrustLeg, bodyToPropagate, centralBody, accelerationModelMap, dependentVariablesToSave );

    throw std::runtime_error( "Error, cannot compute mas yet for low-thrust" );
    double massHalfOfTimeOfFlight;// = computeCurrentMass( timeOfFlight_ / 2.0, specificImpulseFunction, integratorSettings );

    // Create settings for propagating the mass of the vehicle.
    std::pair< std::shared_ptr< propagators::MassPropagatorSettings< double > >,
            std::shared_ptr< propagators::MassPropagatorSettings< double > > > massPropagatorSettings;

    // Create mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ bodyToPropagate ] = simulation_setup::createMassRateModel(
            bodyToPropagate, std::make_shared< simulation_setup::FromThrustMassRateSettings >(1 ), bodies, accelerationModelMap );

    // Define backward mass propagation settings.
    massPropagatorSettings.first = std::make_shared< propagators::MassPropagatorSettings< double > >(
                std::vector< std::string >{ bodyToPropagate }, massRateModels,
                ( Eigen::Matrix< double, 1, 1 >( ) << massHalfOfTimeOfFlight ).finished( ), terminationConditions.first );

    // Define forward mass propagation settings.
    massPropagatorSettings.second = std::make_shared< propagators::MassPropagatorSettings< double > >(
                std::vector< std::string >{ bodyToPropagate }, massRateModels,
                ( Eigen::Matrix< double, 1, 1 >( ) << massHalfOfTimeOfFlight ).finished( ), terminationConditions.second );

    // Create list of propagation settings.
    std::pair< std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > >,
            std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > > > propagatorSettingsVector;

    // Backward propagator settings vector.
    propagatorSettingsVector.first.push_back( translationalStatePropagatorSettings.first );
    propagatorSettingsVector.first.push_back( massPropagatorSettings.first );

    // Forward propagator settings vector.
    propagatorSettingsVector.second.push_back( translationalStatePropagatorSettings.second );
    propagatorSettingsVector.second.push_back( massPropagatorSettings.second );

    // Backward hybrid propagation settings.
    propagatorSettings.first = std::make_shared< propagators::MultiTypePropagatorSettings< double > >(
                propagatorSettingsVector.first, terminationConditions.first, dependentVariablesToSave );

    // Forward hybrid propagation settings.
    propagatorSettings.second = std::make_shared< propagators::MultiTypePropagatorSettings< double > >(
                propagatorSettingsVector.second, terminationConditions.second, dependentVariablesToSave );


    return propagatorSettings;
}



} // namespace simulation_setup

} // namespace tudat
