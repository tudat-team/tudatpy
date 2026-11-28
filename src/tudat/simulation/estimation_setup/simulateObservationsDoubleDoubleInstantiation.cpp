#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"

namespace tudat
{

namespace observation_models
{

template std::vector< std::shared_ptr< ObservationSimulatorBase< double, double > > > createObservationSimulators< double, double >(
        const std::vector< std::shared_ptr< ObservationModelSettings > >& observationSettingsList,
        const simulation_setup::SystemOfBodies& bodies );

}  // namespace observation_models

namespace simulation_setup
{

template std::shared_ptr< observation_models::ObservationCollection< double, double > > simulateObservations< double, double >(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< double > > >& observationsToSimulate,
        const std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< double, double > > >& observationSimulators,
        const SystemOfBodies bodies );

}  // namespace simulation_setup

}  // namespace tudat
