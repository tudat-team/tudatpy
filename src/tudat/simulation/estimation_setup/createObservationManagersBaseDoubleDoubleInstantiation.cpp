#include "tudat/simulation/estimation_setup/createObservationManager.h"

namespace tudat
{

namespace observation_models
{

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS
template std::map< ObservableType, std::shared_ptr< ObservationManagerBase< double, double > > >
createObservationManagersBase< double, double >(
        const std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >& observationSettingsList,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > fullParameters,
        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface,
        const std::shared_ptr< propagators::DependentVariablesInterface< double > > dependentVariablesInterface );
#endif

}  // namespace observation_models

}  // namespace tudat
