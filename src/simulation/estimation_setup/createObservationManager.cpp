#include "tudat/simulation/estimation_setup/createObservationManager.h"

namespace tudat
{

namespace observation_models
{

//template std::shared_ptr< ObservationManagerBase< double, double > > createObservationManagerBase< double, double >(
//        const ObservableType observableType,
//        const std::map< LinkEnds, std::shared_ptr< ObservationModelSettings  > > settingsPerLinkEnds,
//        const simulation_setup::SystemOfBodies &bodies,
//        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate,
//        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface );
//template std::shared_ptr< ObservationManagerBase< double, Time > > createObservationManagerBase< double, Time >(
//        const ObservableType observableType,
//        const std::map< LinkEnds, std::shared_ptr< ObservationModelSettings  > > settingsPerLinkEnds,
//        const simulation_setup::SystemOfBodies &bodies,
//        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate,
//        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface );
//template std::shared_ptr< ObservationManagerBase< long double, double > > createObservationManagerBase< long double, double >(
//        const ObservableType observableType,
//        const std::map< LinkEnds, std::shared_ptr< ObservationModelSettings  > > settingsPerLinkEnds,
//        const simulation_setup::SystemOfBodies &bodies,
//        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate,
//        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface );
//template std::shared_ptr< ObservationManagerBase< long double, Time > > createObservationManagerBase< long double, Time >(
//        const ObservableType observableType,
//        const std::map< LinkEnds, std::shared_ptr< ObservationModelSettings  > > settingsPerLinkEnds,
//        const simulation_setup::SystemOfBodies &bodies,
//        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate,
//        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface );

//template std::shared_ptr< ObservationManagerBase< double, double > > createObservationManager< 1, double, double >(
//        const ObservableType observableType,
//        const std::map< LinkEnds, std::shared_ptr< ObservationModelSettings  > > settingsPerLinkEnds,
//        const simulation_setup::SystemOfBodies &bodies,
//        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate,
//        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface );
//template std::shared_ptr< ObservationManagerBase< double, Time > > createObservationManager< 1, double, Time >(
//        const ObservableType observableType,
//        const std::map< LinkEnds, std::shared_ptr< ObservationModelSettings  > > settingsPerLinkEnds,
//        const simulation_setup::SystemOfBodies &bodies,
//        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate,
//        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface );
//template std::shared_ptr< ObservationManagerBase< long double, double > > createObservationManager< 1, long double, double >(
//        const ObservableType observableType,
//        const std::map< LinkEnds, std::shared_ptr< ObservationModelSettings  > > settingsPerLinkEnds,
//        const simulation_setup::SystemOfBodies &bodies,
//        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate,
//        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface );
//template std::shared_ptr< ObservationManagerBase< long double, Time > > createObservationManager< 1, long double, Time >(
//        const ObservableType observableType,
//        const std::map< LinkEnds, std::shared_ptr< ObservationModelSettings  > > settingsPerLinkEnds,
//        const simulation_setup::SystemOfBodies &bodies,
//        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate,
//        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface );

//template std::shared_ptr< ObservationManagerBase< double, double > > createObservationManager< 2, double, double >(
//        const ObservableType observableType,
//        const std::map< LinkEnds, std::shared_ptr< ObservationModelSettings  > > settingsPerLinkEnds,
//        const simulation_setup::SystemOfBodies &bodies,
//        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate,
//        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface );
//template std::shared_ptr< ObservationManagerBase< double, Time > > createObservationManager< 2, double, Time >(
//        const ObservableType observableType,
//        const std::map< LinkEnds, std::shared_ptr< ObservationModelSettings  > > settingsPerLinkEnds,
//        const simulation_setup::SystemOfBodies &bodies,
//        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate,
//        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface );
//template std::shared_ptr< ObservationManagerBase< long double, double > > createObservationManager< 2, long double, double >(
//        const ObservableType observableType,
//        const std::map< LinkEnds, std::shared_ptr< ObservationModelSettings  > > settingsPerLinkEnds,
//        const simulation_setup::SystemOfBodies &bodies,
//        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate,
//        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface );
//template std::shared_ptr< ObservationManagerBase< long double, Time > > createObservationManager< 2, long double, Time >(
//        const ObservableType observableType,
//        const std::map< LinkEnds, std::shared_ptr< ObservationModelSettings  > > settingsPerLinkEnds,
//        const simulation_setup::SystemOfBodies &bodies,
//        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate,
//        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface );


//template std::shared_ptr< ObservationManagerBase< double, double > > createObservationManager< 3, double, double >(
//        const ObservableType observableType,
//        const std::map< LinkEnds, std::shared_ptr< ObservationModelSettings  > > settingsPerLinkEnds,
//        const simulation_setup::SystemOfBodies &bodies,
//        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate,
//        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface );
//template std::shared_ptr< ObservationManagerBase< double, Time > > createObservationManager< 3, double, Time >(
//        const ObservableType observableType,
//        const std::map< LinkEnds, std::shared_ptr< ObservationModelSettings  > > settingsPerLinkEnds,
//        const simulation_setup::SystemOfBodies &bodies,
//        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate,
//        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface );
//template std::shared_ptr< ObservationManagerBase< long double, double > > createObservationManager< 3, long double, double >(
//        const ObservableType observableType,
//        const std::map< LinkEnds, std::shared_ptr< ObservationModelSettings  > > settingsPerLinkEnds,
//        const simulation_setup::SystemOfBodies &bodies,
//        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate,
//        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface );
//template std::shared_ptr< ObservationManagerBase< long double, Time > > createObservationManager< 3, long double, Time >(
//        const ObservableType observableType,
//        const std::map< LinkEnds, std::shared_ptr< ObservationModelSettings  > > settingsPerLinkEnds,
//        const simulation_setup::SystemOfBodies &bodies,
//        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate,
//        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface );


}


}
