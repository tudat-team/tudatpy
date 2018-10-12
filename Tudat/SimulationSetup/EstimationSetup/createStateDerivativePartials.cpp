#include "Tudat/SimulationSetup/EstimationSetup/createStateDerivativePartials.h"

namespace tudat
{

namespace simulation_setup
{

template std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap > createStateDerivativePartials< double, double >(
        const std::unordered_map< propagators::IntegratedStateType,
        std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< double, double > > > >
        stateDerivativeModels,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > >
        parametersToEstimate );

#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap > createStateDerivativePartials< long double, double >(
        const std::unordered_map< propagators::IntegratedStateType,
        std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< long double, double > > > >
        stateDerivativeModels,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > >
        parametersToEstimate );
template std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap > createStateDerivativePartials< double, Time >(
        const std::unordered_map< propagators::IntegratedStateType,
        std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< double, Time > > > >
        stateDerivativeModels,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > >
        parametersToEstimate );
template std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap > createStateDerivativePartials< long double, Time >(
        const std::unordered_map< propagators::IntegratedStateType,
        std::vector< std::shared_ptr< propagators::SingleStateTypeDerivative< long double, Time > > > >
        stateDerivativeModels,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > >
        parametersToEstimate );
#endif

} // namespace simulation_setup

} // namespace tudat

