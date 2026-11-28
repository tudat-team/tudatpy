#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"

namespace tudat
{

namespace simulation_setup
{

template std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > >
getInitialStateParameterSettings< double, double >( const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
                                                    const SystemOfBodies& bodies,
                                                    const std::vector< double > arcStartTimes );

template std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > >
getInitialMultiArcParameterSettings< double, double >(
        const std::shared_ptr< propagators::MultiArcPropagatorSettings< double, double > > propagatorSettings,
        const SystemOfBodies& bodies,
        const std::vector< double > arcStartTimes );

template std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > createParametersToEstimate< double, double >(
        const std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > >& parameterNames,
        const SystemOfBodies& bodies,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > >& considerParameterNames,
        bool printParameterOrderWarning );

}  // namespace simulation_setup

}  // namespace tudat
