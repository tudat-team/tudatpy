#include "Tudat/SimulationSetup/EstimationSetup/orbitDeterminationManager.h"

namespace tudat
{

namespace simulation_setup
{


template class OrbitDeterminationManager< double, double >;

#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template class OrbitDeterminationManager< double, Time >;
template class OrbitDeterminationManager< long double, double >;
template class OrbitDeterminationManager< long double, Time >;
#endif

template Eigen::Matrix< double, Eigen::Dynamic, 1 > getConcatenatedWeightsVector< double >(
        const typename std::map< observation_models::ObservableType, std::map<
        observation_models::LinkEnds, Eigen::VectorXd > >& weightsData );
template Eigen::Matrix< long double, Eigen::Dynamic, 1 > getConcatenatedWeightsVector< long double >(
        const typename std::map< observation_models::ObservableType, std::map<
        observation_models::LinkEnds, Eigen::VectorXd > >& weightsData );


}

}
