#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/propagation_setup/propagationResults.h"

namespace tudat
{

namespace simulation_setup
{

template class OrbitDeterminationManager< double, double >;

template struct CovarianceAnalysisOutput< double, double >;
template struct EstimationOutput< double, double >;

}  // namespace simulation_setup

}  // namespace tudat
