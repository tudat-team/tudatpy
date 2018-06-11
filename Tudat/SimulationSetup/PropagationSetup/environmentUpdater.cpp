#include "Tudat/SimulationSetup/PropagationSetup/environmentUpdater.h"

namespace tudat
{

namespace propagators
{

template class EnvironmentUpdater< double, double >;
template class EnvironmentUpdater< double, Time >;
template class EnvironmentUpdater< long double, double >;
template class EnvironmentUpdater< long double, Time >;

} // namespace propagators

} // namespace tudat

