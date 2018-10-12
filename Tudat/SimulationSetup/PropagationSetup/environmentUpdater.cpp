#include "Tudat/SimulationSetup/PropagationSetup/environmentUpdater.h"

namespace tudat
{

namespace propagators
{

template class EnvironmentUpdater< double, double >;

#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template class EnvironmentUpdater< double, Time >;
template class EnvironmentUpdater< long double, double >;
template class EnvironmentUpdater< long double, Time >;
#endif

} // namespace propagators

} // namespace tudat

