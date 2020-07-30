#include "tudat/simulation/propagation_setup/environmentUpdater.h"

namespace tudat
{

namespace propagators
{

template class EnvironmentUpdater< double, double >;

#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template class EnvironmentUpdater< double, Time >;
template class EnvironmentUpdater< long double, double >;
template class EnvironmentUpdater< long double, Time >;
#endif

} // namespace propagators

} // namespace tudat

