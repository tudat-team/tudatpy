#include "tudat/simulation/propagation_setup/dynamicsSimulatorBase.h"

namespace tudat
{

namespace propagators
{

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS
template class DynamicsSimulator< double, double >;
#endif

}  // namespace propagators

}  // namespace tudat
