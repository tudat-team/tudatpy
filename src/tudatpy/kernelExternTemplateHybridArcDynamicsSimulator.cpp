#include "tudatpy/scalarTypes.h"

#if TUDAT_BUILD_KERNEL_EXTERN_TEMPLATES

#include "tudat/simulation/propagation_setup/hybridArcDynamicsSimulator.h"

namespace tudat
{

namespace propagators
{

template class HybridArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >;

}  // namespace propagators

}  // namespace tudat

#endif
