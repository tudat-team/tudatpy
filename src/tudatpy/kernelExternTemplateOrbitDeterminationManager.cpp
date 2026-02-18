#include "tudatpy/scalarTypes.h"

#if TUDAT_BUILD_KERNEL_EXTERN_TEMPLATES

#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"

namespace tudat
{

namespace simulation_setup
{

template class OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE, 0 >;

}  // namespace simulation_setup

}  // namespace tudat

#endif
