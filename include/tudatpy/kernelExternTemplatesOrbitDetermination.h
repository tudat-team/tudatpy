#ifndef TUDATPY_KERNEL_EXTERN_TEMPLATES_ORBIT_DETERMINATION_H
#define TUDATPY_KERNEL_EXTERN_TEMPLATES_ORBIT_DETERMINATION_H

#include "tudatpy/scalarTypes.h"

#if TUDAT_BUILD_KERNEL_EXTERN_TEMPLATES
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#endif

#if TUDAT_BUILD_KERNEL_EXTERN_TEMPLATES

namespace tudat
{
namespace simulation_setup
{

extern template class OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE, 0 >;

}  // namespace simulation_setup
}  // namespace tudat

#endif

#endif  // TUDATPY_KERNEL_EXTERN_TEMPLATES_ORBIT_DETERMINATION_H
