#include "tudatpy/scalarTypes.h"

#if TUDAT_BUILD_KERNEL_EXTERN_TEMPLATES

#include "tudat/simulation/estimation_setup/multiArcVariationalEquationsSolver.h"

namespace tudat
{

namespace propagators
{

template class MultiArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >;

}  // namespace propagators

}  // namespace tudat

#endif
