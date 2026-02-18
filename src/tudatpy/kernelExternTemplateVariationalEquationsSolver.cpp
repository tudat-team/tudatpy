#include "tudatpy/scalarTypes.h"

#if TUDAT_BUILD_KERNEL_EXTERN_TEMPLATES

#include "tudat/simulation/estimation_setup/variationalEquationsSolverBase.h"

namespace tudat
{

namespace propagators
{

template class VariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE, 0 >;

}  // namespace propagators

}  // namespace tudat

#endif
