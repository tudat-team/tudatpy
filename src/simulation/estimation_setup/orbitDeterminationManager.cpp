#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"

namespace tudat
{

namespace simulation_setup
{


template class OrbitDeterminationManager< double, double >;

#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template class OrbitDeterminationManager< double, Time >;
template class OrbitDeterminationManager< long double, double >;
template class OrbitDeterminationManager< long double, Time >;
#endif


}

}
