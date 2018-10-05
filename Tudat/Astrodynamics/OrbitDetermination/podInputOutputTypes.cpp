#include "Tudat/Astrodynamics/OrbitDetermination/podInputOutputTypes.h"

namespace tudat
{

namespace simulation_setup
{

template class PodInput< double, double >;
template struct PodOutput< double >;

#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template class PodInput< long double, double >;
template class PodInput< double, Time >;
template class PodInput< long double, Time >;
template struct PodOutput< long double >;
#endif

}

}

