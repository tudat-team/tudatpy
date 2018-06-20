#include "Tudat/Astrodynamics/OrbitDetermination/podInputOutputTypes.h"

namespace tudat
{

namespace simulation_setup
{

template class PodInput< double, double >;
template class PodInput< long double, double >;
template class PodInput< double, Time >;
template class PodInput< long double, Time >;

template struct PodOutput< double >;
template struct PodOutput< long double >;


}

}

