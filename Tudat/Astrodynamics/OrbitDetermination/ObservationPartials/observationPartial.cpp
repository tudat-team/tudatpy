#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/observationPartial.h"

namespace tudat
{

namespace observation_partials
{

template class ObservationPartial< 1 >;
template class ObservationPartial< 2 >;
template class ObservationPartial< 3 >;

template class ObservationPartialWrtConstantAbsoluteBias< 1 >;
template class ObservationPartialWrtConstantAbsoluteBias< 2 >;
template class ObservationPartialWrtConstantAbsoluteBias< 3 >;

template class ObservationPartialWrtArcWiseAbsoluteBias< 1 >;
template class ObservationPartialWrtArcWiseAbsoluteBias< 2 >;
template class ObservationPartialWrtArcWiseAbsoluteBias< 3 >;

template class ObservationPartialWrtConstantRelativeBias< 1 >;
template class ObservationPartialWrtConstantRelativeBias< 2 >;
template class ObservationPartialWrtConstantRelativeBias< 3 >;

template class ObservationPartialWrtArcWiseRelativeBias< 1 >;
template class ObservationPartialWrtArcWiseRelativeBias< 2 >;
template class ObservationPartialWrtArcWiseRelativeBias< 3 >;

}

}
