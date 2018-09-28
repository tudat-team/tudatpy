#include "Tudat/Astrodynamics/ObservationModels/observationManager.h"

namespace tudat
{

namespace observation_models
{

template class ObservationManagerBase< double, double >;
template class ObservationManager< 1, double, double >;
template class ObservationManager< 2, double, double >;
template class ObservationManager< 3, double, double >;
template class ObservationManager< 6, double, double >;

#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template class ObservationManagerBase< double, Time >;
template class ObservationManagerBase< long double, double >;
template class ObservationManagerBase< long double, Time >;

template class ObservationManager< 1, double, Time >;
template class ObservationManager< 1, long double, double >;
template class ObservationManager< 1, long double, Time >;

template class ObservationManager< 2, double, Time >;
template class ObservationManager< 2, long double, double >;
template class ObservationManager< 2, long double, Time >;

template class ObservationManager< 3, double, Time >;
template class ObservationManager< 3, long double, double >;
template class ObservationManager< 3, long double, Time >;

template class ObservationManager< 6, double, Time >;
template class ObservationManager< 6, long double, double >;
template class ObservationManager< 6, long double, Time >;

#endif
}

}
