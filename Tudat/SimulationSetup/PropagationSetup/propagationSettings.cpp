#include "Tudat/SimulationSetup/PropagationSetup/propagationSettings.h"

namespace tudat
{

namespace propagators
{

template class PropagatorSettings< double >;
template class PropagatorSettings< long double >;

template class SingleArcPropagatorSettings< double >;
template class SingleArcPropagatorSettings< long double >;

template class MultiArcPropagatorSettings< double >;
template class MultiArcPropagatorSettings< long double >;

template class TranslationalStatePropagatorSettings< double >;
template class TranslationalStatePropagatorSettings< long double >;

template class RotationalStatePropagatorSettings< double >;
template class RotationalStatePropagatorSettings< long double >;

template class MassPropagatorSettings< double >;
template class MassPropagatorSettings< long double >;

template class CustomStatePropagatorSettings< double >;
template class CustomStatePropagatorSettings< long double >;
template class CustomStatePropagatorSettings< double, Time >;
template class CustomStatePropagatorSettings< long double, Time >;

template class MultiTypePropagatorSettings< double >;
template class MultiTypePropagatorSettings< long double >;

} // namespace propagators

} // namespace tudat

