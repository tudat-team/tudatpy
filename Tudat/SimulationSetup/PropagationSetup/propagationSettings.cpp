#include "Tudat/SimulationSetup/PropagationSetup/propagationSettings.h"

namespace tudat
{

namespace propagators
{

template class PropagatorSettings< double >;
template class SingleArcPropagatorSettings< double >;
template class MultiArcPropagatorSettings< double >;
template class TranslationalStatePropagatorSettings< double >;
template class RotationalStatePropagatorSettings< double >;
template class MassPropagatorSettings< double >;
template class CustomStatePropagatorSettings< double >;
template class MultiTypePropagatorSettings< double >;

#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template class PropagatorSettings< long double >;
template class SingleArcPropagatorSettings< long double >;
template class MultiArcPropagatorSettings< long double >;
template class TranslationalStatePropagatorSettings< long double >;
template class RotationalStatePropagatorSettings< long double >;
template class MassPropagatorSettings< long double >;
template class MultiTypePropagatorSettings< long double >;
template class CustomStatePropagatorSettings< long double >;
template class CustomStatePropagatorSettings< double, Time >;
template class CustomStatePropagatorSettings< long double, Time >;
#endif

template std::map< IntegratedStateType, std::vector< std::pair< std::string, std::string > > > getIntegratedTypeAndBodyList< double >(
        const std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings );

} // namespace propagators

} // namespace tudat

