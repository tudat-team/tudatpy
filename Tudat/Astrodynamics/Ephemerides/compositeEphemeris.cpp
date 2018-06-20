#include "Tudat/Astrodynamics/Ephemerides/compositeEphemeris.h"

namespace tudat
{

namespace ephemerides
{

template class CompositeEphemeris< double, double >;
template class CompositeEphemeris< Time, double >;
template class CompositeEphemeris< double, long double >;
template class CompositeEphemeris< Time, long double >;

} // namespace ephemerides

} // namespace tudat
