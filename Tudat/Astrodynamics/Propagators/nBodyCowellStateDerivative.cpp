#include "Tudat/Astrodynamics/Propagators/nBodyCowellStateDerivative.h"

namespace tudat
{

namespace propagators
{

template class NBodyCowellStateDerivative< double, double >;
template class NBodyCowellStateDerivative< long double, double >;
template class NBodyCowellStateDerivative< double, Time >;
template class NBodyCowellStateDerivative< long double, Time >;



} // namespace propagators

} // namespace tudat

