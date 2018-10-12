#include "Tudat/Astrodynamics/Propagators/nBodyCowellStateDerivative.h"

namespace tudat
{

namespace propagators
{

template class NBodyCowellStateDerivative< double, double >;

#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template class NBodyCowellStateDerivative< long double, double >;
template class NBodyCowellStateDerivative< double, Time >;
template class NBodyCowellStateDerivative< long double, Time >;
#endif


} // namespace propagators

} // namespace tudat

