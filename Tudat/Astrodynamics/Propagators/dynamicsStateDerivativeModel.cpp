#include "Tudat/Astrodynamics/Propagators/dynamicsStateDerivativeModel.h"

namespace tudat
{
namespace propagators
{

template class DynamicsStateDerivativeModel< double, double >;
template class DynamicsStateDerivativeModel< double, long double >;
template class DynamicsStateDerivativeModel< Time, double >;
template class DynamicsStateDerivativeModel< Time, long double >;

} // namespace propagators

} // namespace tudat
