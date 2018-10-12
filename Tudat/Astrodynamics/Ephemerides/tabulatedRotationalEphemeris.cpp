#include "Tudat/Astrodynamics/Ephemerides/tabulatedRotationalEphemeris.h"

namespace tudat
{

namespace ephemerides
{

template class TabulatedRotationalEphemeris< double, double >;

#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template class TabulatedRotationalEphemeris< long double, double >;
template class TabulatedRotationalEphemeris< double, Time >;
template class TabulatedRotationalEphemeris< long double, Time >;
#endif

}

}

