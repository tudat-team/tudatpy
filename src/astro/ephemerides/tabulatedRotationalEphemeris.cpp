#include "tudat/astro/ephemerides/tabulatedRotationalEphemeris.h"

namespace tudat
{

namespace ephemerides
{

template class TabulatedRotationalEphemeris< double, double >;

#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template class TabulatedRotationalEphemeris< long double, double >;
template class TabulatedRotationalEphemeris< double, Time >;
template class TabulatedRotationalEphemeris< long double, Time >;
#endif

}

}

