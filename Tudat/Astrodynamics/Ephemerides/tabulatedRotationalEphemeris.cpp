#include "Tudat/Astrodynamics/Ephemerides/tabulatedRotationalEphemeris.h"

namespace tudat
{

namespace ephemerides
{

template class TabulatedRotationalEphemeris< double, double >;
template class TabulatedRotationalEphemeris< long double, double >;
template class TabulatedRotationalEphemeris< double, Time >;
template class TabulatedRotationalEphemeris< long double, Time >;

}

}

