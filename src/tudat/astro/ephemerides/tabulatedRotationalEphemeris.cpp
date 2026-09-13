#include "tudat/astro/ephemerides/tabulatedRotationalEphemeris.h"

namespace tudat
{

namespace ephemerides
{

// template class TabulatedRotationalEphemeris< double, double >;

//! Function to check whether an ephemeris is a (type of) tabulated ephemeris
bool isTabulatedRotationalEphemeris( const std::shared_ptr< RotationalEphemeris > rotationalEphemeris )
{
    bool objectIsTabulated = 0;
    if( ( std::dynamic_pointer_cast< TabulatedRotationalEphemeris< double, double > >( rotationalEphemeris ) != nullptr ) ||
#if TUDAT_BUILD_WITH_HIGH_PRECISION_STATE_SCALAR
        ( std::dynamic_pointer_cast< TabulatedRotationalEphemeris< HighPrecisionStateScalar, double > >( rotationalEphemeris ) !=
          nullptr ) ||
        ( std::dynamic_pointer_cast< TabulatedRotationalEphemeris< HighPrecisionStateScalar, Time > >( rotationalEphemeris ) != nullptr ) ||
#endif
        ( std::dynamic_pointer_cast< TabulatedRotationalEphemeris< double, Time > >( rotationalEphemeris ) != nullptr ) )
    {
        objectIsTabulated = 1;
    }
    return objectIsTabulated;
}

}  // namespace ephemerides

}  // namespace tudat
