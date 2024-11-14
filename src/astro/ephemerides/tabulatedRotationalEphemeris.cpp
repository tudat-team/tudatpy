#include "tudat/astro/ephemerides/tabulatedRotationalEphemeris.h"

namespace tudat
{

namespace ephemerides
{

//template class TabulatedRotationalEphemeris< double, double >;


//! Function to check whether an ephemeris is a (type of) tabulated ephemeris
bool isTabulatedRotationalEphemeris( const std::shared_ptr< RotationalEphemeris > rotationalEphemeris )
{
    bool objectIsTabulated = 0;
    if( ( std::dynamic_pointer_cast< TabulatedRotationalEphemeris< double, double > >( rotationalEphemeris ) != nullptr ) ||
            ( std::dynamic_pointer_cast< TabulatedRotationalEphemeris< long double, double > >( rotationalEphemeris ) != nullptr ) ||
            ( std::dynamic_pointer_cast< TabulatedRotationalEphemeris< long double, Time > >( rotationalEphemeris ) != nullptr ) ||
            ( std::dynamic_pointer_cast< TabulatedRotationalEphemeris< double, Time > >( rotationalEphemeris ) != nullptr ) )
    {
        objectIsTabulated = 1;
    }
    return objectIsTabulated;
}

}

}

