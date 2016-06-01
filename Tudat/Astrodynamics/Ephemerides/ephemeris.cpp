/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"

namespace tudat
{
namespace ephemerides
{

//! Get state from ephemeris, with state scalar as template type (double specialization).
template<  >
Eigen::Matrix< double, 6, 1 > Ephemeris::getTemplatedStateFromEphemeris( const double& time )
{
    return getCartesianStateFromEphemeris( time, basic_astrodynamics::JULIAN_DAY_ON_J2000 );
}

//! Get state from ephemeris, with state scalar as template type (long double specialization).
template<  >
Eigen::Matrix< long double, 6, 1 > Ephemeris::getTemplatedStateFromEphemeris( const double& time )
{
    return getCartesianLongStateFromEphemeris( time, basic_astrodynamics::JULIAN_DAY_ON_J2000 );
}


} // namespace ephemerides

} // namespace tudat

