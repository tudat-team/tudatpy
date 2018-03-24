/*    Copyright (c) 2010-2018, Delft University of Technology
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
    return getCartesianState( time );
}

//! Get state from ephemeris, with state scalar as template type (long double specialization).
template<  >
Eigen::Matrix< long double, 6, 1 > Ephemeris::getTemplatedStateFromEphemeris( const double& time )
{
    return getCartesianLongState( time );
}

//! Get state from ephemeris, with state scalar as template type (double specialization with Time input).
template<  >
Eigen::Matrix< double, 6, 1 > Ephemeris::getTemplatedStateFromEphemeris( const Time& time )
{
    return getCartesianStateFromExtendedTime( time );
}

//! Get state from ephemeris, with state scalar as template type (long double specialization with Time input).
template<  >
Eigen::Matrix< long double, 6, 1 > Ephemeris::getTemplatedStateFromEphemeris( const Time& time )
{
    return getCartesianLongStateFromExtendedTime( time );
}

//! Function to compute the relative state from two state functions.
void getRelativeState(
        Eigen::Vector6d& relativeState,
        const boost::function< Eigen::Vector6d( ) > stateFunctionOfBody,
        const boost::function< Eigen::Vector6d( ) > stateFunctionOfCentralBody )
{
    relativeState = stateFunctionOfBody( ) - stateFunctionOfCentralBody( );
}

} // namespace ephemerides

} // namespace tudat

