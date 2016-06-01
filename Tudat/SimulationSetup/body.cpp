/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/body.h"

namespace tudat
{

namespace simulation_setup
{

//! Templated function to set the state manually.
template< >
void Body::setTemplatedState( const Eigen::Matrix< double, 6, 1 >& state )
{
    setState( state );
}

//! Templated function to set the state manually.
template< >
void Body::setTemplatedState( const Eigen::Matrix< long double, 6, 1 >& state )
{
    setLongState( state );
}

//! Templated function to get the current state of the body from its ephemeris and
//! global-to-ephemeris-frame function.
template<  >
Eigen::Matrix< double, 6, 1 > Body::getTemplatedStateInBaseFrameFromEphemeris( const double& time )
{
    return getStateInBaseFrameFromEphemeris( time );
}


//! Templated function to get the current state of the body from its ephemeris and
//! global-to-ephemeris-frame function.
template<  >
Eigen::Matrix< long double, 6, 1 >
    Body::getTemplatedStateInBaseFrameFromEphemeris( const double& time )
{
    return getLongStateInBaseFrameFromEphemeris( time );
}

//! Templated function to set the current state of the body from its ephemeris and
//! global-to-ephemeris-frame function.
template< >
void Body::setTemplatedStateFromEphemeris< double, double >( const double& time )
{
    currentState_ = getTemplatedStateInBaseFrameFromEphemeris< double, double >( time );
}

//! Templated function to set the current state of the body from its ephemeris and
//! global-to-ephemeris-frame function.
template< >
void Body::setTemplatedStateFromEphemeris< long double, double >( const double& time )
{
    currentLongState_ = getTemplatedStateInBaseFrameFromEphemeris< long double, double >( time );
    currentState_ = currentLongState_.cast< double >( );
}

} // namespace simulation_setup

} // namespace tudat

