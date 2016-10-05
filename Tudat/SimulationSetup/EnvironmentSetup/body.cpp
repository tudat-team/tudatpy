/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"

namespace tudat
{

namespace simulation_setup
{

template< >
Eigen::Matrix< double, 6, 1 > BaseStateInterface::getBaseFrameState( const double time )
{
    return getBaseFrameDoubleState( time );
}

template< >
Eigen::Matrix< long double, 6, 1 > BaseStateInterface::getBaseFrameState( const double time )
{
    return getBaseFrameLongDoubleState( time );
}

template< >
Eigen::Matrix< double, 6, 1 > BaseStateInterface::getBaseFrameState( const Time time )
{
    return getBaseFrameDoubleState( time );
}

template< >
Eigen::Matrix< long double, 6, 1 > BaseStateInterface::getBaseFrameState( const Time time )
{
    return getBaseFrameLongDoubleState( time );
}






template< >
Eigen::Matrix< double, 6, 1 > Body::getTemplatedState( )
{
    return getState( );
}

template< >
Eigen::Matrix< long double, 6, 1 > Body::getTemplatedState( )
{
    return getLongState( );
}

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


} // namespace simulation_setup

} // namespace tudat

