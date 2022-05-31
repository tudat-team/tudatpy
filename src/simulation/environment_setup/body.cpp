/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include "tudat/astro/ephemerides/synchronousRotationalEphemeris.h"
#include "tudat/simulation/environment_setup/body.h"

namespace tudat
{

namespace simulation_setup
{

void Body::getPositionByReference( Eigen::Vector3d& position ) { position = currentState_.segment( 0, 3 ); }


//template void Body::setStateFromEphemeris< double, double >( const double& time );


//! Function through which the state of baseFrameId_ in the inertial frame can be determined
template< >
Eigen::Matrix< double, 6, 1 > BaseStateInterface::getBaseFrameState( const double time )
{
    return getBaseFrameDoubleState( time );
}

//! Function through which the state of baseFrameId_ in the inertial frame can be determined
template< >
Eigen::Matrix< long double, 6, 1 > BaseStateInterface::getBaseFrameState( const double time )
{
    return getBaseFrameLongDoubleState( time );
}

//! Function through which the state of baseFrameId_ in the inertial frame can be determined
template< >
Eigen::Matrix< double, 6, 1 > BaseStateInterface::getBaseFrameState( const Time time )
{
    return getBaseFrameDoubleState( time );
}

//! Function through which the state of baseFrameId_ in the inertial frame can be determined
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


//! Function to define whether the body is currently being propagated, or not
void Body::setIsBodyInPropagation( const bool isBodyInPropagation )
{
    isBodyInPropagation_ = isBodyInPropagation;

    if( rotationalEphemeris_ != nullptr )
    {
        rotationalEphemeris_->setIsBodyInPropagation( isBodyInPropagation );
    }

    if( !isBodyInPropagation )
    {
        isStateSet_ = false;
        isRotationSet_ = false;
    }
}

double getBodyGravitationalParameter( const SystemOfBodies& bodies, const std::string bodyName )
{
    if( bodies.count( bodyName ) == 0 )
    {
        throw std::runtime_error( "Error when getting gravitational parameter of body " + bodyName + ", no such body is found" );
    }
    else if( bodies.at( bodyName )->getGravityFieldModel( ) == nullptr )
    {
        throw std::runtime_error( "Error when getting gravitational parameter of body " + bodyName + ", body has not gravity field" );
    }
    return  bodies.at( bodyName )->getGravityFieldModel( )->getGravitationalParameter( );
}


//! Function ot retrieve the common global translational state origin of the environment
std::string getGlobalFrameOrigin( const SystemOfBodies& bodies )
{
    std::string globalFrameOrigin = "SSB";

    for( auto bodyIterator : bodies.getMap( ) )
    {
        if( bodyIterator.second->getIsBodyGlobalFrameOrigin( ) == -1 )
        {
            throw std::runtime_error( "Error, body " + bodyIterator.first + " does not have global frame origin set" );
        }
        else if( bodyIterator.second->getIsBodyGlobalFrameOrigin( ) == 1 )
        {
            if( globalFrameOrigin != "SSB" )
            {
                throw std::runtime_error( "Error, body " + bodyIterator.first + " found as global frame origin, but body " +
                                          globalFrameOrigin + " has already been detected as global frame origin." );
            }
            else
            {
               globalFrameOrigin = bodyIterator.first;
            }
        }
    }
    return globalFrameOrigin;
}

std::shared_ptr< ephemerides::ReferenceFrameManager > createFrameManager(
        const std::unordered_map< std::string, std::shared_ptr< Body > > bodies )
{
    // Get ephemerides from bodies
    std::map< std::string, std::shared_ptr< ephemerides::Ephemeris > > ephemerides;
    for( auto bodyIterator : bodies  )
    {
        if( bodyIterator.second->getEphemeris( ) != nullptr )
        {
            ephemerides[ bodyIterator.first ] = bodyIterator.second->getEphemeris( );
        }
    }
    return std::make_shared< ephemerides::ReferenceFrameManager >(
                ephemerides );
}

//! Function to set whether the bodies are currently being propagated, or not
void setAreBodiesInPropagation( const SystemOfBodies& bodies,
                                const bool areBodiesInPropagation )
{
    for( auto bodyIterator : bodies.getMap( )  )
    {
        bodyIterator.second->setIsBodyInPropagation( areBodiesInPropagation );
    }
}


} // namespace simulation_setup

} // namespace tudat

