/*    Copyright (c) 2010-2017, Delft University of Technology
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



void updateBodyInertiaTensor(
        const boost::shared_ptr< Body > body,
        const Eigen::MatrixXd& unnormalizedCosineCoefficients,
        const Eigen::MatrixXd& unnormalizedSineCoefficients,
        const double bodyMass,
        const double referenceRadius )
{
    std::cerr<<"Update turned off A "<<std::endl;

//    body->setBodyInertiaTensor( getInertiaTensor(
//            unnormalizedCosineCoefficients,
//            unnormalizedSineCoefficients,
//            bodyMass / ( bodyMass * referenceRadius * referenceRadius ),
//            bodyMass, referenceRadius ) );

}

void updateBodyInertiaTensor(
        const boost::shared_ptr< Body > body,
        const boost::shared_ptr< gravitation::SphericalHarmonicsGravityField > gravityField )

{
    std::cerr<<"Update turned off A "<<std::endl;

//    updateBodyInertiaTensor( body, gravityField->getCosineCoefficients( 2, 2 ), gravityField->getSineCoefficients( 2, 2 ),
//                             gravityField->getGravitationalParameter( ) / physical_constants::GRAVITATIONAL_CONSTANT,
//                             gravityField->getReferenceRadius( ) );
}

} // namespace simulation_setup

} // namespace tudat

