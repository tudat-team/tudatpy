/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/propagators/rotationalMotionStateDerivative.h"

namespace tudat
{

namespace propagators
{


std::string getRotationalPropagatorName( const RotationalPropagatorType propagatorType )
{
    std::string propagatorName;
    switch( propagatorType )
    {
        case quaternions:
            propagatorName = "Quaternions";
            break;
        case modified_rodrigues_parameters:
            propagatorName = "Modified Rodrigues parameters";
            break;
        case exponential_map:
            propagatorName = "Exponential map";
            break;
        default:
            throw std::runtime_error( "Error when getting rotational propagator name, propagator" + std::to_string( static_cast< int >( propagatorType ) ) + " not found" );
            break;
    }
    return propagatorName;
}


int getRotationalStateSize( const RotationalPropagatorType propagatorType )
{
    int stateSize;
    switch( propagatorType )
    {
        case quaternions:
            stateSize = 7;
            break;
        case modified_rodrigues_parameters:
            stateSize = 7;
            break;
        case exponential_map:
            stateSize = 7;
            break;
        default:
            throw std::runtime_error( "Error when getting rotational propagator size, " + std::to_string( static_cast< int >( propagatorType ) ) + " not found" );
            break;
    }
    return stateSize;
}

//! Function to evaluated the classical rotational equations of motion (Euler equations)
Eigen::Vector3d evaluateRotationalEquationsOfMotion(
        const Eigen::Matrix3d& inertiaTensor, const Eigen::Vector3d& totalTorque,
        const Eigen::Vector3d& angularVelocityVector,
        const Eigen::Matrix3d& inertiaTensorTimeDerivative )
{
    Eigen::Matrix3d inverseInertiaTensor = inertiaTensor.inverse( );
    if( inverseInertiaTensor.hasNaN( ) )
    {
        std::cout<<inertiaTensor<<std::endl<<std::endl;
        std::cout<<inverseInertiaTensor<<std::endl<<std::endl;

        throw std::runtime_error( "Error when evaluating rotational equations of motion, inverse inertia tensor contains NaN. ");
    }
    return inverseInertiaTensor * ( totalTorque );
}

//template class RotationalMotionStateDerivative< double, double >;

}

}
