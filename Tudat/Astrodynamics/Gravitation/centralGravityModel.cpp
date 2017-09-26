/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <iomanip>
#include <cmath>

#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"

namespace tudat
{
namespace gravitation
{

 //! Compute gravitational acceleration.
Eigen::Vector3d computeGravitationalAcceleration(
        const double universalGravitationalConstant,
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double massOfBodyExertingAcceleration,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration )
{
    return computeGravitationalAcceleration(
                positionOfBodySubjectToAcceleration,
                universalGravitationalConstant * massOfBodyExertingAcceleration,
                positionOfBodyExertingAcceleration );
}

//! Compute gravitational acceleration.
Eigen::Vector3d computeGravitationalAcceleration(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBodyExertingAcceleration,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration )
{
//    std::cout<<"Calculating gravity "<<std::endl<<std::setprecision( 8 )<<
//              positionOfBodySubjectToAcceleration.transpose( )<<" "<<std::endl<<
//                positionOfBodyExertingAcceleration.transpose( )<<" "<<std::endl<<
//                ( positionOfBodySubjectToAcceleration - positionOfBodyExertingAcceleration ).transpose( )<<" "<<std::endl<<
//               gravitationalParameterOfBodyExertingAcceleration<<std::endl<<std::endl;
    double distance = ( positionOfBodySubjectToAcceleration - positionOfBodyExertingAcceleration ).norm( );
    return -gravitationalParameterOfBodyExertingAcceleration
            * ( positionOfBodySubjectToAcceleration - positionOfBodyExertingAcceleration )
            / ( distance * distance * distance );
}

//! Compute gravitational force.
Eigen::Vector3d computeGravitationalForce(
        const double universalGravitationalParameter,
        const double massOfBodySubjectToForce,
        const Eigen::Vector3d& positionOfBodySubjectToForce,
        const double massOfBodyExertingForce,
        const Eigen::Vector3d& positionOfBodyExertingForce )
{
    return massOfBodySubjectToForce * computeGravitationalAcceleration(
                universalGravitationalParameter, positionOfBodySubjectToForce,
                massOfBodyExertingForce, positionOfBodyExertingForce );
}

//! Compute gravitational force.
Eigen::Vector3d computeGravitationalForce(
        const double massOfBodySubjectToForce,
        const Eigen::Vector3d& positionOfBodySubjectToForce,
        const double gravitationalParameterOfBodyExertingForce,
        const Eigen::Vector3d& positionOfBodyExertingForce )
{
    return massOfBodySubjectToForce * computeGravitationalAcceleration(
                positionOfBodySubjectToForce, gravitationalParameterOfBodyExertingForce,
                positionOfBodyExertingForce );
}

} // namespace gravitation
} // namespace tudat
