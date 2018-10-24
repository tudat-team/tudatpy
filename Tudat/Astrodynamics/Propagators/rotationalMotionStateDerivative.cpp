/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/Propagators/rotationalMotionStateDerivative.h"

namespace tudat
{

namespace propagators
{

//! Function to evaluated the classical rotational equations of motion (Euler equations)
Eigen::Vector3d evaluateRotationalEquationsOfMotion(
        const Eigen::Matrix3d& inertiaTensor, const Eigen::Vector3d& totalTorque,
        const Eigen::Vector3d& angularVelocityVector,
        const Eigen::Matrix3d& inertiaTensorTimeDerivative )
{
    return inertiaTensor.inverse( ) * ( totalTorque );
}

template class RotationalMotionStateDerivative< double, double >;

#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template class RotationalMotionStateDerivative< long double, double >;
template class RotationalMotionStateDerivative< double, Time >;
template class RotationalMotionStateDerivative< long double, Time >;
#endif

}

}
