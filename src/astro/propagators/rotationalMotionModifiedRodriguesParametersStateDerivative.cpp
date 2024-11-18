/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/astro/propagators/rotationalMotionModifiedRodriguesParametersStateDerivative.h"

namespace tudat
{

namespace propagators
{

//! Function to obtain the time derivative of modified Rodrigues parameters of body-fixed to inertial frame.
Eigen::Vector4d calculateModifiedRodriguesParametersDerivative(
        const Eigen::Vector4d& currentModifiedRodriguesParametersToBaseFrame,
        const Eigen::Vector3d& angularVelocityVectorInBodyFixedFrame )
{
    // Declare eventual output vector
    Eigen::Vector4d modifiedRodriguesParametersDerivative = Eigen::Vector4d::Zero( );

    // Get intermediate variables
    Eigen::Vector3d modifiedRodriguesParametersVector = currentModifiedRodriguesParametersToBaseFrame.segment( 0, 3 );

    // Compute kinematic equation, i.e., derivative of modified Rodrigues parameters (also valid for SMRP)
    Eigen::Matrix3d skewModifiedRodriguesParametersVectorMatrix =
            linear_algebra::getCrossProductMatrix( modifiedRodriguesParametersVector );
    modifiedRodriguesParametersDerivative.segment( 0, 3 ) = 0.5 * (
                0.5 * ( 1.0 - std::pow( modifiedRodriguesParametersVector.norm( ), 2 ) ) * angularVelocityVectorInBodyFixedFrame +
                ( skewModifiedRodriguesParametersVectorMatrix + modifiedRodriguesParametersVector *
                  modifiedRodriguesParametersVector.transpose( ) ) * angularVelocityVectorInBodyFixedFrame );

    // Give output
    return modifiedRodriguesParametersDerivative;
}

//template class RotationalMotionModifiedRodriguesParametersStateDerivative< double, double >;

} // namespace propagators

} // namespace tudat
