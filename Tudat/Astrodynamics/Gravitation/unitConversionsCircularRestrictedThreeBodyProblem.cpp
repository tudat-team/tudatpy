/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *        Wakker, K.F.,"Astrodynamics I, AE4-874", Delft University of Technology, 2007.
 *
 *    Notes
 *      Position dimension-scale is distance between the two massive bodies in the CRTBP.
 *      Time dimension-scale is based on orbital period of 2*pi.
 *
 */

#include <cmath>

#include "Tudat/Astrodynamics/Gravitation/unitConversionsCircularRestrictedThreeBodyProblem.h"

namespace tudat
{

namespace circular_restricted_three_body_problem
{

//! Convert dimensionless Cartesian state to dimensional units.
Eigen::VectorXd convertDimensionlessCartesianStateToDimensionalUnits(
        const Eigen::Vector6d &dimensionlessCartesianState,
        const double gravitationalParameterOfPrimaryBody,
        const double gravitationalParameterOfSecondaryBody,
        const double distanceBetweenPrimaries )
{
    // Declare dimensional Cartesian state.
    Eigen::Vector6d dimensionalCartesianState;

    // Convert position to dimensional units.
    dimensionalCartesianState.segment( 0, 3 )
            = dimensionlessCartesianState.segment( 0, 3 ) * distanceBetweenPrimaries;

    // Convert velocity to dimensional units.
    dimensionalCartesianState.segment( 3, 3 )
            = dimensionlessCartesianState.segment( 3, 3 )
            * std::sqrt( ( gravitationalParameterOfPrimaryBody
                           + gravitationalParameterOfSecondaryBody ) / distanceBetweenPrimaries );

    // Return state in dimensional units.
    return dimensionalCartesianState;
}

//! Convert dimensionless time to dimensional units.
double convertDimensionlessTimeToDimensionalTime(
        const double timeInDimensionlessUnits, const double gravitationalParameterOfPrimaryBody,
        const double gravitationalParameterOfSecondaryBody,
        const double distanceBetweenPrimaries )
{
    return timeInDimensionlessUnits * std::sqrt( std::pow( distanceBetweenPrimaries, 3.0 )
                                                 / ( gravitationalParameterOfPrimaryBody
                                                     + gravitationalParameterOfSecondaryBody ) );
}

//! Convert dimensional Cartesian state to dimensionless state.
Eigen::VectorXd convertDimensionalCartesianStateToDimensionlessState(
        const Eigen::Vector6d &dimensionalCartesianState,
        const double gravitationalParameterOfPrimaryBody,
        const double gravitationalParameterOfSecondaryBody,
        const double distanceBetweenPrimaries )
{
    // Declare dimensional Cartesian state.
    Eigen::Vector6d dimensionlessCartesianState;

    // Convert position to dimensional units.
    dimensionlessCartesianState.segment( 0, 3 )
            = dimensionalCartesianState.segment( 0, 3 ) / distanceBetweenPrimaries;

    // Convert velocity to dimensional units.
    dimensionlessCartesianState.segment( 3, 3 )
            = dimensionalCartesianState.segment( 3, 3 )
            * std::sqrt( distanceBetweenPrimaries /
                         ( gravitationalParameterOfPrimaryBody + gravitationalParameterOfSecondaryBody ) );

    // Return state in dimensional units.
    return dimensionlessCartesianState;
}



//! Convert dimensional time to dimensionless time.
double convertDimensionalTimeToDimensionlessTime(
        const double dimensionalTime,
        const double gravitationalParameterOfPrimaryBody,
        const double gravitationalParameterOfSecondaryBody,
        const double distanceBetweenPrimaries){

    return dimensionalTime * std::sqrt((gravitationalParameterOfPrimaryBody + gravitationalParameterOfSecondaryBody)
                                       / std::pow( distanceBetweenPrimaries, 3.0 ));

}


//! Convert corotating normalized state to inertial cartesian state.
Eigen::Vector6d convertCorotatingNormalizedToCartesianCoordinates(
        const double gravitationalParameterPrimary,
        const double gravitationalParameterSecondary,
        const double distancePrimarySecondary,
        const Eigen::Vector6d& normalizedState,
        const double normalizedTime )
{
    Eigen::Vector3d normalizedPosition = normalizedState.segment( 0, 3 );
    Eigen::Vector3d normalizedVelocity = normalizedState.segment( 3, 3 );

    Eigen::Matrix3d rotationMatrix;
    rotationMatrix.setZero( );
    Eigen::Matrix3d derivativeRotationMatrix;
    derivativeRotationMatrix.setZero( );

    rotationMatrix( 0, 0 ) = std::cos( normalizedTime );
    rotationMatrix( 0, 1 ) = - std::sin( normalizedTime );
    rotationMatrix( 1, 0 ) = std::sin( normalizedTime );
    rotationMatrix( 1, 1 ) = std::cos( normalizedTime );
    rotationMatrix( 2, 2 ) = 1.0;

    derivativeRotationMatrix( 0, 0 ) = - std::sin( normalizedTime );
    derivativeRotationMatrix( 0, 1 ) = - std::cos( normalizedTime );
    derivativeRotationMatrix( 1, 0 ) = std::cos( normalizedTime );
    derivativeRotationMatrix( 1, 1 ) = - std::sin( normalizedTime );

    Eigen::Vector6d inertialNormalizedState;
    inertialNormalizedState.segment( 0, 3 ) = rotationMatrix * normalizedPosition;
    inertialNormalizedState.segment( 3, 3 ) = derivativeRotationMatrix * normalizedPosition + rotationMatrix * normalizedVelocity;
    Eigen::Vector6d cartesianState = circular_restricted_three_body_problem::convertDimensionlessCartesianStateToDimensionalUnits(
                inertialNormalizedState, gravitationalParameterPrimary, gravitationalParameterSecondary, distancePrimarySecondary);

    return cartesianState;
}



//! Convert inertial cartesian state to co-rotating normalized state.
Eigen::Vector6d convertCartesianToCorotatingNormalizedCoordinates(
        const double gravitationalParameterPrimary,
        const double gravitationalParameterSecondary,
        const double distancePrimarySecondary,
        const Eigen::Vector6d& cartesianState,
        const double time )
{
    Eigen::Vector3d cartesianPosition = cartesianState.segment( 0, 3 );
    Eigen::Vector3d cartesianVelocity = cartesianState.segment( 3, 3 );

    double meanMotion = std::sqrt( ( gravitationalParameterPrimary + gravitationalParameterSecondary ) /
                                   std::pow( distancePrimarySecondary, 3 ) );

    Eigen::Matrix3d rotationMatrix;
    rotationMatrix.setZero( );
    rotationMatrix( 0, 0 ) = std::cos( meanMotion * time );
    rotationMatrix( 0, 1 ) = std::sin( meanMotion * time );
    rotationMatrix( 1, 0 ) = -std::sin( meanMotion * time );
    rotationMatrix( 1, 1 ) = std::cos( meanMotion * time );
    rotationMatrix( 2, 2 ) = 1.0;

    Eigen::Matrix3d derivativeRotationMatrix;
    derivativeRotationMatrix.setZero( );
    derivativeRotationMatrix( 0, 0 ) = -std::sin( meanMotion * time );
    derivativeRotationMatrix( 0, 1 ) = std::cos( meanMotion * time );
    derivativeRotationMatrix( 1, 0 ) = -std::cos( meanMotion * time );
    derivativeRotationMatrix( 1, 1 ) = -std::sin( meanMotion * time );
    derivativeRotationMatrix = meanMotion * derivativeRotationMatrix;

    Eigen::Vector6d corotatingDimensionalState;
    corotatingDimensionalState.segment( 0, 3 ) = rotationMatrix * cartesianPosition;
    corotatingDimensionalState.segment( 3, 3 ) = derivativeRotationMatrix * cartesianPosition + rotationMatrix * cartesianVelocity;

    Eigen::Vector6d normalizedState = circular_restricted_three_body_problem::convertDimensionalCartesianStateToDimensionlessState(
                corotatingDimensionalState, gravitationalParameterPrimary, gravitationalParameterSecondary, distancePrimarySecondary );

    return normalizedState;
}



} // namespace circular_restricted_three_body_problem
} // namespace tudat
