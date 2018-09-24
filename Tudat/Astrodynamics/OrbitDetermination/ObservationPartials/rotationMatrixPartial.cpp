/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/rotationMatrixPartial.h"

namespace tudat
{

namespace observation_partials
{

//! Function to calculate a partial of rotation matrix from a body-fixed to inertial frame w.r.t. a constant rotation rate.
Eigen::Matrix3d calculatePartialOfRotationMatrixFromLocalFrameWrtConstantRotationRate(
        const Eigen::Quaterniond initialBodyFixedToIntegrationFrame,
        const double rotationRate, const double timeSinceEpoch )
{
    double currentRotationAngle = rotationRate * timeSinceEpoch;

    // Compute partial of rotation term containing rotation rate
    Eigen::Matrix3d rotationMatrixDerivative;
    double sineOfAngle = sin( currentRotationAngle );
    double cosineOfAngle = cos( currentRotationAngle );
    rotationMatrixDerivative << -sineOfAngle, -cosineOfAngle, 0.0, cosineOfAngle, - sineOfAngle, 0.0, 0.0, 0.0, 0.0;

    return timeSinceEpoch * ( initialBodyFixedToIntegrationFrame.toRotationMatrix( ) ) * rotationMatrixDerivative;

}

//! Function to calculate partial of a rotation matrix derivative (body-fixed to inertial) w.r.t. a constant rotation rate.
Eigen::Matrix3d calculatePartialOfRotationMatrixFromLocalFrameDerivativeWrtConstantRotationRate(
        const Eigen::Matrix3d currentRotationFromLocalToGlobalFrame,
        const double rotationRate, const double timeSinceEpoch )
{
    return currentRotationFromLocalToGlobalFrame * reference_frames::Z_AXIS_ROTATION_MATRIX_DERIVATIVE_PREMULTIPLIER *
            ( rotationRate * timeSinceEpoch * reference_frames::Z_AXIS_ROTATION_MATRIX_DERIVATIVE_PREMULTIPLIER -
              Eigen::Matrix3d::Identity( ) );

}
//! Function to calculate a partial of rotation matrix from a body-fixed to inertial frame w.r.t. a constant pole right
//! ascension and declination
std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixFromLocalFrameWrtPoleOrientation(
        const Eigen::Vector3d initialOrientationAngles,
        const double rotationRate, const double timeSinceEpoch )
{
    Eigen::Matrix3d commonTerm = Eigen::AngleAxisd(
                -1.0 * ( -initialOrientationAngles.z( ) - rotationRate * timeSinceEpoch ),
                Eigen::Vector3d::UnitZ( ) ).toRotationMatrix( );


    double rightAscension = initialOrientationAngles.x( );
    double declination = initialOrientationAngles.y( );

    // Compute partial of rotation term containing right ascension.
    Eigen::Matrix3d rightAscensionPartial = -reference_frames::getDerivativeOfZAxisRotationWrtAngle(
                -( rightAscension + mathematical_constants::PI / 2.0 ) );


    // Compute partial of rotation term containing declination.
    Eigen::Matrix3d declinationPartial = reference_frames::getDerivativeOfXAxisRotationWrtAngle(
                - ( mathematical_constants::PI / 2.0 - declination ) );

    // Compute partials.
    std::vector< Eigen::Matrix3d > rotationMatrixPartials;
    rotationMatrixPartials.push_back(
                rightAscensionPartial * Eigen::Matrix3d(
                    Eigen::AngleAxisd( -( declination - mathematical_constants::PI / 2.0 ),
                                       Eigen::Vector3d::UnitX( ) ) ) * commonTerm );
    rotationMatrixPartials.push_back(
                Eigen::AngleAxisd( ( rightAscension + mathematical_constants::PI / 2.0 ),
                                   Eigen::Vector3d::UnitZ( ) ).toRotationMatrix( ) *
                declinationPartial * commonTerm );

    return rotationMatrixPartials;
}

//! Function to calculate a partial of rotation matrix derivative from a body-fixed to inertial frame w.r.t. a constant
//! pole right  ascension and declination.
std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixFromLocalFrameDerivativeWrtPoleOrientation(
        const Eigen::Vector3d initialOrientationAngles,
        const double rotationRate, const double timeSinceEpoch )
{
    std::vector< Eigen::Matrix3d > partialsOfRotationMatrix =
            calculatePartialOfRotationMatrixFromLocalFrameWrtPoleOrientation(
                initialOrientationAngles, rotationRate, timeSinceEpoch );
    partialsOfRotationMatrix[ 0 ] = -rotationRate *
            partialsOfRotationMatrix.at( 0 ) * reference_frames::Z_AXIS_ROTATION_MATRIX_DERIVATIVE_PREMULTIPLIER ;
    partialsOfRotationMatrix[ 1 ] = -rotationRate *
            partialsOfRotationMatrix.at( 1 )* reference_frames::Z_AXIS_ROTATION_MATRIX_DERIVATIVE_PREMULTIPLIER ;

    return partialsOfRotationMatrix;
}


//! Function to calculate the partial of the position of a vector, which is given in a body-fixed frame, in the inertial
//! frame wrt a parameter.
Eigen::Matrix< double, 3, Eigen::Dynamic > RotationMatrixPartial::calculatePartialOfInertialPositionWrtParameter(
        const double time,
        const Eigen::Vector3d vectorInLocalFrame )
{
    std::vector< Eigen::Matrix3d > rotationMatrixPartials = calculatePartialOfRotationMatrixToBaseFrameWrParameter( time );
    Eigen::Matrix< double, 3, Eigen::Dynamic > rotatedVectorPartial = Eigen::Matrix< double, 3, Eigen::Dynamic >::Zero(
                3, rotationMatrixPartials.size( ) );

    for( unsigned int i = 0; i < rotationMatrixPartials.size( ); i++ )
    {
        rotatedVectorPartial.block( 0, i, 3, 1 ) = rotationMatrixPartials[ i ] * vectorInLocalFrame;
    }
    return rotatedVectorPartial;
}

//! Function to calculate the partial of the velocity of a vector, which is given in a body-fixed frame, in the inertial
//! frame wrt a parameter.
Eigen::Matrix< double, 3, Eigen::Dynamic > RotationMatrixPartial::calculatePartialOfInertialVelocityWrtParameter(
        const double time,
        const Eigen::Vector3d vectorInLocalFrame )
{
    if( rotationModel_ == nullptr )
    {
        throw std::runtime_error( "Error when caling RotationMatrixPartial::calculatePartialOfInertialVelocityWrtParameter, rotation model is nullptr" );
    }

    // Compute rotation matrix (derivative) partials
    std::vector< Eigen::Matrix3d > rotationMatrixPartials =
            calculatePartialOfRotationMatrixToBaseFrameWrParameter( time );
    std::vector< Eigen::Matrix3d > rotationMatrixDerivativePartials =
            calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter( time );

    // Compute current rotation matrix and derivative
    Eigen::Matrix3d currentRotationToBaseFrame = rotationModel_->getRotationToBaseFrame( time ).toRotationMatrix( );
    Eigen::Matrix3d currentRotationToTargetFrameDerivative = rotationModel_->getDerivativeOfRotationToTargetFrame(
                time );

    Eigen::Matrix< double, 3, Eigen::Dynamic > rotatedVectorPartial = Eigen::Matrix< double, 3, Eigen::Dynamic >::Zero(
                3, rotationMatrixPartials.size( ) );

    // Compute inertial velocity partial
    for( unsigned int i = 0; i < rotationMatrixPartials.size( ); i++ )
    {
        rotatedVectorPartial.block( 0, i, 3, 1 ) =
                -( rotationMatrixPartials.at( i ) * currentRotationToTargetFrameDerivative * currentRotationToBaseFrame +
                 currentRotationToBaseFrame * rotationMatrixDerivativePartials.at( i ).transpose( ) * currentRotationToBaseFrame +
                  currentRotationToBaseFrame * currentRotationToTargetFrameDerivative * rotationMatrixPartials.at( i ) ) * vectorInLocalFrame;
    }
    return rotatedVectorPartial;
}


}

}
