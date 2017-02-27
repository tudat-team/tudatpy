/*    Copyright (c) 2010-2017, Delft University of Technology
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

//! Function to calculate a rotation matrix from a body-fixed to inertial frame w.r.t. a constant rotation rate.
Eigen::Matrix3d calculatePartialOfRotationMatrixFromLocalFrameWrtConstantRotationRate(
        const Eigen::Quaterniond intertialBodyFixedToIntegrationFrame,
        const double rotationRate, const double timeSinceEpoch )
{
    double currentRotationAngle = rotationRate * timeSinceEpoch;

    // Compute partial of rotation term containing rotation rate
    Eigen::Matrix3d rotationMatrixDerivative;
    double sineOfAngle = sin( currentRotationAngle );
    double cosineOfAngle = cos( currentRotationAngle );
    rotationMatrixDerivative << -sineOfAngle, -cosineOfAngle, 0.0, cosineOfAngle, - sineOfAngle, 0.0, 0.0, 0.0, 0.0;

    return timeSinceEpoch * ( intertialBodyFixedToIntegrationFrame.toRotationMatrix( ) ) * rotationMatrixDerivative;

}
//! Function to calculate a rotation matrix from a body-fixed to inertial frame w.r.t. a constant pole right ascension
//! and declination
std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixFromLocalFrameWrtPoleOrientation(
        const Eigen::Vector3d initialOrientationAngles,//right ascension, declination, longitude of prime meridian.
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
                - ( mathematical_constants::PI / 2.0 - declination ) );/*
    declinationPartial<< 0.0, 0.0, 0.0,
            0.0, std::sin( mathematical_constants::PI / 2.0 - declination ),
            std::cos( mathematical_constants::PI / 2.0 - declination ),
            0.0, -std::cos( mathematical_constants::PI / 2.0 - declination ),
            std::sin( mathematical_constants::PI / 2.0 - declination );*/

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

//Eigen::Matrix< double, 3, Eigen::Dynamic > RotationMatrixPartial::calculatePartialOfInertialVelocityWrtParameter(
//        const double time,
//        const Eigen::Vector6d stateInLocalFrame )
//{
//    std::vector< Eigen::Matrix3d > rotationMatrixPartials =
//            calculatePartialOfRotationMatrixToBaseFrameWrParameter( time );
//    std::vector< Eigen::Matrix3d > rotationMatrixPartialPartials =
//            calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter( time );
//    //Eigen::Matrix3d currentRotationToBaseFrame =

//    Eigen::Matrix< double, 3, Eigen::Dynamic > rotatedVectorPartial = Eigen::Matrix< double, 3, Eigen::Dynamic >::Zero(
//                3, rotationMatrixPartials.size( ) );

//    for( unsigned int i = 0; i < rotationMatrixPartials.size( ); i++ )
//    {
//        rotatedVectorPartial.block( 0, i, 3, 1 ) = rotationMatrixPartials[ i ] * vectorInLocalFrame;
//    }
//    return rotatedVectorPartial;
//}


}

}
