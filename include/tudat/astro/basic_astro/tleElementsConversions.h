/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TLE_ELEMENTS_CONVERSIONS_H
#define TUDAT_TLE_ELEMENTS_CONVERSIONS_H

#include <functional>
#include <boost/math/special_functions/atanh.hpp>

#include <cmath>
#include <limits>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/interface/sofa/earthOrientation.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace basic_astrodynamics
{

Eigen::Vector6d convertStateFromTEMEtoJ2000( const double epochSinceJ2000,
                                             const Eigen::Vector6d stateTEME )
{    // First, rotate to the True Of Date (TOD) frame.
    double equationOfEquinoxes = sofa_interface::calculateEquationOfEquinoxes( epochSinceJ2000 );

    Eigen::Vector3d positionTEME = stateTEME.segment( 0, 3 );
    Eigen::Vector3d velocityTEME = stateTEME.segment( 3, 3 );

    // Rotate around pole (z-axis)
    Eigen::AngleAxisd rotationObject = Eigen::AngleAxisd( equationOfEquinoxes, Eigen::Vector3d::UnitZ( ) );
    Eigen::Matrix3d rotationMatrix1 = rotationObject.toRotationMatrix( );
    Eigen::Vector3d positionTOD = rotationMatrix1 * positionTEME;
    Eigen::Vector3d velocityTOD = rotationMatrix1 * velocityTEME;

    // These angles (zeta, z, and theta) do not really have descriptive names. For a description of the precession geometry and these angles,
    // see pages 226-228 and figure 3-31 in Vallado (2013).
    double precessionAngleModToGcrfZeta;
    double precessionAngleModToGcrfZ;
    double precessionAngleModToGcrfTheta;
    sofa_interface::getPrecessionAngles( precessionAngleModToGcrfZeta, precessionAngleModToGcrfZ,
                                         precessionAngleModToGcrfTheta,epochSinceJ2000 );

    // Now that we have our state vector in the TOD frame, we need to obtain the combined precession + nutation matrix from Sofa
    // (according to the 1976/1980 model)
    Eigen::Matrix3d precessionNutationMatrix = sofa_interface::getPrecessionNutationMatrix( epochSinceJ2000 );
    Eigen::Matrix3d rotationMatrix2 = precessionNutationMatrix.transpose( );
    // Multiply by inverted matrix to get to J2000
    Eigen::Vector3d positionJ2000 = rotationMatrix2 * positionTOD;
    Eigen::Vector3d velocityJ2000 = rotationMatrix2 * velocityTOD;

    Eigen::Vector6d stateJ2000;
    stateJ2000.segment( 0, 3 ) = positionJ2000;
    stateJ2000.segment( 3, 3 ) = velocityJ2000;

    return stateJ2000;
}

Eigen::Vector6d convertStateFromTEMEtoEclipJ2000( const double epochSinceJ2000,
                                                  const Eigen::Vector6d stateTEME )
{
    Eigen::Vector6d stateJ2000 = convertStateFromTEMEtoJ2000( epochSinceJ2000, stateTEME );
    Eigen::Vector3d positionJ2000 = stateJ2000.segment( 0, 3 );
    Eigen::Vector3d velocityJ2000 = stateJ2000.segment( 3, 3 );

    Eigen::Quaterniond itrsToEclipticQuaternion = spice_interface::computeRotationQuaternionBetweenFrames(
            "J2000", "ECLIPJ2000", epochSinceJ2000 );
    Eigen::Matrix3d rotationMatrix = itrsToEclipticQuaternion.toRotationMatrix( );
    Eigen::Vector3d positionEclipJ2000 = rotationMatrix * positionJ2000;
    Eigen::Vector3d velocityEclipJ2000 = rotationMatrix * velocityJ2000;

    Eigen::Vector6d stateEclipJ2000;
    stateEclipJ2000.segment( 0, 3 ) = positionEclipJ2000;
    stateEclipJ2000.segment( 3, 3 ) = velocityEclipJ2000;

    return stateEclipJ2000;
}

Eigen::Vector6d convertStateFromJ2000ToTEME( const double epochSinceJ2000,
                                             const Eigen::Vector6d stateJ2000 )
{
    // First, rotate to the True Of Date (TOD) frame.
    double equationOfEquinoxes = sofa_interface::calculateEquationOfEquinoxes( epochSinceJ2000 );

    Eigen::Vector3d positionJ2000 = stateJ2000.segment( 0, 3 );
    Eigen::Vector3d velocityJ2000 = stateJ2000.segment( 3, 3 );

    // Rotate around pole (z-axis)
    Eigen::AngleAxisd rotationObject = Eigen::AngleAxisd( equationOfEquinoxes, Eigen::Vector3d::UnitZ( ) );
    Eigen::Matrix3d rotationMatrix1 = rotationObject.toRotationMatrix( );

    // These angles (zeta, z, and theta) do not really have descriptive names. For a description of the precession geometry and these angles,
    // see pages 226-228 and figure 3-31 in Vallado (2013).
    double precessionAngleModToGcrfZeta;
    double precessionAngleModToGcrfZ;
    double precessionAngleModToGcrfTheta;
    sofa_interface::getPrecessionAngles( precessionAngleModToGcrfZeta, precessionAngleModToGcrfZ,
                                         precessionAngleModToGcrfTheta,epochSinceJ2000 );

    // Now that we have our state vector in the TOD frame, we need to obtain the combined precession + nutation matrix from Sofa
    // (according to the 1976/1980 model)
    Eigen::Matrix3d precessionNutationMatrix = sofa_interface::getPrecessionNutationMatrix( epochSinceJ2000 );
    Eigen::Matrix3d rotationMatrix2 = precessionNutationMatrix.transpose( );

    Eigen::Vector3d positionTEME = rotationMatrix1.transpose( ) * rotationMatrix2.transpose( ) * positionJ2000;
    Eigen::Vector3d velocityTEME = rotationMatrix1.transpose( ) * rotationMatrix2.transpose( ) * velocityJ2000;
    Eigen::Vector6d stateTEME;
    stateTEME.segment( 0, 3 ) = positionTEME;
    stateTEME.segment( 3, 3 ) = velocityTEME;

    return stateTEME;
}

Eigen::Vector6d convertStateFromEclipJ2000ToTEME( const double epochSinceJ2000,
                                                  const Eigen::Vector6d stateEclipJ2000 )
{
    // First, rotate to the True Of Date (TOD) frame.
    double equationOfEquinoxes = sofa_interface::calculateEquationOfEquinoxes( epochSinceJ2000 );

    Eigen::Vector3d positionEclipJ2000 = stateEclipJ2000.segment( 0, 3 );
    Eigen::Vector3d velocityEclipJ2000 = stateEclipJ2000.segment( 3, 3 );

    // Rotate around pole (z-axis)
    Eigen::AngleAxisd rotationObject = Eigen::AngleAxisd( equationOfEquinoxes, Eigen::Vector3d::UnitZ( ) );
    Eigen::Matrix3d rotationMatrix1 = rotationObject.toRotationMatrix( );

    // These angles (zeta, z, and theta) do not really have descriptive names. For a description of the precession geometry and these angles,
    // see pages 226-228 and figure 3-31 in Vallado (2013).
    double precessionAngleModToGcrfZeta;
    double precessionAngleModToGcrfZ;
    double precessionAngleModToGcrfTheta;
    sofa_interface::getPrecessionAngles( precessionAngleModToGcrfZeta, precessionAngleModToGcrfZ,
                                         precessionAngleModToGcrfTheta,epochSinceJ2000 );

    // Now that we have our state vector in the TOD frame, we need to obtain the combined precession + nutation matrix from Sofa
    // (according to the 1976/1980 model)
    Eigen::Matrix3d precessionNutationMatrix = sofa_interface::getPrecessionNutationMatrix( epochSinceJ2000 );
    Eigen::Matrix3d rotationMatrix2 = precessionNutationMatrix.transpose( );

    Eigen::Quaterniond itrsToEclipticQuaternion = spice_interface::computeRotationQuaternionBetweenFrames(
            "J2000", "ECLIPJ2000", epochSinceJ2000 );
    Eigen::Matrix3d rotationMatrix3 = itrsToEclipticQuaternion.toRotationMatrix( );

    Eigen::Vector3d positionTEME = rotationMatrix1.transpose( ) * rotationMatrix2.transpose( ) *
                                   rotationMatrix3.transpose( ) * positionEclipJ2000;
    Eigen::Vector3d velocityTEME = rotationMatrix1.transpose( ) * rotationMatrix2.transpose( ) *
                                   rotationMatrix3.transpose( ) * velocityEclipJ2000;
    Eigen::Vector6d stateTEME;
    stateTEME.segment( 0, 3 ) = positionTEME;
    stateTEME.segment( 3, 3 ) = velocityTEME;

    return stateTEME;

}

} // namespace basic_astrodynamics

} // namespace tudat

#endif // TUDAT_TLE_ELEMENTS_CONVERSIONS_H