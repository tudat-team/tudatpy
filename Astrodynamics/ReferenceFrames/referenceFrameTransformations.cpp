/*! \file referenceFrameTransformations.cpp
 *    This file contains the implementation of the reference frame transformations
 *    namespace included in Tudat.
 *
 *    Path              : /Astrodynamics/ReferenceFrames/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Checker           : F. M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 19 May, 2011
 *    Last modified     : 1 July, 2011
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110519    F.M. Engelen      First creation of code.
 *      110628    K. Kumar          Minor comment and layout changes; changed
 *                                  input arguments to pass-by-reference.
 *      110701    K. Kumar          Updated file path.
 */

// Include statements.
#include "referenceFrameTransformations.h"

//! Get rotating planetocentric (R) to intertial (I) reference frame transformation matrix.
Matrix3d reference_frame_transformations::getRotatingPlanetocentricToInertialFrameTransformation(
        const double& angleBetweenXIandXR )
{
    // Declare local variables.
    // Declare local matrix.
    Matrix3d localMatrix_;

    // Set local matrix.
    localMatrix_ = reference_frame_transformations::getInertialToPlanetocentricFrameTransformation(
                angleBetweenXIandXR );

    // Return transformation matrix.
    return localMatrix_.transpose();
}

//! Get inertial (I) to rotating planetocentric (R) reference frame transformtion matrix.
Matrix3d reference_frame_transformations::getInertialToPlanetocentricFrameTransformation(
        const double& angleBetweenXIandXR )
{
    // Compute rotation about Z-Axis.
    AngleAxisd angleAxisd = AngleAxisd( -angleBetweenXIandXR, Vector3d::UnitZ( ) );

    // Return transformation matrix.
    return angleAxisd.toRotationMatrix( );
}

//! Get Aerodynamic (airspeed based) (AA) to body reference frame (B) tranformation matrix.
Matrix3d reference_frame_transformations::getAirspeedBasedAerodynamicToBodyFrameTransformation(
        const double& angleOfAttack, const double& angleOfSideslip )
{
    // Compute rotation by side-slip angle about Z-Axis, followed by rotation negative angle of
    // attack angle about Y-Axis.
    AngleAxisd angleAxisd = AngleAxisd( angleOfSideslip, Vector3d::UnitZ( ) )*
                            AngleAxisd( -angleOfAttack, Vector3d::UnitY( ) );

    // Return transformation matrix.
    return angleAxisd.toRotationMatrix( );
}

// End of file.
