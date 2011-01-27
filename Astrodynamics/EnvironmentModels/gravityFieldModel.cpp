/*! \file gravityFieldModel.h
 *    Header file that defines the gravity field model included in Tudat.
 *
 *    Path              : /Astrodynamics/EnvironmentModel/
 *    Version           : 2
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J.C.P Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 8 December, 2010
 *    Last modified     : 15 December, 2010
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
 *      YYMMDD    author              comment
 *      101208    K. Kumar            File created.
 *      101215    K. Kumar            Removed consts from end of get functions.
 */

// Include statements.
#include "gravityFieldModel.h"

//! Default constructor.
GravityFieldModel::GravityFieldModel( )
{
    // Initialize variables.
    gravitationalParameter_ = -0.0;
    positionVectorOfOrigin_.setConstant( -0.0 );
}

//! Default destructor.
GravityFieldModel::~GravityFieldModel( )
{
}

//! Set the gravitational parameter.
void GravityFieldModel::setGravitationalParameter(
        const double& gravitationalParameter )
{
    gravitationalParameter_ = gravitationalParameter;
}

//! Set origin of gravity field.
void GravityFieldModel::setOrigin( Vector3d& positionVectorOfOrigin )
{
    positionVectorOfOrigin_ = positionVectorOfOrigin;
}

//! Get the gravitational parameter.
double GravityFieldModel::getGravitationalParameter( )
{
    return gravitationalParameter_;
}

//! Get origin of gravity field.
Vector3d GravityFieldModel::getOrigin( )
{
    return positionVectorOfOrigin_;
}

// End of file.
