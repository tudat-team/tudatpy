/*! \file forceModel.cpp
 *    This source file contains the force model base class included in Tudat.
 *
 *    Path              : /Astrodynamics/MomentModels/
 *    Version           : 2
 *    Check status      : Unchecked
 *
 *    Checker           : F.M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 10 May, 2011
 *    Last modified     : 10 August, 2011
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
 *      110510    F.M. Engelen      First creation of code.
 *      110810    K. Kumar          Minor corrections.
 */

// Include statements.
#include "forceModel.h"

//! Default constructor.
ForceModel::ForceModel( )
{
    force_.setZero( );
}

//! Default destructor.
ForceModel::~ForceModel( )
{
}

//! Set force.
void ForceModel::setForce( const Vector3d& force )
{
    force_ = force;
}

//! Get force.
Vector3d& ForceModel::getForce( )
{
   return force_;
}

// End of file.
