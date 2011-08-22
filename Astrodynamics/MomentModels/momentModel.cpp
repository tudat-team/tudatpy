/*! \file momentModel.cpp
 *    This source file contains the moment model base class included in Tudat.
 *
 *    Path              : /Astrodynamics/ForceAndMomentModels/
 *    Version           : 2
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
 *    Date created      : 10 May, 2011
 *    Last modified     : 10 August, 2011
 *
 *    References
 *
 *    Notes
 *    See note in header file.
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
 *      110810    K. Kumar          Minor corrections; changed function names
 *                                  and removed redundant functions.
 */

// Include statements.
#include "momentModel.h"

//! Default constructor.
MomentModel::MomentModel( ) : pointerToForceModel_( NULL )
{
    forceApplicationPoint_.setZero( );
}

//! Default destructor.
MomentModel::~MomentModel( )
{
}

//! Set force application point.
void MomentModel::setForceApplicationPoint( Vector3d& forceApplicationPoint )
{
    forceApplicationPoint_ = forceApplicationPoint;
}

//! Get force application point.
Vector3d& MomentModel::getForceApplicationPoint( )
{
    return forceApplicationPoint_;
}

//! Set force model.
void MomentModel::setForceModel( ForceModel* pointerToForceModel )
{
    pointerToForceModel_ = pointerToForceModel;
}

//! Get force model.
ForceModel* MomentModel::getForceModel( )
{
    return pointerToForceModel_;
}

//! Get moment.
Vector3d& MomentModel::getMoment( )
{
    return moment_;
}

//! Set moment.
void MomentModel::setMoment( Vector3d& moment )
{
    moment_  = moment;
}

// End of file.
