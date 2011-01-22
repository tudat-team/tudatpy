/*! \file cartesianElements.cpp
 *    This source file contains the Cartesian elements class included in Tudat.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 22 Checked, 2010
 *    Last modified     : 02 December, 2010
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
 *      YYMMDD    author        comment
 *      101022    K. Kumar      First creation of code.
 *      101026    K. Kumar      Added position and velocity vector functions.
 *      101201    E. Iorfida    Modified punctuation.
 *      101202    J. Melman     Large modification: only set state_ and use
 *                              state_ to get components. Does not use private
 *                              variables, so those are obsolete.
 */

// Include statements.
#include "cartesianElements.h"

//! Default constructor.
CartesianElements::CartesianElements( )
{
}

//! Default destructor.
CartesianElements::~CartesianElements( )
{
}

//! Set Cartesian element: x.
void CartesianElements::setCartesianElementX( const double& cartesianElementX )
{
    state_( 0 ) = cartesianElementX;
}

//! Set Cartesian element: y.
void CartesianElements::setCartesianElementY( const double& cartesianElementY )
{
    state_( 1 ) = cartesianElementY;
}

//! Set Cartesian element: z.
void CartesianElements::setCartesianElementZ( const double& cartesianElementZ )
{
    state_( 2 ) = cartesianElementZ;
}

//! Set Cartesian element: xDot.
void CartesianElements::setCartesianElementXDot( const double&
                                                 cartesianElementXDot )
{
     state_( 3 ) = cartesianElementXDot;
}

//! Set Cartesian element: yDot.
void CartesianElements::setCartesianElementYDot( const double&
                                                 cartesianElementYDot )
{
    state_( 4 ) = cartesianElementYDot;
}

//! Set Cartesian element: zDot.
void CartesianElements::setCartesianElementZDot( const double&
                                                 cartesianElementZDot )
{
    state_( 5 ) = cartesianElementZDot;
}

//! Set position vector.
void CartesianElements::setPositionVector (Vector3d &positionVector)
{
    state_.segment( 0, 3 ) = positionVector;
}

//! Set velocity vector.
void CartesianElements::setVelocityVector (Vector3d &velocityVector)
{
    state_.segment( 3, 3 ) = velocityVector;
}

//! Get Cartesian element: x.
double& CartesianElements::getCartesianElementX( )
{
    return state_( 0 );
}

//! Get Cartesian element: y.
double& CartesianElements::getCartesianElementY( )
{
    return state_( 1 );
}

//! Get Cartesian element: z.
double& CartesianElements::getCartesianElementZ( )
{
    return state_( 2 );
}

//! Get Cartesian element: xDot.
double& CartesianElements::getCartesianElementXDot( )
{
    return state_( 3 );
}

//! Get Cartesian element: yDot.
double& CartesianElements::getCartesianElementYDot( )
{
    return state_( 4 );
}

//! Get Cartesian element: zDot.
double& CartesianElements::getCartesianElementZDot( )
{
    return state_( 5 );
}

//! Get position vector.
Vector3d& CartesianElements::getPositionVector( )
{
    positionVector_ = state_.segment( 0, 3 );
    return positionVector_;
}

//! Get velocity vector.
Vector3d& CartesianElements::getVelocityVector( )
{
    velocityVector_ = state_.segment( 3, 3 );
    return velocityVector_;
}

// End of file.
