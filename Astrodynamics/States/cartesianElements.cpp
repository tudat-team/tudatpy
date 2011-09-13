/*! \file cartesianElements.cpp
 *    This source file contains the Cartesian elements class included in Tudat.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 6
 *    Check status      : Checked
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 22 October, 2010
 *    Last modified     : 4 February, 2010
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      101022    K. Kumar          First creation of code.
 *      101026    K. Kumar          Added position, velocity vector functions.
 *      101201    E. Iorfida        Modified punctuation.
 *      101202    J. Melman         Large modification: only set state_ and use
 *                                  state_ to get components. Does not use
 *                                  private variables, so those are obsolete.
 *      110110    K. Kumar          Minor modifications; initialised values
 *                                  using constructor.
 *      110204    K. Kumar          Removed "vector" from naming.
 */

// Include statements.
#include "cartesianElements.h"

// Using declarations.
using std::endl;

//! Default constructor.
CartesianElements::CartesianElements( )
{
    // Initialize values.
    position_.setConstant( -0.0 );

    // Initialize variables.
    state.setZero( 6 );
}

//! Default destructor.
CartesianElements::~CartesianElements( )
{
}

//! Set Cartesian element: xDot.
void CartesianElements::setCartesianElementXDot( const double&
                                                 cartesianElementXDot )
{
    state( 3 ) = cartesianElementXDot;
}

//! Set Cartesian element: yDot.
void CartesianElements::setCartesianElementYDot( const double&
                                                 cartesianElementYDot )
{
    state( 4 ) = cartesianElementYDot;
}

//! Set Cartesian element: zDot.
void CartesianElements::setCartesianElementZDot( const double&
                                                 cartesianElementZDot )
{
    state( 5 ) = cartesianElementZDot;
}

//! Set position.
void CartesianElements::setPosition( Vector3d& position )
{
    state.segment( 0, 3 ) = position;
}

//! Set velocity.
void CartesianElements::setVelocity( Vector3d& velocity )
{
    state.segment( 3, 3 ) = velocity;
}

//! Get Cartesian element: xDot.
double& CartesianElements::getCartesianElementXDot( )
{
    // Return Cartesian element: xDot.
    return state( 3 );
}

//! Get Cartesian element: yDot.
double& CartesianElements::getCartesianElementYDot( )
{
    // Return Cartesian element: yDot.
    return state( 4 );
}

//! Get Cartesian element: zDot.
double& CartesianElements::getCartesianElementZDot( )
{
    // Return Cartesian element: zDot.
    return state( 5 );
}

//! Get position.
Vector3d& CartesianElements::getPosition( )
{
    // Return position.
    position_ = state.segment( 0, 3 );
    return position_;
}

//! Get velocity.
Vector3d& CartesianElements::getVelocity( )
{
    // Return velocity.
    velocity_ = state.segment( 3, 3 );
    return velocity_;
}

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream,
                          CartesianElements& cartesianElements )
{
    stream << "The state is set to: " << cartesianElements.state << endl;

    // Return stream.
    return stream;
}

// End of file.
