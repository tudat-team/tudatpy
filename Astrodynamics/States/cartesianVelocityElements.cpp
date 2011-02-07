/*! \file cartesianVelocityElements.cpp
 *    This source file contains the Cartesian velocity elements class included
 *    in Tudat.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@tudelft.nl
 *
 *    Date created      : 7 February, 2011
 *    Last modified     : 7 February, 2011
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
 *      110207    K. Kumar          First creation of code.
 */

// Include statements.
#include "cartesianVelocityElements.h"

//! Default constructor.
CartesianVelocityElements::CartesianVelocityElements( )
{
    // Initialize values.
    state.setZero( 3 );
}

//! Default destructor.
CartesianVelocityElements::~CartesianVelocityElements( )
{
}

//! Set Cartesian element: xDot.
void CartesianVelocityElements::setCartesianElementXDot( const double&
                                                         cartesianElementXDot )
{
    state( 0 ) = cartesianElementXDot;
}

//! Set Cartesian element: yDot.
void CartesianVelocityElements::setCartesianElementYDot( const double&
                                                         cartesianElementYDot )
{
    state( 1 ) = cartesianElementYDot;
}

//! Set Cartesian element: zDot.
void CartesianVelocityElements::setCartesianElementZDot( const double&
                                                         cartesianElementZDot )
{
    state( 2 ) = cartesianElementZDot;
}

//! Get Cartesian element: xDot.
double& CartesianVelocityElements::getCartesianElementXDot( )
{
    // Return Cartesian element: xDot.
    return state( 0 );
}

//! Get Cartesian element: yDot.
double& CartesianVelocityElements::getCartesianElementYDot( )
{
    // Return Cartesian element: yDot.
    return state( 1 );
}

//! Get Cartesian element: zDot.
double& CartesianVelocityElements::getCartesianElementZDot( )
{
    // Return Cartesian element: zDot.
    return state( 2 );
}


// End of file.
