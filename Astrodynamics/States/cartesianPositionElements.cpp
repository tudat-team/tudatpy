/*! \file cartesianPositionElements.cpp
 *    This source file contains the Cartesian position elements class included
 *    in Tudat.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 2
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
 *    Date created      : 31 January, 2011
 *    Last modified     : 4 February, 2011
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
 *      110131    K. Kumar          First creation of code.
 *      110204    K. Kumar          Removed positionVector_, velocityVector_.
 */

// Include statements.
#include "cartesianPositionElements.h"

// Using declarations.
using std::endl;

//! Default constructor.
CartesianPositionElements::CartesianPositionElements( )
{
    // Initialize values.
    state.setZero( 3 );
}

//! Default destructor.
CartesianPositionElements::~CartesianPositionElements( )
{
}

//! Set Cartesian element: x.
void CartesianPositionElements::setCartesianElementX( const double&
                                                      cartesianElementX )
{
    state( 0 ) = cartesianElementX;
}

//! Set Cartesian element: y.
void CartesianPositionElements::setCartesianElementY( const double&
                                                      cartesianElementY )
{
    state( 1 ) = cartesianElementY;
}

//! Set Cartesian element: z.
void CartesianPositionElements::setCartesianElementZ( const double&
                                                      cartesianElementZ )
{
    state( 2 ) = cartesianElementZ;
}

//! Get Cartesian element: x.
double& CartesianPositionElements::getCartesianElementX( )
{
    // Return Cartesian element: x.
    return state( 0 );
}

//! Get Cartesian element: y.
double& CartesianPositionElements::getCartesianElementY( )
{
    // Return Cartesian element: y.
    return state( 1 );
}

//! Get Cartesian element: z.
double& CartesianPositionElements::getCartesianElementZ( )
{
    // Return Cartesian element: z.
    return state( 2 );
}

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream,
                          CartesianPositionElements&
                          cartesianPositionElements )
{
    stream << "The state is set to: " << cartesianPositionElements.state
           << endl;

    // Return stream.
    return stream;
}

// End of file.
