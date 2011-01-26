/*! \file body.cpp
 *    This file contains a class that describes planets
 *    with all their (derived) characteristics.
 *
 *    Path              : /Astrodynamics/Bodies/
 *    Version           : 2
 *    Check status      : Checked
 *
 *    Author            : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 10 September, 2010
 *    Last modified     : 15 January, 2011
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
 *      100910    J. Melman     First creation of code.
 *      100929    K. Kumar      Minor comment changes and Body scope for
 *                              setShapeModel( ) function added.
 *      110115    J. Melman     Added set and get shape model functions.
 */

// Include statements.
#include "body.h"

//! Default constructor.
Body::Body( )
{
}

//! Sets the shape model.
void Body::setShapeModel( GeometricShape* bodyGeometricShape )
{
    bodyGeometricShape_ = bodyGeometricShape;
}

//! Get the shape model.
GeometricShape* Body::getShapeModel( )
{
    return bodyGeometricShape_;
}

// End of file.

