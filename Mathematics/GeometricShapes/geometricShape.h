/*! \file geometricShape.h
 *    This file contains the definitions of the GeometricShape base class.
 *
 *    Path              : /Mathematics/GeometricShapes/
 *    Version           : 7
 *    Check status      : Checked
 *
 *    Author            : Dominic Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 29 September, 2010
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
 *      100910    D. Dirkx          First version of file.
 *      100915    D. Dirkx          Modified to correct comments, 80-lines
 *                                  rule, etc.
 *      100928    D. Dirkx          Modifications following first checking
 *                                  iteration.
 *      100929    D. Dirkx          Creation of separate file for class.
 *      110124    K. Kumar          Updated path; removed "Notes"; minor
 *                                  changes to comments and layout; added
 *                                  "End of line" comment.
 *      110204    K. Kumar          Minor comment and layout modifications.
 *      110207    D. Dirkx          Removed overloaded ostream operator.
 */

#ifndef GEOMETRICSHAPE_H
#define GEOMETRICSHAPE_H

// Include statements.
#include <cmath>
#include <cstdio>
#include <iostream>
#include "basicMathematicsFunctions.h"
#include "linearAlgebra.h"
#include "unitConversions.h"

// Using declarations.
using unit_conversions::convertRadiansToDegrees;

//! Geometry base class.
/*!
 *  Base class for geometric shapes.
 */
class GeometricShape
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    GeometricShape( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~GeometricShape( );

protected:

private:
};

#endif // GEOMETRICSHAPE_H

// End of file.
