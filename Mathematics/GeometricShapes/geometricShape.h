/*! \file geometricShape.h
 *  This file contains the definitions of the GeometricShape
 *  base class.
 *
 *  Path              : Mathematics/Geometry/
 *  Version           : 1
 *  Check status      : Checked
 *
 *  Author            : Dominic Dirkx
 *  Affiliation       : TU Delft
 *  E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *  Checker           : J. Melman
 *  Affiliation       : Delft University of Technology
 *  E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *  Date created      : 29 September, 2010
 *  Last modified     : 29 September, 2010
 *
 *  References
 *
 *  Notes
 *    The setParameter and getParameter functions could benefit from
 *    taking an enum as input instead of an int.
 *
 *    The original file was split into several parts, so that each class is
 *    defined separately in a file. Hence, the this file was created after the
 *    first entries in the changelog and the file version is lower than the
 *    number of entries in changelog.
 *
 *  Copyright (c) 2010 Delft University of Technology.
 *
 *  This software is protected by national and international copyright.
 *  Any unauthorized use, reproduction or modification is unlawful and
 *  will be prosecuted. Commercial and non-private application of the
 *  software in any form is strictly prohibited unless otherwise granted
 *  by the authors.
 *
 *  The code is provided without any warranty; without even the implied
 *  warranty of merchantibility or fitness for a particular purpose.
 *
 *  Changelog
 *    100910   D. Dirkx                    First version of file
 *    100915   D. Dirkx                    Modified to correct comments, 80
 *                                         lines rule, etc.
 *    100928   D. Dirkx                    Modifications following first
 *                                         checking iteration.
 *    100929   D. Dirkx                    Creation of separate file for class
 *
 */

#ifndef GEOMETRICSHAPE_H
#define GEOMETRICSHAPE_H

#include "linearAlgebra.h"
#include "basicMathematicsFunctions.h"
#include "unitConversions.h"
#include <iostream>
#include <cstdio>
#include <cmath>

using unit_conversions::convertRadiansToDegrees;

//! Geometric shape base class.
/*!
 *  Base class for geometric shapes.
 */
class GeometricShape
{
public:

    //! Default destructor
    /*!
     *  Default destructor
     */
    GeometricShape( ) ;

    //! Default constructor
    /*!
     *  Default constructor
     */
    virtual ~GeometricShape( ) ;
};

#endif // GEOMETRICSHAPE_H
