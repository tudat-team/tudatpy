/*! \file surfaceGeometry.h
 *    This file contains the definition of the Surface Geometry base class.
 *
 *    Path              : /Mathematics/GeometricShapes/
 *    Version           : 7
 *    Check status      : Checked
 *
 *    Author            : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
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
 *    Last modified     : 5 September, 2011
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
 *      100910    D. Dirkx          First version of file.
 *      100915    D. Dirkx          Modified to correct comments, 80-lines
 *                                  rule, etc.
 *      100928    D. Dirkx          Modifications following first checking
 *                                  iteration.
 *      100929    D. Dirkx          Creation of separate file for class.
 *      110124    K. Kumar          Minor comment and layout changes; removed
 *                                  comment from "Notes"; updated Doxygen
 *`                                 comments.
 *      110204    K. Kumar          Minor comment and layout modifications.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 */

#ifndef SURFACEGEOMETRY_H
#define SURFACEGEOMETRY_H

// Include statements.
#include "Mathematics/GeometricShapes/geometricShape.h"

//! Surface geometry base class.
/*!
 * Base class for surface geometry representations in terms of two
 * parameterizing variables 1 and 2.
 */
class SurfaceGeometry : public GeometricShape
{
public:

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~SurfaceGeometry( ) { }

protected:

private:
};

#endif // SURFACEGEOMETRY_H

// End of file.
