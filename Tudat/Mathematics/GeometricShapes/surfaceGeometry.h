/*!   Copyright (c) 2010-2012 Delft University of Technology.
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
 *
 *    References
 *
 */

// Temporary notes (move to class/function doxygen):
// 
// `                                 comments.
// 110204    K. Kumar          Minor comment and layout modifications.
// 110905    S. Billemont      Reorganized includes.
// Moved (con/de)structors and getter/setters to header.

#ifndef SURFACEGEOMETRY_H
#define SURFACEGEOMETRY_H

#include <iostream>

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
namespace tudat
{

//! Surface geometry base class.
/*!
 * Base class for surface geometry representations in terms of two
 * parameterizing variables 1 and 2.
 */
class SurfaceGeometry
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

}

#endif // SURFACEGEOMETRY_H

// End of file.
