/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
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

#ifndef TUDAT_SURFACEGEOMETRY_H
#define TUDAT_SURFACEGEOMETRY_H

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

} // namespace tudat

#endif // TUDAT_SURFACEGEOMETRY_H
