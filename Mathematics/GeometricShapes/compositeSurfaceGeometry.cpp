/*! \file compositeSurfaceGeometry.cpp
 *    This file contains the definition of the CompositeSurfaceGeometry class
 *
 *    Path              : /Mathematics/GeometricShapes/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author            : Dominic Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 25 November, 2010
 *    Last modified     : 4 February, 2011
 *
 *    References
 *
 *    Notes
 *      The code uses pointers to pointers to represent an array of pointers,
 *      instead of the <vector> or <map> from the STL for heritage reasons. In
 *      future, it might be fruitful to consider the use of standard STL
 *      containers instead.
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
 *      102511    D. Dirkx          First version of file.
 *      110119    K. Kumar          Minor comments changes; path updated;
 *                                  Doxygen comments updated; updated function
 *                                  arguments to use references and unsigned
 *                                  ints; added "End of file" comment; minor
 *                                  changes to layout.
 *      110204    K. Kumar          Minor comment and layout modifications;
 *                                  corrected Doxygen comments.
 */

// Include statements.
#include "compositeSurfaceGeometry.h"

//! Default constructor.
CompositeSurfaceGeometry::CompositeSurfaceGeometry( )
    : numberOfSingleSurfaceGeometries_( 0 ),
      numberOfCompositeSurfaceGeometries_( 0 ),
      singleSurfaceGeometryList_( NULL ),
      compositeSurfaceGeometryList_( NULL )
{
}

//! Default destructor.
CompositeSurfaceGeometry::~CompositeSurfaceGeometry( )
{
    // Delete single surface geometry list.
    delete [ ] singleSurfaceGeometryList_;
    singleSurfaceGeometryList_ = NULL;

    // Delete composite surface geometry list.
    delete [ ] compositeSurfaceGeometryList_ ;
    compositeSurfaceGeometryList_ = NULL;
}

//! Set pointer to SingleSurfaceGeometry object.
void CompositeSurfaceGeometry::setSingleSurfaceGeometry(
    SingleSurfaceGeometry* pointerToSingleSurfaceGeometry,
    const unsigned int& index )
{
    // Put surface in given index in list.
    singleSurfaceGeometryList_[ index ] = pointerToSingleSurfaceGeometry;
}

//! Set pointer to a CompositeSurfaceGeometry object.
void CompositeSurfaceGeometry::setCompositeSurfaceGeometry(
    CompositeSurfaceGeometry* pointerToCompositeSurfaceGeometry,
    const unsigned int& index )
{
    // Put surface in given index in list.
    compositeSurfaceGeometryList_[ index ] = pointerToCompositeSurfaceGeometry;
}

//! Set number of single surface geometries.
void CompositeSurfaceGeometry::setNumberOfSingleSurfaceGeometries(
    const unsigned int& numberOfSingleSurfaceGeometries )
{
    numberOfSingleSurfaceGeometries_ = numberOfSingleSurfaceGeometries;

    // Allocate memory for singleSurfaceGeometryList_.
    singleSurfaceGeometryList_
            = new SingleSurfaceGeometry*[ numberOfSingleSurfaceGeometries_ ];
}

//! Set number of composite surface geometries.
void CompositeSurfaceGeometry::setNumberOfCompositeSurfaceGeometries(
    const unsigned int& numberOfCompositeSurfaceGeometries )
{
    numberOfCompositeSurfaceGeometries_ = numberOfCompositeSurfaceGeometries;

    // Allocate memory for compositeSurfaceGeometryList_.
    compositeSurfaceGeometryList_
            = new CompositeSurfaceGeometry*[
                    numberOfCompositeSurfaceGeometries_ ];
}

//! Get pointer to stored SingleSurfaceGeometry object.
SingleSurfaceGeometry* CompositeSurfaceGeometry::getSingleSurfaceGeometry(
    const unsigned int& index )
{
    // Return surface from given index in list.
    return singleSurfaceGeometryList_[ index ];
}

//! Get pointer to a stored CompositeSurfaceGeometry object.
CompositeSurfaceGeometry* CompositeSurfaceGeometry::getCompositeSurfaceGeometry(
    const unsigned int& index )
{
    // Return surface from given index in list.
    return compositeSurfaceGeometryList_[index];
}

//! Get number of single surface geometries.
unsigned int& CompositeSurfaceGeometry::getNumberOfSingleSurfaceGeometries( )
{
    return numberOfSingleSurfaceGeometries_;
}

//! Get number of composite surface geometries.
unsigned int& CompositeSurfaceGeometry::getNumberOfCompositeSurfaceGeometries( )
{
    return numberOfCompositeSurfaceGeometries_;
}

//! Overload ostream to print class information.
std::ostream &operator<<( std::ostream &stream,
                          CompositeSurfaceGeometry& compositeSurfaceGeometry )
{
    stream << "This is a composite surface geometry." << std::endl;
    stream << "The number of SingleSurfaceGeometries is: "
           << compositeSurfaceGeometry.numberOfSingleSurfaceGeometries_
           << std::endl;
    stream << "The number of CompositeSurfaceGeometries is: "
           << compositeSurfaceGeometry.numberOfCompositeSurfaceGeometries_
           << std::endl;

    // Return stream.
    return stream;
}

// End of file.
