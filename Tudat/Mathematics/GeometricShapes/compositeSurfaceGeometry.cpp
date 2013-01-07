/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
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
 *      110810    J. Leloux         Corrected doxygen documentation.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120323    D. Dirkx          Removed set functions; moved functionality to constructor
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *
 *    References
 *
 */

#include "Tudat/Mathematics/GeometricShapes/compositeSurfaceGeometry.h"

namespace tudat
{
namespace geometric_shapes
{

CompositeSurfaceGeometry::CompositeSurfaceGeometry(
                          std::vector< boost::shared_ptr< SingleSurfaceGeometry > >
                          singleSurfaceGeometryList,
                          std::vector< boost::shared_ptr< CompositeSurfaceGeometry > >
                          compositeSurfaceGeometryList )
{
    setNumberOfSingleSurfaceGeometries( singleSurfaceGeometryList.size( ) );
    setNumberOfCompositeSurfaceGeometries( compositeSurfaceGeometryList.size( ) );

    for ( unsigned int i = 0; i < numberOfSingleSurfaceGeometries_; i++ )
    {
        setSingleSurfaceGeometry( singleSurfaceGeometryList[ i ], i );
    }

    for ( unsigned int i = 0; i < numberOfCompositeSurfaceGeometries_; i++ )
    {
        setCompositeSurfaceGeometry( compositeSurfaceGeometryList[ i ], i );
    }
}

//! Overload ostream to print class information.
std::ostream &operator<<( std::ostream &stream,
                          CompositeSurfaceGeometry& compositeSurfaceGeometry )
{
    stream << "This is a composite surface geometry." << std::endl;
    stream << "The number of SingleSurfaceGeometries is: "
           << compositeSurfaceGeometry.numberOfSingleSurfaceGeometries_ << std::endl;
    stream << "The number of CompositeSurfaceGeometries is: "
           << compositeSurfaceGeometry.numberOfCompositeSurfaceGeometries_ << std::endl;

    // Return stream.
    return stream;
}


} // namespace geometric_shapes
} // namespace tudat

