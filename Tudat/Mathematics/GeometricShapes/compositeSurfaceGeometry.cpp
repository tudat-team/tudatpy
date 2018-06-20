/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "Tudat/Mathematics/GeometricShapes/compositeSurfaceGeometry.h"

namespace tudat
{
namespace geometric_shapes
{

CompositeSurfaceGeometry::CompositeSurfaceGeometry(
                          std::vector< std::shared_ptr< SingleSurfaceGeometry > >
                          singleSurfaceGeometryList,
                          std::vector< std::shared_ptr< CompositeSurfaceGeometry > >
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
std::ostream &operator << ( std::ostream &stream,
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

