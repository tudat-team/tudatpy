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

#ifndef TUDAT_GEOMETRIC_SHAPES_TO_FILE_H
#define TUDAT_GEOMETRIC_SHAPES_TO_FILE_H

#include <vector>

#include "Tudat/Mathematics/GeometricShapes/compositeSurfaceGeometry.h"
#include "Tudat/Mathematics/GeometricShapes/singleSurfaceGeometry.h"

namespace tudat
{
namespace geometric_shapes
{

//! Write single surface geometry to a file.
/*!
 * Writes the points on a SingleSurfaceGeometry object to a file.
 * Each row contains the x, y and z coordinate of the point, each next row
 * defines a single new point.
 * \param singleSurfaceGeometry Geometry which is to be written.
 * \param numberOfLines Defines how many points are taken over the 1st
 *          independent variable.
 * \param numberOfPoints Defines how many points are taken over the 2nd
 *          independent variable.
 * \param filename Name of the file to which the points are written.
 * \param writeType Defines whether to append or write to the file given by
 *          filename, should be "a" for append and "w" for write.
 * \param isIndependentVariableInverted Boolean flag which if set to true
 *          inverts which independent variable is treated as 1st and which
 *          as 2nd.
 */
void writeSingleSurfaceGeometryPointsToFile(
        geometric_shapes::SingleSurfaceGeometryPointer singleSurfaceGeometry,
        const int numberOfLines, const int numberOfPoints,
        const std::string& filename, const int writeType,
        const bool isIndependentVariableInverted );

//! Write composite surface geometry to a file.
/*!
 *  Writes the single surface geometries in a composite surface geometry to
 *  a file. The writeSingleGeometryPointsToFile( ) function is called for
 *  each surface geometry.
 *  \param compositeSurfaceGeometry Geometry from which there
 *          is to be written.
 *  \param arrayOfNumberOfLines Array of how many points to take over the 1st
 *          independent variables of single surface geometries.
 *  \param arrayOfNumberOfPoints Array of how many points to take over the 2nd
 *          independent variables of single surface geometries.
 *  \param filename Name of the file to which the points are written.
 *  \param writeType Defines whether to append or write to the file given
 *          by filename,  should be "a" for append and "w" for write.
 *  \param isIndependentVariableInvertedArray Array of booleans which if
 *          set to true invert which independent variable is treated as 1st
 *          and which as 2nd for each single surface geometry.
 */
void writeCompositeSurfaceGeometryPointsToFile(
        geometric_shapes::CompositeSurfaceGeometryPointer compositeSurfaceGeometry,
        const std::vector< int >& arrayOfNumberOfLines,
        const std::vector< int >& arrayOfNumberOfPoints,
        const std::string& filename, const int writeType,
        const std::vector< bool >& isIndependentVariableInvertedArray );

} // namespace geometric_shapes
} // namespace tudat

#endif // TUDAT_GEOMETRIC_SHAPES_TO_FILE_H
