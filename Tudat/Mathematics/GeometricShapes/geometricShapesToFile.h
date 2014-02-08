/*    Copyright (c) 2010-2014, Delft University of Technology
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
 *      YYMMDD    Author              Comment
 *      100914    K. Kumar            File created.
 *      100928    K. Kumar            Completed missing comments, changed
 *                                    writeIntegrationHistoryToFile( ) to
 *                                    writePropagationHistoryToFile( ).
 *      100929    B. Romgens          Spelling mistakes corrected and output to
 *                                    file corrected.
 *      100929    K. Kumar            Added checked code written by D. Dirkx.
 *      110202    K. Kumar            Updated writePropagationHistoryToFile( )
 *                                    to work with State*.
 *      120207    D. Dirkx            Split writingOutputToFile to separate free functions.
 *      120323    D. Dirkx            Changed raw pointers to shared pointers.
 *
 *    References
 *
 *    Notes
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
