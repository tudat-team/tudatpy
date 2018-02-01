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

#include <fstream>

#include "Tudat/Mathematics/GeometricShapes/geometricShapesToFile.h"

namespace tudat
{
namespace geometric_shapes
{

//! Write single surface geometry to a file.
void writeSingleSurfaceGeometryPointsToFile(
        geometric_shapes::SingleSurfaceGeometryPointer singleSurfaceGeometry,
        const int numberOfLines, const int numberOfPoints,
        const std::string& filename, const int writeType,
        const bool isIndependentVariableInverted )
{
    std::ofstream outputFile_;

    // Declare local variables.
    // Declare grid sizes.
    double gridSize1_;
    double gridSize2_;

    // Declaration of vector which will be iteratively retrieved from geometry
    // and written to file.
    Eigen::VectorXd point = Eigen::VectorXd( 3 );

    // Open the file to which writing will take place if it is to overwrite any
    // existing content.
    if ( writeType == 0 )
    {
        // Open output file.
        outputFile_.open( filename.c_str( ) );
    }

    // Open the file to which writing will take place if it is to append to
    // existing content.
    else if ( writeType == 1 )
    {
        // Open output file with append option.
        outputFile_.open( filename.c_str( ), std::ios::app );
    }


    // Sets the grid size in both directions.
    if ( isIndependentVariableInverted  == false )
    {
        // Set grid size 1.
        gridSize1_ = ( singleSurfaceGeometry->getMaximumIndependentVariable( 1 )
                      - singleSurfaceGeometry->getMinimumIndependentVariable( 1 ) )
                / ( numberOfLines - 1 );

        // Set grid size 2.
        gridSize2_ = ( singleSurfaceGeometry->getMaximumIndependentVariable( 2 )
                      - singleSurfaceGeometry->getMinimumIndependentVariable( 2 ) )
                / ( numberOfPoints - 1 );
    }

    // Sets the grid size in both directions, inverted
    // ( i.e., numberOfPoints corresponds to number of points in independent
    // variable 2 direction and numberOfLines corresponds to number of points
    // in indepedent variable 1 direction )
    else
    {
        // Set grid size 1.
        gridSize1_ = ( singleSurfaceGeometry->getMaximumIndependentVariable( 1 )
                      - singleSurfaceGeometry->getMinimumIndependentVariable( 1 ) )
                / ( numberOfPoints - 1 );

        // Set grid size 1.
        gridSize2_ = ( singleSurfaceGeometry->getMaximumIndependentVariable( 2 )
                      - singleSurfaceGeometry->getMinimumIndependentVariable( 2 ) )
                / ( numberOfLines - 1 );
    }

    // Iterate over all points, and retrieve the surface point from the
    // surfaceGeometry object, and write it to the file.
    for ( int i = 0; i < numberOfLines ; i++ )
    {
        for ( int j = 0; j < numberOfPoints ; j++ )
        {
            // If not inverted, all points from a single value of independent
            // variable 2 are written first.
            if ( isIndependentVariableInverted  == false )
            {
                // Set point.
                point = singleSurfaceGeometry->getSurfacePoint(
                            singleSurfaceGeometry->getMinimumIndependentVariable( 1 )
                            + i * gridSize1_,
                            singleSurfaceGeometry->getMinimumIndependentVariable( 2 )
                            + j * gridSize2_ );
            }

            // If not inverted, all points from a single value of independent
            // variable 2 are written first.
            else
            {
                // Set point.
                point = singleSurfaceGeometry->getSurfacePoint(
                            singleSurfaceGeometry->getMinimumIndependentVariable( 1 )
                            + j * gridSize1_,
                            singleSurfaceGeometry->getMinimumIndependentVariable( 2 )
                            + i * gridSize2_);
            }

            // Write the x-, y- and z-coordinates, followed by a next line command.
            outputFile_ << i+1 << " " << j+1 << " " << point( 0 ) << " " << point( 1 ) << " ";
            outputFile_ << point( 2 ) << " " << std::endl;
        }
    }

    // Closing output file.
    outputFile_.close( );
}

//! Write composite surface geometry to a file.
void writeCompositeSurfaceGeometryPointsToFile(
        geometric_shapes::CompositeSurfaceGeometryPointer compositeSurfaceGeometry,
        const std::vector< int >& arrayOfNumberOfLines,
        const std::vector< int >& arrayOfNumberOfPoints,
        const std::string& filename, const int writeType,
        const std::vector< bool >& isIndependentVariableInvertedArray )
{
    // Remove file of same name, if writeType is not append.
    if ( writeType == 0 )
    {
        std::remove( filename.c_str( ) );
    }

    // Iterate over all single geometries in composite geometry.
    for ( unsigned int i = 0; i < compositeSurfaceGeometry->
          getNumberOfSingleSurfaceGeometries( ); i++ )
    {
        // Append single geometry to file.
        writeSingleSurfaceGeometryPointsToFile(
                    compositeSurfaceGeometry->getSingleSurfaceGeometry( i ),
                    arrayOfNumberOfLines[ i ], arrayOfNumberOfPoints[ i ],
                    filename, 1, isIndependentVariableInvertedArray[ i ] );
    }
}

} // namespace geometric_shapes
} // namespace tudat
