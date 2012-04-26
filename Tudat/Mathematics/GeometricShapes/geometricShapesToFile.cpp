/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantability or fitness for a particular purpose.
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
 *      120207    D. Dirkx            Split writingOutputToFile to separate free functions
 *      120326    D. Dirkx            Changed raw pointers to shared pointers.
 *
 *    References
 *
 */

#include <fstream>

#include "Tudat/Mathematics/GeometricShapes/geometricShapesToFile.h"

namespace tudat
{
namespace output
{

using namespace tudat::mathematics::geometric_shapes;

//! Write single surface geometry to a file.
void writeSingleSurfaceGeometryPointsToFile(
    boost::shared_ptr< SingleSurfaceGeometry > pointerToSingleSurfaceGeometry,
    int numberOfLines, int numberOfPoints,
    const std::string& filename, int writeType, const bool& isIndependentVariableInverted )
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
        gridSize1_ = ( pointerToSingleSurfaceGeometry->getMaximumIndependentVariable( 1 )
                      - pointerToSingleSurfaceGeometry->getMinimumIndependentVariable( 1 ) )
                / ( numberOfLines - 1 );

        // Set grid size 2.
        gridSize2_ = ( pointerToSingleSurfaceGeometry->getMaximumIndependentVariable( 2 )
                      - pointerToSingleSurfaceGeometry->getMinimumIndependentVariable( 2 ) )
                / ( numberOfPoints - 1 );
    }

    // Sets the grid size in both directions, inverted
    // ( i.e., numberOfPoints corresponds to number of points in independent
    // variable 2 direction and numberOfLines corresponds to number of points
    // in indepedent variable 1 direction )
    else
    {
        // Set grid size 1.
        gridSize1_ = ( pointerToSingleSurfaceGeometry->getMaximumIndependentVariable( 1 )
                      - pointerToSingleSurfaceGeometry->getMinimumIndependentVariable( 1 ) )
                / ( numberOfPoints - 1 );

        // Set grid size 1.
        gridSize2_ = ( pointerToSingleSurfaceGeometry->getMaximumIndependentVariable( 2 )
                      - pointerToSingleSurfaceGeometry->getMinimumIndependentVariable( 2 ) )
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
                point = pointerToSingleSurfaceGeometry->getSurfacePoint(
                            pointerToSingleSurfaceGeometry->
                            getMinimumIndependentVariable( 1 ) + i * gridSize1_,
                            pointerToSingleSurfaceGeometry->
                            getMinimumIndependentVariable( 2 ) + j * gridSize2_ );
            }

            // If not inverted, all points from a single value of independent
            // variable 2 are written first.
            else
            {
                // Set point.
                point = pointerToSingleSurfaceGeometry->getSurfacePoint(
                            pointerToSingleSurfaceGeometry->
                            getMinimumIndependentVariable( 1 ) + j * gridSize1_,
                            pointerToSingleSurfaceGeometry->
                            getMinimumIndependentVariable( 2 ) + i * gridSize2_);
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
        boost::shared_ptr< CompositeSurfaceGeometry > pointerToCompositeSurfaceGeometry,
        std::vector< int > arrayOfNumberOfLines, std::vector< int > arrayOfNumberOfPoints,
        const std::string& filename, int writeType,
        std::vector< bool > isIndependentVariableInvertedArray )
{
    // Remove file of same name, if writeType is not append.
    if ( writeType == 0 )
    {
        std::remove( filename.c_str( ) );
    }

    // Iterate over all single geometries in composite geometry.
    for ( unsigned int i = 0; i < pointerToCompositeSurfaceGeometry->
          getNumberOfSingleSurfaceGeometries( ); i++ )
    {
        // Append single geometry to file.
        writeSingleSurfaceGeometryPointsToFile(
                    pointerToCompositeSurfaceGeometry->getSingleSurfaceGeometry( i ),
                    arrayOfNumberOfLines[ i ], arrayOfNumberOfPoints[ i ],
                    filename, 1, isIndependentVariableInvertedArray[ i ] );
    }
}

} // namespace output
} // namespace tudat
