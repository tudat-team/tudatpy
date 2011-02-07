/*! \file writingOutputToFile.cpp
 *    Header file that defines the class containing all funtionality pertaining
 *    to writing output to file included in Tudat.
 *
 *    Path              : /Astrodynamics/Output/
 *    Version           : 5
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author            : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Checker           : B. Romgens
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : bart.romgens@gmail.com
 *
 *    Date created      : 12 August, 2010
 *    Last modified     : 2 February, 2011
 *
 *    References
 *
 *    Notes
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
 *      110202    K. Kumar            Updated writePropagationHistoryToFile()
 *                                    to work with State*.
 */

// Include statements.
#include "writingOutputToFile.h"

// Ofstream object.
std::ofstream WritingOutputToFile::outputFile_;

//! Default constructor.
WritingOutputToFile::WritingOutputToFile( )
{
}

//! Default destructor.
WritingOutputToFile::~WritingOutputToFile( )
{
}

//! Write propagation history to file.
void WritingOutputToFile::writePropagationHistoryToFile(
        std::map < double, State* >& propagationHistory,
        const std::string& outputFilename )
{
    // Declare local variables.
    // Declare iterator for propagation history.
    std::map < double, State* >::iterator iteratorPropagationHistory_;

    // Open output file.
    outputFile_.open( outputFilename.c_str( ) );

    // Loop over map of propagation history.
    for ( iteratorPropagationHistory_ = propagationHistory.begin( );
          iteratorPropagationHistory_ != propagationHistory.end( );
          iteratorPropagationHistory_++ )
    {
        // Print map key to output file.
        outputFile_ << iteratorPropagationHistory_->first;

        // Loop over map data.
        for ( int i = 0;
              i < iteratorPropagationHistory_->second->state.rows( ); i++ )
        {
            // Print map data to file.
            outputFile_ << ", "
                        << iteratorPropagationHistory_->second->state[ i ];
        }

        // End line of output file.
        outputFile_ << std::endl;
    }

    // Close output file.
    outputFile_.close( );
}

//! Write single surface geometry to a file.
void WritingOutputToFile:: writeSingleSurfaceGeometryPointsToFile(
        SingleSurfaceGeometry* pointerToSingleSurfaceGeometry,
        const int& numberOfLines, const int& numberOfPoints,
        const std::string& filename, const int& writeType,
        const bool& isIndependentVariableInverted )
{
    // Declare local variables.
    // Declare grid sizes.
    double gridSize1_;
    double gridSize2_;

    // Declaration of vector which will be iteratively retrieved from geometry
    // and written to file.
    VectorXd point = VectorXd( 3 );

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
        gridSize1_ = ( pointerToSingleSurfaceGeometry->
                      getMaximumIndependentVariable( 1 )
                      - pointerToSingleSurfaceGeometry->
                      getMinimumIndependentVariable( 1 ) )
                    / ( numberOfLines - 1 );

        // Set grid size 2.
        gridSize2_ = ( pointerToSingleSurfaceGeometry->
                      getMaximumIndependentVariable( 2 )
                      - pointerToSingleSurfaceGeometry->
                      getMinimumIndependentVariable( 2 ) )
                    / ( numberOfPoints - 1 );
    }

    // Sets the grid size in both directions, inverted
    // ( i.e., numberOfPoints corresponds to number of points in independent
    // variable 2 direction and numberOfLines corresponds to number of points
    // in indepedent variable 1 direction )
    else
    {
        // Set grid size 1.
        gridSize1_ = ( pointerToSingleSurfaceGeometry->
                      getMaximumIndependentVariable( 1 )
                      - pointerToSingleSurfaceGeometry->
                      getMinimumIndependentVariable( 1 ) ) /
                      ( numberOfPoints - 1 );

        // Set grid size 1.
        gridSize2_ = ( pointerToSingleSurfaceGeometry->
                      getMaximumIndependentVariable( 2 )
                      - pointerToSingleSurfaceGeometry->
                      getMinimumIndependentVariable( 2 ) ) /
                     ( numberOfLines - 1 );
    }

    // Iterate over all points, and retrieve the surface point from the
    // surfaceGeometry object, and write it to the file.
    for( int i = 0; i < numberOfLines ; i++ )
    {
        for( int j = 0; j < numberOfPoints ; j++ )
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
            outputFile_ << i+1 << " ";
            outputFile_ << j+1 << " ";
            outputFile_ << point( 0 ) << " ";
            outputFile_ << point( 1 ) << " ";
            outputFile_ << point( 2 ) << " ";
            outputFile_ << std::endl;
        }
    }

    // Closing output file.
    outputFile_.close( );
}

//! Write composite surface geometry to a file.
void WritingOutputToFile:: writeCompositeSurfaceGeometryPointsToFile(
        CompositeSurfaceGeometry* pointerToCompositeSurfaceGeometry,
        int* arrayOfNumberOfLines, int* arrayOfNumberOfPoints,
        const std::string& filename, const int& writeType,
        bool* isIndependentVariableInvertedArray )
{
    // Remove file of same name, if writeType is not append.
    if ( writeType == 0 )
    {
        remove( filename.c_str( ) );
    }

    // Iterate over all single geometries in composite geometry.
    for ( unsigned int i = 0; i < pointerToCompositeSurfaceGeometry->
          getNumberOfSingleSurfaceGeometries( ); i++ )
    {
        // Append single geometry to file.
        writeSingleSurfaceGeometryPointsToFile(
                pointerToCompositeSurfaceGeometry
                ->getSingleSurfaceGeometry( i ),
                arrayOfNumberOfLines[ i ], arrayOfNumberOfPoints[ i ],
                filename, 1, isIndependentVariableInvertedArray[ i ] );
    }
}

// End of file.
