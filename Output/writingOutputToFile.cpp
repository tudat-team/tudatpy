/*! \file writingOutputToFile.h
 *    Header file that defines the class containing all funtionality pertaining
 *    to writing output to file included in Tudat.
 *
 *    Path              : /Astrodynamics/Output/
 *    Version           : 4
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
 *    Date created      : 12 august, 2010
 *    Last modified     : 29 september, 2010
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
 *      YYMMDD    author              comment
 *      100914    K. Kumar            File created.
 *      100928    K. Kumar            Completed missing comments, changed
 *                                    writeIntegrationHistoryToFile( ) to
 *                                    writePropagationHistoryToFile( ).
 *      100929    B. Romgens          Spelling mistakes corrected and output to
 *                                    file corrected.
 *      100929    K. Kumar            Added checked code written by D. Dirkx.
 *
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
        std::map < double, VectorXd >& propagationHistory,
        const std::string& outputFilename )
{
    // Declare local variables.
    // Declare iterator for propagation history.
    std::map < double, VectorXd >::iterator iteratorPropagationHistory_;

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
        for ( int i = 0; i < iteratorPropagationHistory_->second.rows( ); i++ )
        {
            // Print map data to file.
            outputFile_ << ", " << iteratorPropagationHistory_->second[ i ];
        }

        // End line of output file.
        outputFile_ << std::endl;
    }

    // Close output file.
    outputFile_.close( );
}

//! Function to write a surface geometry to a file.
void WritingOutputToFile::
        writeGeometryPointsToFile( SurfaceGeometry* geometry,
                                   int numberOfLines, int numberOfPoints,
                                   const std::string& filename, int writeType,
                                   bool isInvertIndependentVariable)
{
    // Declare local variables.
    // Declare grid sizes.
    double gridSize1;
    double gridSize2;

    // Declaration of vector which will be iteratively retrieved from geometry
    // and written to file.
    VectorXd point = VectorXd(3);

    // Open the file to which writing will take place if it is to overwrite any
    // existing content.
    if ( writeType== 0 )
    {
        // Open output file.
        outputFile_.open(filename.c_str( ));
    }

    // Open the file to which writing will take place if it is to append to
    // existing content
    else if ( writeType == 1 )
    {
        // Open output file with append option.
        outputFile_.open(filename.c_str( ), std::ios::app );
    }


    // Sets the grid size in both directions.
    if ( isInvertIndependentVariable  == false )
    {
        // Set grid size 1.
        gridSize1 = ( geometry->getMaximumIndependentVariable( 1 )
                      - geometry->getMinimumIndependentVariable( 1 ) )
                    / ( numberOfLines - 1 );

        // Set grid size 2.
        gridSize2 = ( geometry->getMaximumIndependentVariable( 2 )
                      - geometry->getMinimumIndependentVariable( 2 ))
                    / ( numberOfPoints - 1 );
    }

    // Sets the grid size in both directions, inverted
    // (i.e., numberOfPoints corresponds to number of points in independent
    // variable 2 direction and numberOfLines corresponds to number of points
    // in indepedent variable 1 direction)
    else
    {
        // Set grid size 1.
        gridSize1 = ( geometry->getMaximumIndependentVariable( 1 )
                      - geometry->getMinimumIndependentVariable( 1 ) )
                    / ( numberOfPoints - 1 );

        // Set grid size 1.
        gridSize2 = ( geometry->getMaximumIndependentVariable( 2 )
                      - geometry->getMinimumIndependentVariable( 2 ) )
                    / ( numberOfLines - 1 );
    }

    // Iterate over all points, and retrieve the surface point from the
    // surfaceGeometry object, and write it to the file.
    for( int i = 0; i < numberOfLines ; i++ )
    {
        for( int j = 0; j < numberOfPoints ; j++ )
        {
            // If not inverted, all points from a single value of independent
            // variable 2 are written first.
            if ( isInvertIndependentVariable  == false )
            {
                // Set point.
                point = geometry->getSurfacePoint(
                        geometry->getMinimumIndependentVariable( 1 )
                        + i * gridSize1,
                        geometry->getMinimumIndependentVariable( 2 )
                        + j * gridSize2 );
            }

            // If not inverted, all points from a single value of independent
            // variable 2 are written first.
            else
            {
                // Set point.
                point = geometry->getSurfacePoint(
                        geometry->getMinimumIndependentVariable( 1 )
                        + j * gridSize1,
                        geometry->getMinimumIndependentVariable( 2 )
                        + i * gridSize2);
            }

            // Write the x-, y- and z-coordinates, followed by a next line command.
            outputFile_ << point( 0 ) << " ";
            outputFile_ << point( 1 ) << " ";
            outputFile_ << point( 2 ) << " ";
            outputFile_ << std::endl;
        }
    }

    // Closing output file.
    outputFile_.close( );
}

// End of file.
