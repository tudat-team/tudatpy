/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      100902    K. Kumar          File header and footer added.
 *      100916    D. Dirkx          Added minor comments and placeholder tag during checking.
 *      100928    K. Kumar          Added reference and adjusted include statements.
 *      100929    K. Kumar          Changed EigenRoutines.h include statement
 *                                  to linearAlgebra.h and removed placeholder comment.
 *                                  Added small comment modifications.
 *      110202    K. Kumar          Added overload for map with State* for
 *                                  computeNearestLeftNeighborUsingBinarySearch( ).
 *      110803    J. Leloux         Added convertStringToTemplate.
 *      110805    J. Leloux         Added outputCurrentRunningTime.
 *      110810    J. Leloux         Minor comment modifications.
 *      110913    K. Kumar          Implemented automatic root-path functions based on
 *                                  suggestions by M. Persson.
 *      111117    K. Kumar          Added listAllFilesInDirectory( ) function.
 *      120127    K. Kumar          Adapted for Tudat Core.
 *      120712    K. Kumar          Updated use of filesystem3 boost-namespace to filesystem.
 *      130110    K. Kumar          Added function that allows formatting of string representation
 *                                  of floating-point number in scientific notation.
 *
 *    References
 *
 *    Notes
 *      The function printStandardScientificNotation() has been implemented to cope with
 *      cross-platform incompatibilities in the printed output of floating-point numbers in
 *      scientific notation.
 *
 */
 
#include <iomanip>
#include <sstream>
#include <string>

#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{
namespace input_output
{

//! Print floating-point number in formatted scientific notation.
std::string printInFormattedScientificNotation( const double floatingPointNumber,
                                                const int basePrecision,
                                                const int exponentWidth )
{
    // Declare string representation of the floating-point number.
    std::string floatingPointNumberString;

    // Extract the decimal part (base) of the floating-point number with specified precision.
    // Default precision is numeric_limits< double >::digits10.
    std::stringstream buffer;
    buffer << std::scientific
           << std::setprecision( basePrecision )
           << std::uppercase
           << floatingPointNumber
           << std::endl;

    // Write base to string.
    buffer >> floatingPointNumberString;

    // Determine the location where the exponent part of the number starts.
    std::string::size_type exponentLocation = floatingPointNumberString.find( 'E' );

    // Set the exponent to zero.
    int exponent = 0;

    // Set flag indicating if exponent is positive or negative to positive.
    bool isExponentPositive = true;

    // Check if the exponent is present.
    if ( exponentLocation != std::string::npos )
    {
        // Extract the exponent to the buffer.
        buffer << floatingPointNumberString.substr( exponentLocation + 1,
                                                    floatingPointNumberString.length( )
                                                    - exponentLocation - 1 )
               << std::endl;

        // Write the exponent part from the buffer.
        buffer >> exponent;

        // Check if the exponent is negative.
        if ( exponent < 0 )
        {
            // If it is negative, set the flag to false, and switch signs of the exponent.
            isExponentPositive = false;
            exponent = -exponent;
        }
    }

    // Extract the exponent with a fixed field width.
    // Default width is 2.
    std::stringstream outputStream;
    outputStream << floatingPointNumberString.substr( 0, exponentLocation + 1 )
                 << ( isExponentPositive ? '+' : '-' )
                 << std::setw( exponentWidth )
                 << std::setfill( '0' )
                 << exponent;

    // Declare output string and write buffer output to the string.
    std::string outputString;
    outputStream >> outputString;

    // Return the formatted string representation of the floating-point number.
    return outputString;
}

//! Lists all files in directory.
std::vector< boost::filesystem::path > listAllFilesInDirectory(
    const boost::filesystem::path& directory, const bool isRecurseIntoSubdirectories )
{
    // Declare local variables.
    std::vector < boost::filesystem::path > listOfFileNamesWithPath_;

    if ( boost::filesystem::exists( directory ) )
    {
        boost::filesystem::directory_iterator iteratorPastEndOfDirectory_;

        for ( boost::filesystem::directory_iterator directoryIterator_( directory );
              directoryIterator_ != iteratorPastEndOfDirectory_ ; ++directoryIterator_ )
        {
            if ( boost::filesystem::is_directory( *directoryIterator_ ) )
            {
                if ( isRecurseIntoSubdirectories )
                {
                    listAllFilesInDirectory( *directoryIterator_ );
                }
            }

            else
            {
                listOfFileNamesWithPath_.push_back( directoryIterator_->path( ).filename( ) );
            }
        }
    }

    // Return container of filenames.
    return listOfFileNamesWithPath_;
}

void writeMatrixToFile( Eigen::MatrixXd matrixToWrite,
                        std::string fileName,
                        const int numberOfDigits  )
{
    std::ofstream outputFile_;

    // Open output file.
    outputFile_.open( fileName.c_str( ) );

    for( int i = 0; i < matrixToWrite.rows( ); i++ )
    {
        for ( int j = 0;   j < matrixToWrite.cols( ); j++ )
        {
            // Print map data to file.
            outputFile_.precision( numberOfDigits );
            outputFile_ << matrixToWrite( i, j );
            if( j != matrixToWrite.cols( ) - 1 )
            {
                outputFile_ <<" ";
            }
        }
        outputFile_<<std::endl;
    }

    // Close output file.
    outputFile_.close( );
}

} // namespace input_output
} // namespace tudat
