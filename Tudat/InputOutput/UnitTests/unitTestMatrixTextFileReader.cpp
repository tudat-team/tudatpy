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
 *      110530    F.M. Engelen      Creation of code.
 *
 *    References
 *
 *    If tabs are used as spaces, it doesn't work. The seperator should also be tabs then.
 *
 */

#include <cmath>
#include <iostream>
#include <limits>

#include <Eigen/Core>

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/matrixTextFileReader.h"

//! Test matrix text file reader.
/*!
 * Tests the matrix text file reader.
 */
int main( )
{
    using std::cerr;
    using std::endl;

    bool isTestMatrixBroken = false;

    Eigen::MatrixXd readMatrix;

    Eigen::MatrixXd expectedMatrix( 4,3 );
    expectedMatrix( 0,0 ) = 1.0;
    expectedMatrix( 0,1 ) = 2.0;
    expectedMatrix( 0,2 ) = 3.0;
    expectedMatrix( 1,0 ) = 4.0;
    expectedMatrix( 1,1 ) = 5.0;
    expectedMatrix( 1,2 ) = 6.0;
    expectedMatrix( 2,0 ) = 7.0;
    expectedMatrix( 2,1 ) = 8.0;
    expectedMatrix( 2,2 ) = 9.0;
    expectedMatrix( 3,0 ) = 10.0;
    expectedMatrix( 3,1 ) = 11.0;
    expectedMatrix( 3,2 ) = 12.0;

    // Test 1: Test for semi-colon-separated files.
    readMatrix = tudat::input_output::readMatrixFromFile( tudat::input_output::getTudatRootPath( )
                                                          + "InputOutput/UnitTests/testMatrix.txt",
                                                          ";"  );

    for ( int i = 0; i < 3; i++ )
    {
        for ( int j = 0; j < 2; j++ )
        {
            if ( std::fabs( expectedMatrix( i,j ) - readMatrix( i,j ) ) >
                 std::numeric_limits< double >::epsilon( ) )
            {
                isTestMatrixBroken = true;
                cerr << "The unit test for the matrix file reader misread testMatrix.txt." << endl;
            }
        }
    }

    if ( readMatrix.rows( ) != 4 || readMatrix.cols ( ) != 3 )
    {
        isTestMatrixBroken = true;
        cerr << "The unit test for the matrix file reader gives an incorrect size matrix.\n";
    }

    // Test 2: Test for space-separated files.
    readMatrix = tudat::input_output::readMatrixFromFile(
                tudat::input_output::getTudatRootPath( ) +
                "InputOutput/UnitTests/testMatrix2.txt", " \t", "#" );

    for ( int i = 0; i < 3; i++ )
    {
        for ( int j = 0; j < 2; j++ )
        {
            if ( std::fabs(expectedMatrix( i,j )- readMatrix( i, j) ) >
                 std::numeric_limits< double >::epsilon( ) )
            {
                isTestMatrixBroken = true;
                cerr << "The unit test for the matrix file reader misread testMatrix2.txt." << endl;
            }
        }
    }

    return isTestMatrixBroken;
}
