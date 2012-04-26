/*    Copyright (c) 2010 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110530    F.M. Engelen      First creation of code.
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
