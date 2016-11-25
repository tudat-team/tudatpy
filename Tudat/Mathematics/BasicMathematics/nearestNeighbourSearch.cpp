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
 *      100916    D. Dirkx          Added minor comments during checking.
 *      100928    K. Kumar          Small comment modifications.
 *      100929    K. Kumar          Small comment modifications.
 *      110202    K. Kumar          Added overload for map with State* for
 *                                  computeNearestLeftNeighborUsingBinarySearch( ).
 *      110803    J. Leloux         Added convertStringToTemplate.
 *      110805    J. Leloux         Added outputCurrentRunningTime( ).
 *      110807    K. Kumar          Minor comment modifications.
 *      110810    J. Leloux         Minor comment modifications.
 *      110913    K. Kumar          Implemented automatic root-path functions based on
 *                                  suggestions by M. Persson.
 *      111117    K. Kumar          Added listAllFilesInDirectory( ) function.
 *      160930    M. Van den Broeck Added computeNearestNeighborUsingBinarySearch function
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2002.
 *
 *    Notes
 *      Need to add bounds checking with respect to targetValue for
 *      computeNearestLeftNeighborUsingBinarySearch( ).
 */

#include <map>
#include <iterator>
#include <Eigen/Core>

namespace tudat
{
namespace basic_mathematics
{

//! Nearest left neighbor binary search.
int computeNearestLeftNeighborUsingBinarySearch(
        const Eigen::VectorXd& vectorOfSortedData,
        const double targetValueInVectorOfSortedData )
{
    // Declare local variables.
    // Declare bounds of vector of sorted data and current position.
    int leftLimitOfVectorOfSortedData = 0;
    int rightLimitOfVectorOfSortedData = vectorOfSortedData.rows( ) - 1;
    int currentPositionInVectorOfSortedData;

    // Check if data is sorted in ascending order.
    // ( true if ascending, else false ).
    bool isVectorOfSortedDataAscending
            = ( vectorOfSortedData[ rightLimitOfVectorOfSortedData ]
                >= vectorOfSortedData[ leftLimitOfVectorOfSortedData ] );

    // Loop through vector of sorted data until left and right limits
    // are neighbours.
    while ( rightLimitOfVectorOfSortedData
            - leftLimitOfVectorOfSortedData > 1 )
    {
        // Compute midpoint ( bitshift is same as division by 2.0 ).
        currentPositionInVectorOfSortedData
                = ( rightLimitOfVectorOfSortedData
                    + leftLimitOfVectorOfSortedData ) >> 1;

        // Check which limit to replace ( if ascending and target datum
        // is in right half, replace left limit ).
        if ( targetValueInVectorOfSortedData
             >= vectorOfSortedData[ currentPositionInVectorOfSortedData ]
             && isVectorOfSortedDataAscending )
        {
            // Set left limit to current position in vector of sorted data.
            leftLimitOfVectorOfSortedData = currentPositionInVectorOfSortedData;
        }

        else
        {
            // Set right limit to current position in vector of sorted data.
            rightLimitOfVectorOfSortedData = currentPositionInVectorOfSortedData;
        }
    }

    // Set current position to left limit.
    currentPositionInVectorOfSortedData = leftLimitOfVectorOfSortedData;

    // Return current position in vector.
    return currentPositionInVectorOfSortedData;
}

//! Nearest neighbor binary search.
int computeNearestNeighborUsingBinarySearch(
        const Eigen::VectorXd& vectorOfSortedData,
        const double targetValueInVectorOfSortedData )
{
    // Declare local variables.
    // Declare bounds of vector of sorted data and current position.
    int leftLimitOfVectorOfSortedData = 0;
    int rightLimitOfVectorOfSortedData = vectorOfSortedData.rows( ) - 1;
    int currentPositionInVectorOfSortedData;

    // Check if data is sorted in ascending order.
    // ( true if ascending, else false ).
    bool isVectorOfSortedDataAscending
            = ( vectorOfSortedData[ rightLimitOfVectorOfSortedData ]
                >= vectorOfSortedData[ leftLimitOfVectorOfSortedData ] );

    // Loop through vector of sorted data until left and right limits
    // are neighbours.
    while ( rightLimitOfVectorOfSortedData
            - leftLimitOfVectorOfSortedData > 1 ) // No rounding off errors because limits are integers
    {
        // Compute midpoint ( bitshift is same as division by 2.0 ).
        currentPositionInVectorOfSortedData
                = ( rightLimitOfVectorOfSortedData
                    + leftLimitOfVectorOfSortedData ) >> 1;

        // Check which limit to replace ( if ascending and target datum
        // is in right half, replace left limit ).
        if ( targetValueInVectorOfSortedData
             >= vectorOfSortedData[ currentPositionInVectorOfSortedData ]
             && isVectorOfSortedDataAscending )
        {
            // Set left limit to current position in vector of sorted data.
            leftLimitOfVectorOfSortedData = currentPositionInVectorOfSortedData;
        }

        else
        {
            // Set right limit to current position in vector of sorted data.
            rightLimitOfVectorOfSortedData = currentPositionInVectorOfSortedData;
        }
    }

    // Set current position to left or right limit, whichever is closer.
    currentPositionInVectorOfSortedData = leftLimitOfVectorOfSortedData;
    if( ( ( vectorOfSortedData[ rightLimitOfVectorOfSortedData ] +
              vectorOfSortedData[ leftLimitOfVectorOfSortedData ] ) / 2.0 ) < targetValueInVectorOfSortedData )
    {
        currentPositionInVectorOfSortedData = rightLimitOfVectorOfSortedData;
    }

    else
    {
        currentPositionInVectorOfSortedData = leftLimitOfVectorOfSortedData;
    }

    // Return current position in vector.
    return currentPositionInVectorOfSortedData;
}

//! Nearest left neighbor binary search.
int computeNearestLeftNeighborUsingBinarySearch(
        const std::map < double, Eigen::VectorXd >& sortedIndepedentAndDependentVariables,
        const double targetValueInMapOfData )
{
    // Declare local variables.
    // Declare bounds of key of map of data and current position.
    int leftLimitOfKeyOfMapOfData = 0;
    int rightLimitOfKeyOfMapOfData = sortedIndepedentAndDependentVariables
                                     .size( ) - 1;
    int currentPositionInKeyOfMapOfData;

    // Declare map iterator
    std::map < double, Eigen::VectorXd >::const_iterator mapIterator;

    // Loop through vector of sorted data until left and right limits
    // are neighbours
    while ( rightLimitOfKeyOfMapOfData - leftLimitOfKeyOfMapOfData > 1 )
    {
        // Compute midpoint ( bitshift is same as division by 2.0 ).
        currentPositionInKeyOfMapOfData
                = ( rightLimitOfKeyOfMapOfData
                   + leftLimitOfKeyOfMapOfData ) >> 1;

        // Set map iterator to begin begin of map of sorted independent and
        // dependent variables.
        mapIterator = sortedIndepedentAndDependentVariables.begin( );

        // Advance iterator to location of current position in key of map of
        // data.
        advance( mapIterator, currentPositionInKeyOfMapOfData );

        // Check that target value lies to the right of lower bound.
        if ( targetValueInMapOfData
             >= mapIterator->first )
        {
            // Set left limit to current position in map of data.
            leftLimitOfKeyOfMapOfData = currentPositionInKeyOfMapOfData;
        }

        else
        {
            // Set right limit to current position in map of data.
            rightLimitOfKeyOfMapOfData = currentPositionInKeyOfMapOfData;
        }
    }

    // Set current position to left limit.
    currentPositionInKeyOfMapOfData = leftLimitOfKeyOfMapOfData;

    // Return current position in map of data.
    return currentPositionInKeyOfMapOfData;
}

} // namespace basic_mathematics
} // namespace tudat
