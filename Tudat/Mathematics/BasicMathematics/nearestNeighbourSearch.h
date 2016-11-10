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
 *      120716    D. Dirkx          Updated with new nearest neighbour search algorithms.
 *      130114    D. Dirkx          Added missing include statements; corrected include guard
 *                                  name.
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2002.
 *
 *    Notes
 *
 */

#ifndef TUDAT_NEAREST_NEIGHBOR_SEARCH_H
#define TUDAT_NEAREST_NEIGHBOR_SEARCH_H 

#include <map>
#include <stdexcept>
#include <vector>

#include <boost/exception/all.hpp>

#include <Eigen/Core>

namespace tudat
{
namespace basic_mathematics
{

//! Nearest left neighbor binary search.
/*!
 * Searches for the nearest left neighbor in a vector of sorted data using a
 * binary algorithm (Press W.H., et al., 2002).
 * \param vectorOfSortedData Vector of data sorted in ascending/descending order.
 * \param targetValueInVectorOfSortedData Target value in vector of sorted data.
 * \return Index of nearest left neighbor to target value.
 */
int computeNearestLeftNeighborUsingBinarySearch( const Eigen::VectorXd& vectorOfSortedData,
                                                 const double targetValueInVectorOfSortedData );

//! Nearest neighbor binary search.
/*!
 * Searches for the nearest neighbor in a vector of sorted data using a
 * binary algorithm (Press W.H., et al., 2002).
 * \param vectorOfSortedData Vector of data sorted in ascending/descending order.
 * \param targetValueInVectorOfSortedData Target value in vector of sorted data.
 * \return Index of nearest neighbor to target value.
 */
int computeNearestNeighborUsingBinarySearch( const Eigen::VectorXd& vectorOfSortedData,
                                             const double targetValueInVectorOfSortedData );

//! Nearest left neighbor binary search.
/*!
 * Searches for the nearest left neighbor in a map of sorted data using a
 * binary algorithm (Press W.H., et al., 2002).
 * \param sortedIndepedentAndDependentVariables Map of independent and
 *           dependent data sorted in ascending/descending order.
 * \param targetValueInMapOfData Target value in map of sorted data.
 * \return Index of nearest left neighbor to target value.
 */
int computeNearestLeftNeighborUsingBinarySearch(
        const std::map < double, Eigen::VectorXd >& sortedIndepedentAndDependentVariables,
        const double targetValueInMapOfData );

//! Templated nearest left neighbor binary search.
/*!
 *  Templated nearest left neighbor binary search.
 *  \tparam IndependentVariableType Type of independent variables in which search is to be done.
 *  \param vectorOfSortedData STL vector, sorted in ascending order,
 *  containing independent variable values.
 *  \param targetValueInVectorOfSortedData Value of independent variable of which the nearest left
 *  neighbour is to be determined.
 */
template< typename IndependentVariableType >
int computeNearestLeftNeighborUsingBinarySearch(
        const std::vector< IndependentVariableType > vectorOfSortedData,
        const IndependentVariableType targetValueInVectorOfSortedData )
{
    // Declare local variables.
    // Declare bounds of vector of sorted data and current position.
    int leftLimitOfVectorOfSortedData = 0;
    int rightLimitOfVectorOfSortedData = vectorOfSortedData.size( ) - 1;
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

//! This function checks whether a value is in a given interval of a sorted vector.
/*!
 * This function checks whether a value is in a given interval of an STL vector, sorted in
 * ascending order of value of the entries. The interval is identified by the lower index of
 * the entries of the vector defining the interval.
 * \tparam IndependentVariableType Type of the independent variable values.
 * \param lowerIndex Index of lower bound of interval under consideration.
 * \param independentVariableValue Value of which it is to be checked whether it is in the
 *          interval.
 * \param independentValues Vector of values, srted in ascending order.
 * \return True if value is in interval, false otherwise.
 */
template< typename IndependentVariableType >
bool isIndependentVariableInInterval( const int lowerIndex,
                                  const IndependentVariableType independentVariableValue,
                                  const std::vector< IndependentVariableType >& independentValues )
{
    bool isInInterval = false;

    //Check if value is in interval.
    if ( independentVariableValue >= independentValues[ lowerIndex ] &&
         independentVariableValue <= independentValues[ lowerIndex + 1 ] )
    {
        isInInterval = true;
    }

    return isInInterval;
}

//! Nearest left leighbour search using hunting algorithm.
/*!
 * Nearest left leighbour search using hunting algorithm, using an initial guess to decrease
 * the look-up time. Especially useful for long data vectors where the value for which the
 * nearest neighbour is to be found changes slowly w.r.t. the data entries. Implementation is
 * taken from (Press W.H., et al., 2002). Algorithm is at worst twice as slow as binary search
 * and can be orders of magnitude faster in suitable cases.
 * \tparam IndependentVariableType Type for entries of vector of in which nearest neighbour is
 *          sought.
 * \param independentVariableValue Value of which the nearest left neighbour is to be calculated.
 * \param previousNearestLowerIndex_ Initial guess of nearest lft neighbour.
 * \param independentValues_ Vector of independent variables, sorted in ascending order, in which
 *          the nearest left (lower) neighbour is to be determined.
 * \return Index of independentValues_ that is the nearest left neighbour
 */
template< typename IndependentVariableType >
int findNearestLeftNeighbourUsingHuntingAlgorithm(
        const IndependentVariableType independentVariableValue,
        const int previousNearestLowerIndex_,
        const std::vector< IndependentVariableType >& independentValues_ )
{
    // Initialize return variable for new nearest left neighbor.
    int newNearestLowerIndex = 0;

    // Initialize boolean denoting whether the new value has been found.
    bool isFound = 0;

    int independentValueVectorSize = static_cast< int >( independentValues_.size( ) );

    // Check whether initial estimate is possible.
    if ( previousNearestLowerIndex_ < 0 ||
         previousNearestLowerIndex_ > static_cast< int >( independentValues_.size( ) - 2 ) )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                                    std::runtime_error(
        "Error, initial guess for nearest neighbour search not within allowable bounds." ) ) );
    }

    // Check if independent variable value falls within region of values provided.
    if ( independentVariableValue <= independentValues_[ 0 ] )
    {
        newNearestLowerIndex = 0;
        isFound = 1;
    }

    else if ( independentVariableValue >= independentValues_[ independentValueVectorSize - 1 ] )
    {
        newNearestLowerIndex = independentValueVectorSize - 2;
        isFound = 1;
    }

    // Perform hunting algorithm if requested value is not outside given range.
    if ( !isFound )
    {
        // Create indices delimiting interval under consideration.
        int upperIndex, lowerIndex;

        // Create size and direction of jumps.
        int jumpSize = 1;
        int jumpDirection;

        // Check direction of jumps that is needed.
        if ( independentVariableValue >= independentValues_[ previousNearestLowerIndex_ ] )
        {
            jumpDirection = 1;
            lowerIndex = previousNearestLowerIndex_;
            upperIndex = previousNearestLowerIndex_ + 1;
        }
        else
        {
            jumpDirection = -1;
            lowerIndex = previousNearestLowerIndex_ - 1;
            upperIndex = previousNearestLowerIndex_;
        }

        // Boolean to denote whether the value is in te interval under consideration.
        bool isRegionReached = 0;

        // While nearest lower neighbour is not found, continue algorithm.
        while ( !isFound )
        {
            // Check if value is in regio under consideration.
            if( independentVariableValue >= independentValues_[ lowerIndex ] &&
                    independentVariableValue <= independentValues_[ upperIndex ] )
            {
                isRegionReached = 1;
            }

            // If not in current region:
            else
            {
                // Double jump size.
                jumpSize *= 2;
                if ( jumpDirection == 1 )
                {
                    // Reset interval under consideration.
                    lowerIndex = upperIndex;
                    upperIndex += jumpSize;

                    // Correct jump size if it would bring interval out of bounds.
                    if ( upperIndex > independentValueVectorSize - 1 )
                    {
                        upperIndex = independentValueVectorSize - 1;
                        jumpSize = upperIndex - lowerIndex;
                    }
                }
                else
                {
                    // Reset interval under consideration.
                    upperIndex = lowerIndex;
                    lowerIndex -= jumpSize;

                    // Correct jump size if it would bring interval out of bounds.
                    if ( lowerIndex < 0 )
                    {
                        lowerIndex = 0;
                        jumpSize = upperIndex;
                    }
                }
            }

            // Independent variable value has been identified to be between upperIndex and
            // lowerIndex. Perform bisection of interval to decrease interval to size 1.
            if ( isRegionReached )
            {
                while( !isFound )
                {

                    int middleIndex;
                    assert( upperIndex - lowerIndex > 0 );

                    // If the upper and lower indices have a difference of exactly one, the
                    // interval has been found.
                    if ( upperIndex - lowerIndex == 1 )
                    {
                        isFound = 1;
                        newNearestLowerIndex = lowerIndex;
                    }
                    else
                    {
                        // Half size of interval.
                        jumpSize /= 2;
                        middleIndex = lowerIndex + jumpSize;

                        // Check which half of interval is to be new interval.
                        if ( independentVariableValue < independentValues_[ middleIndex ] )
                        {
                            upperIndex = middleIndex;
                        }
                        else
                        {
                            lowerIndex = middleIndex;
                        }

                        jumpSize = upperIndex - lowerIndex;
                    }
                }
            }
        }
    }

    return newNearestLowerIndex;
}

} // namespace basic_mathematics
} // namespace tudat

#endif // TUDAT_NEAREST_NEIGHBOR_SEARCH_H 
