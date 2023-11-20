/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2002.
 *
 */

#ifndef TUDAT_NEAREST_NEIGHBOR_SEARCH_H
#define TUDAT_NEAREST_NEIGHBOR_SEARCH_H 

#include <map>
#include <stdexcept>
#include <vector>

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
    if( vectorOfSortedData.size( ) <= 0 )
    {
        throw std::runtime_error( "Error when getting nearest neighbour, no input data." );
    }
    int rightLimitOfVectorOfSortedData = static_cast< int >( vectorOfSortedData.size( ) ) - 1;
    int currentPositionInVectorOfSortedData;

    // Check if data is sorted in ascending order.
    // ( true if ascending, else false ).
    bool isVectorOfSortedDataAscending
            = ( vectorOfSortedData[ rightLimitOfVectorOfSortedData ]
                >= vectorOfSortedData[ leftLimitOfVectorOfSortedData ] );

    // Check if target value is to the right of largest limit
    if ( std::max( { vectorOfSortedData[ rightLimitOfVectorOfSortedData ], vectorOfSortedData[ leftLimitOfVectorOfSortedData ] } )
            <= targetValueInVectorOfSortedData )
    {
        if ( isVectorOfSortedDataAscending )
        {
            leftLimitOfVectorOfSortedData = rightLimitOfVectorOfSortedData;
        }
        else
        {
            rightLimitOfVectorOfSortedData = leftLimitOfVectorOfSortedData;
        }
    }

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
    // Check if independent variable coincide with found upper index. In that case the upper index is the true nearest lower index.
    if ( targetValueInVectorOfSortedData == vectorOfSortedData[ rightLimitOfVectorOfSortedData ] )
    {
        currentPositionInVectorOfSortedData = rightLimitOfVectorOfSortedData;
    }

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
         independentVariableValue < independentValues[ lowerIndex + 1 ] )
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

    if( !( independentVariableValue == independentVariableValue ) )
    {
        throw std::runtime_error( "Error in nearest left neighbour search, input is NaN" );
    }

    if( independentValueVectorSize < 2 )
    {
        throw std::runtime_error( "Error in nearest neighbour search, size of input vector is " +
                                  std::to_string( independentValueVectorSize ) );
    }

    // Check whether initial estimate is possible.
    if ( previousNearestLowerIndex_ < 0 ||
         previousNearestLowerIndex_ > static_cast< int >( independentValues_.size( ) - 2 ) )
    {
        // throw std::runtime_error( "Error, initial guess for nearest neighbour search not within allowable bounds." );
        newNearestLowerIndex = computeNearestLeftNeighborUsingBinarySearch< IndependentVariableType >(
                independentValues_, independentVariableValue );
        isFound = 1;
    }

    // Check if independent variable value falls within region of values provided.
    if ( independentVariableValue <= independentValues_[ 0 ] )
    {
        newNearestLowerIndex = 0;
        isFound = 1;
    }

    else if ( independentVariableValue >= independentValues_[ independentValueVectorSize - 1 ] )
    {
        newNearestLowerIndex = independentValueVectorSize - 1;
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

                    if( !( upperIndex - lowerIndex  > 0 ) )
                    {
                        throw std::runtime_error( "Error, upper and lower indices are inconsistent in nearest neighbour search" +
                                                  std::to_string( upperIndex ) + " " +
                                                  std::to_string( lowerIndex ) );
                    }

                    // If the upper and lower indices have a difference of exactly one, the
                    // interval has been found.
                    if ( upperIndex - lowerIndex == 1 )
                    {
                        isFound = 1;
                        newNearestLowerIndex = lowerIndex;

                        // Check if independent variable coincides with upper index. In that case the upper index is
                        // true lower index
                        if ( independentVariableValue == independentValues_[ upperIndex ] )
                        {
                            newNearestLowerIndex = upperIndex;
                        }
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
