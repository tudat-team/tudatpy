/*    Copyright (c) 2010-2018, Delft University of Technology
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
        const std::map< double, Eigen::VectorXd >& sortedIndepedentAndDependentVariables,
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
