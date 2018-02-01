/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Basics/utilities.h"

namespace tudat
{

namespace utilities
{

//! Get indices of pointer to single entry in multi-array (size 1) of doubles
boost::array< boost::multi_array< double, 1 >::index, 1 > getMultiArrayIndexArray(
        const boost::multi_array< double, 1 >& multiArray, const double* requestedElement )
{
    typedef boost::multi_array< double, 1 > NMultiArray;
    boost::array< NMultiArray::index, 1 >  currentIndices;

    for ( unsigned int dir = 0; dir < 1; dir++ )
    {
        currentIndices[ dir ] = getMultiArrayIndex< 1 >( multiArray, requestedElement, dir );
    }

    return currentIndices;
}

//! Get indices of pointer to single entry in multi-array (size 2) of doubles
boost::array< boost::multi_array< double, 2 >::index, 2 > getMultiArrayIndexArray(
        const boost::multi_array< double, 2 >& multiArray, const double* requestedElement )
{
    typedef boost::multi_array< double, 2 > NMultiArray;
    boost::array< NMultiArray::index, 2 >  currentIndices;

    for ( unsigned int dir = 0; dir < 2; dir++ )
    {
        currentIndices[ dir ] = getMultiArrayIndex< 2 >( multiArray, requestedElement, dir );
    }

    return currentIndices;
}

//! Get indices of pointer to single entry in multi-array (size 3) of doubles
boost::array< boost::multi_array< double, 3 >::index, 3 > getMultiArrayIndexArray(
        const boost::multi_array< double, 3 >& multiArray, const double* requestedElement )
{
    typedef boost::multi_array< double, 3 > NMultiArray;
    boost::array< NMultiArray::index, 3 >  currentIndices;

    for ( unsigned int dir = 0; dir < 3; dir++ )
    {
        currentIndices[ dir ] = getMultiArrayIndex< 3 >( multiArray, requestedElement, dir );
    }

    return currentIndices;
}

}

}
