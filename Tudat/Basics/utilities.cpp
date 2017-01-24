/*    Copyright (c) 2010-2016, Delft University of Technology
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

boost::array< boost::multi_array< double, 1 >::index, 1 > getMultiArrayIndexArray(
        const boost::multi_array< double, 1 >& m, const double* requestedElement )
{
    typedef boost::multi_array< double, 1 > NMultiArray;
    boost::array< NMultiArray::index, 1 >  currentIndices;
    for ( unsigned int dir = 0; dir < 1; dir++ )
    {
        currentIndices[ dir ] = getMultiArrayIndex< 1 >( m, requestedElement, dir );
    }

    return currentIndices;
}

boost::array< boost::multi_array< double, 2 >::index, 2 > getMultiArrayIndexArray(
        const boost::multi_array< double, 2 >& m, const double* requestedElement )
{
    typedef boost::multi_array< double, 2 > NMultiArray;
    boost::array< NMultiArray::index, 2 >  currentIndices;
    for ( unsigned int dir = 0; dir < 2; dir++ )
    {
        currentIndices[ dir ] = getMultiArrayIndex< 2 >( m, requestedElement, dir );
    }

    return currentIndices;
}

boost::array< boost::multi_array< double, 3 >::index, 3 > getMultiArrayIndexArray(
        const boost::multi_array< double, 3 >& m, const double* requestedElement )
{
    typedef boost::multi_array< double, 3 > NMultiArray;
    boost::array< NMultiArray::index, 3 >  currentIndices;
    for ( unsigned int dir = 0; dir < 3; dir++ )
    {
        currentIndices[ dir ] = getMultiArrayIndex< 3 >( m, requestedElement, dir );
    }

    return currentIndices;
}

}

}
