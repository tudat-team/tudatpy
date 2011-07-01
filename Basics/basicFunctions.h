/*! \file basicFunctions.h
 *    Header file that defines the basicFunctions namespace, containing all
 *    basic functions contained in Tudat.
 *
 *    Path              : /Basics/
 *    Version           : 5
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Date created      : 10 August, 2010
 *    Last modified     : 2 February, 2010
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
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
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author              Comment
 *      100902    K. Kumar            File header and footer added.
 *      100916    D. Dirkx            Added minor comments and placeholder
 *                                    tag during checking.
 *      100928    K. Kumar            Added reference and adjusted include
 *                                    statements.
 *      100929    K. Kumar            Changed EigenRoutines.h include statement
 *                                    to linearAlgebra.h and removed
 *                                    placeholder comment. Added small comment
 *                                    modifications.
 *      110202    K. Kumar            Added overload for map with State* for
 *                                    computeNearestLeftNeighborUsingBinarySearch().
 */

#ifndef BASICOPERATIONS_H
#define BASICOPERATIONS_H

// Include statements.
#include <iostream>
#include <map>
#include <iterator>
#include <string>
#include "linearAlgebra.h"
#include "state.h"

//! Basic functions namespace.
/*!
 * Basic functions namespace.
 */
namespace basic_functions
{

//! Root path of Tudat library.
/*!
 * Root path of Tudat library must be set for file-reading/file-writing purposes.
 * Path must include a trailing slash.
 */
static std::string ROOT_PATH;

//! Nearest left neighbor binary search.
/*!
 * Searches for the nearest left neighbor in a vector of sorted data using a
 * binary algorithm (Press W.H., et al., 2002).
 * \param vectorOfSortedData Vector of data sorted in ascending/descending order.
 * \param targetValueInVectorOfSortedData Target value in vector of sorted data.
 * \return Index of nearest left neighbor to target value.
 */
int computeNearestLeftNeighborUsingBinarySearch(
        VectorXd& vectorOfSortedData,
        double& targetValueInVectorOfSortedData );

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
        std::map < double, VectorXd >& sortedIndepedentAndDependentVariables,
        double& targetValueInMapOfData );

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
        std::map < double, State* >& sortedIndepedentAndDependentVariables,
        double& targetValueInMapOfData );

}

#endif // BASICOPERATIONS_H

// End of file.
