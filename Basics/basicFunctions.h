/*! \file basicFunctions.h
 *    Header file that defines the basicFunctions namespace, containing all
 *    basic functions contained in Tudat.
 *
 *    Path              : /Basics/
 *    Version           : 8
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author            : J. Leloux
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.Leloux@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 10 August, 2010
 *    Last modified     : 10 August, 2011
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      110803    J. Leloux           Added convertStringToTemplate.
 *      110805    J. Leloux           Added outputCurrentRunningTime.
 *      110810    J. Leloux           Minor comment modifications.
 */

#ifndef BASICOPERATIONS_H
#define BASICOPERATIONS_H

// Include statements.
#include <ctime>
#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <sstream>
#include <vector>
#include "linearAlgebra.h"
#include "state.h"

// Using declarations.
using std::string;
using std::vector;

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

//! Convert string to variable type.
/*!
 * Definition of a template function which converts a string to any variable
 * type which it can be converted to, and is used as input.
 * \param inputString String to be converted.
 * \param outputTemplate Template type to which the string will be converted.
 * \return Boolean stating the success or failure of the conversion.
 */
template < class T >
bool convertStringToTemplate( const std::string& inputString, T& outputTemplate )
{
    // Create stringstream containing input string.
    std::istringstream inputStringStream( inputString );

    // Return output template containing contents of stringstream.
    return inputStringStream >> outputTemplate;
}

//! Write the current running time and status to vector.
/*!
 * Definition of a function which determines the running time of a program with respect to a start
 * clock and outputs this time together with a string containing the current state of the program.
 * \param start_clock Input starting clock.
 * \param status Current status of executed application.
 * \return Container of current status of executed application and current running time statement.
 *          The first element contains the current status and the second contains the current
 *          running time statement.
 */
vector< string > outputCurrentRunningTime( clock_t start_clock, const string& status );

}

#endif // BASICOPERATIONS_H

// End of file.
