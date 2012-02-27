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
 */

// Temporary notes (move to class/function doxygen):
// If tabs are used as spaces, it doesn't work. The seperator should also be tabs then.
// 

#ifndef TUDAT_MATRIX_TEXT_FILEREADER_H
#define TUDAT_MATRIX_TEXT_FILEREADER_H

#include <Eigen/Core>
#include <string>

namespace tudat
{
namespace input_output
{

//! Read the file and return the data matrix.
/*!
 * Read a textfile whith separated (space, tab, comma etc...) numbers. The class returns these
 * numbers as a matrixXd. The first line with numbers is used to define the number of columns.
 * \param relativePath Relative path to file.
 * \param separators Separators used, every character in the string will be used as separators.
 *         (multiple seperators possible).
 * \param skipLinesCharacter Skip lines starting with this character.
 * \return The datamatrix.
 */
Eigen::MatrixXd readMatrixFromFile( const std::string& relativePath,
                                    const std::string& separators = "\t ;,",
                                    const std::string& skipLinesCharacter = "%" );

} // namespace input_output
} // namespace tudat

#endif // TUDAT_MATRIX_TEXT_FILEREADER_H
