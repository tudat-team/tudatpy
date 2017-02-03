/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      If tabs are used as spaces, it doesn't work. The seperator should also be tabs then.
 *
 */

#ifndef TUDAT_MATRIX_TEXT_FILEREADER_H
#define TUDAT_MATRIX_TEXT_FILEREADER_H

#include <string>

#include <Eigen/Core>

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
