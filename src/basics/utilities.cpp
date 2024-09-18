/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/basics/utilities.h"

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

/*!
 * Function to extract a map from string to 3d vector from a file. The first 4 columns are used and the rest is ignored if present
 * @param fileName path to file of interest
 * @param commentSymbol Lines starting with this character are ignored
 * @param separators String of characters that mark a new column
 * @return map with string to 3d vector
 */
template< >
std::map<std::string, Eigen::Vector3d> getMapFromFile<std::string, Eigen::Vector3d>(std::string fileName,
                                                                                    char commentSymbol,
                                                                                    std::string separators)
{
  std::ifstream file(fileName);
  if (!file.good()) {
    throw std::runtime_error("Error when opening file: " + fileName + " could not be opened.");
  }
  std::string currentLine;
  std::vector<std::string> currentSplitLine;
  std::map<std::string, Eigen::Vector3d> namesAndPositions;
  while (std::getline(file, currentLine)) {
    if (!currentLine.empty() && currentLine.at(0) != commentSymbol) {
      boost::algorithm::trim(currentLine);
      boost::algorithm::split(currentSplitLine,
                              currentLine,
                              boost::is_any_of(separators),
                              boost::algorithm::token_compress_on);
      try {
        namesAndPositions[currentSplitLine.at(0)] = Eigen::Vector3d(std::stod(currentSplitLine.at(1)),
                                                                    std::stod(currentSplitLine.at(2)),
                                                                    std::stod(currentSplitLine.at(3)));
      } catch (...) {
        continue;
      }// Ignore lines that cannot be read

      if (namesAndPositions.empty()) {
        throw std::runtime_error("Error when reading file: " + fileName + "has no lines in acceptable format.");
      }
    }
  }

  return namesAndPositions;
}

/*!
 * Function to extract a map from string to 3d vector from a file. The first 4 columns are used and the rest is ignored if present
 * @param fileName path to file of interest
 * @param commentSymbol Lines starting with this character are ignored
 * @param separators String of characters that mark a new column
 * @return map with string to 3d vector
 */
template< >
std::map<std::string, std::string> getMapFromFile<std::string, std::string>(std::string fileName, char commentSymbol, std::string separators)
{
  std::ifstream file(fileName);
  if (!file.good()) {
    throw std::runtime_error("Error when opening file: " + fileName + " could not be opened.");
  }
  std::string currentLine;
  std::vector<std::string> currentSplitLine;
  std::map<std::string, std::string> namesAndPositions;
  while (std::getline(file, currentLine)) {
    if (!currentLine.empty() && currentLine.at(0) != commentSymbol) {
      boost::algorithm::trim(currentLine);
      boost::algorithm::split(currentSplitLine,
                              currentLine,
                              boost::is_any_of(separators),
                              boost::algorithm::token_compress_on);
      try {
        namesAndPositions[currentSplitLine.at(0)] = currentSplitLine.at(1);
      } catch (...) {
        continue;
      }// Ignore lines that cannot be read

      if (namesAndPositions.empty()) {
        throw std::runtime_error("Error when reading file: " + fileName + "has no lines in acceptable format.");
      }
    }
  }

  return namesAndPositions;
}

}

}
