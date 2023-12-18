/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/readTrackingTxtFile.h"

namespace tudat
{
namespace input_output
{

void TrackingTxtFileContents::parseData()
{
  std::ifstream dataFile(fileName_);
  if (!dataFile.good()) {
    throw std::runtime_error("Error when opening Jpl Range file, file " + fileName_ + " could not be opened.");
  }
  readRawDataMap(dataFile);
  convertDataMap();
}

} // namespace input_output

} // namespace tudat
