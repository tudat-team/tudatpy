/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_FRAMES_H
#define TUDATPY_FRAMES_H

#include <string>

class ReferenceFrame {
 public:
  ReferenceFrame();
  ReferenceFrame(std::string origin,
                 std::string orientation);

  bool isSet();
  std::string getString();
  std::string getOrigin();
  std::string getOrientation();

 private:
  std::string origin_;
  std::string orientation_;
  bool isSet_;
};

const ReferenceFrame SSB_J2000();
const ReferenceFrame SSB_ECLIPJ2000();
const ReferenceFrame SUN_ECLIPJ2000();
const ReferenceFrame MARS_ECLIPJ2000();
const ReferenceFrame MARS_J2000();

#endif//TUDATPY_FRAMES_H
