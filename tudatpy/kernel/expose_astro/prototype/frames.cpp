/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "frames.h"
#include <sstream>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <tudat/astro/reference_frames.h>

ReferenceFrame::ReferenceFrame(
    std::string origin,
    std::string orientation)
    : isSet_(true),
      origin_(origin),
      orientation_(orientation) {
}

ReferenceFrame::ReferenceFrame()
    : isSet_(false) {
}

bool ReferenceFrame::isSet() { return isSet_; };

std::string ReferenceFrame::getString() {
  std::ostringstream out;
  if (isSet()) {
    out << "ReferenceFrame(" << getOrigin() << ", " << getOrientation() << ")";
  } else {
    out << "ReferenceFrame(None)";
  }
  return out.str();
}

std::string ReferenceFrame::getOrigin() {
  if (isSet()) {
    return origin_;
  } else {
    throw std::logic_error("getOrigin() was invoked without the reference frame being defined.");
  }
}
std::string ReferenceFrame::getOrientation() {
  if (isSet()) {
    return orientation_;
    ;
  } else {
    throw std::logic_error("getOrientation() was invoked without the reference frame being defined.");
  }
}

const ReferenceFrame SSB_J2000() { return ReferenceFrame("SSB", "J2000"); };
const ReferenceFrame SSB_ECLIPJ2000() { return ReferenceFrame("SSB", "ECLIPJ2000"); };
const ReferenceFrame SUN_ECLIPJ2000() { return ReferenceFrame("Sun", "ECLIPJ2000"); };
const ReferenceFrame SUN_J2000() { return ReferenceFrame("Sun", "J2000"); };

#define SSB_J2000 SSB_J2000()
#define SSB_ECLIPJ2000 SSB_ECLIPJ2000()
#define SUN_ECLIPJ2000 SUN_ECLIPJ2000()
#define SUN_J2000 SUN_J2000()