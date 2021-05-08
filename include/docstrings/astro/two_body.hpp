/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "string"

#ifndef TUDATPY_TWO_BODY_H
#define TUDATPY_TWO_BODY_H

std::string lambert_targeter_docstring() {
  return R"(//! Lambert targeting algorithm class.
)";
}

std::string lambert_targeter_ctor_docstring() {
  return R"(
Parameters
----------
departure_position : np.ndarray[float64[3,1]]
    The position at departure in Cartesian coordinates. [m]
arrival_position : np.ndarray[float64[3,1]]
    The position at arrival in Cartesian coordinates. [m]
time_of_flight : float
    The time-of-flight between departure and arrival. [s]
gravitational_param : float
    The gravitational parameter of the main body. [m^3 s^-2]
)";
}

#endif//TUDATPY_TWO_BODY_H
