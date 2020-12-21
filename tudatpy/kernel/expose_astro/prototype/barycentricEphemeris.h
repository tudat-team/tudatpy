/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_BARYCENTRICEPHEMERIS_H
#define TUDATPY_BARYCENTRICEPHEMERIS_H

#include <tudat/astro/ephemerides.h>

namespace tudat {

namespace ephemerides {
class BarycentricEphemeris : public Ephemeris {
 public:
  BarycentricEphemeris(std::vector<std::pair<std::function<double()>, std::shared_ptr<Ephemeris>>> gravitationalParameterEphemerisPairs);
  Eigen::Vector6d getCartesianState(const double secondsSinceEpoch);

 private:
  std::vector<std::pair<std::function<double()>, std::shared_ptr<Ephemeris>>> gravitationalParameterEphemerisPairs_;
};
}// namespace ephemerides
}// namespace tudat

#endif//TUDATPY_BARYCENTRICEPHEMERIS_H
