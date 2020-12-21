/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "barycentricEphemeris.h"

namespace tudat {

namespace ephemerides {

BarycentricEphemeris::BarycentricEphemeris(
    std::vector<std::pair<std::function<double()>, std::shared_ptr<Ephemeris>>> gravitationalParameterEphemerisPairs)
    : gravitationalParameterEphemerisPairs_(gravitationalParameterEphemerisPairs){};

Eigen::Vector6d BarycentricEphemeris::getCartesianState(const double secondsSinceEpoch) {
  Eigen::Vector6d weightedStates;
  weightedStates << 0, 0, 0, 0, 0, 0;
  double sumGravitationalParameter = 0;
  for (auto pair : gravitationalParameterEphemerisPairs_) {
    double gravitationalParameter = pair.first();
    weightedStates += gravitationalParameter * pair.second->getCartesianState(secondsSinceEpoch);
    sumGravitationalParameter += gravitationalParameter;
  }
};

}// namespace ephemerides
}// namespace tudat
