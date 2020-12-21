/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "prototype/auxilliary.h"
#include "prototype/gravityAssist.h"
#include "prototype/planetaryDeparture.h"
#include "prototype/planetaryRendezvous.h"
#include <pybind11/pybind11.h>

using namespace tudat::prototype;

namespace py = pybind11;

namespace tudatpy {

void expose_iod(py::module &m) {
  py::class_<PlanetaryDeparture, std::shared_ptr<PlanetaryDeparture>>(m, "PlanetaryDeparture")
      .def(py::init<Eigen::Vector3d, std::shared_ptr<BaseState>, double, double, bool>(),
           py::arg("outgoing_velocity"),
           py::arg("central_body_state"),
           py::arg("periapsis_distance"),
           py::arg("gravitational_parameter"),
           py::arg("prograde_orbit") = true)
      .def_property_readonly("hyperbolic_excess_velocity", &PlanetaryDeparture::getHyperbolicExcessVelocity)
      .def_property_readonly("hyperbolic_excess_speed", &PlanetaryDeparture::getHyperbolicExcessSpeed)
      .def_property_readonly("eccentricity", &PlanetaryDeparture::getEccentricity)
      .def_property_readonly("angular_momentum", &PlanetaryDeparture::getAngularMomentum)
      .def_property_readonly("periapsis_speed", &PlanetaryDeparture::getPeriapsisSpeed)
      .def_property_readonly("parking_orbit_speed", &PlanetaryDeparture::getParkingOrbitSpeed)
      .def_property_readonly("required_delta_v", &PlanetaryDeparture::getRequiredDeltaV)
      .def_property_readonly("beta_angle", &PlanetaryDeparture::getBetaAngle)
      .def_property_readonly("orbital_plane_vector", &PlanetaryDeparture::getOrbitalPlaneVector)
      .def_property_readonly("periapsis_velocity_unit_vector", &PlanetaryDeparture::getPeriapsisVelocityUnitVector)
      .def_property_readonly("periapsis_outgoing_velocity", &PlanetaryDeparture::getPeriapsisOutgoingVelocity)
      .def_property_readonly("periapsis_incoming_velocity", &PlanetaryDeparture::getPeriapsisIncomingVelocity)
      .def_property_readonly("periapsis_position", &PlanetaryDeparture::getPeriapsisPosition)
      .def_property_readonly("sphere_of_influence_distance_at_exit", &PlanetaryDeparture::getSphereOfInfluenceDistanceAtExit)
      .def_property_readonly("true_anomaly_exit", &PlanetaryDeparture::getTrueAnomalyExitAsymptote)
      .def_property_readonly("semi_major_axis_outgoing", &PlanetaryDeparture::getSemiMajorAxisOutgoing)
      .def_property_readonly("hyperbolic_anomaly_at_exit", &PlanetaryDeparture::getHyperbolicAnomalyAtExit)
      .def_property_readonly("time_of_flight_to_exit", &PlanetaryDeparture::getTimeOfFlightToExit)
      .def_property_readonly("outgoing_cartesian_state", &PlanetaryDeparture::getOutgoingCartesianState)
      .def_property_readonly("incoming_cartesian_state", &PlanetaryDeparture::getIncomingCartesianState)
      .def_property_readonly("central_body_state", &PlanetaryDeparture::getCentralBodyState)
      .def_property_readonly("outgoing_velocity", &PlanetaryDeparture::getOutgoingVelocity);

  //////////////////////////////////////////////////////////////////////
  // prototype / planetaryArrival
  //////////////////////////////////////////////////////////////////////
  py::class_<PlanetaryRendezvous, std::shared_ptr<PlanetaryRendezvous>>(m, "PlanetaryRendezvous")
      .def(py::init<Eigen::Vector3d, std::shared_ptr<BaseState>, double, double, bool>(),
           py::arg("incoming_velocity"),
           py::arg("central_body_state"),
           py::arg("periapsis_distance"),
           py::arg("gravitational_parameter"),
           py::arg("prograde_orbit") = true)
      .def_property_readonly("hyperbolic_excess_velocity", &PlanetaryRendezvous::getHyperbolicExcessVelocity)
      .def_property_readonly("hyperbolic_excess_speed", &PlanetaryRendezvous::getHyperbolicExcessSpeed)
      .def_property_readonly("eccentricity", &PlanetaryRendezvous::getEccentricity)
      .def_property_readonly("angular_momentum", &PlanetaryRendezvous::getAngularMomentum)
      .def_property_readonly("periapsis_speed", &PlanetaryRendezvous::getPeriapsisSpeed)
      .def_property_readonly("parking_orbit_speed", &PlanetaryRendezvous::getParkingOrbitSpeed)
      .def_property_readonly("required_delta_v", &PlanetaryRendezvous::getRequiredDeltaV)
      .def_property_readonly("beta_angle", &PlanetaryRendezvous::getBetaAngle)
      .def_property_readonly("turning_angle", &PlanetaryRendezvous::getTurningAngle)
      .def_property_readonly("orbital_plane_vector", &PlanetaryRendezvous::getOrbitalPlaneVector)
      .def_property_readonly("periapsis_velocity_unit_vector", &PlanetaryRendezvous::getPeriapsisVelocityUnitVector)
      .def_property_readonly("periapsis_outgoing_velocity", &PlanetaryRendezvous::getPeriapsisOutgoingVelocity)
      .def_property_readonly("periapsis_incoming_velocity", &PlanetaryRendezvous::getPeriapsisIncomingVelocity)
      .def_property_readonly("periapsis_position", &PlanetaryRendezvous::getPeriapsisPosition)
      .def_property_readonly("sphere_of_influence_distance_at_entrance", &PlanetaryRendezvous::getSphereOfInfluenceDistanceAtEntrance)
      .def_property_readonly("true_anomaly_entrance", &PlanetaryRendezvous::getTrueAnomalyEntranceAsymptote)
      .def_property_readonly("semi_major_axis_incoming", &PlanetaryRendezvous::getSemiMajorAxisIncoming)
      .def_property_readonly("hyperbolic_anomaly_at_entrance", &PlanetaryRendezvous::getHyperbolicAnomalyAtEntrance)
      .def_property_readonly("time_of_flight_from_entrance", &PlanetaryRendezvous::getTimeOfFlightFromEntrance)
      .def_property_readonly("outgoing_cartesian_state", &PlanetaryRendezvous::getOutgoingCartesianState)
      .def_property_readonly("incoming_cartesian_state", &PlanetaryRendezvous::getIncomingCartesianState)
      .def_property_readonly("central_body_state", &PlanetaryRendezvous::getCentralBodyState)
      .def_property_readonly("incoming_velocity", &PlanetaryRendezvous::getIncomingVelocity);
};

}// namespace tudatpy
