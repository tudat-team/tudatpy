/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_fundamentals.h"

#include <tudat/astro/basic_astro.h>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;

namespace tudatpy {

void expose_fundamentals(py::module &m) {

  py::class_<
      tba::AccelerationModel<Eigen::Vector3d>,
      std::shared_ptr<tba::AccelerationModel<Eigen::Vector3d>>>
      acceleration_model(m, "AccelerationModel");

  // TODO: This should be moved in the tudat source to propagation_setup.
//  py::enum_<tba::AvailableAcceleration>(m, "AvailableAcceleration")
//      .value("undefined_acceleration", tba::AvailableAcceleration::undefined_acceleration)
//      .value("point_mass_gravity", tba::AvailableAcceleration::point_mass_gravity)
//      .value("central_gravity", tba::AvailableAcceleration::central_gravity)
//      .value("aerodynamic", tba::AvailableAcceleration::aerodynamic)
//      .value("cannon_ball_radiation_pressure", tba::AvailableAcceleration::cannon_ball_radiation_pressure)
//      .value("spherical_harmonic_gravity", tba::AvailableAcceleration::spherical_harmonic_gravity)
//      .value("mutual_spherical_harmonic_gravity", tba::AvailableAcceleration::mutual_spherical_harmonic_gravity)
//      .value("third_body_point_mass_gravity", tba::AvailableAcceleration::third_body_point_mass_gravity)
//      .value("third_body_central_gravity", tba::AvailableAcceleration::third_body_central_gravity)
//      .value("third_body_spherical_harmonic_gravity", tba::AvailableAcceleration::third_body_spherical_harmonic_gravity)
//      .value("third_body_mutual_spherical_harmonic_gravity", tba::AvailableAcceleration::third_body_mutual_spherical_harmonic_gravity)
//      .value("thrust_acceleration", tba::AvailableAcceleration::thrust_acceleration)
//      .value("relativistic_correction_acceleration", tba::AvailableAcceleration::relativistic_correction_acceleration)
//      .value("empirical_acceleration", tba::AvailableAcceleration::empirical_acceleration)
//      .value("direct_tidal_dissipation_in_central_body_acceleration", tba::AvailableAcceleration::direct_tidal_dissipation_in_central_body_acceleration)
//      .value("direct_tidal_dissipation_in_orbiting_body_acceleration", tba::AvailableAcceleration::direct_tidal_dissipation_in_orbiting_body_acceleration)
//      .value("panelled_radiation_pressure_acceleration", tba::AvailableAcceleration::panelled_radiation_pressure_acceleration)
//      .value("momentum_wheel_desaturation_acceleration", tba::AvailableAcceleration::momentum_wheel_desaturation_acceleration)
//      .value("solar_sail_acceleration", tba::AvailableAcceleration::solar_sail_acceleration)
//      .export_values();

}

}// namespace tudatpy
