/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_reference_frames.h"

#include <tudat/astro/reference_frames.h>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

namespace trf = tudat::reference_frames;

namespace py = pybind11;

namespace tudatpy {

void expose_frames(py::module &m) {

    py::enum_<trf::AerodynamicsReferenceFrameAngles>(m, "AerodynamicsReferenceFrameAngles" )
            .value("latitude_angle", trf::AerodynamicsReferenceFrameAngles::latitude_angle)
            .value("longitude_angle", trf::AerodynamicsReferenceFrameAngles::longitude_angle)
            .value("heading_angle", trf::AerodynamicsReferenceFrameAngles::heading_angle)
            .value("flight_path_angle", trf::AerodynamicsReferenceFrameAngles::flight_path_angle)
            .value("angle_of_attack", trf::AerodynamicsReferenceFrameAngles::angle_of_attack)
            .value("angle_of_sideslip", trf::AerodynamicsReferenceFrameAngles::angle_of_sideslip)
            .value("bank_angle", trf::AerodynamicsReferenceFrameAngles::bank_angle)
            .export_values( );

    py::enum_<trf::AerodynamicsReferenceFrames>(m, "AerodynamicsReferenceFrames" )
            .value("inertial_frame", trf::AerodynamicsReferenceFrames::inertial_frame)
            .value("corotating_frame", trf::AerodynamicsReferenceFrames::corotating_frame)
            .value("vertical_frame", trf::AerodynamicsReferenceFrames::vertical_frame)
            .value("trajectory_frame", trf::AerodynamicsReferenceFrames::trajectory_frame)
            .value("aerodynamic_frame", trf::AerodynamicsReferenceFrames::aerodynamic_frame)
            .value("body_frame", trf::AerodynamicsReferenceFrames::body_frame)
            .export_values( );

    py::class_<trf::AerodynamicAngleCalculator,
            std::shared_ptr<trf::AerodynamicAngleCalculator>>(m, "AerodynamicAngleCalculator")
            .def("set_orientation_angle_functions",
                 py::overload_cast<
                 const std::function<double()>,
                 const std::function<double()>,
                 const std::function<double()>,
                 const std::function<void(const double)>>(&trf::AerodynamicAngleCalculator::setOrientationAngleFunctions),
                 py::arg("angle_of_attack_function") = std::function<double()>(),       // <pybind11/functional.h>
                 py::arg("angle_of_sideslip_function") = std::function<double()>(),     // <pybind11/functional.h>
                 py::arg("bank_angle_function") = std::function<double()>(),            // <pybind11/functional.h>
                 py::arg("angle_update_function") = std::function<void(const double)>(),// <pybind11/functional.h>
                 "<no_doc>")
            .def("set_orientation_angle_functions",
                 py::overload_cast<
                 const double,
                 const double,
                 const double>(&trf::AerodynamicAngleCalculator::setOrientationAngleFunctions),
                 py::arg("angle_of_attack") = TUDAT_NAN,
                 py::arg("angle_of_sideslip") = TUDAT_NAN,
                 py::arg("bank_angle") = TUDAT_NAN,
                 "<no_doc>")
            .def("get_rotation_quaternion_between_frames",
                 &trf::AerodynamicAngleCalculator::getRotationQuaternionBetweenFrames,
                 py::arg("original_frame"),
                 py::arg("target_frame" ) )
            .def("get_rotation_matrix_between_frames",
                 &trf::AerodynamicAngleCalculator::getRotationMatrixBetweenFrames,
                 py::arg("original_frame"),
                 py::arg("target_frame" ) )
            .def("get_angle",
                 &trf::AerodynamicAngleCalculator::getAerodynamicAngle,
                 py::arg("angle_type") );


    m.def("inertial_to_rsw_rotation_matrix",
          &trf::getInertialToRswSatelliteCenteredFrameRotationMatrix,
          py::arg("inertial_cartesian_state") );


    m.def("lvlh_to_inertial_rotation_matrix",
          &trf::getVelocityBasedLvlhToInertialRotation,
          py::arg("inertial_cartesian_state"),
          py::arg("central_body_cartesian_state") = Eigen::Vector6d::Zero( ),
          py::arg("local_y_points_away_from_central_body" ) = true);


    m.def("inertial_to_body_fixed_rotation_matrix",
         py::overload_cast< const double, const double, const double >(
              &trf::getInertialToPlanetocentricFrameTransformationMatrix ),
          py::arg("pole_declination"),
          py::arg("pole_right_ascension"),
          py::arg("prime_meridian_longitude") );

    m.def("corotating_to_inertial",
          py::overload_cast< const double, const double, const double >(
              &trf::getRotatingPlanetocentricToInertialFrameTransformationQuaternion ),
          py::arg("pole_declination"),
          py::arg("pole_right_ascension"),
          py::arg("pole_meridian") );


}

}// namespace tudatpy
