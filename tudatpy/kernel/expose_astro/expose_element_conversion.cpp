/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudatpy/docstrings.h"

#include "expose_element_conversion.h"

#include <tudat/math/basic.h>
#include <tudat/astro/conversions.h>
#include <tudat/astro/basic_astro/stateRepresentationConversions.h>
#include <tudat/astro/ephemerides/rotationalEphemeris.h>

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
namespace toec = tudat::orbital_element_conversions;
namespace tcc = tudat::coordinate_conversions;
namespace tla = tudat::linear_algebra;
namespace te = tudat::ephemerides;
namespace tba = tudat::basic_astrodynamics;
namespace tmg = tudat::mission_geometry;

namespace tudatpy {

namespace astro {
namespace element_conversion {

void expose_element_conversion(py::module &m) {


    py::enum_<toec::KeplerianElementIndices>(m, "KeplerianElementIndices")
            .value("semi_major_axis_index", toec::KeplerianElementIndices::semiMajorAxisIndex)
            .value("eccentricity_index", toec::KeplerianElementIndices::eccentricityIndex)
            .value("inclination_index", toec::KeplerianElementIndices::inclinationIndex)
            .value("argument_of_periapsis_index", toec::KeplerianElementIndices::argumentOfPeriapsisIndex)
            .value("longitude_of_ascending_node_index", toec::KeplerianElementIndices::longitudeOfAscendingNodeIndex)
            .value("true_anomaly_index", toec::KeplerianElementIndices::trueAnomalyIndex)
            .value("semi_latus_rectum_index", toec::KeplerianElementIndices::semiLatusRectumIndex)
            .export_values();

    py::enum_<toec::SphericalOrbitalStateElementIndices>(m, "SphericalOrbitalStateElementIndices")
            .value("radius_index", toec::SphericalOrbitalStateElementIndices::radiusIndex)
            .value("latitude_index", toec::SphericalOrbitalStateElementIndices::latitudeIndex)
            .value("longitude_index", toec::SphericalOrbitalStateElementIndices::longitudeIndex)
            .value("speed_index", toec::SphericalOrbitalStateElementIndices::speedIndex)
            .value("flight_path_index", toec::SphericalOrbitalStateElementIndices::flightPathIndex)
            .value("heading_angle_index", toec::SphericalOrbitalStateElementIndices::headingAngleIndex)
            .export_values();

    py::enum_<tcc::PositionElementTypes>(m, "PositionElementTypes")
            .value("cartesian_position_type", tcc::PositionElementTypes::cartesian_position)
            .value("spherical_position_type", tcc::PositionElementTypes::spherical_position)
            .value("geodetic_position_type", tcc::PositionElementTypes::geodetic_position)
            .export_values();
    /*!
     **************   KEPLER ELEMENTS  ******************
     */

    m.def("cartesian_to_keplerian",
          &toec::convertCartesianToKeplerianElements< double >,
          py::arg("cartesian_elements"),
          py::arg("gravitational_parameter"),
           get_docstring("cartesian_to_keplerian").c_str());


    m.def("keplerian_to_cartesian",
          py::overload_cast< const Eigen::Vector6d&, double >(
              &toec::convertKeplerianToCartesianElements< double > ),
          py::arg("keplerian_elements"),
          py::arg("gravitational_parameter"),
           get_docstring("keplerian_to_cartesian").c_str());

    m.def("keplerian_to_cartesian_elementwise",
          py::overload_cast<
          double, double, double, double, double, double, double >(
              &toec::convertKeplerianToCartesianElements< double > ),
          py::arg("semi_major_axis"),
          py::arg("eccentricity"),
          py::arg("inclination"),
          py::arg("argument_of_periapsis"),
          py::arg("longitude_of_ascending_node"),
          py::arg("true_anomaly"),
          py::arg("gravitational_parameter"),
           get_docstring("keplerian_to_cartesian_elementwise").c_str());

    m.def("mean_to_true_anomaly",
          &toec::convertMeanAnomalyToTrueAnomaly< double >,
          py::arg("eccentricity"),
          py::arg("mean_anomaly"),
          py::arg("use_default_initial_guess") = true,
          py::arg("non_default_initial_guess") = TUDAT_NAN,
          py::arg("root_finder") = nullptr ,
           get_docstring("mean_to_true_anomaly").c_str());

    m.def("true_to_mean_anomaly",
          &toec::convertTrueAnomalyToMeanAnomaly< double >,
          py::arg("eccentricity"),
          py::arg("true_anomaly") ,
           get_docstring("true_to_mean_anomaly").c_str());

    m.def("true_to_eccentric_anomaly",
          &toec::convertTrueAnomalyToEccentricAnomaly< double >,
          py::arg("true_anomaly"),
          py::arg("eccentricity") ,
           get_docstring("true_to_eccentric_anomaly").c_str());

    m.def("eccentric_to_true_anomaly",
          &toec::convertEccentricAnomalyToTrueAnomaly< double >,
          py::arg("eccentric_anomaly"),
          py::arg("eccentricity") ,
           get_docstring("eccentric_to_true_anomaly").c_str());


    m.def("eccentric_to_mean_anomaly",
          &toec::convertEccentricAnomalyToMeanAnomaly< double >,
          py::arg("eccentric_anomaly"),
          py::arg("eccentricity") ,
           get_docstring("eccentric_to_mean_anomaly").c_str());


    m.def("mean_to_eccentric_anomaly",
          &toec::convertMeanAnomalyToEccentricAnomaly< double >,
          py::arg("eccentricity"),
          py::arg("mean_anomaly"),
          py::arg("use_default_initial_guess") = true,
          py::arg("non_default_initial_guess") = TUDAT_NAN,
          py::arg("root_finder") = nullptr ,
           get_docstring("mean_to_eccentric_anomaly").c_str());


    m.def("elapsed_time_to_delta_mean_anomaly",
          &toec::convertElapsedTimeToMeanAnomalyChange< double >,
          py::arg("elapsed_time"),
          py::arg("gravitational_parameter"),
          py::arg("semi_major_axis") ,
           get_docstring("elapsed_time_to_delta_mean_anomaly").c_str());

    m.def("delta_mean_anomaly_to_elapsed_time",
          &toec::convertMeanAnomalyChangeToElapsedTime< double >,
          py::arg("mean_anomaly_change"),
          py::arg("gravitational_parameter"),
          py::arg("semi_major_axis") ,
           get_docstring("delta_mean_anomaly_to_elapsed_time").c_str());

    m.def("mean_motion_to_semi_major_axis",
          &toec::convertEllipticalMeanMotionToSemiMajorAxis< double >,
          py::arg("mean_motion"),
          py::arg("gravitational_parameter") ,
           get_docstring("mean_motion_to_semi_major_axis").c_str());

    m.def("semi_major_axis_to_mean_motion",
          &toec::convertSemiMajorAxisToEllipticalMeanMotion< double >,
          py::arg("semi_major_axis"),
          py::arg("gravitational_parameter") ,
           get_docstring("semi_major_axis_to_mean_motion").c_str());


    /*!
     **************   MODIFIED EQUIONOCTIAL ELEMENTS  ******************
     */

    m.def("keplerian_to_mee_manual_singularity",
          py::overload_cast< const Eigen::Vector6d&, const bool >(
              &toec::convertKeplerianToModifiedEquinoctialElements< double > ),
          py::arg("keplerian_elements"),
          py::arg("singularity_at_zero_inclination"),
          get_docstring("keplerian_to_mee").c_str());

    m.def("keplerian_to_mee",
          py::overload_cast< const Eigen::Vector6d& >(
              &toec::convertKeplerianToModifiedEquinoctialElements< double > ),
          py::arg("keplerian_elements"),
          get_docstring("keplerian_to_mee").c_str());

//    m.def("flip_mee_singularity",
//          &tmg::isOrbitRetrograde,
//          py::arg("keplerian_elements"),
//          get_docstring("flip_mee_singularity").c_str());

    m.def("mee_to_keplerian",
          &toec::convertModifiedEquinoctialToKeplerianElements< double >,
          py::arg("modified_equinoctial_elements"),
          py::arg("singularity_at_zero_inclination"),
          get_docstring("mee_to_keplerian").c_str());

    m.def("cartesian_to_mee",
          py::overload_cast< const Eigen::Vector6d&, const double >(
              &toec::convertCartesianToModifiedEquinoctialElements< double > ),
          py::arg("cartesian_elements"),
          py::arg("gravitational_parameter"),
          get_docstring("cartesian_to_mee").c_str());

    m.def("cartesian_to_mee_manual_singularity",
          py::overload_cast< const Eigen::Vector6d&, const double, const bool >(
              &toec::convertCartesianToModifiedEquinoctialElements< double > ),
          py::arg("cartesian_elements"),
          py::arg("gravitational_parameter"),
          py::arg("singularity_at_zero_inclination"),
          get_docstring("cartesian_to_mee_manual_singularity").c_str());

    m.def("mee_to_cartesian",
          py::overload_cast< const Eigen::Vector6d&, const double, const bool >(
              &toec::convertModifiedEquinoctialToCartesianElements< double > ),
          py::arg("modified_equinoctial_elements"),
          py::arg("gravitational_parameter"),
          py::arg("singularity_at_zero_inclination"),
          get_docstring("mee_to_cartesian").c_str());

    /*!
     **************   SPHERICAL ELEMENTS  ******************
     */

    m.def("spherical_to_cartesian_elementwise",
          py::overload_cast<
          double, double, double, double, double, double >(
              &toec::convertSphericalOrbitalToCartesianState< double > ),
          py::arg("radial_distance"),
          py::arg("latitude"),
          py::arg("longitude"),
          py::arg("speed"),
          py::arg("flight_path_angle"),
          py::arg("heading_angle"),
          get_docstring("spherical_to_cartesian_elementwise").c_str());

    m.def("spherical_to_cartesian",
          py::overload_cast<
          const Eigen::Vector6d& >(
              &toec::convertSphericalOrbitalToCartesianState< double > ),
          py::arg("spherical_elements"),
          get_docstring("spherical_to_cartesian").c_str());

    m.def("cartesian_to_spherical",
          &toec::convertCartesianToSphericalOrbitalState,
          py::arg("cartesian_elements"),
          get_docstring("cartesian_to_spherical").c_str());


    /*!
     **************   QUATERNIONS  ******************
     */

    m.def("quaternion_entries_to_rotation_matrix",
          &tla::convertVectorQuaternionToMatrixFormat,
          py::arg( "quaternion_entries" ) ,
          get_docstring("quaternion_entries_to_rotation_matrix").c_str());

    m.def("rotation_matrix_to_quaternion_entries",
          &tla::convertMatrixToVectorQuaternionFormat,
          py::arg( "rotation_matrix" ) ,
          get_docstring("rotation_matrix_to_quaternion_entries").c_str());

    /*!
     **************   TLE  ******************

     */
    m.def("teme_state_to_j2000",
          py::overload_cast< double, Eigen::Vector6d >(
                  &tba::convertStateFromTEMEtoJ2000 ),
          py::arg("epoch"),
          py::arg("teme_state") );

    m.def("teme_state_to_eclipj2000",
          py::overload_cast< double, Eigen::Vector6d >(
                  &tba::convertStateFromTEMEtoEclipJ2000 ),
          py::arg("epoch"),
          py::arg("teme_state") );

    m.def("j2000_state_to_teme",
          py::overload_cast< double, Eigen::Vector6d >(
                  &tba::convertStateFromJ2000ToTEME ),
          py::arg("epoch"),
          py::arg("j2000_state") );

    m.def("eclipj2000_state_to_teme",
          py::overload_cast< double, Eigen::Vector6d >(
                  &tba::convertStateFromEclipJ2000ToTEME ),
          py::arg("epoch"),
          py::arg("eclipj2000_state") );


}
} // namespace element_conversion
} // namespace astro
}// namespace tudatpy
