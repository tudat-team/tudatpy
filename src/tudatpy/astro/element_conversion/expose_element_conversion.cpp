/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "tudat/astro/basic_astro/stateRepresentationConversions.h"
#include "tudat/astro/conversions.h"
#include "tudat/astro/ephemerides/rotationalEphemeris.h"
#include "tudat/math/basic.h"
#include "tudatpy/docstrings.h"

namespace py = pybind11;
namespace toec = tudat::orbital_element_conversions;
namespace tcc = tudat::coordinate_conversions;
namespace tla = tudat::linear_algebra;
namespace te = tudat::ephemerides;
namespace tba = tudat::basic_astrodynamics;
namespace tmg = tudat::mission_geometry;


PYBIND11_MODULE(expose_element_conversion, m) {
    py::enum_<toec::KeplerianElementIndices>(m, "KeplerianElementIndices")
        .value("semi_major_axis_index",
               toec::KeplerianElementIndices::semiMajorAxisIndex)
        .value("eccentricity_index",
               toec::KeplerianElementIndices::eccentricityIndex)
        .value("inclination_index",
               toec::KeplerianElementIndices::inclinationIndex)
        .value("argument_of_periapsis_index",
               toec::KeplerianElementIndices::argumentOfPeriapsisIndex)
        .value("longitude_of_ascending_node_index",
               toec::KeplerianElementIndices::longitudeOfAscendingNodeIndex)
        .value("true_anomaly_index",
               toec::KeplerianElementIndices::trueAnomalyIndex)
        .value("semi_latus_rectum_index",
               toec::KeplerianElementIndices::semiLatusRectumIndex)
        .export_values();

    py::enum_<toec::SphericalOrbitalStateElementIndices>(
        m, "SphericalOrbitalStateElementIndices")
        .value("radius_index",
               toec::SphericalOrbitalStateElementIndices::radiusIndex)
        .value("latitude_index",
               toec::SphericalOrbitalStateElementIndices::latitudeIndex)
        .value("longitude_index",
               toec::SphericalOrbitalStateElementIndices::longitudeIndex)
        .value("speed_index",
               toec::SphericalOrbitalStateElementIndices::speedIndex)
        .value("flight_path_index",
               toec::SphericalOrbitalStateElementIndices::flightPathIndex)
        .value("heading_angle_index",
               toec::SphericalOrbitalStateElementIndices::headingAngleIndex)
        .export_values();

    py::enum_<tcc::PositionElementTypes>(m, "PositionElementTypes")
        .value("cartesian_position_type",
               tcc::PositionElementTypes::cartesian_position)
        .value("spherical_position_type",
               tcc::PositionElementTypes::spherical_position)
        .value("geodetic_position_type",
               tcc::PositionElementTypes::geodetic_position)
        .export_values();
    /*!
     **************   KEPLER ELEMENTS  ******************
     */
    m.def("convert_position_elements", &tcc::convertPositionElements,
          py::arg("originalElements"), py::arg("original_elemet_types"),
          py::arg("new_element_types"), py::arg("shape_model"),
          py::arg("tolerance"),
          tudatpy::get_docstring("convert_position_elements").c_str());


    m.def("cartesian_to_keplerian",
          &toec::convertCartesianToKeplerianElements<double>,
          py::arg("cartesian_elements"), py::arg("gravitational_parameter"),
          R"doc(Convert Cartesian to Keplerian elements.

	.. note:: See module level documentation for the standard ordering
	          convention of Keplerian elements used.


	:param cartesian_elements:
		Cartesian state that is to be converted to Keplerian elements
	:param gravitational_parameter:
		Gravitational parameter of central body used for conversion
	:return:
		Keplerian elements, as computed from Cartesian element input.
)doc");


    m.def("keplerian_to_cartesian",
          py::overload_cast<const Eigen::Vector6d&, double>(
              &toec::convertKeplerianToCartesianElements<double>),
          py::arg("keplerian_elements"), py::arg("gravitational_parameter"),
          R"doc(Convert Keplerian elements to Cartesian.

	.. note:: See module level documentation for the standard ordering
	          convention of Keplerian elements used.


	:param keplerian_elements:
		Keplerian state that is to be converted to Cartesian elements
	:param gravitational_parameter:
		Gravitational parameter of central body used for conversion
	:return:
		Cartesian elements, as computed from Keplerian element input.
)doc");

    m.def("keplerian_to_cartesian_elementwise",
          py::overload_cast<double, double, double, double, double, double,
                            double>(
              &toec::convertKeplerianToCartesianElements<double>),
          py::arg("semi_major_axis"), py::arg("eccentricity"),
          py::arg("inclination"), py::arg("argument_of_periapsis"),
          py::arg("longitude_of_ascending_node"), py::arg("true_anomaly"),
          py::arg("gravitational_parameter"),
          R"doc(Convert Keplerian elements to Cartesian, with elementwise input.

	.. note:: The final Keplerian element is always the true anomaly.


	:param semi_major_axis:
		Semi-major axis (except if eccentricity = 1.0, then represents semi-latus rectum)
	:param eccentricity:
		Eccentricity
	:param inclination:
		Inclination
	:param argument_of_periapsis:
		Argument of periapsis
	:param longitude_of_ascending_node:
		Longitude of ascending node
	:param true_anomaly:
		True anomaly
	:param gravitational_parameter:
		Gravitational parameter of central body used for conversion
	:return:
		Cartesian elements, as computed from Keplerian element input.
)doc");

    m.def("mean_to_true_anomaly",
          &toec::convertMeanAnomalyToTrueAnomaly<double>,
          py::arg("eccentricity"), py::arg("mean_anomaly"),
          py::arg("use_default_initial_guess") = true,
          py::arg("non_default_initial_guess") = TUDAT_NAN,
          py::arg("root_finder") = nullptr,
          R"doc(Convert mean to true anomaly.

	Convert the mean anomaly of the orbit to its true anomaly. This conversion first converts mean to eccentric anomaly
	(hyperbolic eccentric anomaly, if eccentriciy is larger than 1, elliptical eccentric anomaly if it is smaller than 1), and subsequently to true anomaly.


	:param eccentricity:
		Value of the orbital eccentricity
	:param mean_anomaly:
		Value of the mean anomaly
	:param use_default_initial_guess:
		Boolean to determine whether the user-defined initial guess (for mean-to-eccentric anomaly conversion) is used, or an automatically generated one.
	:param non_default_initial_guess:
		User-defined initial guess for mean-to-eccentric anomaly conversion, to be used only if ``use_default_initial_guess`` is set to ``True``.
	:param root_finder:
		User-defined root finder, overriding default root-finding algorithm for mean-to-eccentric anomaly conversion (default is used if this input is left empty)
	:return:
		Value of the true anomaly
)doc");

    m.def("true_to_mean_anomaly",
          &toec::convertTrueAnomalyToMeanAnomaly<double>,
          py::arg("eccentricity"), py::arg("true_anomaly"),
          R"doc(Convert true to mean anomaly.

	Convert the true anomaly of the orbit to its mean anomaly. This conversion first converts true to eccentric anomaly
	(hyperbolic eccentric anomaly, if eccentriciy is larger than 1, elliptical eccentric anomaly if it is smaller than 1),
	and subsequently to mean anomaly.


	:param eccentricity:
		Value of the orbital eccentricity
	:param true_anomaly:
		Value of the true anomaly
	:return:
		Value of the mean anomaly
)doc");

    m.def("true_to_eccentric_anomaly",
          &toec::convertTrueAnomalyToEccentricAnomaly<double>,
          py::arg("true_anomaly"), py::arg("eccentricity"),
          R"doc(Convert true to eccentric anomaly.

	:param eccentricity:
		Value of the orbital eccentricity
	:param true_anomaly:
		Value of the true anomaly
	:return:
		Hyperbolic eccentric anomaly, if eccentriciy is larger than 1, elliptical eccentric anomaly if it is smaller than 1
)doc");

    m.def("eccentric_to_true_anomaly",
          &toec::convertEccentricAnomalyToTrueAnomaly<double>,
          py::arg("eccentric_anomaly"), py::arg("eccentricity"),
          R"doc(Convert eccentric to true anomaly.

	:param eccentricity:
		Value of the orbital eccentricity
	:param eccentric_anomaly:
		Hyperbolic eccentric anomaly, if eccentriciy is larger than 1, elliptical eccentric anomaly if it is smaller than 1
	:return:
		Value of the true anomaly
)doc");


    m.def("eccentric_to_mean_anomaly",
          &toec::convertEccentricAnomalyToMeanAnomaly<double>,
          py::arg("eccentric_anomaly"), py::arg("eccentricity"),
          R"doc(Convert eccentric to mean anomaly.

	:param eccentricity:
		Value of the orbital eccentricity
	:param eccentric_anomaly:
		Hyperbolic eccentric anomaly, if eccentriciy is larger than 1, elliptical eccentric anomaly if it is smaller than 1
	:return:
		Mean of the true anomaly
)doc");


    m.def("mean_to_eccentric_anomaly",
          &toec::convertMeanAnomalyToEccentricAnomaly<double>,
          py::arg("eccentricity"), py::arg("mean_anomaly"),
          py::arg("use_default_initial_guess") = true,
          py::arg("non_default_initial_guess") = TUDAT_NAN,
          py::arg("root_finder") = nullptr,
          R"doc(Convert mean to eccentric anomaly.

	:param eccentricity:
		Value of the orbital eccentricity
	:param mean_anomaly:
		Value of the mean anomaly
	:param use_default_initial_guess:
		Boolean to determine whether the user-defined initial guess is used for conversion, or an automatically generated one.
	:param non_default_initial_guess:
		User-defined initial guess for conversion, to be used only if ``use_default_initial_guess`` is set to ``True``.
	:param root_finder:
		User-defined root finder, overriding default root-finding algorithm for conversion (default is used if this input is left empty)
	:return:
		Value of the eccentric anomaly
)doc");


    m.def(
        "elapsed_time_to_delta_mean_anomaly",
        &toec::convertElapsedTimeToMeanAnomalyChange<double>,
        py::arg("elapsed_time"), py::arg("gravitational_parameter"),
        py::arg("semi_major_axis"),
        R"doc(Convert elapsed time to the corresponding change in mean anomaly along a Keplerian orbit.

	:param elapsed_time:
		Elapsed time (in seconds)
	:param gravitational_parameter:
		Gravitational parameter of central body
	:param semi_major_axis:
		Semi-major axis of orbit
	:return:
		Total change in mean anomaly along the Kepler orbit, accumulated in the provided time.
)doc");

    m.def(
        "delta_mean_anomaly_to_elapsed_time",
        &toec::convertMeanAnomalyChangeToElapsedTime<double>,
        py::arg("mean_anomaly_change"), py::arg("gravitational_parameter"),
        py::arg("semi_major_axis"),
        R"doc(Convert change in mean anomaly along a Keplerian orbit to the corresponding elapsed time.

	:param mean_anomaly_change:
		Total change in mean anomaly along the Kepler orbit
	:param gravitational_parameter:
		Gravitational parameter of central body
	:param semi_major_axis:
		Semi-major axis of orbit
	:return:
		Time required for the provided mean anomaly change to be accumulated
)doc");

    m.def(
        "mean_motion_to_semi_major_axis",
        &toec::convertEllipticalMeanMotionToSemiMajorAxis<double>,
        py::arg("mean_motion"), py::arg("gravitational_parameter"),
        R"doc(Convert mean motion to corresponding semi-major axis (in a Keplerian orbit).

	:param mean_motion:
		Orbital mean motion
	:param gravitational_parameter:
		Gravitational parameter of central body
	:return:
		Semi-major axis corresponding to mean motion
)doc");

    m.def(
        "semi_major_axis_to_mean_motion",
        &toec::convertSemiMajorAxisToEllipticalMeanMotion<double>,
        py::arg("semi_major_axis"), py::arg("gravitational_parameter"),
        R"doc(Convert semi-major axis to corresponding mean motion (along a Keplerian orbit).

	:param semi_major_axis:
		Semi-major axis of orbit
	:param gravitational_parameter:
		Gravitational parameter of central body
	:return:
		Semi-major axis corresponding to mean motion
)doc");


    /*!
     **************   MODIFIED EQUIONOCTIAL ELEMENTS  ******************
     */

    m.def("keplerian_to_mee_manual_singularity",
          py::overload_cast<const Eigen::Vector6d&, const bool>(
              &toec::convertKeplerianToModifiedEquinoctialElements<double>),
          py::arg("keplerian_elements"),
          py::arg("singularity_at_zero_inclination"),
          R"doc(Convert Keplerian to Modified equinoctial elements.

	Convert Keplerian to Modified equinoctial elements (without intermediate step to Cartesian elements). The singularity-flipping
	element :math:`I` is to be provided manually for this function
	.. note:: See module level documentation for the standard ordering
	          convention of Modified Equinoctial elements used.


	:param keplerian_elements:
		Keplerian elements that are to be converted to Modified equinoctial elements
	:param singularity_at_zero_inclination:
		Singularity at 0 degrees inclination if false, 180 degrees if true
	:return:
		Modified equinoctial elements, as computed from Keplerian element input.
)doc");

    m.def("keplerian_to_mee",
          py::overload_cast<const Eigen::Vector6d&>(
              &toec::convertKeplerianToModifiedEquinoctialElements<double>),
          py::arg("keplerian_elements"),
          R"doc(Convert Keplerian to Modified equinoctial elements.

	Convert Keplerian to Modified equinoctial elements (without intermediate step to Cartesian elements). The singularity-flipping
	element :math:`I` is computed automatically by this function (using :func:`flip_mee_singularity`)
	.. note:: See module level documentation for the standard ordering
	          convention of Modified Equinoctial elements used.


	:param keplerian_elements:
		Keplerian elements that are to be converted to Modified equinoctial elements
	:return:
		Modified equinoctial elements, as computed from Keplerian element input (with element :math:`I` defined by :func:`flip_mee_singularity`).
)doc");

    m.def(
        "flip_mee_singularity",
        py::overload_cast<const Eigen::Vector6d&>(&tmg::isOrbitRetrograde),
        py::arg("keplerian_elements"),
        R"doc(Function to determine 'optimal' location of the singularity-flipping modified equinoctial element.

	Function to determine 'optimal' location of the singularity-flipping modified equinoctial element :math:`I`, if orbit inclination is less than
	90 degrees, it puts the singularity at 180 degrees, if it is larger than 90 degrees, it puts it at 0 degrees.


	:param keplerian_elements:
		Keplerian elements that are to be converted to Modified equinoctial elements
	:return:
		Singularity at 0 degrees inclination if false, 180 degrees if true
)doc");

    m.def("mee_to_keplerian",
          &toec::convertModifiedEquinoctialToKeplerianElements<double>,
          py::arg("modified_equinoctial_elements"),
          py::arg("singularity_at_zero_inclination"),
          R"doc(Convert Modified equinoctial to Keplerian elements.

	Modified equinoctial elements to Keplerian (without intermediate step to Cartesian elements).
	.. note:: See module level documentation for the standard ordering
	          convention of Modified Equinoctial elements used.


	:param modified_equinoctial_elements:
		Modified equinoctial elements that are to be converted to Keplerian elements
	:param singularity_at_zero_inclination:
		Singularity at 0 degrees inclination if false, 180 degrees if true
	:return:
		Keplerian elements, as computed from Modified equinoctial element input.
)doc");

    m.def("cartesian_to_mee",
          py::overload_cast<const Eigen::Vector6d&, const double>(
              &toec::convertCartesianToModifiedEquinoctialElements<double>),
          py::arg("cartesian_elements"), py::arg("gravitational_parameter"),
          R"doc(Convert Cartesian to Modified equinoctial elements.

	Convery cartesian to Modified equinoctial elements. The singularity-flipping
	element :math:`I` is computed automatically by this function (using :func:`flip_mee_singularity`)
	.. note:: See module level documentation for the standard ordering
	          convention of Modified Equinoctial elements used.


	:param cartesian_elements:
		Cartesian elements that are to be converted to Modified equinoctial elements
	:param gravitational_parameter:
		Gravitational parameter of central body
	:return:
		Modified equinoctial elements, as computed from Cartesian element input.
)doc");

    m.def("cartesian_to_mee_manual_singularity",
          py::overload_cast<const Eigen::Vector6d&, const double, const bool>(
              &toec::convertCartesianToModifiedEquinoctialElements<double>),
          py::arg("cartesian_elements"), py::arg("gravitational_parameter"),
          py::arg("singularity_at_zero_inclination"),
          R"doc(Convert Cartesian to Modified equinoctial elements.

	Convery cartesian to Modified equinoctial elements. The singularity-flipping
	element :math:`I` is to be provided manually for this function
	.. note:: See module level documentation for the standard ordering
	          convention of Modified Equinoctial elements used.


	:param cartesian_elements:
		Cartesian elements that are to be converted to Modified equinoctial elements
	:param gravitational_parameter:
		Gravitational parameter of central body
	:param singularity_at_zero_inclination:
		Singularity at 0 degrees inclination if false, 180 degrees if true
	:return:
		Modified equinoctial elements, as computed from Cartesian element input.
)doc");

    m.def("mee_to_cartesian",
          py::overload_cast<const Eigen::Vector6d&, const double, const bool>(
              &toec::convertModifiedEquinoctialToCartesianElements<double>),
          py::arg("modified_equinoctial_elements"),
          py::arg("gravitational_parameter"),
          py::arg("singularity_at_zero_inclination"),
          R"doc(Convert Modified equinoctial to Cartesian elements.

	Convert Modified equinoctial to Cartesian elements
	.. note:: See module level documentation for the standard ordering
	          convention of Modified Equinoctial elements used.


	:param modified_equinoctial_elements:
		Modified equinoctial elements that are to be converted to Cartesian elements
	:param gravitational_parameter:
		Gravitational parameter of central body
	:param singularity_at_zero_inclination:
		Singularity at 0 degrees inclination if false, 180 degrees if true
	:return:
		Cartesian elements, as computed from Modified equinoctial element input.
)doc");

    /*!
     **************   SPHERICAL ELEMENTS  ******************
     */

    m.def("spherical_to_cartesian_elementwise",
          py::overload_cast<double, double, double, double, double, double>(
              &toec::convertSphericalOrbitalToCartesianState<double>),
          py::arg("radial_distance"), py::arg("latitude"), py::arg("longitude"),
          py::arg("speed"), py::arg("flight_path_angle"),
          py::arg("heading_angle"),
          R"doc(Convert Spherical elements to Cartesian, with elementwise input.

	:param radial_distance:
		Distance from origin of central body
	:param latitude:
		Central body-fixed latitude
	:param longitude:
		Central body-fixed longitude
	:param speed:
		Central body-fixed speed (norm of velocity vector). Note that this is *not* the norm of the inertial velocity
	:param flight_path_angle:
		Flight-path angle (of central body-fixed velocity vector)
	:param heading_angle:
		Heading angle (of central body-fixed velocity vector)
	:return:
		Cartesian elements, as computed from spherical element input.
)doc");

    m.def("spherical_to_cartesian",
          py::overload_cast<const Eigen::Vector6d&>(
              &toec::convertSphericalOrbitalToCartesianState<double>),
          py::arg("spherical_elements"),
          R"doc(Convert spherical elements to Cartesian.

	.. note:: See module level documentation for the standard ordering
	          convention of spherical state elements used.


	:param spherical_elements:
		Spherical state that is to be converted to Cartesian elements
	:return:
		Cartesian elements, as computed from spherical element input.
)doc");

    m.def("cartesian_to_spherical",
          &toec::convertCartesianToSphericalOrbitalState,
          py::arg("cartesian_elements"),
          R"doc(Convert Cartesian to spherical elements.

	.. note:: See module level documentation for the standard ordering
	          convention of spherical state elements used.


	:param cartesian_elements:
		Cartesian state that is to be converted to spherical elements
	:return:
		Spherial elements, as computed from Cartesian element input.
)doc");


    /*!
     **************   QUATERNIONS  ******************
     */

    m.def(
        "quaternion_entries_to_rotation_matrix",
        &tla::convertVectorQuaternionToMatrixFormat,
        py::arg("quaternion_entries"),
        R"doc(Converts an array of four quaternion elements to the equivalent rotation matrix.

	Function to convert an array of four quaternion elements to the equivalent rotation matrix. These quaternion elements
	are for instance used when propagating rotational dynamics in Tudat, and this function can be used to convert the
	numerical results to a usable rotation matrix. See `our user guide <https://tudat-space.readthedocs.io/en/latest/_src_user_guide/state_propagation/environment_setup/use_of_reference_frames.html#rotational-states>`_ for more details.


	:param quaternion_entries:
		Quaternion elements, as per the convention used in the `Eigen library <https://eigen.tuxfamily.org/dox/classEigen_1_1Quaternion.html>`_
	:return:
		Rotation matrix defining the equivalent rotation.
)doc");

    m.def(
        "rotation_matrix_to_quaternion_entries",
        &tla::convertMatrixToVectorQuaternionFormat, py::arg("rotation_matrix"),
        R"doc(Converts a rotation matrix to the equivalent array of four quaternion elements.

	Inverse function of :func:`quaternion_entries_to_rotation_matrix`.


	:param rotation_matrix:
		Rotation matrix
	:return:
		Equivalent quaternion elements, as per the convention used in the `Eigen library <https://eigen.tuxfamily.org/dox/classEigen_1_1Quaternion.html>`_
)doc");

    /*!
     **************   TLE  ******************

     */
    m.def("teme_state_to_j2000",
          py::overload_cast<double, Eigen::Vector6d>(
              &tba::convertStateFromTEMEtoJ2000),
          py::arg("epoch"), py::arg("teme_state"));

    m.def("teme_state_to_eclipj2000",
          py::overload_cast<double, Eigen::Vector6d>(
              &tba::convertStateFromTEMEtoEclipJ2000),
          py::arg("epoch"), py::arg("teme_state"));

    m.def("j2000_state_to_teme",
          py::overload_cast<double, Eigen::Vector6d>(
              &tba::convertStateFromJ2000ToTEME),
          py::arg("epoch"), py::arg("j2000_state"));

    m.def("eclipj2000_state_to_teme",
          py::overload_cast<double, Eigen::Vector6d>(
              &tba::convertStateFromEclipJ2000ToTEME),
          py::arg("epoch"), py::arg("eclipj2000_state"));
}
