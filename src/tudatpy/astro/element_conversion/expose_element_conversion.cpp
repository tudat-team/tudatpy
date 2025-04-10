/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_element_conversion.h"

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <tudat/astro/basic_astro/attitudeElementConversions.h>
#include <tudat/astro/basic_astro/stateRepresentationConversions.h>
#include <tudat/astro/conversions.h>
#include <tudat/astro/ephemerides/rotationalEphemeris.h>
#include <tudat/interface/spice/spiceInterface.h>
#include <tudat/math/basic.h>

namespace py = pybind11;
namespace toec = tudat::orbital_element_conversions;
namespace tcc = tudat::coordinate_conversions;
namespace tla = tudat::linear_algebra;
namespace te = tudat::ephemerides;
namespace tba = tudat::basic_astrodynamics;
namespace tmg = tudat::mission_geometry;
namespace tsi = tudat::spice_interface;

namespace tudatpy
{

namespace astro
{
namespace element_conversion
{

void expose_element_conversion( py::module& m )
{
    py::module_::import( "tudatpy.math.root_finders" );
    py::enum_< toec::KeplerianElementIndices >( m,
                                                "KeplerianElementIndices",
                                                R"doc(
         Enumeration for indices of Keplerian elements"
      )doc" )
            .value( "semi_major_axis_index",
                    toec::KeplerianElementIndices::semiMajorAxisIndex,
                    R"doc(
 Element 0 in vector of Keplerian elements (for eccentricity not equal to 1.0)
      )doc" )
            .value( "semi_latus_rectum_index",
                    toec::KeplerianElementIndices::semiLatusRectumIndex,
                    R"doc(
 Element 0 in vector of Keplerian elements (for eccentricity equal to 1.0)
      )doc" )
            .value( "eccentricity_index",
                    toec::KeplerianElementIndices::eccentricityIndex,
                    R"doc(
 Element 1 in vector of Keplerian elements
      )doc" )
            .value( "inclination_index",
                    toec::KeplerianElementIndices::inclinationIndex,
                    R"doc(
 Element 2 in vector of Keplerian elements
      )doc" )
            .value( "argument_of_periapsis_index",
                    toec::KeplerianElementIndices::argumentOfPeriapsisIndex,
                    R"doc(
 Element 3 in vector of Keplerian elements
      )doc" )
            .value( "longitude_of_ascending_node_index",
                    toec::KeplerianElementIndices::longitudeOfAscendingNodeIndex,
                    R"doc(
 Element 4 in vector of Keplerian elements
      )doc" )
            .value( "true_anomaly_index",
                    toec::KeplerianElementIndices::trueAnomalyIndex,
                    R"doc(
 Element 5 in vector of Keplerian elements
      )doc" )
            .export_values( );

    py::enum_< toec::SphericalOrbitalStateElementIndices >( m,
                                                            "SphericalOrbitalStateElementIndices",
                                                            R"doc(
         Enumeration for indices of spherical orbital state elements"
      )doc" )
            .value( "radius_index",
                    toec::SphericalOrbitalStateElementIndices::radiusIndex,
                    R"doc(
 Element 0 in vector of spherical orbital state elements
      )doc" )
            .value( "latitude_index",
                    toec::SphericalOrbitalStateElementIndices::latitudeIndex,
                    R"doc(
 Element 1 in vector of spherical orbital state elements
      )doc" )
            .value( "longitude_index",
                    toec::SphericalOrbitalStateElementIndices::longitudeIndex,
                    R"doc(
 Element 2 in vector of spherical orbital state elements
      )doc" )
            .value( "speed_index",
                    toec::SphericalOrbitalStateElementIndices::speedIndex,
                    R"doc(
 Element 3 in vector of spherical orbital state elements
      )doc" )
            .value( "flight_path_index",
                    toec::SphericalOrbitalStateElementIndices::flightPathIndex,
                    R"doc(
 Element 4 in vector of spherical orbital state elements
      )doc" )
            .value( "heading_angle_index",
                    toec::SphericalOrbitalStateElementIndices::headingAngleIndex,
                    R"doc(
 Element 5 in vector of spherical orbital state elements
      )doc" )
            .export_values( );

    py::enum_< tcc::PositionElementTypes >( m,
                                            "PositionElementTypes",
                                            R"doc(
 Enumeration describing different types of position element types (typically used for body-centered, body-0fixed position)
      )doc" )
            .value( "cartesian_position_type", tcc::PositionElementTypes::cartesian_position )
            .value( "spherical_position_type", tcc::PositionElementTypes::spherical_position )
            .value( "geodetic_position_type", tcc::PositionElementTypes::geodetic_position )
            .export_values( );

    m.def( "convert_cartesian_to_geodetic_coordinates",
           &tcc::convertCartesianToGeodeticCoordinates,
           py::arg( "cartesian_coordinates" ),
           py::arg( "equatorial_radius" ),
           py::arg( "flattening" ),
           py::arg( "tolerance" ),
           R"doc(

 Convert cartesian position to geodetic position, given equatorial radius and flattening (based on the Body Shape Model).

 Parameters
 ----------
 cartesian_coordinates : numpy.ndarray
     Position from which the conversion is to be performed
 equatorial_radius : float
 flattening : float
 tolerance : float
     Tolerance (in meters) used as convergence criterion for converting to/from geodetic altitude
 Returns
 -------
 numpy.ndarray
     Geodetic coordinates, as computed from Cartesian element input.



     )doc" );

    m.def( "convert_position_elements",
           &tcc::convertPositionElements,
           py::arg( "original_elements" ),
           py::arg( "original_element_type" ),
           py::arg( "new_element_type" ),
           py::arg( "shape_model" ),
           py::arg( "tolerance" ),
           R"doc(

 Convert position from one element set to another

 Parameters
 ----------
 original_elements : numpy.ndarray
     Position from which the conversion is to be performed
 original_element_type : PositionElementTypes
     Element type in which ``original_elements`` is provided
 new_element_type : PositionElementTypes
     Element type to which the ``original_elements`` are to be converted
 shape_model : BodyShapeModel
     Shape model of body sed to transform altitudes (w.r.t. this shape)
     to/from distances from the center of the body (can be set to ``None`` if requested conversion does not require this model)
 tolerance : float
     Tolerance (in meters) used as convergence criterion for converting to/from geodetic altitude
 Returns
 -------
 numpy.ndarray
     Keplerian elements, as computed from Cartesian element input.






     )doc" );

    /*!
     **************   KEPLER ELEMENTS  ******************
     */

    m.def( "cartesian_to_keplerian",
           &toec::convertCartesianToKeplerianElements< double >,
           py::arg( "cartesian_elements" ),
           py::arg( "gravitational_parameter" ),
           R"doc(

 Convert Cartesian to Keplerian elements.

 .. note:: See module level documentation for the standard ordering
           convention of Keplerian elements used.


 Parameters
 ----------
 cartesian_elements : numpy.ndarray
     Cartesian state that is to be converted to Keplerian elements
 gravitational_parameter : float
     Gravitational parameter of central body used for conversion
 Returns
 -------
 numpy.ndarray
     Keplerian elements, as computed from Cartesian element input.






     )doc" );

    m.def( "keplerian_to_cartesian",
           py::overload_cast< const Eigen::Vector6d&, double >(
                   &toec::convertKeplerianToCartesianElements< double > ),
           py::arg( "keplerian_elements" ),
           py::arg( "gravitational_parameter" ),
           R"doc(

 Convert Keplerian elements to Cartesian.

 .. note:: See module level documentation for the standard ordering
           convention of Keplerian elements used.


 Parameters
 ----------
 keplerian_elements : numpy.ndarray
     Keplerian state that is to be converted to Cartesian elements
 gravitational_parameter : float
     Gravitational parameter of central body used for conversion
 Returns
 -------
 numpy.ndarray
     Cartesian elements, as computed from Keplerian element input.






     )doc" );

    m.def( "keplerian_to_cartesian_elementwise",
           py::overload_cast< double, double, double, double, double, double, double >(
                   &toec::convertKeplerianToCartesianElements< double > ),
           py::arg( "semi_major_axis" ),
           py::arg( "eccentricity" ),
           py::arg( "inclination" ),
           py::arg( "argument_of_periapsis" ),
           py::arg( "longitude_of_ascending_node" ),
           py::arg( "true_anomaly" ),
           py::arg( "gravitational_parameter" ),
           R"doc(

 Convert Keplerian elements to Cartesian, with elementwise input.

 .. note:: The final Keplerian element is always the true anomaly.


 Parameters
 ----------
 semi_major_axis : float
     Semi-major axis (except if eccentricity = 1.0, then represents semi-latus rectum)
 eccentricity : float
     Eccentricity
 inclination : float
     Inclination
 argument_of_periapsis : float
     Argument of periapsis
 longitude_of_ascending_node : float
     Longitude of ascending node
 true_anomaly : float
     True anomaly
 gravitational_parameter : float
     Gravitational parameter of central body used for conversion
 Returns
 -------
 numpy.ndarray
     Cartesian elements, as computed from Keplerian element input.






     )doc" );

    m.def( "cartesian_to_usm_em",
           &toec::convertCartesianToUnifiedStateModelExponentialMapElements,
           py::arg( "cartesian_elements" ),
           py::arg( "gravitational_parameter" ),
           R"doc(

 Convert Cartesian elements to Unified State Model (USM) elements with Exponential map for rotational coordinates.

 .. note:: See `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/available_state_definitions_conversions.html#unified-state-model-elements>`_ for details on Unified State Model elements


 Parameters
 ----------
 cartesian_elements : numpy.ndarray
     Cartesian state that is to be converted to USM elements
 gravitational_parameter : float
     Gravitational parameter of central body used for conversion
 Returns
 -------
 numpy.ndarray
     USM elements using exponential map, as computed from Cartesian element input.






     )doc" );

    m.def( "cartesian_to_usm_7",
           &toec::convertCartesianToUnifiedStateModelQuaternionsElements,
           py::arg( "cartesian_elements" ),
           py::arg( "gravitational_parameter" ),
           R"doc(

 Convert Cartesian elements to Unified State Model (USM) elements with quaternion for rotational coordinates.

 .. note:: See `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/available_state_definitions_conversions.html#unified-state-model-elements>`_ for details on Unified State Model elements


 Parameters
 ----------
 cartesian_elements : numpy.ndarray
     Cartesian state that is to be converted to USM elements
 gravitational_parameter : float
     Gravitational parameter of central body used for conversion
 Returns
 -------
 numpy.ndarray
     USM elements using quaternion, as computed from Cartesian element input.






     )doc" );

    m.def( "cartesian_to_usm_6",
           &toec::convertCartesianToUnifiedStateModelModifiedRodriguesParameterElements,
           py::arg( "cartesian_elements" ),
           py::arg( "gravitational_parameter" ),
           R"doc(

 Convert Cartesian elements to Unified State Model (USM) elements with Modified Rodrigues parameters map for rotational coordinates.

 .. note:: See `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/available_state_definitions_conversions.html#unified-state-model-elements>`_ for details on Unified State Model elements


 Parameters
 ----------
 cartesian_elements : numpy.ndarray
     Cartesian state that is to be converted to USM elements
 gravitational_parameter : float
     Gravitational parameter of central body used for conversion
 Returns
 -------
 numpy.ndarray
     USM elements using Modified Rodrigues parameters, as computed from Cartesian element input.






     )doc" );

    m.def( "usm_em_to_cartesian",
           &toec::convertUnifiedStateModelExponentialMapToCartesianElements,
           py::arg( "usm_em_elements" ),
           py::arg( "gravitational_parameter" ),
           R"doc(

 Convert Unified State Model (USM) elements with Exponential map for rotational coordinates to Cartesian elements.

 .. note:: See `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/available_state_definitions_conversions.html#unified-state-model-elements>`_ for details on Unified State Model elements


 Parameters
 ----------
 usm_em_elements : numpy.ndarray
     USM elements using exponential map that is to be converted to Cartesian elements
 gravitational_parameter : float
     Gravitational parameter of central body used for conversion
 Returns
 -------
 numpy.ndarray
     Cartesian elements, as computed from USM element input.






     )doc" );

    m.def( "usm_7_to_cartesian",
           &toec::convertUnifiedStateModelQuaternionsToCartesianElements,
           py::arg( "usm_7_elements" ),
           py::arg( "gravitational_parameter" ),
           py::arg( "normalize_quaternion" ) = true,
           R"doc(

 Convert Unified State Model (USM) elements with quaternion for rotational coordinates to Cartesian elements.

 .. note:: See `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/available_state_definitions_conversions.html#unified-state-model-elements>`_ for details on Unified State Model elements


 Parameters
 ----------
 usm_7_elements : numpy.ndarray
     USM elements using quaternion that is to be converted to Cartesian elements
 gravitational_parameter : float
     Gravitational parameter of central body used for conversion
 Returns
 -------
 numpy.ndarray
     Cartesian elements, as computed from USM element input.






     )doc" );

    m.def( "usm_6_to_cartesian",
           &toec::convertUnifiedStateModelModifiedRodriguesParametersToCartesianElements,
           py::arg( "usm_6_elements" ),
           py::arg( "gravitational_parameter" ),
           R"doc(

 Convert Unified State Model (USM) elements with Modified Rodrigues parameters for rotational coordinates to Cartesian elements.

 .. note:: See `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/available_state_definitions_conversions.html#unified-state-model-elements>`_ for details on Unified State Model elements


 Parameters
 ----------
 usm_6_elements : numpy.ndarray
     USM elements using Modified Rodrigues parameters that is to be converted to Cartesian elements
 gravitational_parameter : float
     Gravitational parameter of central body used for conversion
 Returns
 -------
 numpy.ndarray
     Cartesian elements, as computed from USM element input.






     )doc" );

    m.def( "mean_to_true_anomaly",
           &toec::convertMeanAnomalyToTrueAnomaly< double >,
           py::arg( "eccentricity" ),
           py::arg( "mean_anomaly" ),
           py::arg( "use_default_initial_guess" ) = true,
           py::arg( "non_default_initial_guess" ) = TUDAT_NAN,
           py::arg( "root_finder" ) = nullptr,
           R"doc(

 Convert mean to true anomaly.

 Convert the mean anomaly of the orbit to its true anomaly. This conversion first converts mean to eccentric anomaly
 (hyperbolic eccentric anomaly, if eccentricity is larger than 1, elliptical eccentric anomaly if it is smaller than 1), and subsequently to true anomaly.


 Parameters
 ----------
 eccentricity : float
     Value of the orbital eccentricity
 mean_anomaly : float
     Value of the mean anomaly
 use_default_initial_guess : bool, default = True
     Boolean to determine whether the user-defined initial guess (for mean-to-eccentric anomaly conversion) is used, or an automatically generated one.
 non_default_initial_guess : float, default = NaN
     User-defined initial guess for mean-to-eccentric anomaly conversion, to be used only if ``use_default_initial_guess`` is set to ``False``.
 root_finder : RootFinder, default = None
     User-defined root finder, overriding default root-finding algorithm for mean-to-eccentric anomaly conversion (default is used if this input is left empty)
 Returns
 -------
 float
     Value of the true anomaly






     )doc" );

    m.def( "true_to_mean_anomaly",
           &toec::convertTrueAnomalyToMeanAnomaly< double >,
           py::arg( "eccentricity" ),
           py::arg( "true_anomaly" ),
           R"doc(

 Convert true to mean anomaly.

 Convert the true anomaly of the orbit to its mean anomaly. This conversion first converts true to eccentric anomaly
 (hyperbolic eccentric anomaly, if eccentricity is larger than 1, elliptical eccentric anomaly if it is smaller than 1),
 and subsequently to mean anomaly.


 Parameters
 ----------
 eccentricity : float
     Value of the orbital eccentricity
 true_anomaly : float
     Value of the true anomaly
 Returns
 -------
 float
     Value of the mean anomaly






     )doc" );

    m.def( "true_to_eccentric_anomaly",
           &toec::convertTrueAnomalyToEccentricAnomaly< double >,
           py::arg( "true_anomaly" ),
           py::arg( "eccentricity" ),
           R"doc(

 Convert true to eccentric anomaly.


 Parameters
 ----------
 eccentricity : float
     Value of the orbital eccentricity
 true_anomaly : float
     Value of the true anomaly
 Returns
 -------
 float
     Hyperbolic eccentric anomaly, if eccentricity is larger than 1, elliptical eccentric anomaly if it is smaller than 1






     )doc" );

    m.def( "eccentric_to_true_anomaly",
           &toec::convertEccentricAnomalyToTrueAnomaly< double >,
           py::arg( "eccentric_anomaly" ),
           py::arg( "eccentricity" ),
           R"doc(

 Convert eccentric to true anomaly.


 Parameters
 ----------
 eccentric_anomaly : float
     Hyperbolic eccentric anomaly, if eccentricity is larger than 1, elliptical eccentric anomaly if it is smaller than 1
 eccentricity : float
     Value of the orbital eccentricity
 Returns
 -------
 float
     Value of the true anomaly






     )doc" );

    m.def( "eccentric_to_mean_anomaly",
           &toec::convertEccentricAnomalyToMeanAnomaly< double >,
           py::arg( "eccentric_anomaly" ),
           py::arg( "eccentricity" ),
           R"doc(

 Convert eccentric to mean anomaly.


 Parameters
 ----------
 eccentric_anomaly : float
     Hyperbolic eccentric anomaly, if eccentricity is larger than 1, elliptical eccentric anomaly if it is smaller than 1
 eccentricity : float
     Value of the orbital eccentricity
 Returns
 -------
 float
     Value of the mean anomaly






     )doc" );

    m.def( "mean_to_eccentric_anomaly",
           &toec::convertMeanAnomalyToEccentricAnomaly< double >,
           py::arg( "eccentricity" ),
           py::arg( "mean_anomaly" ),
           py::arg( "use_default_initial_guess" ) = true,
           py::arg( "non_default_initial_guess" ) = TUDAT_NAN,
           py::arg( "root_finder" ) = nullptr,
           R"doc(

 Convert mean to eccentric anomaly.


 Parameters
 ----------
 eccentricity : float
     Value of the orbital eccentricity
 mean_anomaly : float
     Value of the mean anomaly
 use_default_initial_guess : bool, default = True
     Boolean to determine whether the user-defined initial guess is used for conversion, or an automatically generated one.
 non_default_initial_guess : float, default = NaN
     User-defined initial guess for conversion, to be used only if ``use_default_initial_guess`` is set to ``False``.
 root_finder : RootFinder, default = None
     User-defined root finder, overriding default root-finding algorithm for conversion (default is used if this input is left empty)
 Returns
 -------
 float
     Value of the eccentric anomaly






     )doc" );

    m.def( "elapsed_time_to_delta_mean_anomaly",
           &toec::convertElapsedTimeToMeanAnomalyChange< double >,
           py::arg( "elapsed_time" ),
           py::arg( "gravitational_parameter" ),
           py::arg( "semi_major_axis" ),
           R"doc(

 Convert elapsed time to the corresponding change in mean anomaly along a Keplerian orbit.


 Parameters
 ----------
 elapsed_time : float
     Elapsed time (in seconds)
 gravitational_parameter : float
     Gravitational parameter of central body
 semi_major_axis : float
     Semi-major axis of orbit
 Returns
 -------
 float
     Total change in mean anomaly along the Kepler orbit, accumulated in the provided time.






     )doc" );

    m.def( "delta_mean_anomaly_to_elapsed_time",
           &toec::convertMeanAnomalyChangeToElapsedTime< double >,
           py::arg( "mean_anomaly_change" ),
           py::arg( "gravitational_parameter" ),
           py::arg( "semi_major_axis" ),
           R"doc(

 Convert change in mean anomaly along a Keplerian orbit to the corresponding elapsed time.


 Parameters
 ----------
 mean_anomaly_change : float
     Total change in mean anomaly along the Kepler orbit
 gravitational_parameter : float
     Gravitational parameter of central body
 semi_major_axis : float
     Semi-major axis of orbit
 Returns
 -------
 float
     Time required for the provided mean anomaly change to be accumulated






     )doc" );

    m.def( "mean_motion_to_semi_major_axis",
           &toec::convertEllipticalMeanMotionToSemiMajorAxis< double >,
           py::arg( "mean_motion" ),
           py::arg( "gravitational_parameter" ),
           R"doc(

 Convert mean motion to corresponding semi-major axis (in a Keplerian orbit).


 Parameters
 ----------
 mean_motion : float
     Orbital mean motion
 gravitational_parameter : float
     Gravitational parameter of central body
 Returns
 -------
 float
     Semi-major axis corresponding to mean motion






     )doc" );

    m.def( "semi_major_axis_to_mean_motion",
           &toec::convertSemiMajorAxisToEllipticalMeanMotion< double >,
           py::arg( "semi_major_axis" ),
           py::arg( "gravitational_parameter" ),
           R"doc(

 Convert semi-major axis to corresponding mean motion (along a Keplerian orbit).


 Parameters
 ----------
 semi_major_axis : float
     Semi-major axis of orbit
 gravitational_parameter : float
     Gravitational parameter of central body
 Returns
 -------
 float
     Semi-major axis corresponding to mean motion






     )doc" );

    /*!
     **************   MODIFIED EQUIONOCTIAL ELEMENTS
     *******************
     */

    m.def( "keplerian_to_mee_manual_singularity",
           py::overload_cast< const Eigen::Vector6d&, const bool >(
                   &toec::convertKeplerianToModifiedEquinoctialElements< double > ),
           py::arg( "keplerian_elements" ),
           py::arg( "singularity_at_zero_inclination" ),
           R"doc(

 Convert Keplerian to Modified equinoctial elements.

 Convert Keplerian to Modified equinoctial elements (without intermediate step to Cartesian elements). The singularity-flipping
 element :math:`I` is to be provided manually for this function

 .. note:: See module level documentation for the standard ordering convention of Modified Equinoctial elements used.


 Parameters
 ----------
 keplerian_elements : numpy.ndarray
     Keplerian elements that are to be converted to Modified equinoctial elements
 singularity_at_zero_inclination : bool
     Singularity at 0 degrees inclination if ``True``, 180 degrees if ``False``
 Returns
 -------
 numpy.ndarray
     Modified equinoctial elements, as computed from Keplerian element input.






     )doc" );

    m.def( "keplerian_to_mee",
           py::overload_cast< const Eigen::Vector6d& >(
                   &toec::convertKeplerianToModifiedEquinoctialElements< double > ),
           py::arg( "keplerian_elements" ),
           R"doc(

 Convert Keplerian to Modified equinoctial elements.

 Convert Keplerian to Modified equinoctial elements (without intermediate step to Cartesian elements). The singularity-flipping
 element :math:`I` is computed automatically by this function (using :func:`flip_mee_singularity`)

 .. note:: See module level documentation for the standard ordering convention of Modified Equinoctial elements used.


 Parameters
 ----------
 keplerian_elements : numpy.ndarray
     Keplerian elements that are to be converted to Modified equinoctial elements
 Returns
 -------
 numpy.ndarray
     Modified equinoctial elements, as computed from Keplerian element input (with element :math:`I` defined by :func:`flip_mee_singularity`).






     )doc" );

    m.def( "flip_mee_singularity",
           py::overload_cast< const Eigen::Vector6d& >( &tmg::isOrbitRetrograde ),
           py::arg( "keplerian_elements" ),
           R"doc(

 Function to determine 'optimal' location of the singularity-flipping modified equinoctial element.

 Function to determine 'optimal' location of the singularity-flipping modified equinoctial element :math:`I`, if orbit inclination is less than
 90 degrees, it puts the singularity at 180 degrees, if it is larger than 90 degrees, it puts it at 0 degrees.


 Parameters
 ----------
 keplerian_elements : numpy.ndarray
     Keplerian elements that are to be converted to Modified equinoctial elements
 Returns
 -------
 bool
     Singularity at 0 degrees inclination if false, 180 degrees if true






     )doc" );

    m.def( "mee_to_keplerian",
           &toec::convertModifiedEquinoctialToKeplerianElements< double >,
           py::arg( "modified_equinoctial_elements" ),
           py::arg( "singularity_at_zero_inclination" ),
           R"doc(

 Convert Modified equinoctial to Keplerian elements.

 Modified equinoctial elements to Keplerian (without intermediate step to Cartesian elements).

 .. note:: See module level documentation for the standard ordering convention of Modified Equinoctial elements used.


 Parameters
 ----------
 modified_equinoctial_elements : numpy.ndarray
     Modified equinoctial elements that are to be converted to Keplerian elements
 singularity_at_zero_inclination : bool
     Singularity at 0 degrees inclination if false, 180 degrees if true
 Returns
 -------
 numpy.ndarray
     Keplerian elements, as computed from Modified equinoctial element input.






     )doc" );

    m.def( "cartesian_to_mee",
           py::overload_cast< const Eigen::Vector6d&, const double >(
                   &toec::convertCartesianToModifiedEquinoctialElements< double > ),
           py::arg( "cartesian_elements" ),
           py::arg( "gravitational_parameter" ),
           R"doc(

 Convert Cartesian to Modified equinoctial elements.

 Convert cartesian to Modified equinoctial elements. The singularity-flipping
 element :math:`I` is computed automatically by this function (using :func:`flip_mee_singularity`)

 .. note:: See module level documentation for the standard ordering convention of Modified Equinoctial elements used.


 Parameters
 ----------
 cartesian_elements : numpy.ndarray
     Cartesian elements that are to be converted to Modified equinoctial elements
 gravitational_parameter : float
     Gravitational parameter of central body
 Returns
 -------
 numpy.ndarray
     Modified equinoctial elements, as computed from Cartesian element input.






     )doc" );

    m.def( "cartesian_to_mee_manual_singularity",
           py::overload_cast< const Eigen::Vector6d&, const double, const bool >(
                   &toec::convertCartesianToModifiedEquinoctialElements< double > ),
           py::arg( "cartesian_elements" ),
           py::arg( "gravitational_parameter" ),
           py::arg( "singularity_at_zero_inclination" ),
           R"doc(

 Convert Cartesian to Modified equinoctial elements.

 Convert cartesian to Modified equinoctial elements. The singularity-flipping
 element :math:`I` is to be provided manually for this function

 .. note:: See module level documentation for the standard ordering convention of Modified Equinoctial elements used.


 Parameters
 ----------
 cartesian_elements : numpy.ndarray
     Cartesian elements that are to be converted to Modified equinoctial elements
 gravitational_parameter : float
     Gravitational parameter of central body
 singularity_at_zero_inclination : bool
     Singularity at 0 degrees inclination if false, 180 degrees if true
 Returns
 -------
 numpy.ndarray
     Modified equinoctial elements, as computed from Cartesian element input.






     )doc" );

    m.def( "mee_to_cartesian",
           py::overload_cast< const Eigen::Vector6d&, const double, const bool >(
                   &toec::convertModifiedEquinoctialToCartesianElements< double > ),
           py::arg( "modified_equinoctial_elements" ),
           py::arg( "gravitational_parameter" ),
           py::arg( "singularity_at_zero_inclination" ),
           R"doc(

 Convert Modified equinoctial to Cartesian elements.


 .. note:: See module level documentation for the standard ordering convention of Modified Equinoctial elements used.


 Parameters
 ----------
 modified_equinoctial_elements : numpy.ndarray
     Modified equinoctial elements that are to be converted to Cartesian elements
 gravitational_parameter : float
     Gravitational parameter of central body
 singularity_at_zero_inclination : bool
     Singularity at 0 degrees inclination if false, 180 degrees if true
 Returns
 -------
 numpy.ndarray
     Cartesian elements, as computed from Modified equinoctial element input.






     )doc" );

    /*!
     **************   SPHERICAL ELEMENTS  ******************
     */

    m.def( "spherical_to_cartesian_elementwise",
           py::overload_cast< double, double, double, double, double, double >(
                   &toec::convertSphericalOrbitalToCartesianState< double > ),
           py::arg( "radial_distance" ),
           py::arg( "latitude" ),
           py::arg( "longitude" ),
           py::arg( "speed" ),
           py::arg( "flight_path_angle" ),
           py::arg( "heading_angle" ),
           R"doc(

 Convert Spherical elements to Cartesian, with elementwise input.


 Parameters
 ----------
 radial_distance : float
     Distance from origin of central body
 latitude : float
     Central body-fixed latitude
 longitude : float
     Central body-fixed longitude
 speed : float
     Central body-fixed speed (norm of velocity vector). Note that this is *not* the norm of the inertial velocity
 flight_path_angle : float
     Flight-path angle (of central body-fixed velocity vector)
 heading_angle : float
     Heading angle (of central body-fixed velocity vector)
 Returns
 -------
 numpy.ndarray
     Cartesian elements, as computed from spherical element input.






     )doc" );

    m.def( "spherical_to_cartesian",
           py::overload_cast< const Eigen::Vector6d& >(
                   &toec::convertSphericalOrbitalToCartesianState< double > ),
           py::arg( "spherical_elements" ),
           R"doc(

 Convert spherical elements to Cartesian.

 .. note:: See module level documentation for the standard ordering convention of spherical state elements used.


 Parameters
 ----------
 spherical_elements : numpy.ndarray
     Spherical state that is to be converted to Cartesian elements
 Returns
 -------
 numpy.ndarray
     Cartesian elements, as computed from spherical element input.






     )doc" );

    m.def( "cartesian_to_spherical",
           &toec::convertCartesianToSphericalOrbitalState,
           py::arg( "cartesian_elements" ),
           R"doc(

 Convert Cartesian to spherical elements.

 .. note:: See module level documentation for the standard ordering  convention of spherical state elements used.


 Parameters
 ----------
 cartesian_elements : numpy.ndarray
     Cartesian state that is to be converted to spherical elements
 Returns
 -------
 numpy.ndarray
     Spherical elements, as computed from Cartesian element input.






     )doc" );

    /*!
     **************   QUATERNIONS  ******************
     */

    m.def( "quaternion_entries_to_rotation_matrix",
           &tla::convertVectorQuaternionToMatrixFormat,
           py::arg( "quaternion_entries" ),
           R"doc(

 Converts an array of four quaternion elements to the equivalent rotation matrix.

 Function to convert an array of four quaternion elements to the equivalent rotation matrix. These quaternion elements
 are for instance used when propagating rotational dynamics in Tudat, and this function can be used to convert the
 numerical results to a usable rotation matrix. See `our user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/frames_in_environment.html?highlight=rotational%20states#rotational-states>`_ for more details.


 Parameters
 ----------
 quaternion_entries : numpy.ndarray
     Quaternion elements, as per the convention used in the `Eigen library <https://eigen.tuxfamily.org/dox/classEigen_1_1Quaternion.html>`_
 Returns
 -------
 numpy.ndarray
     Rotation matrix defining the equivalent rotation.






     )doc" );

    m.def( "rotation_matrix_to_quaternion_entries",
           &tla::convertMatrixToVectorQuaternionFormat,
           py::arg( "rotation_matrix" ),
           R"doc(

 Converts a rotation matrix to the equivalent array of four quaternion elements.

 Inverse function of :func:`quaternion_entries_to_rotation_matrix`.


 Parameters
 ----------
 rotation_matrix : numpy.ndarray
     Rotation matrix
 Returns
 -------
 numpy.ndarray
     Equivalent quaternion elements, as per the convention used in the `Eigen library <https://eigen.tuxfamily.org/dox/classEigen_1_1Quaternion.html>`_






     )doc" );

    m.def( "quaternion_to_modified_rodrigues_parameters",
           &toec::convertQuaternionsToModifiedRodriguesParameterElements,
           py::arg( "quaternion_entries" ),
           R"doc(

 Converts quaternion elements to the equivalent modified Rodrigues parameters (for rotation representation).

 Parameters
 ----------
 quaternion_entries : numpy.ndarray
     Quaternion elements
 Returns
 -------
 numpy.ndarray
     Equivalent modified Rodrigues parameters






     )doc" );

    m.def( "modified_rodrigues_parameters_to_quaternion",
           &toec::convertModifiedRodriguesParametersToQuaternionElements,
           py::arg( "modified_rodrigues_parameters" ),
           R"doc(

 Converts modified Rodrigues parameters to the equivalent array of four quaternion elements (for rotation representation).

 Parameters
 ----------
 modified_rodrigues_parameters : numpy.ndarray
     Modified Rodrigues parameters
 Returns
 -------
 numpy.ndarray
     Equivalent quaternion elements


     )doc" );

    m.def( "quaternion_to_exponential_map",
           &toec::convertQuaternionsToExponentialMapElements,
           py::arg( "quaternion_entries" ),
           R"doc(

 Converts quaternion elements to the equivalent exponential map (for rotation representation).

 Parameters
 ----------
 quaternion_entries : numpy.ndarray
     Quaternion elements
 Returns
 -------
 numpy.ndarray
     Equivalent exponential map rotation elements






     )doc" );

    m.def( "exponential_map_to_quaternion",
           &toec::convertExponentialMapToQuaternionElements,
           py::arg( "exponential_map" ),
           R"doc(

 Converts modified Rodrigues parameters to the equivalent exponential map (for rotation representation).

 Parameters
 ----------
 exponential_map : numpy.ndarray
     Exponential map rotation elements
 Returns
 -------
 numpy.ndarray
     Equivalent quaternion elements


     )doc" );
    /*!
     **************   TLE  ******************

     */
    m.def( "teme_to_j2000",
           &te::getRotationMatrixFromTemeToJ2000,
           py::arg( "epoch" ),
           R"doc(

 Computes the rotation matrix from the TEME (True Equator Mean Equinox) frame to the J2000 frame, using the following:

 .. math::
     \mathbf{R}^{(\text{J2000}/\text{TEME})}=\mathbf{PN}(t)\mathbf{R}_{z}(-\theta(t)))

 where :math:`\theta` is the difference between the actual and mean position of the first point of Aries (or 'equation of the equinoxes`), computes using the ``iauEe00b`` function of the SOFA
 library (which computes this angle compatible with IAU 2000 resolutions but using the truncated nutation model IAU 2000B), and the precession-nutation matrix :math:`\mathbf{PN}` is computes using
 the function ``iauPnm80`` of the Sofa library, which uses th IAU 1976 precession model, and the IAU 1980 nutation model. The choice of slightly inconsistent IAU conventions is made for computational
 efficiency, as the 1976/1980 model includes fewer terms that the newer resolutions, and is combined with the truncated nutation model. Since the definition of the TEME frame is slightly ambiguous,
 the possible computational error that is incurred is insignificant for the application of the TEME frame: the transformation of SGP4-propagated TLEs to other epochs.

 Parameters
 ----------
 epoch : float
     Time as TDB seconds since J2000
 Returns
 -------
 numpy.ndarray
     Rotation matrix from TEME to J2000 frame


     )doc" );

    m.def( "j2000_to_teme",
           &te::getRotationMatrixFromTemeToJ2000,
           py::arg( "epoch" ),
           R"doc(

 Computes the rotation matrix from the J2000 to the TEME (True Equator Mean Equinox) frame, which is the inverse of the :func:`~teme_to_j2000` function.

 Parameters
 ----------
 epoch : float
     Time as TDB seconds since J2000
 Returns
 -------
 numpy.ndarray
     Rotation matrix from J2000 to TEME frame


     )doc" );

    /*!
 **************   STANDARD FRAMES  ******************

 */
    m.def( "j2000_to_eclipj2000",
           &tsi::getRotationFromJ2000ToEclipJ2000,
           R"doc(

 Provides the (constant) rotation matrix from the J2000 to the ECLIPJ2000 frame, as defined in the SPICE library (see :ref:`\`\`spice\`\`` for more details on our interface with this library).

 Returns
 -------
 numpy.ndarray
     Rotation matrix from J2000 to ECLIPJ2000 frame


     )doc" );

    m.def( "eclipj2000_to_j2000",
           &tsi::getRotationFromEclipJ2000ToJ2000,
           R"doc(

 Provides the (constant) rotation matrix from the ECLIPJ2000 to the J2000 frame, as defined in the SPICE library (see :ref:`\`\`spice\`\`` for more details on our interface with this library).

 Returns
 -------
 numpy.ndarray
     Rotation matrix from ECLIPJ2000 to J2000 frame


     )doc" );
}
}  // namespace element_conversion
}  // namespace astro
}  // namespace tudatpy
