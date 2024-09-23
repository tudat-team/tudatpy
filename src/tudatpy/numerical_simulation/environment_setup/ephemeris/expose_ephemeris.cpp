/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include <tudat/astro/reference_frames/referenceFrameTransformations.h>
#include <tudat/basics/deprecationWarnings.h>
#include <tudat/simulation/environment_setup.h>

#include "tudatpy/scalarTypes.h"

// #include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
// #include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace te = tudat::ephemerides;
namespace ti = tudat::interpolators;


namespace tudat {
    namespace simulation_setup {

        inline std::shared_ptr<EphemerisSettings>
        customEphemerisSettingsDeprecated(
            const std::function<Eigen::Vector6d(const double)>
                customStateFunction,
            const std::string& frameOrigin = "SSB",
            const std::string& frameOrientation = "ECLIPJ2000") {
            static bool isWarningPrinted = false;
            if(isWarningPrinted == false) {
                tudat::utilities::printDeprecationWarning(
                    "tudatpy.numerical_simulation.environment_setup.ephemeris."
                    "custom",
                    "tudatpy.numerical_simulation.environment_setup.ephemeris."
                    "custom_ephemeris");
                isWarningPrinted = true;
            }

            return customEphemerisSettings(customStateFunction, frameOrigin,
                                           frameOrientation);
        }
    }  // namespace simulation_setup
}  // namespace tudat

namespace tudatpy {
    namespace numerical_simulation {
        namespace environment_setup {
            namespace ephemeris {

                PYBIND11_MODULE(expose_ephemeris, m) {
                    py::module_::import("tudatpy.math.interpolators");

                    /////////////////////////////////////////////////////////////////////////////
                    // createEphemeris.h (complete, unverified)
                    /////////////////////////////////////////////////////////////////////////////
                    py::class_<tss::EphemerisSettings,
                               std::shared_ptr<tss::EphemerisSettings>>(
                        m, "EphemerisSettings",
                        R"doc(Base class for providing settings for ephemeris model.

	Functional (base) class for settings of ephemeris models that require no information in addition to their type (and frame origin and orientation).
	Ephemeris model classes requiring additional information must be created using an object derived from this class.

)doc")
                        //            .def(py::init<const tss::EphemerisType,
                        //                 const std::string &,
                        //                 const std::string &>(),
                        //                 py::arg("ephemeris_type"),
                        //                 py::arg("frame_origin") = "SSB",
                        //                 py::arg("frame_orientation") =
                        //                 "ECLIPJ2000")
                        .def_property(
                            "frame_origin",
                            &tss::EphemerisSettings::getFrameOrigin,
                            &tss::EphemerisSettings::resetFrameOrigin,
                            R"doc(Origin of frame in which ephemeris data is to be defined.
	)doc")
                        .def_property(
                            "frame_orientation",
                            &tss::EphemerisSettings::getFrameOrientation,
                            &tss::EphemerisSettings::resetFrameOrientation,
                            R"doc(Orientation of frame in which ephemeris data is to be defined.
	)doc")
                        .def_property(
                            "make_multi_arc_ephemeris",
                            &tss::EphemerisSettings::getMakeMultiArcEphemeris,
                            &tss::EphemerisSettings::resetMakeMultiArcEphemeris,
                            R"doc(Boolean denoting whether the ephemeris that is to be created is a multi-arc ephemeris.
	)doc")
                        .def_property_readonly(
                            "ephemeris_type",
                            &tss::EphemerisSettings::getEphemerisType,
                            R"doc(Type of ephemeris that is to be created.
	)doc");

                    py::class_<
                        tss::DirectSpiceEphemerisSettings,
                        std::shared_ptr<tss::DirectSpiceEphemerisSettings>,
                        tss::EphemerisSettings>(
                        m, "DirectSpiceEphemerisSettings",
                        R"doc(Class for defining settings of an ephemeris linked directly to Spice.

	`EphemerisSettings` derived class for ephemeris which are directly linked to Spice.
)doc")
                        //           .def(py::init<const std::string, const
                        //           std::string, const bool,
                        //                const bool, const bool, const
                        //                tss::EphemerisType>(),
                        //                py::arg("frame_origin") = "SSB",
                        //                py::arg("frame_orientation") =
                        //                "ECLIPJ2000",
                        //                py::arg("correct_for_stellar_aberration")
                        //                = false,
                        //                py::arg("correct_for_light_time_aberration")
                        //                = false,
                        //                py::arg("converge_light_time_aberration")
                        //                = false, py::arg("ephemeris_type") =
                        //                tss::direct_spice_ephemeris)
                        .def_property_readonly(
                            "correct_for_stellar_aberration",
                            &tss::DirectSpiceEphemerisSettings::
                                getCorrectForStellarAberration,
                            R"doc(Boolean defining whether to correct for stellar aberrations in retrieved values (of observed state).
	)doc")
                        .def_property_readonly(
                            "correct_for_light_time_aberration",
                            &tss::DirectSpiceEphemerisSettings::
                                getCorrectForLightTimeAberration,
                            R"doc(Boolean defining whether to correct for light time in retrieved values (of observed state).
	)doc")
                        .def_property_readonly(
                            "converge_light_time_aberration",
                            // TODO : Fix getConvergeLighTimeAberration typo in
                            // Tudat.
                            &tss::DirectSpiceEphemerisSettings::
                                getConvergeLighTimeAberration,
                            R"doc(Boolean defining whether to use single iteration or max. 3 iterations for calculating light time correction.
	)doc");


                    py::class_<tss::InterpolatedSpiceEphemerisSettings,
                               std::shared_ptr<
                                   tss::InterpolatedSpiceEphemerisSettings>,
                               tss::DirectSpiceEphemerisSettings>(
                        m, "InterpolatedSpiceEphemerisSettings",
                        R"doc(Class for defining settings of an ephemeris interpolated from Spice data.

	`DirectSpiceEphemerisSettings` derived class for setting ephemerides to be created from interpolated Spice ephemeris data.
)doc")
                        //            .def(py::init<
                        //                 double, double, double, std::string,
                        //                 std::string,
                        //                 std::shared_ptr<tudat::interpolators::InterpolatorSettings>>(),
                        //                 py::arg("initial_time"),
                        //                 py::arg("final_time"),
                        //                 py::arg("time_step"),
                        //                 py::arg("frame_origin") = "SSB",
                        //                 py::arg("frame_orientation") =
                        //                 "ECLIPJ2000",
                        //                 py::arg("interpolator_settings") =
                        //                 std::make_shared<
                        //            tudat::interpolators::LagrangeInterpolatorSettings>(6))
                        .def_property_readonly(
                            "initial_time",
                            &tss::InterpolatedSpiceEphemerisSettings::
                                getInitialTime,
                            R"doc(Initial time from which interpolated data from Spice should be created.
	)doc")
                        .def_property_readonly(
                            "final_time",
                            &tss::InterpolatedSpiceEphemerisSettings::
                                getFinalTime,
                            R"doc(Final time from which interpolated data from Spice should be created.
	)doc")
                        .def_property_readonly(
                            "time_step",
                            &tss::InterpolatedSpiceEphemerisSettings::
                                getTimeStep,
                            R"doc(Time step setting to be used for the state interpolation.
	)doc");


                    py::class_<
                        tss::ApproximateJplEphemerisSettings,
                        std::shared_ptr<tss::ApproximateJplEphemerisSettings>,
                        tss::EphemerisSettings>(
                        m, "ApproximateJplEphemerisSettings",
                        R"doc(Class for creating settings of approximate ephemeris for major planets.

	`EphemerisSettings` derived class for approximate ephemeris for major planets as implemented in ApproximateJplEphemerisSettings class and derived class (described in `this document <https://ssd.jpl.nasa.gov/planets/approx_pos.html>`_).
)doc")
                        .def_property_readonly(
                            "body_name",
                            &tss::ApproximateJplEphemerisSettings::getBodyName,
                            "");


                    py::class_<tss::ScaledEphemerisSettings,
                               std::shared_ptr<tss::ScaledEphemerisSettings>,
                               tss::EphemerisSettings>(
                        m, "ScaledEphemerisSettings",
                        R"doc(Class for defining settings from scaling existing ephemeris settings.

	`EphemerisSettings` derived class for a new ephemeris created from scaling an existing ephemeris settings object. It allows the user to apply a scaling factor to the resulting Cartesian states (for instance for an uncertainty analysis).
)doc");


                    py::class_<tss::ConstantEphemerisSettings,
                               std::shared_ptr<tss::ConstantEphemerisSettings>,
                               tss::EphemerisSettings>(
                        m, "ConstantEphemerisSettings",
                        R"doc(Class for defining settings of constant ephemerides.

	`EphemerisSettings` derived class for ephemerides producing a constant (time-independent) state.
)doc");
                    //            .def(py::init<const Eigen::Vector6d &,
                    //                 const std::string &,
                    //                 const std::string &>(),
                    //                 py::arg("constant_state"),
                    //                 py::arg("frame_origin") = "SSB",
                    //                 py::arg("frame_orientation") =
                    //                 "ECLIPJ2000");

                    py::class_<tss::CustomEphemerisSettings,
                               std::shared_ptr<tss::CustomEphemerisSettings>,
                               tss::EphemerisSettings>(
                        m, "CustomEphemerisSettings",
                        R"doc(Class for defining settings of a custom ephemeris.

	`EphemerisSettings` derived class for ephemerides which represent an ideal Kepler orbit.
)doc")
                        //           .def(py::init<const
                        //           std::function<Eigen::Vector6d(const
                        //           double)>,
                        //                const std::string &,
                        //                const std::string &>(),
                        //                py::arg("custom_state_function"),
                        //                py::arg("frame_origin") = "SSB",
                        //                py::arg("frame_orientation") =
                        //                "ECLIPJ2000")
                        .def_property_readonly("get_custom_state_function",
                                               &tss::CustomEphemerisSettings::
                                                   getCustomStateFunction,
                                               "");


                    py::class_<tss::KeplerEphemerisSettings,
                               std::shared_ptr<tss::KeplerEphemerisSettings>,
                               tss::EphemerisSettings>(
                        m, "KeplerEphemerisSettings", "")
                        //            .def(py::init<const Eigen::Vector6d &,
                        //            const double, const double,
                        //                 const std::string &, const
                        //                 std::string &, const double, const
                        //                 double>(),
                        //                 py::arg("initial_state_in_keplerian_elements"),
                        //                 py::arg("epoch_of_initial_state"),
                        //                 py::arg("central_body_gravitational_parameter"),
                        //                 py::arg("frame_origin") = "SSB",
                        //                 py::arg("frame_orientation") =
                        //                 "ECLIPJ2000",
                        //                 py::arg("root_finder_absolute_tolerance")
                        //                 =
                        //            200.0 *
                        //            std::numeric_limits<double>::epsilon(),
                        //                 py::arg("root_finder_maximum_number_of_iterations")
                        //                 = 1000.0)
                        .def_property_readonly(
                            "initial_state_in_keplerian_elements",
                            &tss::KeplerEphemerisSettings::
                                getInitialStateInKeplerianElements,
                            "")
                        .def_property_readonly("epoch_of_initial_state",
                                               &tss::KeplerEphemerisSettings::
                                                   getEpochOfInitialState,
                                               "")
                        .def_property_readonly(
                            "central_body_gravitational_parameter",
                            &tss::KeplerEphemerisSettings::
                                getCentralBodyGravitationalParameter,
                            "")
                        .def_property_readonly(
                            "root_finder_absolute_tolerance",
                            &tss::KeplerEphemerisSettings::
                                getRootFinderAbsoluteTolerance,
                            "")
                        .def_property_readonly(
                            "root_finder_maximum_number_of_iterations",
                            &tss::KeplerEphemerisSettings::
                                getRootFinderMaximumNumberOfIterations,
                            "");


                    py::class_<tss::TabulatedEphemerisSettings,
                               std::shared_ptr<tss::TabulatedEphemerisSettings>,
                               tss::EphemerisSettings>(
                        m, "TabulatedEphemerisSettings",
                        R"doc(Class for defining settings of ephemeris to be created from tabulated data.

	`EphemerisSettings` derived class for ephemeris created from tabulated data. The provided data is interpolated into ephemerides.
)doc")
                        //            .def(py::init<const std::map<double,
                        //            Eigen::Vector6d> &, std::string,
                        //                 std::string>())
                        .def_property_readonly(
                            "body_state_history",
                            &tss::TabulatedEphemerisSettings::
                                getBodyStateHistory,
                            R"doc(Dictionary of the discrete state history data from which ephemeris is to be created.
	)doc");


                    m.def("create_ephemeris",
                          &tss::createBodyEphemeris<double, TIME_TYPE>,
                          py::arg("ephemeris_settings"), py::arg("body_name"),
                          "");


                    m.def(
                        "keplerian", &tss::keplerEphemerisSettings,
                        py::arg("initial_keplerian_state"),
                        py::arg("initial_state_epoch"),
                        py::arg("central_body_gravitational_parameter"),
                        py::arg("frame_origin") = "SSB",
                        py::arg("frame_orientation") = "ECLIPJ2000",
                        py::arg("root_finder_absolute_tolerance") =
                            200.0 * std::numeric_limits<double>::epsilon(),
                        py::arg("root_finder_maximum_iterations") = 1000.0,
                        R"doc(Factory function for creating Keplerian ephemeris model settings.

	Factory function for settings object, defining ephemeris model which represents an ideal Kepler orbit from the given Kepler elements.
	These are taken as the elements at the ``initial_state_epoch`` and propagated to any other time using the provided ``central_body_gravitational_parameter``.
	See `Frame/State Transformations <https://tudat-space.readthedocs.io/en/latest/_src_user_guide/astrodynamics/transformations.html#frame-state-transformations>`_ and the :ref:`\`\`astro\`\`` module for more details on orbital elements in tudat.


	:param initial_state_in_keplerian_elements:
		Kepler elements at epoch given by ``initial_state_epoch``.

	:param initial_state_epoch:
		Epoch at which ``initial_state_epoch`` represents the Keplerian state.

	:param central_body_gravitational_parameter:
		Effective gravitational parameter of the central body that is used in the computations. Note that when
		the Keplerian orbit is to represent the relative state of two massive bodies, with one of these bodies as the origin
		this values should be the *sum* of the two bodies' gravitational parameters

	:param frame_origin:
		Origin of frame in which ephemeris data is defined.
	:param frame_orientation:
		Orientation of frame in which ephemeris data is defined.
	:param root_finder_absolute_tolerance:
		Convergence tolerance on iterative conversion from mean to eccentric anomaly;
		applies every time a cartesian state is requested from the kepler ephemeris, such as during propagation.

	:param root_finder_maximum_number_of_iterations:
		Maximum iteration on iterative conversion from mean to eccentric anomaly;
		applies every time a cartesian state is requested from the kepler ephemeris, such as during propagation.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.KeplerEphemerisSettings` class
)doc");

                    m.def(
                        "keplerian_from_spice",
                        &tss::keplerEphemerisFromSpiceSettings, py::arg("body"),
                        py::arg("initial_state_epoch"),
                        py::arg("central_body_gravitational_parameter"),
                        py::arg("frame_origin") = "SSB",
                        py::arg("frame_orientation") = "ECLIPJ2000",
                        py::arg("root_finder_absolute_tolerance") =
                            200.0 * std::numeric_limits<double>::epsilon(),
                        py::arg("root_finder_maximum_iterations") = 1000.0,
                        R"doc(Factory function for creating Keplerian ephemeris model settings with initial state from Spice.

	Factory function for settings object, defining ephemeris model which represents an ideal Kepler orbit from an initial state from Spice.
	The Kepler elements inferred from the initial state are propagated to any other time using the provided ``central_body_gravitational_parameter``.
	See `Frame/State Transformations <https://tudat-space.readthedocs.io/en/latest/_src_user_guide/astrodynamics/transformations.html#frame-state-transformations>`_ and the :ref:`\`\`astro\`\`` module for more details on orbital elements in tudat.


	:param body:
		Name of body for which to create ephemeris settings and infer initial state from Spice.
	:param initial_state_epoch:
		Epoch at which ``initial_state_epoch`` represents the Keplerian state.

	:param central_body_gravitational_parameter:
		Gravitational parameter of the central body that is used in the computations.
	:param frame_origin:
		Origin of frame in which ephemeris data is defined.
	:param frame_orientation:
		Orientation of frame in which ephemeris data is defined.
	:param root_finder_absolute_tolerance:
		Convergence tolerance on iterative conversion from mean to eccentric anomaly;
		applies every time a cartesian state is requested from the kepler ephemeris, such as during propagation.

	:param root_finder_maximum_number_of_iterations:
		Maximum iteration on iterative conversion from mean to eccentric anomaly;
		applies every time a cartesian state is requested from the kepler ephemeris, such as during propagation.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.KeplerEphemerisSettings` class
)doc");


                    m.def(
                        "approximate_jpl_model",
                        py::overload_cast<const std::string>(
                            &tss::approximateJplEphemerisSettings),
                        py::arg("body_name"),
                        R"doc(Factory function for creating approximate ephemeris model settings for major planets.

	Factory function for settings object, defining approximate ephemeris model for major planets.
	In this highly simplified ephemeris model, Keplerian elements of the major solar system bodies are modelled as linear functions of time and several sinusoidal variations (described `this document <https://ssd.jpl.nasa.gov/planets/approx_pos.html>`_).
	Note that this option is only available for solar system planets. For the case of the Earth the approximate ephemeris of the Earth-Moon barycenter is returned.


	:param body_name:
		String that is attempted to be matched to an identifier for the body that the ephemeris is to be created for.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.ApproximateJPLEphemerisSettings` class
)doc");

                    m.def(
                        "direct_spice",
                        py::overload_cast<const std::string, const std::string,
                                          const std::string>(
                            &tss::directSpiceEphemerisSettings),
                        py::arg("frame_origin") = "SSB",
                        py::arg("frame_orientation") = "ECLIPJ2000",
                        py::arg("body_name_to_use") = "",
                        R"doc(Factory function for creating ephemeris model settings entirely from Spice.

	Factory function for settings object, defining ephemeris model directly and entirely from Spice.
	Requires an appropriate Spice kernel to be loaded.


	:param frame_origin:
		Origin of frame in which ephemeris data is defined.
	:param frame_orientation:
		Orientation of frame in which ephemeris data is defined.
	:param body_name_to_use:
		Body from which Spice ephemeris is to be created.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.DirectSpiceEphemerisSettings` class
)doc");

                    m.def(
                        "interpolated_spice",
                        &tss::interpolatedSpiceEphemerisSettings,
                        py::arg("initial_time"), py::arg("final_time"),
                        py::arg("time_step"), py::arg("frame_origin") = "SSB",
                        py::arg("frame_orientation") = "ECLIPJ2000",
                        py::arg("interpolator_settings") =
                            std::make_shared<ti::LagrangeInterpolatorSettings>(
                                6),
                        py::arg("body_name_to_use") = "",
                        R"doc(Factory function for creating ephemeris model settings using interpolated Spice data.

	Factory function for settings object defining an ephemeris model from interpolated Spice data.
	Using this option the state of the body is retrieved from Spice at regular intervals `before` the environment propagation (as opposed to during the propagation).
	These data are then used to create an interpolator, which is put into the environment, and called during the propagation.
	This option has the downside of being applicable only during a limited time interval and requiring the tabulated data to be stored in RAM,
	but may for `some special cases <https://tudat-space.readthedocs.io/en/latest/_src_user_guide/environment_setup/valid_time_range.html>`_
	offer an advantage over a direct Spice ephemeris (:func:`~tudatpy.numerical_simulation.environment_setup.ephemeris.direct_spice`).


	:param initial_time:
		Initial time from which interpolated data from Spice should be created.
	:param final_time:
		Final time from which interpolated data from Spice should be created.
	:param time_step:
		Time step with which interpolated data from Spice should be created.
	:param frame_origin:
		Origin of frame in which ephemeris data is defined.
	:param frame_orientation:
		Orientation of frame in which ephemeris data is defined.
	:param interpolator_settings:
		Settings to be used for the state interpolation.
	:param body_name_to_use:
		Body from which Spice ephemeris is to be created.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.DirectSpiceEphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.InterpolatedSpiceEphemerisSettings` class
)doc");

                    m.def(
                        "tabulated",
                        py::overload_cast<
                            const std::map<double, Eigen::Vector6d>&,
                            std::string, std::string>(
                            &tss::tabulatedEphemerisSettings),
                        py::arg("body_state_history"),
                        py::arg("frame_origin") = "SSB",
                        py::arg("frame_orientation") = "ECLIPJ2000",
                        R"doc(Factory function for creating ephemeris model settings from tabulated data.

	Factory function for settings object, defining ephemeris model to be created from tabulated data.
	Currently the data that is provided gets interpolated by a 6th order Lagrange interpolator (hardcoded).
	At the edges of the interpolation interval a cubic spline interpolator is used to suppress the influence of Runge's phenomenon.


	:param body_state_history:
		Dictionary of the discrete state history data from which ephemeris is to be created. Keys representing the time (float) and values representing Cartesian states (numpy.ndarray).
	:param frame_origin:
		Origin of frame in which ephemeris data is defined.
	:param frame_orientation:
		Orientation of frame in which ephemeris data is defined.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.TabulatedEphemerisSettings` class
)doc");


                    m.def(
                        "tabulated_from_existing",
                        py::overload_cast<
                            const std::shared_ptr<tss::EphemerisSettings>,
                            const double, const double, const double,
                            const std::shared_ptr<ti::InterpolatorSettings>>(
                            &tss::tabulatedEphemerisSettings),
                        py::arg("ephemeris_settings"), py::arg("start_time"),
                        py::arg("end_time"), py::arg("time_step"),
                        py::arg("interpolator_settings") =
                            std::make_shared<ti::LagrangeInterpolatorSettings>(
                                8),
                        R"doc(Factory function for creating tabulated ephemeris model settings from existing ephemeris.

	Factory function for creating tabulated ephemeris model settings from existing ephemeris.
	The ephemeris that is provided gets tabulated in a given time frame, for a given time step.
	When called, this tabulated ephemeris will use interpolation, when needed, from the specified interpolator.

	.. note:: Creating tabulated ephemeris from existing ephemeris can for instance be used when combined with estimation.
	          This is because estimation needs the ephemeris to be tabulated to work.


	:param ephemeris_settings:
		Existing ephemeris settings that have to be tabulated.
	:param start_time:
		Initial time for which to create the tabulated ephemeris.
	:param end_time:
		Final time for which to create the tabulated ephemeris.
	:param time_step:
		Time step to use to tabulate the existing ephemeris.
	:param interpolator_settings:
		Interpolator settings to use when interpolating between two tabulated ephemeris.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.TabulatedEphemerisSettings` class
)doc");

                    m.def(
                        "constant", &tss::constantEphemerisSettings,
                        py::arg("constant_state"),
                        py::arg("frame_origin") = "SSB",
                        py::arg("frame_orientation") = "ECLIPJ2000",
                        R"doc(Factory function for creating constant ephemeris model settings.

	Factory function for settings object, defining ephemeris model with a constant, time-independent state.


	:param constant_state:
		Constant state that will be provided as output of the ephemeris at all times.
	:param frame_origin:
		Origin of frame in which ephemeris data is defined.
	:param frame_orientation:
		Orientation of frame in which ephemeris data is defined.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.constantEphemerisSettings` class
)doc");

                    m.def(
                        "scaled_by_constant",
                        py::overload_cast<
                            const std::shared_ptr<tss::EphemerisSettings>,
                            const double, const bool>(
                            &tss::scaledEphemerisSettings),
                        py::arg("unscaled_ephemeris_settings"),
                        py::arg("scaling_constant"),
                        py::arg("is_scaling_absolute") = false,
                        R"doc(Factory function for creating scaled ephemeris model settings.

	Factory function for settings object, defining ephemeris model based on an scaling of an existing ephemeris settings object.
	The user can apply a scaling factor (or an absolute value) to the resulting Cartesian states (for instance for an uncertainty analysis).


	:param unscaled_ephemeris_settings:
		Sets base settings of ephemeris to be scaled.
	:param scaling_constant:
		Constant scaling factor to be applied to all elements of the Cartesian state.
	:param is_scaling_absolute:
		Boolean indicating whether ephemeris scaling is absolute. Setting this boolean to true will add the scaling value to the state, instead of the default behaviour of multiplying the state by the scaling value.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.ScaledEphemerisSettings` class
)doc");

                    m.def(
                        "scaled_by_vector",
                        py::overload_cast<
                            const std::shared_ptr<tss::EphemerisSettings>,
                            const Eigen::Vector6d, const bool>(
                            &tss::scaledEphemerisSettings),
                        py::arg("unscaled_ephemeris_settings"),
                        py::arg("scaling_vector"),
                        py::arg("is_scaling_absolute") = false,
                        R"doc(Factory function for creating scaled ephemeris model settings.

	Factory function for settings object, defining ephemeris model based on an scaling of an existing ephemeris settings object.
	The user can apply a scaling factor (or an absolute value) to the resulting Cartesian states (for instance for an uncertainty analysis).


	:param unscaled_ephemeris_settings:
		Sets base settings of ephemeris to be scaled.
	:param scaling_vector:
		Vector containing scaling factors to be applied to each element of the Cartesian state.
	:param is_scaling_absolute:
		Boolean indicating whether ephemeris scaling is absolute. Setting this boolean to true will add the scaling value to the state, instead of the default behaviour of multiplying the state by the scaling value.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.ScaledEphemerisSettings` class
)doc");

                    m.def(
                        "scaled_by_vector_function",
                        py::overload_cast<
                            const std::shared_ptr<tss::EphemerisSettings>,
                            const std::function<Eigen::Vector6d(const double)>,
                            const bool>(&tss::scaledEphemerisSettings),
                        py::arg("unscaled_ephemeris_settings"),
                        py::arg("scaling_vector_function"),
                        py::arg("is_scaling_absolute") = false,
                        R"doc(Factory function for creating scaled ephemeris model settings.

	Factory function for settings object, defining ephemeris model based on an scaling of an existing ephemeris settings object.
	The user can apply a scaling factor (or an absolute value) to the resulting Cartesian states (for instance for an uncertainty analysis).


	:param unscaled_ephemeris_settings:
		Sets base settings of ephemeris to be scaled.
	:param scaling_vector_function:
		Function returning a vector with the scaling factors to be applied to each element of the Cartesian state.
	:param is_scaling_absolute:
		Boolean indicating whether ephemeris scaling is absolute. Setting this boolean to true will add the scaling value to the state, instead of the default behaviour of multiplying the state by the scaling value.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.ScaledEphemerisSettings` class
)doc");

                    m.def(
                        "custom_ephemeris", &tss::customEphemerisSettings,
                        py::arg("custom_state_function"),
                        py::arg("frame_origin") = "SSB",
                        py::arg("frame_orientation") = "ECLIPJ2000",
                        R"doc(Factory function for creating custom ephemeris model settings.

	Factory function for settings object, defining ephemeris model with a custom state.
	This allows the user to provide an custom state function as ephemeris model.


	:param custom_state_function:
		Function returning the state as a function of time.
	:param frame_origin:
		Origin of frame in which ephemeris data is defined.
	:param frame_orientation:
		Orientation of frame in which ephemeris data is defined.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.EphemerisSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.ephemeris.CustomEphemerisSettings` class
)doc");

                    m.def("custom", &tss::customEphemerisSettingsDeprecated,
                          py::arg("custom_state_function"),
                          py::arg("frame_origin") = "SSB",
                          py::arg("frame_orientation") = "ECLIPJ2000");
                }


            }  // namespace ephemeris
        }  // namespace environment_setup
    }  // namespace numerical_simulation
}  // namespace tudatpy
