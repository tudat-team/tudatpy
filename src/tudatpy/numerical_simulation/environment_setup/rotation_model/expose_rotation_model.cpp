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
#include <tudat/simulation/environment_setup.h>


// #include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
// #include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tba = tudat::basic_astrodynamics;

namespace tudatpy {
    namespace numerical_simulation {
        namespace environment_setup {
            namespace rotation_model {

                PYBIND11_MODULE(expose_rotation_model, m) {
                    /////////////////////////////////////////////////////////////////////////////
                    // createRotationalModel.h
                    /////////////////////////////////////////////////////////////////////////////
                    py::enum_<tss::RotationModelType>(
                        m, "RotationModelType",
R"doc(Enumeration of rotation model types.

	Enumeration of rotation model types supported by tudat.


	:member simple_rotation_model:
	:member spice_rotation_model:
	:member gcrs_to_itrs_rotation_model:
	:member synchronous_rotation_model:
	:member planetary_rotation_model:
)doc")
                        .value("simple_rotational_model",
                               tss::RotationModelType::simple_rotation_model,
"")
                        .value("spice_rotation_model",
                               tss::RotationModelType::spice_rotation_model,
"")
                        .value(
                            "gcrs_to_itrs_rotation_model",
                            tss::RotationModelType::gcrs_to_itrs_rotation_model,
"")
                        .value(
                            "synchronous_rotation_model",
                            tss::RotationModelType::synchronous_rotation_model,
"")
                        .value("planetary_rotation_model",
                               tss::RotationModelType::planetary_rotation_model,
"")
                        .export_values();

                    py::enum_<tba::IAUConventions>(
                        m, "IAUConventions",
R"doc(Enumeration of IAU conventions for Earth rotation.

	Enumeration of IAU conventions for Earth rotation supported by tudat.


	:member iau_2000_a:
	:member iau_2000_b:
	:member iau_2006:
)doc")
                        .value(
                            "iau_2000_a", tba::IAUConventions::iau_2000_a,
"")
                        .value(
                            "iau_2000_b", tba::IAUConventions::iau_2000_b,
"")
                        .value("iau_2006", tba::IAUConventions::iau_2006,
"")
                        .export_values();

                    py::class_<tss::RotationModelSettings,
                               std::shared_ptr<tss::RotationModelSettings>>(
                        m, "RotationModelSettings",
R"doc(Base class for providing settings for automatic rotation model creation.

	This class is a functional base class for settings of rotation models that require no information in addition to their type.
	Basic rotation model has constant orientation of the rotation axis (body-fixed z-axis) and constant rotation rate about this axis.
	Rotation models requiring additional information must be created using the factory functions which create the specific object derived from this base class.

)doc")
                        //            .def(py::init<const
                        //            tss::RotationModelType, const std::string
                        //            &,
                        //                 const std::string &>(),
                        //                 py::arg("rotation_type"),
                        //                 py::arg("base_frame"),
                        //                 py::arg("target_frame"))
                        .def_property_readonly(
                            "rotation_type",
                            &tss::RotationModelSettings::getRotationType,
R"doc(Type of rotation model that is to be created.
	)doc")
                        .def_property(
                            "base_frame",
                            &tss::RotationModelSettings::getOriginalFrame,
                            &tss::RotationModelSettings::resetOriginalFrame,
R"doc(Name of the base frame of rotation model.
	)doc")
                        .def_property_readonly(
                            "target_frame",
                            &tss::RotationModelSettings::getTargetFrame,
R"doc(Name of the target frame of rotation model.
	)doc");

                    py::class_<
                        tss::SimpleRotationModelSettings,
                        std::shared_ptr<tss::SimpleRotationModelSettings>,
                        tss::RotationModelSettings>(
                        m, "SimpleRotationModelSettings",
"");

                    py::class_<
                        tss::PlanetaryRotationModelSettings,
                        std::shared_ptr<tss::PlanetaryRotationModelSettings>,
                        tss::RotationModelSettings>(
                        m, "PlanetaryRotationModelSettings",
"");


                    m.def("simple",
                          py::overload_cast<const std::string &,
                                            const std::string &,
                                            const Eigen::Matrix3d &,
                                            const double, const double>(
                              &tss::simpleRotationModelSettings),
                          py::arg("base_frame"), py::arg("target_frame"),
                          py::arg("initial_orientation"),
                          py::arg("initial_time"), py::arg("rotation_rate"),
R"doc(Factory function for creating simple rotation model settings.

	Factory function for settings object, defining a basic rotation model with constant orientation of the rotation axis and constant rotation rate about this axis.
	Rotation from original (inertial) to target (body-fixed) frame at some reference time ``initial_time`` (:math:`t_{0}`) is defined by the ``initial_orientation`` (:math:`\mathbf{R}^{(B/I)}(t_{0})`) rotation matrix.
	Rotation about the body-fixed z-axis is defined by the ``rotation_rate`` (:math:`\omega`) float variable (in rad/s). The rotation matrix is computed from:

	.. math::
	   \mathbf{R}^{(B/I)}(t)=\mathbf{R}_{z}(\omega(t-t_{0}))(t_{0})\mathbf{R}^{(B/I)}(t_{0})

	where :math:`\mathbf{R}^{(B/I)}` denotes the rotation matrix from inertial to body-fixed frame, and :math:`\mathbf{R}_{z}` denotes a rotaion matrix about the z-axis.

	The matrix :math:`\mathbf{R}^{(B/I)}(t_{0})` is sometimes parameterized by pole right ascension and declination (:math:`\alpha` and :math:`\delta`), as well as the meridian of date :math:`W_{0}` with

	.. math::
	   \mathbf{R}^{(B/I)}(t_{0})=\mathbf{R}_{z}(W_{0})\mathbf{R}_{x}(\pi/2-\delta)\mathbf{R}_{z}(\pi/2+\alpha)


	:param base_frame:
		Name of the base frame of rotation model.
	:param target_frame:
		Name of the target frame of rotation model.
	:param initial_orientation:
		Orientation of target frame in base frame at initial time.
	:param initial_time:
		Initial time (reference epoch for rotation matrices).
	:param rotation_rate:
		Constant rotation rate [rad/s] about rotational axis.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.SimpleRotationModelSettings` class
)doc");

                    m.def("simple_from_spice",
                          &tss::simpleRotationModelFromSpiceSettings,
                          py::arg("base_frame"), py::arg("target_frame"),
                          py::arg("target_frame_spice"),
                          py::arg("initial_time"),
R"doc(Factory function for creating simple rotation model settings using initial orientation and rotation rates from Spice.

	Factory function for settings object, defining a :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.simple` rotation model with the added functionality that the initial orientation and rotation rate are extracted from Spice, as opposed to provided manually.
	Note that `only` the initial orientation and rotation rate ( at the time defined by ``initial_time`` ) are extracted from Spice - for
	the full Spice rotation model see :func:`~tudatpy.numerical_simulation.environment_setup.rotation_model.spice`.
	Also note the distinction between the ``target_frame`` and ``target_frame_spice`` parameters.


	:param base_frame:
		Name of the base frame of rotation model.
	:param target_frame:
		Target frame of rotation model - name of frame that Tudat assigns to the body-fixed frame
	:param target_frame_spice:
		Spice reference of target frame - name of the frame in Spice for which the initial orientation and rotation rate are extracted.
	:param initial_time:
		Initial time (reference epoch for rotation matrices).
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.SimpleRotationModelSettings` class
)doc");

                    m.def("synchronous", &tss::synchronousRotationModelSettings,
                          py::arg("central_body_name"), py::arg("base_frame"),
                          py::arg("target_frame"),
R"doc(Factory function for creating synchronous rotational ephemeris settings.

	Factory function for settings object, defining a synchronous rotation model where rotation of a body is defined from its relative orbit w.r.t. some central body. Specifically
	- the body-fixed x-axis is *always* pointing towards the central body
	- the body-fixed z-axis is *always* perpendicular to the orbital plane (along the direction of :math:`\mathbf{x}\times\mathbf{v}` )
	- the body-fixed y-axis completes the right-handed reference frame

	Such a model can be useful for, for instance, approximate rotation of tidally locked natural satellites or nadir-pointing spacecraft.


	:param central_body_name:
		Name of the base frame of rotation model.
	:param base_frame:
		Name of the base frame of rotation model.
	:param target_frame:
		Spice reference of target frame.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.SynchronousRotationModelSettings` class
)doc");

                    m.def("spice", &tss::spiceRotationModelSettings,
                          py::arg("base_frame"), py::arg("target_frame"),
                          py::arg("spice_frame_name") = "",
R"doc(Factory function for creating rotation model settings from the Spice interface.

	Factory function for settings object, defining a rotation model directly (and entirely) from Spice interface.


	:param base_frame:
		Name of the base frame of rotation model.
	:param target_frame:
		Name of the target frame of rotation model.
	:return:
		Instance of :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` class.
)doc");

                    m.def("gcrs_to_itrs", &tss::gcrsToItrsRotationModelSettings,
                          py::arg("precession_nutation_theory") = tba::iau_2006,
                          py::arg("base_frame") = "GCRS",
                          py::arg("cio_interpolation_settings") = nullptr,
                          py::arg("tdb_to_tt_interpolation_settings") = nullptr,
                          py::arg("short_term_eop_interpolation_settings") =
                              nullptr,
R"doc(Factory function for creating high-accuracy Earth rotation model settings.

	Factory function for settings object, defining high-accuracy Earth rotation model according to the IERS 2010 Conventions.
	This settings class has various options to deviate from the default settings, typical applications will use default.
	Note that for this model the original frame must be J2000 or GCRS (in the case of the former, the frame bias between GCRS and J2000 is automatically corrected for). The target frame (e.g. body-fixed frame) name is ITRS.
	The precession-nutation theory may be any member of :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.IAUConventions` (``iau_2000a`` / ``iau_2000b`` or ``iau_2006``).
	Alternative options to modify the input (not shown here) include the EOP correction file, input time scale, short period UT1 and polar motion variations.
	The target frame (e.g. body-fixed frame) name is ITRS.


	:param precession_nutation_theory:
		Setting theory for modelling Earth nutation.
	:param base_frame:
		Base frame of rotation model
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.GcrsToItrsRotationModelSettings` class
)doc");

                    m.def("aerodynamic_angle_based",
                          &tss::aerodynamicAngleRotationSettings,
                          py::arg("central_body"), py::arg("base_frame"),
                          py::arg("target_frame"),
                          py::arg("angle_funcion") = nullptr,
R"doc(Factory function for creating rotation model settings based on custom aerodynamic angles (attack, sideslip, bank).

	Factory function for creating rotation model settings based on custom aerodynamic angles:
	angle of attack :math:`\alpha`, sideslip angle :math:`\beta` and bank angle :math:`\sigma`. The use of this function is typical for
	simulating the dynamics of a (guided) re-entry vehicle. It calculates the rotation matrix from inertial frame to the body-fixed frame
	of the current body B (typically a vehicle) w.r.t. the body-fixed frame of a central body C (e.g., the body at which the re-entry is taking place.
	The full algorithm for :math:`R^{(I/B)}` is described by Mooij (1994), and is composed of:

	*  The rotation from inertial frame to the body fixed frame of body C, using the existing rotation model of body C
	*  The rotation from body-fixed frame of body C to the vehicle's vertical frame V. This rotation uses the current latitude and longitude angles.
	*  The rotation of the vehicle's vertical frame V to its trajectory frame T. This rotation uses the current heading and flight path angles.
	*  The rotation of the vehicle's trajectory frame T to its aerodynamic frame A. This rotation uses the current bank angle
	*  The rotation of the vehicle's aerodynamic frame A to its body-fixed frame. This rotation uses the current angle of attack and sideslip angles

	In the above algorithm, the latitude, longitude, heading and flight-path angles are computed from the vehicle's current translational state, in the body-fixed
	frame of body C. The angle of attack, sideslip angle and bank angle are to be defined by the user, through a single custom function that is passed to
	the ``angle_function`` argument of this functions


	:param central_body:
		Name of the central body C that is to be used.
	:param base_frame:
		Name of the base frame of rotation model.
	:param target_frame:
		Name of the target frame of rotation model.
	:param angle_function:
		Custom function provided by the user, which returns an array of three values as a function of time. The output of this function *must* be ordered as :math:`[\alpha,\beta,\sigma]`. If this input is left empty, these angles are both fixed to 0.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.CustomRotationModelSettings` class, which defines the required settings for the rotation model.
)doc");

                    m.def("zero_pitch_moment_aerodynamic_angle_based",
                          &tss::pitchTrimRotationSettings,
                          py::arg("central_body"), py::arg("base_frame"),
                          py::arg("target_frame"),
                          py::arg("angle_funcion") = nullptr,
R"doc(Factory function for creating rotation model settings based on an angle of attack calculated from pitch-trim, and custom aerodynamic angles sideslip, bank.

	Factory function for creating rotation model settings based on an angle of attack calculated from pitch-trim, and custom aerodynamic angles sideslip, bank. This function is
	largely identical to the :func:`~aerodynamic_angle_based`, with the difference that the angle of attack :math:`\alpha` is not provided as a custom value by the user, but is
	calculated from the body's aerodynamic moment coefficients, such that we have :math:`C_{m}=0`. This requires aerodynamic moment coefficients to be defined for the vehicle that
	depend on (among others) the body's angle of attack


	:param central_body:
		Name of the central body C that is to be used.
	:param base_frame:
		Name of the base frame of rotation model.
	:param target_frame:
		Name of the target frame of rotation model.
	:param angle_funcion:
		Custom function provided by the user, which returns an array of three values as a function of time. The output of this function *must* be ordered as :math:`[\beta,\sigma]`. If this input is left empty, these angles are both fixed to 0.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.CustomRotationModelSettings` class, which defines the required settings for the rotation model.
)doc");

                    m.def("custom_inertial_direction_based",
                          &tss::bodyFixedDirectionBasedRotationSettings,
                          py::arg("inertial_body_axis_direction"),
                          py::arg("base_frame"), py::arg("target_frame"),
                          py::arg("free_rotation_angle_function") = nullptr,
R"doc(Factory function for creating rotation model settings where the body-fixed x-axis is imposed to lie in a user-defined inertial direction

	Factory function for creating rotation model settings where the body-fixed x-axis is imposed to lie in a user-defined inertial direction :math:`\hat{\mathbf{T}}_{I}`. Specifically, it ensures
	that the rotation matrix from body-fixed to inertial frame is set up such that :math:`\hat{\mathbf{T}}_{I}=R^{(I/B)}\hat{\mathbf{i}}` (where :math:`\mathbf{i}` is the unit-vector in local x-direction).
	The complete rotation matrix requires an additional angle :math:`\phi` (rotation of the body about its body-fixed x-axis), which is set to 0 by default.

	The full rotation matrix is computed from a 3-2-1 Euler angle rotation
	:math:`R^{(I/B)}=R_{z}(\psi)R_{y}(\theta)R_{x}(\phi)`, with :math:`\psi` and :math:`\theta` computed from the suitable decomposition of :math:`\hat{\mathbf{T}}_{I}`.
	This function is typically used for simulating the (guided) dynamics of a spacecraft under thrust, where the thrust is provided in the x-direction of the body-fixed frame. By providing a suitable
	``inertial_body_axis_direction``, this thrust can be defined to point in an arbitrary direction (typically defined by a guidance algorithm) in the inertial frame as a function of time.

	NOTE: this function may be extended in the future to allow an arbitrary body-fixed direction to align with an arbitrary inertial direction. At present, its functionality is limited to imposing the inertial direction of the body-fixed x-axis.


	:param inertial_body_axis_direction:
		Custom function defined by the user, which imposes the inertial orientation of the body-fixed x-axis, by providing :math:`\hat{\mathbf{T}}_{I}(t)`.
	:param base_frame:
		Name of the base frame of rotation model.
	:param target_frame:
		Name of the target frame of rotation model.
	:param free_rotation_angle_function:
		Custom function provided by the user, which returns a value for the free rotation angle :math:`\phi` about the body-fixed x-axis as a function of time. If this input is left empty, this angle is fixed to 0.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.BodyFixedDirectionBasedRotationSettings` class, which defines the required settings for the rotation model.
)doc");

                    m.def(
                        "orbital_state_direction_based",
                        &tss::orbitalStateBasedRotationSettings,
                        py::arg("central_body"),
                        py::arg("is_colinear_with_velocity"),
                        py::arg("direction_is_opposite_to_vector"),
                        py::arg("base_frame"), py::arg("target_frame") = "",
                        py::arg("free_rotation_angle_function") = nullptr,
R"doc(Factory function for creating rotation model settings where the body-fixed x-axis is imposed to lie in the direction of a relative position or velocity vector.

	Factory function for creating rotation model settings where the body-fixed x-axis is imposed to lie in the direction of a relative position or velocity vector. This function is
	similar to the :func:`~custom_inertial_direction_based` function, with the exception that the :math:`\hat{\mathbf{T}}_{I}` vector is not defined by thee user, but is defined by the
	relative position vector :math:`\mathbf{r}_{C}` or velocity vector :math:`\mathbf{r}_{C}` of the vehicle w.r.t. some body C. The inputs to this function allow :math:`\hat{\mathbf{T}}_{I}` to
	be set to :math:`\pm\mathbf{r}_{C}` or :math:`\pm\mathbf{v}_{C}`, for any body C. It is typically used for simplified or preliminary thrust analyses.


	:param central_body:
		Name of central body w.r.t. which the position/velocity vector is to be computed
	:param is_colinear_with_velocity:
		Boolean defining whether :math:`\hat{\mathbf{T}}_{I}` is to be aligned with velocity (if true) or position (if false)
	:param direction_is_opposite_to_vector:
		Boolean defining whether :math:`\hat{\mathbf{T}}_{I}` is to be in the same direction as position/velocity (if false), or in the opposite direction (if true).
	:param base_frame:
		Name of the base frame of rotation model.
	:param target_frame:
		Name of the target frame of rotation model.
	:param free_rotation_angle_function:
		Custom function provided by the user, which returns a value for the free rotation angle :math:`\phi` about the body-fixed x-axis as a function of time. If this input is left empty, this angle is fixed to 0.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.BodyFixedDirectionBasedRotationSettings` class, which defines the required settings for the rotation model.
)doc");


                    m.def("constant_rotation_model",
                          py::overload_cast<const std::string &,
                                            const std::string &,
                                            const Eigen::Matrix3d &>(
                              &tss::constantRotationModelSettings),
                          py::arg("base_frame"), py::arg("target_frame"),
                          py::arg("initial_orientation"),
R"doc(Factory function for creating simple rotation model settings for target-frames with constant orientation.

	Factory function for settings object, defining simple rotation model setting objects with constant rotation matrix.
	These model settings are for target frames which do not have a rotational rate in the base frame and are fully defined by their initial orientation.


	:param base_frame:
		Name of the base frame of rotation model.
	:param target_frame:
		Name of the target frame of rotation model.
	:param initial_orientation:
		Rotation matrix from inertial to body-fixed (base to target) frame at initial time (constant throughout).
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.SimpleRotationModelSettings` class.
)doc");

                    m.def("custom_rotation_model",
                          &tss::customRotationModelSettings,
                          py::arg("base_frame"), py::arg("target_frame"),
                          py::arg("custom_rotation_matrix_function"),
                          py::arg("finite_difference_time_step"),
"");


                    m.def("mars_high_accuracy",
                          &tss::getHighAccuracyMarsRotationModel,
                          py::arg("base_frame") = "ECLIPJ2000",
                          py::arg("target_frame") = "Mars_Fixed",
R"doc(Factory function for creating a high-accuracy Mars rotation model.

	Factory function for creating a high-accuracy Mars rotation model, using the default parameters of `Konopliv et al. (2016) <https://www.sciencedirect.com/science/article/abs/pii/S0019103516001305>`_
	and the mathematical model of ` Konopliv et al. (2006) <https://www.sciencedirect.com/science/article/pii/S0019103506000297>`_. The rotation matrix formulation is given in Eq. (13)-(19) of that paper.
	Note that, at the moment, all parameters in this rotation model are hard coded, and cannot be adapted by the user (except by estimating a number of its constituent parameters, see :ref:`\`\`parameter\`\`` module )
	As such, this model is at present applicable to Mars rotation only. If you require more fine-grained control of the parameters, please contact the Tudat support team

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.RotationModelSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.rotation_model.PlanetaryRotationModelSettings` class, which defines the required settings for the rotation model.
)doc");
                }

            }  // namespace rotation_model
        }  // namespace environment_setup
    }  // namespace numerical_simulation
}  // namespace tudatpy
