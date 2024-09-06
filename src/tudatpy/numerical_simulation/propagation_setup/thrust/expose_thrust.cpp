/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

// #include "kernel/expose_numerical_simulation/deprecation_support.h"

#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/simulation/propagation_setup.h>


namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tinterp = tudat::interpolators;
namespace te = tudat::ephemerides;
namespace tni = tudat::numerical_integrators;
namespace trf = tudat::reference_frames;
namespace tmrf = tudat::root_finders;

namespace tudat {
    namespace simulation_setup {
        // Deprecated (We still rely on it tho)
        enum ThrustFrames {
            unspecified_thrust_frame = -1,
            inertial_thrust_frame = 0,
            tnw_thrust_frame = 1
        };
    }  // namespace simulation_setup
}  // namespace tudat

namespace tudatpy {
    namespace numerical_simulation {
        namespace propagation_setup {
            namespace thrust {

                PYBIND11_MODULE(expose_thrust, m) {
                    py::module_::import("tudatpy.math.interpolators");
                    py::module_::import("tudatpy.math.root_finders");
                    py::module_::import(
                        "tudatpy.numerical_simulation.environment");
                    py::enum_<tss::ThrustMagnitudeTypes>(m,
                                                         "ThrustMagnitudeTypes")
                        .value("constant_thrust_magnitude",
                               tss::ThrustMagnitudeTypes::
                                   constant_thrust_magnitude)
                        //            .value("from_engine_properties_thrust_magnitude",
                        //            tss::ThrustMagnitudeTypes::from_engine_properties_thrust_magnitude)
                        .value("thrust_magnitude_from_time_function",
                               tss::ThrustMagnitudeTypes::
                                   thrust_magnitude_from_time_function)
                        .value("thrust_magnitude_from_dependent_variables",
                               tss::ThrustMagnitudeTypes::
                                   thrust_magnitude_from_dependent_variables);
                    //            .value("bang_bang_thrust_magnitude_from_mee_costates",
                    //            tss::ThrustMagnitudeTypes::bang_bang_thrust_magnitude_from_mee_costates);

                    py::class_<tss::ThrustMagnitudeSettings,
                               std::shared_ptr<tss::ThrustMagnitudeSettings>>(
                        m, "ThrustMagnitudeSettings",
R"doc(Functional base class to define settings for the thrust magnitude.

)doc")
                        .def_readonly(
                            "thrust_magnitude_type",
                            &tss::ThrustMagnitudeSettings::thrustMagnitudeType_)
                        .def_readonly(
                            "thrust_origin_id",
                            &tss::ThrustMagnitudeSettings::thrustOriginId_);

                    py::class_<
                        tss::ConstantThrustMagnitudeSettings,
                        std::shared_ptr<tss::ConstantThrustMagnitudeSettings>,
                        tss::ThrustMagnitudeSettings>(
                        m, "ConstantThrustMagnitudeSettings",
R"doc(`ThrustMagnitudeSettings`-derived class to define settings for constant thrust magnitude.

	Derived class to provide settings for the thrust magnitude. This class should be used to define a constant thrust
	magnitude.

)doc")
                        .def_readonly("thrust_magnitude",
                                      &tss::ConstantThrustMagnitudeSettings::
                                          thrustMagnitude_)
                        .def_readonly("specific_impulse",
                                      &tss::ConstantThrustMagnitudeSettings::
                                          specificImpulse_);

                    py::class_<
                        tss::CustomThrustMagnitudeSettings,
                        std::shared_ptr<tss::CustomThrustMagnitudeSettings>,
                        tss::ThrustMagnitudeSettings>(
                        m, "CustomThrustMagnitudeSettings",
"");

                    m.def("get_propulsion_input_variables",
                          &tss::getPropulsionInputVariables,
                          py::arg("body_with_guidance") =
                              std::shared_ptr<tss::Body>(),
                          py::arg("independent_variables") = std::vector<
                              tudat::propulsion::ThrustIndependentVariables>(),
                          py::arg("guidance_input_functions") =
                              std::vector<std::function<double()>>());  //,


                    // Thrust orientation factory functions

                    m.def("constant_thrust_magnitude",
                          &tss::constantThrustMagnitudeSettings,
                          py::arg("thrust_magnitude"),
                          py::arg("specific_impulse"),
R"doc(Create thrust magnitude settings from a constant thrust magnitude and Isp.

	Factory function that creates constant thrust magnitude settings. The specific impulse to use for the thrust is
	also supplied when applying a mass rate model in the propagation of the vehicle dynamics, relating the thrust
	to the mass decrease of the vehicle.


	:param thrust_magnitude:
		Value of the constant thrust magnitude.
	:param specific_impulse:
		Value of the constant specific impulse, used to link the thrust model to the mass propagation.
	:return:
		Constant thrust magnitude settings object.
)doc");

                    m.def("custom_thrust_magnitude",
                          &tss::fromFunctionThrustMagnitudeSettings,
                          py::arg("thrust_magnitude_function"),
                          py::arg("specific_impulse_function"),
R"doc(Create thrust magnitude settings from a custom thrust force magnitude function.

	Factory function that creates thrust magnitude from a custom thrust force magnitude function.
	This model defines a thrust force and specific impulse that can vary with time. The thrust acceleration
	is computed during the propagation by dividing the thrust force by the current vehicle mass.
	The specific impulse can be used to apply a mass rate model in the propagation the vehicle dynamics, relating the thrust to the mass
	decrease of the vehicle.


	:param thrust_magnitude_function:
		Function of time returning the value of the thrust force magnitude.
	:param specific_impulse_function:
		Function of time returning the value of the specific impulse, useful to link the mass propagation to the thrust model.
	:return:
		From function thrust magnitude settings object.
)doc");

                    m.def("custom_thrust_magnitude_fixed_isp",
                          &tss::fromFunctionThrustMagnitudeFixedIspSettings,
                          py::arg("thrust_magnitude_function"),
                          py::arg("specific_impulse"),
R"doc(Same as :func:`~custom_thrust_magnitude`, but with a fixed value for the specific impulse.

	:param thrust_magnitude_function:
		Function of time returning the value of the thrust force magnitude.
	:param specific_impulse:
		Constant value for specific impulse, useful to link the mass propagation to the thrust model.
	:return:
		From function thrust magnitude settings object.
)doc");


                    m.def("custom_thrust_acceleration_magnitude",
                          &tss::customThrustAccelerationMagnitudeSettings,
                          py::arg("thrust_acceleration_magnitude_function"),
                          py::arg("specific_impulse_function"),
R"doc(Create thrust magnitude settings from a custom thrust acceleration magnitude function.

	Factory function that creates thrust magnitude from a custom thrust acceleration magnitude function.
	This model is similar to the :func:`~custom_thrust_magnitude`, with the difference being that this function
	directly provides the thrust *acceleration*, not the thrust *force*.


	:param thrust_acceleration_magnitude_function:
		Function of time returning the value of the thrust acceleration magnitude.
	:param specific_impulse_function:
		Function of time returning the value of the specific impulse, useful to link the mass propagation to the thrust model.
	:return:
		From function thrust magnitude settings object.
)doc");

                    m.def(
                        "custom_thrust_acceleration_magnitude_fixed_isp",
                        &tss::customThrustAccelerationMagnitudeFixedIspSettings,
                        py::arg("thrust_acceleration_magnitude_function"),
                        py::arg("specific_impulse"),
R"doc(Same as :func:`~custom_thrust_acceleration_magnitude`, but with a fixed value for the specific impulse.

	:param thrust_acceleration_magnitude_function:
		Function of time returning the value of the thrust acceleration magnitude.
	:param specific_impulse:
		Constant value for specific impulse, useful to link the mass propagation to the thrust model.
	:return:
		From function thrust magnitude settings object.
)doc");


                    /*!
                     *  To be removed, no longer used, kept only to support
                     * printing of useful error output when 'legacy' code is
                     * used
                     */

                    py::enum_<tss::ThrustDirectionTypes>(
                        m, "ThrustDirectionGuidanceTypes",
"")
                        .value(
                            "colinear_with_state_segment_thrust_direction_type",
                            tss::ThrustDirectionTypes::
                                colinear_with_state_segment_thrust_direction)
                        .value(
                            "thrust_direction_from_existing_body_orientation_"
                            "type",
                            tss::ThrustDirectionTypes::
                                thrust_direction_from_existing_body_orientation)
                        .value(
                            "custom_thrust_direction_type",
                            tss::ThrustDirectionTypes::custom_thrust_direction)
                        .value("custom_thrust_orientation_type",
                               tss::ThrustDirectionTypes::
                                   custom_thrust_orientation)
                        .value("mee_costate_based_thrust_direction_type",
                               tss::ThrustDirectionTypes::
                                   mee_costate_based_thrust_direction);


                    py::enum_<tss::ThrustFrames>(m, "ThrustFrames")
                        .value("unspecified_thrust_frame_type",
                               tss::ThrustFrames::unspecified_thrust_frame)
                        .value("inertial_thrust_frame_type",
                               tss::ThrustFrames::inertial_thrust_frame)
                        .value("tnw_thrust_frame_type",
                               tss::ThrustFrames::tnw_thrust_frame)
                        .export_values();

                    py::class_<tss::ThrustDirectionSettings,
                               std::shared_ptr<tss::ThrustDirectionSettings>>(
                        m, "ThrustDirectionSettings",
"")
                        .def_readonly(
                            "thrust_direction_type",
                            &tss::ThrustDirectionSettings::thrustDirectionType_)
                        .def_readonly(
                            "relative_body",
                            &tss::ThrustDirectionSettings::relativeBody_);


                    py::class_<
                        tss::ThrustDirectionFromStateGuidanceSettings,
                        std::shared_ptr<
                            tss::ThrustDirectionFromStateGuidanceSettings>,
                        tss::ThrustDirectionSettings>(
                        m, "ThrustDirectionFromStateGuidanceSettings",
"")
                        .def_readonly(
                            "is_colinear_with_velocity",
                            &tss::ThrustDirectionFromStateGuidanceSettings::
                                isColinearWithVelocity_)
                        .def_readonly(
                            "direction_is_opposite_to_vector",
                            &tss::ThrustDirectionFromStateGuidanceSettings::
                                directionIsOppositeToVector_);

                    py::class_<
                        tss::CustomThrustDirectionSettings,
                        std::shared_ptr<tss::CustomThrustDirectionSettings>,
                        tss::ThrustDirectionSettings>(
                        m, "CustomThrustDirectionSettings",
"")
                        .def_readonly("thrust_direction_function",
                                      &tss::CustomThrustDirectionSettings::
                                          thrustDirectionFunction_);

                    py::class_<
                        tss::CustomThrustOrientationSettings,
                        std::shared_ptr<tss::CustomThrustOrientationSettings>,
                        tss::ThrustDirectionSettings>(
                        m, "CustomThrustOrientationSettings",
"")
                        .def_readonly("thrust_orientation_function",
                                      &tss::CustomThrustOrientationSettings::
                                          thrustOrientationFunction_);


                    m.def("thrust_direction_from_state_guidance",
                          &tss::thrustDirectionFromStateGuidanceSettings,
                          py::arg("central_body"),
                          py::arg("is_colinear_with_velocity"),
                          py::arg("direction_is_opposite_to_vector"));

                    m.def("thrust_from_existing_body_orientation",
                          &tss::thrustFromExistingBodyOrientation);

                    m.def("custom_thrust_orientation",
                          py::overload_cast<
                              std::function<Eigen::Matrix3d(const double)>>(
                              &tss::customThrustOrientationSettings),
                          py::arg("thrust_orientation_function"));

                    m.def("custom_thrust_direction",
                          &tss::customThrustDirectionSettings,
                          py::arg("thrust_direction_function"));
                }

            }  // namespace thrust
        }  // namespace propagation_setup
    }  // namespace numerical_simulation
}  // namespace tudatpy
