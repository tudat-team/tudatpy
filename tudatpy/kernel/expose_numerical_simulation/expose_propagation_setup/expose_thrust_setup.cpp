/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_acceleration_setup.h"
//#include "kernel/expose_numerical_simulation/deprecation_support.h"

#include "tudatpy/docstrings.h"
#include <tudat/simulation/propagation_setup.h>

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tinterp = tudat::interpolators;
namespace te = tudat::ephemerides;
namespace tni = tudat::numerical_integrators;
namespace trf = tudat::reference_frames;
namespace tmrf = tudat::root_finders;

namespace tudatpy {
namespace numerical_simulation {
namespace propagation_setup {
namespace thrust {

void expose_thrust_setup(py::module &m) {

    py::enum_<tss::ThrustMagnitudeTypes>(m, "ThrustMagnitudeTypes" )
//                                         get_docstring("ThrustMagnitudeTypes").c_str())
            .value("constant_thrust_magnitude", tss::ThrustMagnitudeTypes::constant_thrust_magnitude)
//            .value("from_engine_properties_thrust_magnitude", tss::ThrustMagnitudeTypes::from_engine_properties_thrust_magnitude)
            .value("thrust_magnitude_from_time_function", tss::ThrustMagnitudeTypes::thrust_magnitude_from_time_function)
            .value("thrust_magnitude_from_dependent_variables", tss::ThrustMagnitudeTypes::thrust_magnitude_from_dependent_variables);
//            .value("bang_bang_thrust_magnitude_from_mee_costates", tss::ThrustMagnitudeTypes::bang_bang_thrust_magnitude_from_mee_costates);

    py::class_<
            tss::ThrustMagnitudeSettings,
            std::shared_ptr<tss::ThrustMagnitudeSettings>>(m, "ThrustMagnitudeSettings",
                                                           get_docstring("ThrustMagnitudeSettings").c_str())
            .def_readonly("thrust_magnitude_type", &tss::ThrustMagnitudeSettings::thrustMagnitudeType_)
            .def_readonly("thrust_origin_id", &tss::ThrustMagnitudeSettings::thrustOriginId_);

    py::class_<
            tss::ConstantThrustMagnitudeSettings,
            std::shared_ptr<tss::ConstantThrustMagnitudeSettings>,
            tss::ThrustMagnitudeSettings>(m, "ConstantThrustMagnitudeSettings",
                                          get_docstring("ConstantThrustMagnitudeSettings").c_str())
            .def_readonly("thrust_magnitude", &tss::ConstantThrustMagnitudeSettings::thrustMagnitude_)
            .def_readonly("specific_impulse", &tss::ConstantThrustMagnitudeSettings::specificImpulse_);

    py::class_<
            tss::CustomThrustMagnitudeSettings,
            std::shared_ptr<tss::CustomThrustMagnitudeSettings>,
            tss::ThrustMagnitudeSettings>(m, "CustomThrustMagnitudeSettings",
                                          get_docstring("CustomThrustMagnitudeSettings").c_str());

    m.def("get_propulsion_input_variables",
          &tss::getPropulsionInputVariables,
          py::arg("body_with_guidance") = std::shared_ptr<tss::Body>(),
          py::arg("independent_variables") = std::vector<tudat::propulsion::ThrustIndependentVariables>(),
          py::arg("guidance_input_functions") = std::vector<std::function<double()>>() );//,
//          get_docstring("get_propulsion_input_variables").c_str());


    // Thrust orientation factory functions

    m.def("constant_thrust_magnitude", &tss::constantThrustMagnitudeSettings,
          py::arg("thrust_magnitude"),
          py::arg("specific_impulse"),
          get_docstring("constant_thrust_magnitude").c_str());

    m.def("custom_thrust_magnitude", &tss::fromFunctionThrustMagnitudeSettings,
          py::arg("thrust_magnitude_function"),
          py::arg("specific_impulse_function"),
          get_docstring("custom_thrust_magnitude").c_str());

    m.def("custom_thrust_magnitude_fixed_isp", &tss::fromFunctionThrustMagnitudeFixedIspSettings,
          py::arg("thrust_magnitude_function"),
          py::arg("specific_impulse"),
          get_docstring("custom_thrust_magnitude_fixed_isp").c_str());


    m.def("custom_thrust_acceleration_magnitude", &tss::customThrustAccelerationMagnitudeSettings,
          py::arg("thrust_acceleration_magnitude_function"),
          py::arg("specific_impulse_function"),
          get_docstring("custom_thrust_acceleration_magnitude").c_str());

    m.def("custom_thrust_acceleration_magnitude_fixed_isp", &tss::customThrustAccelerationMagnitudeFixedIspSettings,
          py::arg("thrust_acceleration_magnitude_function"),
          py::arg("specific_impulse"),
          get_docstring("custom_thrust_acceleration_magnitude_fixed_isp").c_str());





    /*!
     *  To be removed, no longer used, kept only to support printing of useful error output when 'legacy' code is used
     */

    py::enum_<tss::ThrustDirectionTypes>(m, "ThrustDirectionGuidanceTypes",
                                         get_docstring("ThrustDirectionGuidanceTypes").c_str())
            .value("colinear_with_state_segment_thrust_direction_type", tss::ThrustDirectionTypes::colinear_with_state_segment_thrust_direction)
            .value("thrust_direction_from_existing_body_orientation_type", tss::ThrustDirectionTypes::thrust_direction_from_existing_body_orientation)
            .value("custom_thrust_direction_type", tss::ThrustDirectionTypes::custom_thrust_direction)
            .value("custom_thrust_orientation_type", tss::ThrustDirectionTypes::custom_thrust_orientation)
            .value("mee_costate_based_thrust_direction_type", tss::ThrustDirectionTypes::mee_costate_based_thrust_direction);


    py::enum_<tss::ThrustFrames>(m, "ThrustFrames" )
//                                 get_docstring("ThrustFrames").c_str())
            .value("unspecified_thrust_frame_type", tss::ThrustFrames::unspecified_thrust_frame)
            .value("inertial_thrust_frame_type", tss::ThrustFrames::inertial_thrust_frame)
            .value("tnw_thrust_frame_type", tss::ThrustFrames::tnw_thrust_frame)
            .export_values();

    py::class_<
            tss::ThrustDirectionSettings,
            std::shared_ptr<tss::ThrustDirectionSettings>>(m, "ThrustDirectionSettings",
                    get_docstring("ThrustDirectionSettings").c_str())
            .def_readonly("thrust_direction_type", &tss::ThrustDirectionSettings::thrustDirectionType_)
            .def_readonly("relative_body", &tss::ThrustDirectionSettings::relativeBody_);


    py::class_<
            tss::ThrustDirectionFromStateGuidanceSettings,
            std::shared_ptr<tss::ThrustDirectionFromStateGuidanceSettings>,
            tss::ThrustDirectionSettings>(m, "ThrustDirectionFromStateGuidanceSettings",
                                          get_docstring("ThrustDirectionFromStateGuidanceSettings").c_str())
            .def_readonly("is_colinear_with_velocity", &tss::ThrustDirectionFromStateGuidanceSettings::isColinearWithVelocity_)
            .def_readonly("direction_is_opposite_to_vector", &tss::ThrustDirectionFromStateGuidanceSettings::directionIsOppositeToVector_);

    py::class_<
            tss::CustomThrustDirectionSettings,
            std::shared_ptr<tss::CustomThrustDirectionSettings>,
            tss::ThrustDirectionSettings>(m, "CustomThrustDirectionSettings",
                                          get_docstring("CustomThrustDirectionSettings").c_str())
            .def_readonly("thrust_direction_function", &tss::CustomThrustDirectionSettings::thrustDirectionFunction_);

    py::class_<
            tss::CustomThrustOrientationSettings,
            std::shared_ptr<tss::CustomThrustOrientationSettings>,
            tss::ThrustDirectionSettings>(m, "CustomThrustOrientationSettings",
                                          get_docstring("CustomThrustOrientationSettings").c_str())
            .def_readonly("thrust_orientation_function", &tss::CustomThrustOrientationSettings::thrustOrientationFunction_);



        m.def("thrust_direction_from_state_guidance", &tss::thrustDirectionFromStateGuidanceSettings,
              py::arg( "central_body"),
              py::arg("is_colinear_with_velocity"),
              py::arg("direction_is_opposite_to_vector") );

        m.def("thrust_from_existing_body_orientation", &tss::thrustFromExistingBodyOrientation );

        m.def("custom_thrust_orientation",
              py::overload_cast< std::function< Eigen::Matrix3d( const double ) > >(
                      &tss::customThrustOrientationSettings ),
              py::arg( "thrust_orientation_function" ) );

        m.def("custom_thrust_direction", &tss::customThrustDirectionSettings,
              py::arg( "thrust_direction_function" ) );



}

}// namespace thrust
}// namespace propagation_setup
}// namespace numerical_simulation
}// namespace tudatpy
