/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include "expose_environment_setup.h"

#include "expose_environment_setup/expose_aerodynamic_coefficient_setup.h"
#include "expose_environment_setup/expose_atmosphere_setup.h"
#include "expose_environment_setup/expose_ephemeris_setup.h"
#include "expose_environment_setup/expose_gravity_field_setup.h"
#include "expose_environment_setup/expose_gravity_field_variation_setup.h"
#include "expose_environment_setup/expose_radiation_pressure_setup.h"
#include "expose_environment_setup/expose_ground_station_setup.h"
#include "expose_environment_setup/expose_rotation_model_setup.h"
#include "expose_environment_setup/expose_shape_setup.h"
#include "expose_environment_setup/expose_shape_deformation_setup.h"
#include "expose_environment_setup/expose_rigid_body_setup.h"
#include "expose_environment_setup/expose_vehicle_systems_setup.h"

#include "tudatpy/docstrings.h"
#include "tudatpy/scalarTypes.h"

#include <tudat/simulation/environment_setup.h>
#include <tudat/astro/reference_frames/referenceFrameTransformations.h>

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace te = tudat::ephemerides;
namespace ti = tudat::interpolators;
namespace tba = tudat::basic_astrodynamics;
namespace ta = tudat::aerodynamics;
namespace trf = tudat::reference_frames;
namespace tg = tudat::gravitation;
namespace tcc = tudat::coordinate_conversions;
namespace tp = tudat::propagators;



namespace tudatpy {
namespace numerical_simulation {
namespace environment_setup {

    void expose_environment_setup(py::module &m) {

//        m.def("get_body_gravitational_parameter",
//              &tss::getBodyGravitationalParameter,
//              py::arg("body_collection"), py::arg("body_name"));


        py::class_<tss::BodySettings, std::shared_ptr<tss::BodySettings>>(
                m, "BodySettings", get_docstring("BodySettings").c_str())
                .def_readwrite("constant_mass", &tss::BodySettings::constantMass, get_docstring("BodySettings.constant_mass").c_str())
                .def_readwrite("atmosphere_settings", &tss::BodySettings::atmosphereSettings, get_docstring("BodySettings.atmosphere_settings").c_str())
                .def_readwrite("ephemeris_settings", &tss::BodySettings::ephemerisSettings, get_docstring("BodySettings.ephemeris_settings").c_str())
                .def_readwrite("gravity_field_settings", &tss::BodySettings::gravityFieldSettings, get_docstring("BodySettings.gravity_field_settings").c_str())
                .def_readwrite("rotation_model_settings", &tss::BodySettings::rotationModelSettings, get_docstring("BodySettings.rotation_model_settings").c_str())
                .def_readwrite("shape_settings", &tss::BodySettings::shapeModelSettings, get_docstring("BodySettings.shape_settings").c_str())
                .def_readwrite("aerodynamic_coefficient_settings", &tss::BodySettings::aerodynamicCoefficientSettings, get_docstring("BodySettings.aerodynamic_coefficient_settings").c_str())
                .def_readwrite("gravity_field_variation_settings", &tss::BodySettings::gravityFieldVariationSettings, get_docstring("BodySettings.gravity_field_variation_settings").c_str())
                .def_readwrite("shape_deformation_settings", &tss::BodySettings::bodyDeformationSettings, get_docstring("BodySettings.shape_deformation_settings").c_str())
                .def_readwrite("ground_station_settings", &tss::BodySettings::groundStationSettings, get_docstring("BodySettings.ground_station_settings").c_str())
                .def_readwrite("rigid_body_settings", &tss::BodySettings::rigidBodyPropertiesSettings, get_docstring("BodySettings.rigidBodyPropertiesSettings").c_str())
                .def_readwrite("radiation_pressure_target_settings", &tss::BodySettings::radiationPressureTargetModelSettings, get_docstring("BodySettings.rigidBodyPropertiesSettings").c_str())
                .def_readwrite("radiation_source_settings", &tss::BodySettings::radiationSourceModelSettings, get_docstring("BodySettings.rigidBodyPropertiesSettings").c_str())
                .def_readwrite("vehicle_shape_settings", &tss::BodySettings::bodyExteriorPanelSettings_, get_docstring("BodySettings.vehicle_shape_settings").c_str())
                .def_readwrite("radiation_pressure_settings", &tss::BodySettings::radiationPressureSettings, get_docstring("BodySettings.radiation_pressure_settings").c_str());


        py::class_<tss::BodyListSettings,
                std::shared_ptr<tss::BodyListSettings> >(m, "BodyListSettings", get_docstring("BodyListSettings").c_str())
                .def(py::init<const std::string, const std::string>(),
                        py::arg("frame_origin"),
                        py::arg("frame_orientation"))
                .def("get", &tss::BodyListSettings::get, get_docstring("BodyListSettings.get").c_str())
                .def("add_settings", py::overload_cast<std::shared_ptr<tss::BodySettings>, const std::string>
                        (&tss::BodyListSettings::addSettings), py::arg("settings_to_add"), py::arg("body_name"))
                .def("add_empty_settings", py::overload_cast<const std::string>(&tss::BodyListSettings::addSettings),
                     py::arg("body_name"))
                .def_property_readonly("frame_origin", &tss::BodyListSettings::getFrameOrigin, get_docstring("BodyListSettings.frame_origin").c_str())
                .def_property_readonly("frame_orientation", &tss::BodyListSettings::getFrameOrientation, get_docstring("BodyListSettings.frame_orientation").c_str());

        m.def("get_default_body_settings",
              py::overload_cast<const std::vector<std::string> &, const std::string,
                      const std::string>(
                      &tss::getDefaultBodySettings),
              py::arg("bodies"),
              py::arg("base_frame_origin") = "SSB",
              py::arg("base_frame_orientation") = "ECLIPJ2000",
              get_docstring("get_default_body_settings").c_str());

        m.def("get_default_body_settings_time_limited",
              py::overload_cast<const std::vector<std::string> &,
                      const double, const double,
                      const std::string, const std::string, const double>(
                      &tss::getDefaultBodySettings),
              py::arg("bodies"),
              py::arg("initial_time"),
              py::arg("final_time"),
              py::arg("base_frame_origin") = "SSB",
              py::arg("base_frame_orientation") = "ECLIPJ2000",
              py::arg("time_step") = 300.0,
              get_docstring("get_default_body_settings_time_limited").c_str());

        m.def("get_default_single_body_settings",
              py::overload_cast<const std::string&, const std::string&>(
                      &tss::getDefaultSingleBodySettings),
              py::arg("body_name"),
              py::arg("base_frame_orientation") = "ECLIPJ2000",
              get_docstring("get_default_single_body_settings").c_str());

        m.def("get_default_single_body_settings_time_limited",
              py::overload_cast< const std::string&, const double, const double, const std::string&, const double >(
                  &tss::getDefaultSingleBodySettings),
              py::arg("body_name"),
              py::arg("initial_time"),
              py::arg("final_time"),
              py::arg("base_frame_orientation") = "ECLIPJ2000",
              py::arg("time_step") = 300.0,
              get_docstring("get_default_single_body_settings_time_limited").c_str());

        m.def("get_default_single_alternate_body_settings",
              py::overload_cast<const std::string&, const std::string&, const std::string&>(
                  &tss::getDefaultSingleAlternateNameBodySettings),
              py::arg("body_name"),
              py::arg("source_body_name"),
              py::arg("base_frame_orientation") = "ECLIPJ2000",
              get_docstring("get_default_single_alternate_body_settings").c_str());

        m.def("get_default_single_alternate_body_settings_time_limited",
              py::overload_cast< const std::string&, const std::string&, const double, const double, const std::string&, const double >(
                  &tss::getDefaultSingleAlternateNameBodySettings),
              py::arg("body_name"),
              py::arg("source_body_name"),
              py::arg("initial_time"),
              py::arg("final_time"),
              py::arg("base_frame_orientation") = "ECLIPJ2000",
              py::arg("time_step") = 300.0,
              get_docstring("get_default_single_alternate_body_settings_time_limited").c_str());

        m.def("create_simplified_system_of_bodies", &tss::createSimplifiedSystemOfBodies,
              py::arg("initial_time") = 0,
              get_docstring("create_simplified_system_of_bodies").c_str());

        m.def("create_system_of_bodies", &tss::createSystemOfBodies< double, TIME_TYPE >,
              py::arg("body_settings"),
              get_docstring("create_system_of_bodies").c_str());

        m.def("add_empty_tabulated_ephemeris", &tp::addEmptyTabulatedEphemeris< double, TIME_TYPE >,
              py::arg("bodies"),
              py::arg("body_name"),
              py::arg("ephemeris_origin") = "",
              get_docstring("add_empty_tabulated_ephemeris").c_str());

        m.def("create_tabulated_ephemeris_from_spice",
                &tss::createTabulatedEphemerisFromSpice<double, TIME_TYPE >, py::arg("body"),
                py::arg("initial_time"), py::arg("end_time"), py::arg("time_step"),
                py::arg("observer_name"), py::arg("reference_frame_name"),
                py::arg("interpolator_settings") =
                        std::make_shared<tudat::interpolators::LagrangeInterpolatorSettings>(
                                8));

        m.def("create_body_ephemeris", &tss::createBodyEphemeris< double, TIME_TYPE >,
              py::arg("ephemeris_settings"), py::arg("body_name"));

        m.def("create_ground_station_ephemeris",
              py::overload_cast< const std::shared_ptr< tss::Body >, const std::string& >(
                  &tss::createReferencePointEphemeris< TIME_TYPE, double > ),
              py::arg("body"),
              py::arg("station_name") );

        m.def("get_safe_interpolation_interval", &tss::getSafeInterpolationInterval,
              py::arg("ephemeris_model"));


        m.def("add_aerodynamic_coefficient_interface",
              &tss::addAerodynamicCoefficientInterface,
              py::arg("bodies"), py::arg("body_name"), py::arg("coefficient_settings"),
              get_docstring("add_aerodynamic_coefficient_interface").c_str());

        m.def("create_aerodynamic_coefficient_interface",
              &tss::createAerodynamicCoefficientInterfaceDeprecated,
              py::arg("coefficient_settings"), py::arg("body") );

        m.def("create_aerodynamic_coefficient_interface",
              &tss::createAerodynamicCoefficientInterface,
              py::arg("coefficient_settings"), py::arg("body"), py::arg("bodies"),
              get_docstring("create_aerodynamic_coefficient_interface").c_str());

        m.def("add_radiation_pressure_interface",
              &tss::addRadiationPressureInterface,
              py::arg("bodies"), py::arg("body_name"), py::arg("radiation_pressure_settings"),
              get_docstring("add_radiation_pressure_interface").c_str());


        m.def("add_radiation_pressure_target_model",
              &tss::addRadiationPressureTargetModel,
              py::arg("bodies"), py::arg("body_name"), py::arg("radiation_pressure_target_settings"),
              get_docstring("add_radiation_pressure_interface").c_str());

        m.def("add_rotation_model",
              &tss::addRotationModel,
              py::arg("bodies"), py::arg("body_name"), py::arg("rotation_model_settings"),
              get_docstring("add_rotation_model").c_str());

        m.def("add_mass_properties_model",
              &tss::addRigidBodyProperties,
              py::arg("bodies"), py::arg("body_name"), py::arg("mass_property_settings"),
              get_docstring("add_mass_properties_model").c_str());

        m.def("add_engine_model",
              &tss::addEngineModel,
              py::arg("body_name"),
              py::arg("engine_name"),
              py::arg("thrust_magnitude_settings"),
              py::arg("bodies"),
              py::arg("body_fixed_thrust_direction") = Eigen::Vector3d::UnitX( ),
              get_docstring("add_engine_model").c_str());

        m.def("add_variable_direction_engine_model",
              &tss::addVariableDirectionEngineModel,
              py::arg("body_name"),
              py::arg("engine_name"),
              py::arg("thrust_magnitude_settings"),
              py::arg("bodies"),
              py::arg("body_fixed_thrust_direction_function"),
              get_docstring("add_variable_direction_engine_model").c_str());


        m.def("add_flight_conditions",
              &tss::addFlightConditions,
              py::arg("bodies"), py::arg("body_name"), py::arg("central_body_name"),
              get_docstring("add_flight_conditions").c_str());

        m.def("add_ground_station",
              py::overload_cast<
                      const std::shared_ptr<tss::Body>,
                      const std::string,
                      const Eigen::Vector3d,
                      const tcc::PositionElementTypes,
                      const std::vector< std::shared_ptr< tss::GroundStationMotionSettings > > >(&tss::createGroundStation),
              py::arg("body"),
              py::arg("ground_station_name"),
              py::arg("ground_station_position"),
              py::arg("position_type") = tcc::cartesian_position,
              py::arg("station_motion_settings") = std::vector< std::shared_ptr< tss::GroundStationMotionSettings > >( ));

        m.def("add_ground_station",
              py::overload_cast<
              const std::shared_ptr< tss::Body >,
              const std::shared_ptr< tss::GroundStationSettings > >(&tss::createGroundStation),
              py::arg("body"),
              py::arg("ground_station_settings"),
              get_docstring("add_ground_station").c_str());

        m.def("create_radiation_pressure_interface",
              &tss::createRadiationPressureInterface,
              py::arg("radiationPressureInterfaceSettings"), py::arg("body_name"),
              py::arg("body_dict"));

        m.def("get_ground_station_list",
              &tss::getGroundStationsLinkEndList,
              py::arg( "body" ) );

        m.def("get_target_elevation_angles",
              &tss::getTargetElevationAngles,
              py::arg( "observing_body" ),
              py::arg( "target_body" ),
              py::arg( "station_name" ),
              py::arg( "times" ) );

        auto aerodynamic_coefficient_setup = m.def_submodule("aerodynamic_coefficients");
        aerodynamic_coefficients::expose_aerodynamic_coefficient_setup(aerodynamic_coefficient_setup);

        auto radiation_pressure_setup = m.def_submodule("radiation_pressure");
        radiation_pressure::expose_radiation_pressure_setup(radiation_pressure_setup);

        auto rotation_model_setup = m.def_submodule("rotation_model");
        rotation_model::expose_rotation_model_setup(rotation_model_setup);

        auto gravity_field_setup = m.def_submodule("gravity_field");
        gravity_field::expose_gravity_field_setup(gravity_field_setup);

        auto ephemeris_setup = m.def_submodule("ephemeris");
        ephemeris::expose_ephemeris_setup(ephemeris_setup);

        auto atmosphere_setup = m.def_submodule("atmosphere");
        atmosphere::expose_atmosphere_setup(atmosphere_setup);

        auto shape_setup = m.def_submodule("shape");
        shape::expose_shape_setup(shape_setup);

        auto gravity_variation_setup = m.def_submodule("gravity_field_variation");
        gravity_field_variation::expose_gravity_field_variation_setup(gravity_variation_setup);

        auto shape_deformation_setup = m.def_submodule("shape_deformation");
        shape_deformation::expose_shape_deformation_setup(shape_deformation_setup);

        auto ground_station_setup = m.def_submodule("ground_station");
        ground_station::expose_ground_station_setup(ground_station_setup);

        auto rigid_body_setup = m.def_submodule("rigid_body");
        rigid_body::expose_rigid_body_setup(rigid_body_setup);

        auto vehicle_systems_setup = m.def_submodule("vehicle_systems");
        vehicle_systems::expose_vehicle_systems_setup(vehicle_systems_setup);


//        auto system_model_setup = m.def_submodule("system_models");
//        gravity_field_variation::expose_system_model_setup(system_model_setup);



        // Function removed; error is shown
        m.def("set_aerodynamic_guidance",
              py::overload_cast<
                      const std::shared_ptr<ta::AerodynamicGuidance>,
                      const std::shared_ptr<tss::Body >,
                      const bool >
                      (&tss::setGuidanceAnglesFunctions),
              py::arg("aerodynamic_guidance"),
              py::arg("body"),
              py::arg("silence_warnings") = false );

        // Function removed; error is shown
        m.def("set_aerodynamic_orientation_functions", &tss::setAerodynamicOrientationFunctions,
              py::arg("body"),
              py::arg("angle_of_attack_function") = std::function<double()>(),
              py::arg("sideslip_angle_function") = std::function<double()>(),
              py::arg("bank_angle_function") = std::function<double()>(),
              py::arg("update_function") = std::function<void(const double)>());

        // Function removed; error is shown
        m.def("set_constant_aerodynamic_orientation", &tss::setConstantAerodynamicOrientation,
              py::arg("body"),
              py::arg("angle_of_attack"),
              py::arg("sideslip_angle"),
              py::arg("bank_angle"),
              py::arg("silence_warnings") = false );
    }

}// namespace environment_setup
}// namespace numerical_simulation
}// namespace tudatpy
