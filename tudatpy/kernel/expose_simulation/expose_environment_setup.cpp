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
#include "expose_environment_setup/expose_rotation_model_setup.h"
#include "expose_environment_setup/expose_shape_setup.h"

#include "tudatpy/docstrings.h"
#include <tudat/simulation/environment_setup.h>
#include <tudat/astro/reference_frames/referenceFrameTransformations.h>

//#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
//#include <pybind11/numpy.h>
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



namespace tudat {
namespace simulation_setup {

inline std::shared_ptr< GravityFieldVariationSettings > fixedSingleDegreeLoveNumberGravityFieldVariationSettingsSimplified(
        const std::string deformingBody,
        const double loveNumber,
        const int degree )
{
    return fixedSingleDegreeLoveNumberGravityFieldVariationSettings(
                deformingBody, loveNumber, degree );
}

} // namespace simulation_setup
} // namespace tudat


namespace tudatpy {
namespace simulation {
namespace environment_setup {

    void expose_environment_setup(py::module &m) {

        /*  This exposition module follows the structure of
       *  tudat/src/simulation/environment_setup/
       *  environment_setup
       *  ├── body.h
       *  ├── createAerodynamicCoefficientInterface.h
       *  ├── createAerodynamicControlSurfaces.h
       *  ├── createAtmosphereModel.h
       *  ├── createBodies.h
       *  ├── createBodyShapeModel.h
       *  ├── createEphemeris.h
       *  ├── createFlightConditions.h
       *  ├── createGravityField.h
       *  ├── createGravityFieldVariations.h
       *  ├── createGroundStations.h
       *  ├── createRadiationPressureInterface.h
       *  ├── createRotationModel.h
       *  └── defaultBodies.h
       *
       *  environment_setup/
       *  ├── body.cpp
       *  ├── CMakeLists.txt
       *  ├── createAerodynamicCoefficientInterface.cpp
       *  ├── createAerodynamicControlSurfaces.cpp
       *  ├── createAtmosphereModel.cpp
       *  ├── createBodies.cpp
       *  ├── createBodyShapeModel.cpp
       *  ├── createEphemeris.cpp
       *  ├── createFlightConditions.cpp
       *  ├── createGravityField.cpp
       *  ├── createGravityFieldVariations.cpp
       *  ├── createGroundStations.cpp
       *  ├── createRadiationPressureInterface.cpp
       *  ├── createRotationModel.cpp
       *  └── defaultBodies.cpp
       */

        /////////////////////////////////////////////////////////////////////////////
        // body.h ///////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////
        py::class_<tss::Body, std::shared_ptr<tss::Body>>(m, "Body")
//            .def(py::init<const Eigen::Vector6d &>(), py::arg("state") = Eigen::Vector6d::Zero())
                .def("get_ephemeris_frame_to_base_frame", &tss::Body::getEphemerisFrameToBaseFrame)
                .def("set_ephemeris_frame_to_base_frame", &tss::Body::setEphemerisFrameToBaseFrame)
                .def_property("ephemeris_frame_to_base_frame", &tss::Body::getEphemerisFrameToBaseFrame,
                              &tss::Body::setEphemerisFrameToBaseFrame)
                .def("get_state", &tss::Body::getState)
                .def("set_state", &tss::Body::setState)
                .def("get_state_in_based_frame_from_ephemeris",
                     &tss::Body::getStateInBaseFrameFromEphemeris<double, double>)
                .def("get_ephemeris", &tss::Body::getEphemeris)
                .def("set_ephemeris", &tss::Body::setEphemeris)
                .def("get_atmosphere_model", &tss::Body::getAtmosphereModel)
                .def("set_atmosphere_model", &tss::Body::setAtmosphereModel)
                .def("get_gravity_field_model", &tss::Body::getGravityFieldModel)
                .def("get_shape_model", &tss::Body::getShapeModel)
                .def("set_shape_model", &tss::Body::setShapeModel)
                .def_property("shape_model", &tss::Body::getShapeModel, &tss::Body::setShapeModel)
                .def("set_gravity_field_model", &tss::Body::setGravityFieldModel)
                .def("get_gravitational_parameter", &tss::Body::getGravitationalParameter)
                .def_property_readonly("gravitational_parameter", &tss::Body::getGravitationalParameter)
                .def_property("gravity_field_model", &tss::Body::getGravityFieldModel, &tss::Body::setGravityFieldModel)
                .def("get_aerodynamic_coefficient_interface", &tss::Body::getAerodynamicCoefficientInterface)
                .def("set_aerodynamic_coefficient_interface", &tss::Body::setAerodynamicCoefficientInterface)
                .def_property("aerodynamic_coefficient_interface", &tss::Body::getAerodynamicCoefficientInterface,
                              &tss::Body::setAerodynamicCoefficientInterface)
                .def("get_body_mass", &tss::Body::getBodyMass)
                .def("set_constant_mass", &tss::Body::setConstantBodyMass)
//            .def("get_radiation_pressure_interfaces", &tss::Body::getRadiationPressureInterfaces)
//            .def("set_radiation_pressure_interface", &tss::Body::setRadiationPressureInterface,
//                 py::arg( "radiating_body" ),
//                 py::arg( "radiation_pressure_interface" ) )
                .def("set_aerodynamic_coefficient_interface", &tss::Body::setAerodynamicCoefficientInterface)
                .def("get_aerodynamic_coefficient_interface", &tss::Body::getAerodynamicCoefficientInterface)
                .def("get_flight_conditions", &tss::Body::getFlightConditions)
                .def("set_flight_conditions", &tss::Body::setFlightConditions, py::arg("aerodynamic_flight_conditions"))
                .def_property("flight_conditions", &tss::Body::getFlightConditions, &tss::Body::setFlightConditions)
                .def("get_rotation_model", &tss::Body::getRotationalEphemeris)
                .def("set_rotation_model", &tss::Body::setRotationalEphemeris, py::arg("rotational_ephemeris"))
                .def_property("rotation_model", &tss::Body::getRotationalEphemeris, &tss::Body::setRotationalEphemeris)
                .def_property("inertia_tensor", &tss::Body::getBodyInertiaTensor,
                              py::overload_cast<const Eigen::Matrix3d &>(
                                      &tss::Body::setBodyInertiaTensor))
                .def_property_readonly("state", &tss::Body::getState)
                .def_property_readonly("position", &tss::Body::getPosition)
                .def_property_readonly("velocity", &tss::Body::getVelocity)
                .def_property_readonly("angular_velocity_body_fixed",
                                       &tss::Body::getCurrentAngularVelocityVectorInLocalFrame)
                .def_property_readonly("mass", &tss::Body::getBodyMass);


        py::class_<tss::SystemOfBodies,
                std::shared_ptr<tss::SystemOfBodies> >(m, "SystemOfBodies")
//            .def(py::init<//ctor 1
//                 const std::string,
//                 const std::string,
//                 const std::unordered_map< std::string, std::shared_ptr< tss::Body > > >(),
//                 py::arg("frame_origin") = "SSB",
//                 py::arg("frame_orientation") = "ECLIPJ2000",
//                 py::arg("body_map") =
//            std::unordered_map< std::string, std::shared_ptr< tss::Body > >( ) )
                .def("get_body", &tss::SystemOfBodies::getBody)
                .def("count", &tss::SystemOfBodies::count)
                .def("create_empty_body", &tss::SystemOfBodies::createEmptyBody,
                     py::arg("body_name"),
                     py::arg("process_body") = 1)
                .def("add_body", &tss::SystemOfBodies::addBody,
                     py::arg("body_to_add"),
                     py::arg("body_name"),
                     py::arg("process_body") = 1)
                .def("delete_body", &tss::SystemOfBodies::deleteBody,
                     py::arg("body_name"))
                .def("processBodyFrameDefinitions", &tss::SystemOfBodies::processBodyFrameDefinitions)
                .def(py::pickle(
                        [](const tss::SystemOfBodies &p) { // __getstate__
                            /* Return a tuple that fully encodes the state of the object */
                            return py::make_tuple(p.getFrameOrigin(), p.getMap());
                        },
                        [](py::tuple t) { // __setstate__
                            if (t.size() != 3)
                                throw std::runtime_error("Invalid state for SystemOfBodies!");

                            /* Create a new C++ instance */
                            tss::SystemOfBodies p(t[0].cast<std::string>());

                            return p;
                        }
                ));


        m.def("get_body_gravitational_parameter",
              &tss::getBodyGravitationalParameter,
              py::arg("body_collection"), py::arg("body_name"));

        /////////////////////////////////////////////////////////////////////////////
        // createBodies.h ///////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////
        py::class_<tss::BodySettings, std::shared_ptr<tss::BodySettings>>(
                m, "BodySettings", get_docstring("BodySettings").c_str())
                .def_readwrite("constant_mass", &tss::BodySettings::constantMass)
                .def_readwrite("atmosphere_settings", &tss::BodySettings::atmosphereSettings)
                .def_readwrite("ephemeris_settings", &tss::BodySettings::ephemerisSettings)
                .def_readwrite("gravity_field_settings", &tss::BodySettings::gravityFieldSettings)
                .def_readwrite("rotation_model_settings", &tss::BodySettings::rotationModelSettings)
                .def_readwrite("shape_settings", &tss::BodySettings::shapeModelSettings)
                .def_readwrite("radiation_pressure_settings", &tss::BodySettings::radiationPressureSettings)
                .def_readwrite("aerodynamic_coefficient_settings", &tss::BodySettings::aerodynamicCoefficientSettings)
                .def_readwrite("gravity_field_variation_settings", &tss::BodySettings::gravityFieldVariationSettings)
                .def_readwrite("ground_station_settings", &tss::BodySettings::groundStationSettings);

        py::class_<tss::BodyListSettings,
                std::shared_ptr<tss::BodyListSettings> >(m, "BodyListSettings")
//            .def(py::init<//ctor 1
//                 const std::string,
//                 const std::string >(),
//                 py::arg("frame_origin") ,
//                 py::arg("frame_orientation"))
//            .def(py::init<//ctor 2
//                 const std::map< std::string, std::shared_ptr< tss::BodySettings > >&,
//                 const std::string,
//                 const std::string >(),
//                 py::arg("body_settings"),
//                 py::arg("frame_origin"),
//                 py::arg("frame_orientation"))
                .def("at", &tss::BodyListSettings::at)
                .def("get", &tss::BodyListSettings::get)
                .def("add_settings", py::overload_cast<const std::string>(&tss::BodyListSettings::addSettings),
                     py::arg("body_name"))
                .def("add_settings", py::overload_cast<std::shared_ptr<tss::BodySettings>, const std::string>
                        (&tss::BodyListSettings::addSettings), py::arg("settings_to_add"), py::arg("body_name"))
                .def("add_empty_settings", py::overload_cast<const std::string>(&tss::BodyListSettings::addSettings),
                     py::arg("body_name"))
                .def("get_frame_origin", &tss::BodyListSettings::getFrameOrigin)
                .def("get_frame_orientation", &tss::BodyListSettings::getFrameOrientation);

        /////////////////////////////////////////////////////////////////////////////
        // defaultBodies.h //////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////

        // getDefaultBodySettings (overload 1)
        m.def("get_default_body_settings",
              py::overload_cast<const std::vector<std::string> &, const std::string,
                      const std::string>(
                      &tss::getDefaultBodySettings),
              py::arg("bodies"),
              py::arg("base_frame_origin") = "SSB",
              py::arg("base_frame_orientation") = "ECLIPJ2000");

        // getDefaultBodySettings (overload 2)
        m.def("get_default_body_settings",
              py::overload_cast<const std::vector<std::string> &,
                      const double, const double,
                      const std::string, const std::string, const double>(
                      &tss::getDefaultBodySettings),
              py::arg("bodies"),
              py::arg("initial_time"),
              py::arg("final_time"),
              py::arg("base_frame_origin") = "SSB",
              py::arg("base_frame_orientation") = "ECLIPJ2000",
              py::arg("time_step") = 300.0);


        m.def("get_global_frame_origin", &tss::getGlobalFrameOrigin);

        m.def("create_system_of_bodies", &tss::createSystemOfBodies);

        m.def("add_empty_tabulate_ephemeris", &tss::addEmptyTabulateEphemeris,
              py::arg("bodies"),
              py::arg("body_name"),
              py::arg("ephemeris_origin") = "");

        // Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.cpp
        m.def("create_tabulated_ephemeris_from_spice",
                &tss::createTabulatedEphemerisFromSpice<>, py::arg("body"),
                py::arg("initial_time"), py::arg("end_time"), py::arg("time_step"),
                py::arg("observer_name"), py::arg("reference_frame_name"),
                py::arg("interpolator_settings") =
                        std::make_shared<tudat::interpolators::LagrangeInterpolatorSettings>(
                                8));

        // Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.cpp
        m.def("create_body_ephemeris", &tss::createBodyEphemeris,
              py::arg("ephemeris_settings"), py::arg("body_name"));

        m.def("get_safe_interpolation_interval", &tss::getSafeInterpolationInterval,
              py::arg("ephemeris_model"));


        m.def("add_aerodynamic_coefficient_interface",
              &tss::addAerodynamicCoefficientInterface,
              py::arg("bodies"), py::arg("body_name"), py::arg("coefficient_settings"));

        m.def("create_aerodynamic_coefficient_interface",
              &tss::createAerodynamicCoefficientInterface,
              py::arg("coefficient_settings"), py::arg("body"));

        m.def("add_radiation_pressure_interface",
              &tss::addRadiationPressureInterface,
              py::arg("bodies"), py::arg("body_name"), py::arg("radiation_pressure_settings"));

        m.def("add_ground_station",
              py::overload_cast<
                      const std::shared_ptr<tss::Body>,
                      const std::string,
                      const Eigen::Vector3d,
                      const tcc::PositionElementTypes>(&tss::createGroundStation),
              py::arg("body"),
              py::arg("ground_station_name"),
              py::arg("ground_station_position"),
              py::arg("position_type") = tcc::cartesian_position);

        m.def("create_radiation_pressure_interface",
              &tss::createRadiationPressureInterface,
              py::arg("radiationPressureInterfaceSettings"), py::arg("body_name"),
              py::arg("body_dict"));

        m.def("set_aerodynamic_guidance",
              py::overload_cast<
                      const std::shared_ptr<ta::AerodynamicGuidance>,
                      const std::shared_ptr<tss::Body> >
                      (&tss::setGuidanceAnglesFunctions),
              py::arg("aerodynamic_guidance"),
              py::arg("body"));

        m.def("set_aerodynamic_orientation_functions", &tss::setAerodynamicOrientationFunctions,
              py::arg("body"),
              py::arg("angle_of_attack_function") = std::function<double()>(),
              py::arg("sideslip_angle_function") = std::function<double()>(),
              py::arg("bank_angle_function") = std::function<double()>(),
              py::arg("update_function") = std::function<void(const double)>());

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


    }

}// namespace environment_setup
}// namespace simulation
}// namespace tudatpy
