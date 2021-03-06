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

#include "../docstrings.h"
#include <tudat/simulation/environment_setup.h>

//#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
//#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace te = tudat::ephemerides;
namespace ti = tudat::interpolators;

namespace tudatpy {

void expose_aerodynamic_coefficient_setup(py::module &m) {

    /////////////////////////////////////////////////////////////////////////////
    // createAerodynamicCoefficientInterface.h
    /////////////////////////////////////////////////////////////////////////////
    py::class_<tss::AerodynamicCoefficientSettings,
            std::shared_ptr<tss::AerodynamicCoefficientSettings>>
            AerodynamicCoefficientSettings_(m, "AerodynamicCoefficientSettings",
                                            "<no doc>");

    py::class_<tss::ConstantAerodynamicCoefficientSettings,
            std::shared_ptr<tss::ConstantAerodynamicCoefficientSettings>,
            tss::AerodynamicCoefficientSettings>(
                m, "ConstantAerodynamicCoefficientSettings", "<no doc>")
            .def(py::init<const double, const double, const double,
                 const Eigen::Vector3d &, const Eigen::Vector3d &,
                 const Eigen::Vector3d &, const bool, const bool,
                 const std::shared_ptr<ti::InterpolatorSettings>>(),
                 py::arg("reference_length"), py::arg("reference_area"),
                 py::arg("lateral_reference_length"),
                 py::arg("moment_reference_point"),
                 py::arg("constant_force_coefficient"),
                 py::arg("constant_moment_coefficient") = Eigen::Vector3d::Zero(),
                 py::arg("are_coefficients_in_aerodynamic_frame") = true,
                 py::arg("are_coefficients_in_negative_axis_direction") = true,
                 py::arg("interpolator_settings") = nullptr)
            .def(py::init<const double, const Eigen::Vector3d &, const bool,
                 const bool>(),
                 py::arg("reference_area"),
                 py::arg("constant_force_coefficient"),
                 py::arg("are_coefficients_in_aerodynamic_frame") = true,
                 py::arg("are_coefficients_in_negative_axis_direction") = true);

    m.def("constant",
          py::overload_cast<const double, const Eigen::Vector3d &, const bool,
          const bool>(
              &tss::constantAerodynamicCoefficientSettings),
          py::arg("reference_area"),
          py::arg("constant_force_coefficient"),
          py::arg("are_coefficients_in_aerodynamic_frame") = true,
          py::arg("are_coefficients_in_negative_axis_direction") = true);
}

void expose_atmosphere_setup(py::module &m) {

    /////////////////////////////////////////////////////////////////////////////
    py::class_<tss::AtmosphereSettings,
            std::shared_ptr<tss::AtmosphereSettings>>
            AtmosphereSettings(m, "AtmosphereSettings",
                                            "<no doc>");
    m.def("exponential",
          py::overload_cast<const double, const double>(
          &tss::exponentialAtmosphereSettings ),
          py::arg("scale_height"),
          py::arg("surface_density") );
}

void expose_radiation_pressure_setup(py::module &m) {

    /////////////////////////////////////////////////////////////////////////////
    // createRadiationPressureInterface.h
    /////////////////////////////////////////////////////////////////////////////
    py::enum_<tss::RadiationPressureType>(m, "RadiationPressureType", "<no_doc>")
            .value(
                "cannon_ball_radiation_pressure_interface",
                tss::RadiationPressureType::cannon_ball_radiation_pressure_interface)
            .value("panelled_radiation_pressure_interface",
                   tss::RadiationPressureType::panelled_radiation_pressure_interface)
            .value("solar_sailing_radiation_pressure_interface",
                   tss::RadiationPressureType::
                   solar_sailing_radiation_pressure_interface)
            .export_values();


    py::class_<tss::RadiationPressureInterfaceSettings,
            std::shared_ptr<tss::RadiationPressureInterfaceSettings>>(
                m, "RadiationPressureInterfaceSettings", "<no_doc>")
            .def(py::init<const tss::RadiationPressureType, const std::string &,
                 const std::vector<std::string>>(),
                 py::arg("radiation_pressure_type"), py::arg("source_body"),
                 py::arg("occulting_bodies") = std::vector<std::string>());

    py::class_<tss::CannonBallRadiationPressureInterfaceSettings,
            std::shared_ptr<tss::CannonBallRadiationPressureInterfaceSettings>,
            tss::RadiationPressureInterfaceSettings>(
                m, "CannonBallRadiationPressureInterfaceSettings", "<no_doc>")
            .def(py::init<const std::string &, const double, const double,
                 const std::vector<std::string> &>(),
                 py::arg("source_body"), py::arg("area"),
                 py::arg("radiation_pressure_coefficient"),
                 py::arg("occulting_bodies") = std::vector<std::string>());

    m.def("cannonball",
          &tss::cannonBallRadiationPressureSettings,
          py::arg("source_body"), py::arg("reference_area"),
          py::arg("radiation_pressure_coefficient"),
          py::arg("occulting_bodies") = std::vector<std::string>());
}

void expose_rotation_model_setup(py::module &m) {
    /////////////////////////////////////////////////////////////////////////////
    // createRotationalModel.h
    /////////////////////////////////////////////////////////////////////////////
    py::enum_<tss::RotationModelType>(m, "RotationModelType", "<no doc>")
            .value("simple_rotational_model",
                   tss::RotationModelType::simple_rotation_model)
            .value("spice_rotation_model",
                   tss::RotationModelType::spice_rotation_model)
            .value("gcrs_to_itrs_rotation_model",
                   tss::RotationModelType::gcrs_to_itrs_rotation_model)
            .value("synchronous_rotation_model",
                   tss::RotationModelType::synchronous_rotation_model)
            .value("planetary_rotation_model",
                   tss::RotationModelType::planetary_rotation_model)
            .export_values();

    py::class_<tss::RotationModelSettings,
            std::shared_ptr<tss::RotationModelSettings>>(
                m, "RotationalModelSettings", "<no doc>")
            .def(py::init<const tss::RotationModelType, const std::string &,
                 const std::string &>(),
                 py::arg("rotation_type"), py::arg("original_frame"),
                 py::arg("target_frame"))
            .def("get_rotation_type", &tss::RotationModelSettings::getRotationType)
            .def("get_original_frame", &tss::RotationModelSettings::getOriginalFrame)
            .def("get_target_frame", &tss::RotationModelSettings::getTargetFrame)
            .def("reset_original_frame",
                 &tss::RotationModelSettings::resetOriginalFrame);
}

void expose_gravity_field_setup(py::module &m) {
    /////////////////////////////////////////////////////////////////////////////
    // createGravityField.h
    /////////////////////////////////////////////////////////////////////////////
    py::enum_<tss::GravityFieldType>(m, "GravityFieldType", "<no doc>")
            .value("central", tss::GravityFieldType::central)
            .value("central_spice", tss::GravityFieldType::central_spice)
            .value("spherical_harmonic", tss::GravityFieldType::spherical_harmonic)
            .export_values();

    py::enum_<tss::SphericalHarmonicsModel>(m, "SphericalHarmonicsModel", "<no doc>")
            .value("custom_model", tss::SphericalHarmonicsModel::customModel)
            .value("egm96", tss::SphericalHarmonicsModel::egm96)
            .value("ggm02c", tss::SphericalHarmonicsModel::ggm02c)
            .value("ggm02s", tss::SphericalHarmonicsModel::ggm02s)
            .value("glgm3150", tss::SphericalHarmonicsModel::glgm3150)
            .value("lpe200", tss::SphericalHarmonicsModel::lpe200)
            .value("jgmro120d", tss::SphericalHarmonicsModel::jgmro120d)
            .export_values();

    py::class_<tss::GravityFieldSettings, std::shared_ptr<tss::GravityFieldSettings>>(
                m, "GravityFieldSettings", "<no doc>")
            .def(py::init<const tss::GravityFieldType>(),
                 py::arg("gravity_field_type"))
            .def("get_gravity_field_type", &tss::GravityFieldSettings::getGravityFieldType);

    py::class_<tss::CentralGravityFieldSettings, std::shared_ptr<tss::CentralGravityFieldSettings>,
            tss::GravityFieldSettings>(m, "CentralGravityFieldSettings", "<no doc>")
            .def(py::init<double>(), py::arg("gravitational_parameter"));

    py::class_<tss::SphericalHarmonicsGravityFieldSettings, std::shared_ptr<tss::SphericalHarmonicsGravityFieldSettings>,
            tss::GravityFieldSettings>(m, "SphericalHarmonicsGravityFieldSettings", "<no doc>")
            .def(py::init<const double, const double, const Eigen::MatrixXd, const Eigen::MatrixXd, const std::string&>(),
                 py::arg("gravitational_parameter"), py::arg("reference_radius"), py::arg("cosine_coefficients"),
                 py::arg("sine_coefficients"), py::arg("associated_reference_frame"))
            .def("get_gravitational_parameter", &tss::SphericalHarmonicsGravityFieldSettings::getGravitationalParameter)
            .def("reset_gravitational_parameter", &tss::SphericalHarmonicsGravityFieldSettings::resetGravitationalParameter)
            .def("get_reference_radius", &tss::SphericalHarmonicsGravityFieldSettings::getReferenceRadius)
            .def("get_cosine_coefficients", &tss::SphericalHarmonicsGravityFieldSettings::getCosineCoefficients)
            .def("get_sine_coefficients", &tss::SphericalHarmonicsGravityFieldSettings::getSineCoefficients)
            .def("get_associated_reference_frame", &tss::SphericalHarmonicsGravityFieldSettings::getAssociatedReferenceFrame)
            .def("reset_associated_reference_frame", &tss::SphericalHarmonicsGravityFieldSettings::resetAssociatedReferenceFrame)
            .def("get_create_time_dependent_field", &tss::SphericalHarmonicsGravityFieldSettings::getCreateTimeDependentField)
            .def("set_create_time_dependent_field", &tss::SphericalHarmonicsGravityFieldSettings::setCreateTimeDependentField);
}

void expose_ephemeris_setup(py::module &m) {

    /////////////////////////////////////////////////////////////////////////////
    // createEphemeris.h (complete, unverified)
    /////////////////////////////////////////////////////////////////////////////
    py::class_<tss::EphemerisSettings,
            std::shared_ptr<tss::EphemerisSettings>>(m, "EphemerisSettings")
            .def(py::init<const tss::EphemerisType,
                 const std::string &,
                 const std::string &>(),
                 py::arg("ephemeris_type"),
                 py::arg("frame_origin") = "SSB",
                 py::arg("frame_orientation") = "ECLIPJ2000")
            .def("get_ephemeris_type", &tss::EphemerisSettings::getEphemerisType)
            .def("get_frame_origin", &tss::EphemerisSettings::getFrameOrigin)
            .def("get_frame_orientation", &tss::EphemerisSettings::getFrameOrientation)
            .def("get_multi_arc_ephemeris", &tss::EphemerisSettings::getMakeMultiArcEphemeris)
            .def("reset_frame_origin", &tss::EphemerisSettings::resetFrameOrigin)
            .def("reset_frame_orientation", &tss::EphemerisSettings::resetFrameOrientation)
            .def("reset_make_multi_arc_ephemeris", &tss::EphemerisSettings::resetMakeMultiArcEphemeris);

    py::class_<tss::DirectSpiceEphemerisSettings,
            std::shared_ptr<tss::DirectSpiceEphemerisSettings>,
            tss::EphemerisSettings>(m, "DirectSpiceEphemerisSettings")
            .def(py::init<const std::string, const std::string, const bool,
                 const bool, const bool, const tss::EphemerisType>(),
                 py::arg("frame_origin") = "SSB",
                 py::arg("frame_orientation") = "ECLIPJ2000",
                 py::arg("correct_for_stellar_aberration") = false,
                 py::arg("correct_for_light_time_aberration") = false,
                 py::arg("converge_light_time_aberration") = false,
                 py::arg("ephemeris_type") = tss::direct_spice_ephemeris)
            .def("get_correct_for_steller_aberration", &tss::DirectSpiceEphemerisSettings::getCorrectForStellarAberration)
            .def("get_correct_for_steller_aberration", &tss::DirectSpiceEphemerisSettings::getCorrectForLightTimeAberration)
            .def("get_converge_light_time_aberration",
                 // TODO : Fix getConvergeLighTimeAberration typo in Tudat.
                 &tss::DirectSpiceEphemerisSettings::getConvergeLighTimeAberration);

    py::class_<tss::InterpolatedSpiceEphemerisSettings,
            std::shared_ptr<tss::InterpolatedSpiceEphemerisSettings>,
            tss::DirectSpiceEphemerisSettings>(m, "InterpolatedSpiceEphemerisSettings")
            .def(py::init<
                 double, double, double, std::string, std::string,
                 std::shared_ptr<tudat::interpolators::InterpolatorSettings>>(),
                 py::arg("initial_time"), py::arg("final_time"), py::arg("time_step"),
                 py::arg("frame_origin") = "SSB",
                 py::arg("frame_orientation") = "ECLIPJ2000",
                 py::arg("interpolator_settings") = std::make_shared<
            tudat::interpolators::LagrangeInterpolatorSettings>(6));

    py::class_<tss::ApproximatePlanetPositionSettings,
            std::shared_ptr<tss::ApproximatePlanetPositionSettings>,
            tss::EphemerisSettings>(m, "ApproximatePlanetPositionSettings")
            .def(py::init<const tudat::ephemerides::ApproximatePlanetPositionsBase::
                 BodiesWithEphemerisData,
                 const bool>(),
                 py::arg("body_identifier"),
                 py::arg("use_circular_coplanar_approximation"))
            .def("get_body_identifier",
                 &tss::ApproximatePlanetPositionSettings::getBodyIdentifier)
            .def("get_use_circular_coplanar_approximation",
                 &tss::ApproximatePlanetPositionSettings::
                 getUseCircularCoplanarApproximation);

    py::class_<tss::ConstantEphemerisSettings,
            std::shared_ptr<tss::ConstantEphemerisSettings>,
            tss::EphemerisSettings>(m, "ConstantEphemerisSettings")
            .def(py::init<const Eigen::Vector6d &,
                 const std::string &,
                 const std::string &>(),
                 py::arg("constant_state"),
                 py::arg("frame_origin") = "SSB",
                 py::arg("frame_orientation") = "ECLIPJ2000");

    py::class_<tss::CustomEphemerisSettings,
            std::shared_ptr<tss::CustomEphemerisSettings>,
            tss::EphemerisSettings>(m, "CustomEphemerisSettings")
            .def(py::init<const std::function<Eigen::Vector6d(const double)>,
                 const std::string &,
                 const std::string &>(),
                 py::arg("custom_state_function"),
                 py::arg("frame_origin") = "SSB",
                 py::arg("frame_orientation") = "ECLIPJ2000")
            .def("get_custom_state_function",
                 &tss::CustomEphemerisSettings::getCustomStateFunction);

    py::class_<tss::KeplerEphemerisSettings,
            std::shared_ptr<tss::KeplerEphemerisSettings>,
            tss::EphemerisSettings>(m, "KeplerEphemerisSettings")
            .def(py::init<const Eigen::Vector6d &, const double, const double,
                 const std::string &, const std::string &, const double,
                 const double>(),
                 py::arg("initial_state_in_keplerian_elements"),
                 py::arg("epoch_of_initial_state"),
                 py::arg("central_body_gravitational_parameter"),
                 py::arg("reference_frame_origin") = "SSB",
                 py::arg("reference_frame_orientation") = "ECLIPJ2000",
                 py::arg("root_finder_absolute_tolerance") =
            200.0 * std::numeric_limits<double>::epsilon(),
                 py::arg("root_finder_maximum_number_of_iterations") = 1000.0)
            .def("get_initial_state_in_keplerian_elements",
                 &tss::KeplerEphemerisSettings::getInitialStateInKeplerianElements)
            .def("get_epoch_of_initial_state",
                 &tss::KeplerEphemerisSettings::getEpochOfInitialState)
            .def("get_central_body_gravitational_parameter",
                 &tss::KeplerEphemerisSettings::getCentralBodyGravitationalParameter)
            .def("get_root_finder_absolute_tolerance",
                 &tss::KeplerEphemerisSettings::getRootFinderAbsoluteTolerance)
            .def("get_root_finder_maximum_number_of_iterations",
                 &tss::KeplerEphemerisSettings::
                 getRootFinderMaximumNumberOfIterations);

    py::class_<tss::TabulatedEphemerisSettings,
            std::shared_ptr<tss::TabulatedEphemerisSettings>,
            tss::EphemerisSettings>(m, "TabulatedEphemerisSettings")
            .def(py::init<const std::map<double, Eigen::Vector6d> &, std::string,
                 std::string>())
            .def("get_body_state_history",
                 &tss::TabulatedEphemerisSettings::getBodyStateHistory)
            .def("get_use_long_double_states",
                 &tss::TabulatedEphemerisSettings::getUseLongDoubleStates)
            .def("set_use_long_double_states",
                 &tss::TabulatedEphemerisSettings::setUseLongDoubleStates);

    m.def("create_ephemeris", &tss::createBodyEphemeris,
          py::arg("ephemeris_settings"), py::arg("body_name"));


    m.def("keplerian",
          &tss::keplerEphemerisSettings,
          py::arg( "initial_keplerian_state" ),
          py::arg( "initial_state_epoch" ),
          py::arg( "central_body_gravitational_parameter" ),
          py::arg( "frame_origin" ) = "SSB" ,
          py::arg( "frame_orientation" ) = "ECLIPJ2000" ,
          py::arg( "root_finder_absolute_tolerance" ) = 200.0 * std::numeric_limits< double >::epsilon( ),
          py::arg( "root_finder_maximum_iterations" ) = 1000.0 );
}

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
            .def(py::init<const Eigen::Vector6d &>(), py::arg("state") = Eigen::Vector6d::Zero())
            .def("get_ephemeris_frame_to_base_frame", &tss::Body::getEphemerisFrameToBaseFrame)
            .def("set_ephemeris_frame_to_base_frame", &tss::Body::setEphemerisFrameToBaseFrame)
            .def_property("ephemeris_frame_to_base_frame", &tss::Body::getEphemerisFrameToBaseFrame, &tss::Body::setEphemerisFrameToBaseFrame)
            .def("get_state", &tss::Body::getState)
            .def("set_state", &tss::Body::setState)
            .def("get_state_in_based_frame_from_ephemeris", &tss::Body::getStateInBaseFrameFromEphemeris<double,double>)
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
            .def_property("aerodynamic_coefficient_interface", &tss::Body::getAerodynamicCoefficientInterface, &tss::Body::setAerodynamicCoefficientInterface)
            .def("get_body_mass", &tss::Body::getBodyMass)
            .def("set_constant_mass", &tss::Body::setConstantBodyMass)
            .def("get_radiation_pressure_interfaces", &tss::Body::getRadiationPressureInterfaces)
            .def("set_radiation_pressure_interface", &tss::Body::setRadiationPressureInterface)
            .def("set_aerodynamic_coefficient_interface", &tss::Body::setAerodynamicCoefficientInterface)
            .def("get_aerodynamic_coefficient_interface", &tss::Body::getAerodynamicCoefficientInterface)
            .def("get_flight_conditions", &tss::Body::getFlightConditions)
            .def("set_flight_conditions", &tss::Body::setFlightConditions, py::arg("aerodynamic_flight_conditions"))
            .def("get_rotation_model", &tss::Body::getRotationalEphemeris)
            .def("set_rotation_model", &tss::Body::setRotationalEphemeris, py::arg("rotational_ephemeris"))
            .def_property("rotation_model", &tss::Body::getRotationalEphemeris, &tss::Body::setRotationalEphemeris);


    py::class_<tss::SystemOfBodies,
            std::shared_ptr<tss::SystemOfBodies> >(m, "SystemOfBodies")
            .def(py::init<//ctor 1
                 const std::string,
                 const std::string,
                 const std::unordered_map< std::string, std::shared_ptr< tss::Body > > >(),
                 py::arg("frame_origin") = "SSB",
                 py::arg("frame_orientation") = "ECLIPJ2000",
                 py::arg("body_map") =
            std::unordered_map< std::string, std::shared_ptr< tss::Body > >( ) )
            .def("get_body", &tss::SystemOfBodies::getBody)
            .def("count", &tss::SystemOfBodies::count)
            .def("create_empty_body", &tss::SystemOfBodies::createEmptyBody,
                 py::arg("body_name"),
                 py::arg("process_body") = 1 )
            .def("add_body", &tss::SystemOfBodies::addBody,
                 py::arg("body_to_add"),
                 py::arg("body_name"),
                 py::arg("process_body") = 1 )
            .def("processBodyFrameDefinitions", &tss::SystemOfBodies::processBodyFrameDefinitions);


    m.def("get_body_gravitational_parameter",
          &tss::getBodyGravitationalParameter,
          py::arg("body_collection"), py::arg("body_name"));

    /////////////////////////////////////////////////////////////////////////////
    // createBodies.h ///////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    py::class_<tss::BodySettings, std::shared_ptr<tss::BodySettings>>(
                m, "BodySettings", tudatpy::body_settings_docstring().c_str())
            .def_readwrite("constant_mass", &tss::BodySettings::constantMass)
            .def_readwrite("atmosphere_settings", &tss::BodySettings::atmosphereSettings)
            .def_readwrite("ephemeris_settings", &tss::BodySettings::ephemerisSettings)
            .def_readwrite("gravity_field_settings", &tss::BodySettings::gravityFieldSettings)
            .def_readwrite("rotation_model_settings", &tss::BodySettings::rotationModelSettings)
            .def_readwrite("shape_model_settings", &tss::BodySettings::shapeModelSettings)
            .def_readwrite("radiation_pressure_settings", &tss::BodySettings::radiationPressureSettings)
            .def_readwrite("aerodynamic_coefficient_settings", &tss::BodySettings::aerodynamicCoefficientSettings)
            .def_readwrite("gravity_field_variation_settings", &tss::BodySettings::gravityFieldVariationSettings)
            .def_readwrite("ground_station_settings", &tss::BodySettings::groundStationSettings);

    py::class_<tss::BodyListSettings,
            std::shared_ptr<tss::BodyListSettings> >(m, "BodyListSettings")
            .def(py::init<//ctor 1
                 const std::string,
                 const std::string >(),
                 py::arg("frame_origin") ,
                 py::arg("frame_orientation"))
            .def(py::init<//ctor 2
                 const std::map< std::string, std::shared_ptr< tss::BodySettings > >&,
                 const std::string,
                 const std::string >(),
                 py::arg("body_settings"),
                 py::arg("frame_origin"),
                 py::arg("frame_orientation"))
            .def("at", &tss::BodyListSettings::at)
            .def("get", &tss::BodyListSettings::get)
            .def("add_settings", py::overload_cast<const std::string>(&tss::BodyListSettings::addSettings),
                 py::arg("body_name"))
            .def("add_settings", py::overload_cast<std::shared_ptr<tss::BodySettings>, const std::string>
                 (&tss::BodyListSettings::addSettings), py::arg("settings_to_add"), py::arg("body_name"))
            .def("get_frame_origin", &tss::BodyListSettings::getFrameOrigin)
            .def("get_frame_orientation", &tss::BodyListSettings::getFrameOrientation);

    /////////////////////////////////////////////////////////////////////////////
    // defaultBodies.h //////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    // getDefaultBodySettings (overload 1)
    m.def("get_default_body_settings",
          py::overload_cast<const std::vector<std::string> &, const std::string,
          const std::string >(
              &tss::getDefaultBodySettings),
          py::arg("bodies"),
          py::arg("base_frame_origin") = "SSB",
          py::arg("base_frame_orientation") = "ECLIPJ2000");

    // getDefaultBodySettings (overload 2)
    m.def("get_default_body_settings",
          py::overload_cast<const std::vector<std::string> &,
          const double, const double,
          const std::string, const std::string, const double >(
              &tss::getDefaultBodySettings),
          py::arg("bodies"),
          py::arg("initial_time"),
          py::arg("final_time"),
          py::arg("base_frame_origin") = "SSB",
          py::arg("base_frame_orientation") = "ECLIPJ2000",
          py::arg("time_step") = 300.0);


    m.def("get_global_frame_origin", &tss::getGlobalFrameOrigin);

    m.def("create_system_of_bodies", &tss::createSystemOfBodies);


    // Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.cpp
    m.def(
                "create_tabulated_ephemeris_from_spice",
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


    m.def("create_radiation_pressure_interface",
          &tss::createRadiationPressureInterface,
          py::arg("radiationPressureInterfaceSettings"), py::arg("body_name"),
          py::arg("body_dict"));


    m.def("set_aerodynamic_orientation_functions", &tss::setAerodynamicOrientationFunctions,
          py::arg("body"),
          py::arg("angle_of_attack_function") = std::function< double( ) >( ),
          py::arg("sideslip_angle_function") = std::function< double( ) >( ),
          py::arg("bank_angle_function") = std::function< double( ) >( ),
          py::arg("update_function") = std::function< void( const double ) >( ) );

    auto aerodynamic_coefficient_setup = m.def_submodule("aerodynamic_coefficients");
    expose_aerodynamic_coefficient_setup(aerodynamic_coefficient_setup);

    auto radiation_pressure_setup = m.def_submodule("radiation_pressure");
    expose_radiation_pressure_setup(radiation_pressure_setup);

    auto rotation_model_setup = m.def_submodule("rotation_model");
    expose_rotation_model_setup(rotation_model_setup);

    auto gravity_field_setup = m.def_submodule("gravity_field");
    expose_gravity_field_setup(gravity_field_setup);

    auto ephemeris_setup = m.def_submodule("ephemeris");
    expose_ephemeris_setup(ephemeris_setup);

    auto atmosphere_setup = m.def_submodule("atmosphere");
    expose_atmosphere_setup(atmosphere_setup);

}

}// namespace tudatpy
