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

namespace tudatpy {
    namespace numerical_simulation {
        namespace environment_setup {
            namespace radiation_pressure {

                PYBIND11_MODULE(expose_radiation_pressure, m) {
                    /////////////////////////////////////////////////////////////////////////////
                    // createRadiationPressureInterface.h
                    /////////////////////////////////////////////////////////////////////////////
                    py::enum_<tss::RadiationPressureType>(
                        m, "RadiationPressureType",
R"doc(Enumeration of available radiation pressure types.


	:member cannonball_radiation_pressure_interface:
	:member panelled_radiation_pressure_interface:
	:member solar_sailing_radiation_pressure_interface:
)doc")
                        .value("cannonball_radiation_pressure_interface",
                               tss::RadiationPressureType::
                                   cannon_ball_radiation_pressure_interface,
"")
                        //                .value("panelled_radiation_pressure_interface",
                        //                       tss::RadiationPressureType::panelled_radiation_pressure_interface,
                        //                .value("solar_sailing_radiation_pressure_interface",
                        //                       tss::RadiationPressureType::solar_sailing_radiation_pressure_interface,
                        .export_values();
                    py::enum_<tss::RadiationPressureTargetModelType>(
                        m, "RadiationPressureTargetModelType",
"")
                        .value("cannonball_target",
                               tss::RadiationPressureTargetModelType::
                                   cannonball_target,
"")
                        .value("paneled_target",
                               tss::RadiationPressureTargetModelType::
                                   paneled_target,
"")
                        .value("multi_type_target",
                               tss::RadiationPressureTargetModelType::
                                   multi_type_target,
"")
                        .value("undefined_target",
                               tss::RadiationPressureTargetModelType::
                                   undefined_target,
"")
                        .export_values();


                    py::class_<tss::RadiationPressureInterfaceSettings,
                               std::shared_ptr<
                                   tss::RadiationPressureInterfaceSettings>>(
                        m, "RadiationPressureInterfaceSettings",
R"doc(Base class for providing settings for radiation pressure interface models.

	Functional (base) class for settings of radiation pressure interface models that require no information in addition to their type.
	Radiation pressure interface model settings requiring additional information must be defined using an object derived from this class.

)doc");
                    //            .def(py::init<const
                    //            tss::RadiationPressureType, const std::string
                    //            &,
                    //                 const std::vector<std::string>>(),
                    //                 py::arg("radiation_pressure_type"),
                    //                 py::arg("source_body"),
                    //                 py::arg("occulting_bodies") =
                    //                 std::vector<std::string>());

                    py::class_<
                        tss::CannonBallRadiationPressureInterfaceSettings,
                        std::shared_ptr<
                            tss::CannonBallRadiationPressureInterfaceSettings>,
                        tss::RadiationPressureInterfaceSettings>(
                        m, "CannonBallRadiationPressureInterfaceSettings",
R"doc(Class for defining model settings of a cannonball radiation pressure interface.

	`RadiationPressureInterfaceSettings` derived class for cannonball radiation pressure interface model settings.
)doc");

                    //            .def(py::init<const std::string &, const
                    //            double, const double,
                    //                 const std::vector<std::string> &>(),
                    //                 py::arg("source_body"), py::arg("area"),
                    //                 py::arg("radiation_pressure_coefficient"),
                    //                 py::arg("occulting_bodies") =
                    //                 std::vector<std::string>())

                    m.def("cannonball",
                          py::overload_cast<const std::string&, const double,
                                            const double,
                                            const std::vector<std::string>&>(
                              &tss::cannonBallRadiationPressureSettings),
                          py::arg("source_body"), py::arg("reference_area"),
                          py::arg("radiation_pressure_coefficient"),
                          py::arg("occulting_bodies") =
                              std::vector<std::string>(),
R"doc(Factory function for creating cannonball radiation pressure interface model settings.

	Factory function for settings object, defining a cannonball radiation pressure interface model,
	In this model the effective force is co-linear with the vector from radiation source to the body experiencing the force.


	:param source_body:
		Name of body emitting the radiation.
	:param reference_area:
		Surface area that undergoes radiation pressure.
	:param radiation_pressure_coefficient:
		Radiation pressure coefficient.
	:param occulting_bodies:
		List of bodies causing (partial) occultation.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.radiation_pressure.RadiationPressureInterfaceSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.radiation_pressure.CannonBallRadiationPressureInterfaceSettings` class
)doc");

                    //        m.def("panelled",
                    //              &tss::panelledRadiationPressureInterfaceSettings,
                    //              py::arg("source_body"),
                    //              py::arg("emissivities"),
                    //              py::arg("areas"),
                    //              py::arg("diffusion_coefficients"),
                    //              py::arg("surface_normals_in_body_fixed_frame"),
                    //              py::arg("occulting_bodies") =
                    //              std::vector<std::string>(),
                    //              );


                    ///////////////////////////////////////////////////////////
                    ///////////   ENUMS
                    ///////////////////////////////////////////////////////////

                    py::enum_<tss::KnockeTypeSurfacePropertyDistributionModel>(
                        m, "KnockeTypeSurfacePropertyDistributionModel")
                        .value("custom",
                               tss::KnockeTypeSurfacePropertyDistributionModel::
                                   custom)
                        .value("albedo_knocke",
                               tss::KnockeTypeSurfacePropertyDistributionModel::
                                   albedo_knocke)
                        .value("emissivity_knocke",
                               tss::KnockeTypeSurfacePropertyDistributionModel::
                                   emissivity_knocke)
                        .export_values();

                    py::enum_<
                        tss::
                            SphericalHarmonicsSurfacePropertyDistributionModel>(
                        m, "SphericalHarmonicsSurfacePropertyDistributionModel")
                        .value(
                            "albedo_dlam1",
                            tss::
                                SphericalHarmonicsSurfacePropertyDistributionModel::
                                    albedo_dlam1)
                        .export_values();

                    enum class
                        SphericalHarmonicsSurfacePropertyDistributionModel {
                            custom,
                            albedo_dlam1 /**< DLAM-1 lunar albedo model:
                                            Floberghagen, R. et al. "Lunar
                                            Albedo Force Modeling and its Effect
                                            on Low Lunar Orbit and Gravity Field
                                            Determination". ASR 23. 4(1999):
                                            733-738. */
                        };

                    ///////////////////////////////////////////////////////////
                    ///////////   LUMINOSITY MODELS
                    ///////////////////////////////////////////////////////////

                    py::class_<tss::LuminosityModelSettings,
                               std::shared_ptr<tss::LuminosityModelSettings>>(
                        m, "LuminosityModelSettings",
"");

                    m.def("constant_luminosity",
                          &tss::constantLuminosityModelSettings,
                          py::arg("luminosity"),
"");

                    m.def("irradiance_based_constant_luminosity",
                          &tss::irradianceBasedLuminosityModelSettings,
                          py::arg("constant_irradiance"),
                          py::arg("reference_distance"),
"");

                    m.def("time_variable_luminosity",
                          &tss::timeVariableLuminosityModelSettings,
                          py::arg("luminosity_function"),
"");

                    m.def(
                        "irradiance_based_time_variable_luminosity",
                        &tss::
                            timeVariableIrradianceBasedLuminosityModelSettings,
                        py::arg("irradiance_function"),
                        py::arg("reference_distance"),
"");


                    ///////////////////////////////////////////////////////////
                    ///////////   SURFACE PROPERTY MODELS
                    ///////////////////////////////////////////////////////////


                    py::class_<tss::SurfacePropertyDistributionSettings,
                               std::shared_ptr<
                                   tss::SurfacePropertyDistributionSettings>>(
                        m, "SurfacePropertyDistributionSettings",
"");


                    m.def(
                        "constant_surface_property_distribution",
                        &tss::constantSurfacePropertyDistributionSettings,
                        py::arg("constant_value"),
"");

                    m.def(
                        "spherical_harmonic_surface_property_distribution",
                        py::overload_cast<const Eigen::MatrixXd&,
                                          const Eigen::MatrixXd&>(
                            &tss::
                                sphericalHarmonicsSurfacePropertyDistributionSettings),
                        py::arg("cosine_coefficients"),
                        py::arg("sine_coefficients"),
"");

                    m.def(
                        "predefined_spherical_harmonic_surface_property_"
                        "distribution",
                        py::overload_cast<
                            tss::
                                SphericalHarmonicsSurfacePropertyDistributionModel>(
                            &tss::
                                sphericalHarmonicsSurfacePropertyDistributionSettings),
                        py::arg("predefined_model"),
"");


                    m.def(
                        "knocke_type_surface_property_distribution",
                        &tss::
                            manualSecondDegreeZonalPeriodicSurfacePropertyDistributionSettings,
                        py::arg("constant_contribution"),
                        py::arg("constant_degree_one_contribution"),
                        py::arg("cosine_periodic_degree_one_contribution"),
                        py::arg("sine_periodic_degree_one_contribution"),
                        py::arg("constant_degree_two_contribution"),
                        py::arg("reference_epoch"), py::arg("period"),
"");


                    m.def(
                        "predefined_knocke_type_surface_property_distribution",
                        py::overload_cast<
                            tss::KnockeTypeSurfacePropertyDistributionModel>(
                            &tss::
                                secondDegreeZonalPeriodicSurfacePropertyDistributionSettings),
                        py::arg("predefined_model"),
"");

                    m.def("custom_surface_property_distribution",
                          &tss::customSurfacePropertyDistributionSettings,
                          py::arg("custom_function"),
"");


                    ///////////////////////////////////////////////////////////
                    ///////////   PANEL RADIOSITY MODELS
                    ///////////////////////////////////////////////////////////


                    py::class_<
                        tss::PanelRadiosityModelSettings,
                        std::shared_ptr<tss::PanelRadiosityModelSettings>>(
                        m, "PanelRadiosityModelSettings",
"");

                    m.def("constant_radiosity",
                          &tss::constantPanelRadiosityModelSettings,
                          py::arg("radiosity"),
"");

                    m.def("constant_albedo_surface_radiosity",
                          py::overload_cast<double, const std::string&>(
                              &tss::albedoPanelRadiosityModelSettings),
                          py::arg("constant_albedo"),
                          py::arg("original_source_name"),
"");

                    m.def("variable_albedo_surface_radiosity",
                          &tss::albedoPanelRadiosityModelSettingsGeneric,
                          py::arg("albedo_distribution_settings"),
                          py::arg("original_source_name"),
"");

                    m.def("thermal_emission_blackbody_constant_emissivity",
                          py::overload_cast<double, const std::string&>(
                              &tss::delayedThermalPanelRadiosityModelSettings),
                          py::arg("constant_emissivity"),
                          py::arg("original_source_name"),
"");

                    m.def(
                        "thermal_emission_blackbody_variable_emissivity",
                        &tss::delayedThermalPanelRadiosityModelSettingsGeneric,
                        py::arg("emissivity_distribution_model"),
                        py::arg("original_source_name"),
"");

                    m.def(
                        "thermal_emission_angle_based_radiosity",
                        &tss::angleBasedThermalPanelRadiosityModelSettings,
                        py::arg("minimum_temperature"),
                        py::arg("maximum_temperature"),
                        py::arg("constant_emissivity"),
                        py::arg("original_source_name"),
"");


                    ///////////////////////////////////////////////////////////
                    ///////////   PANEL REFLECTION MODELS
                    ///////////////////////////////////////////////////////////

                    py::class_<
                        tss::BodyPanelReflectionLawSettings,
                        std::shared_ptr<tss::BodyPanelReflectionLawSettings>>(
                        m, "BodyPanelReflectionLawSettings",
"");


                    m.def(
                        "specular_diffuse_body_panel_reflection",
                        &tss::specularDiffuseBodyPanelReflectionLawSettings,
                        py::arg("specular_reflectivity"),
                        py::arg("diffuse_reflectivity"),
                        py::arg("with_instantaneous_reradiation"),
"");

                    m.def("lambertian_body_panel_reflection",
                          &tss::lambertainBodyPanelReflectionLawSettings,
                          py::arg("reflectivity"),
"");

                    ///////////////////////////////////////////////////////////
                    ///////////   RADIATION SOURCE MODELS
                    ///////////////////////////////////////////////////////////


                    py::class_<
                        tss::RadiationSourceModelSettings,
                        std::shared_ptr<tss::RadiationSourceModelSettings>>(
                        m, "RadiationSourceModelSettings",
"");


                    m.def("isotropic_radiation_source",
                          &tss::isotropicPointRadiationSourceModelSettings,
                          py::arg("luminosity_model"),
"");

                    m.def(
                        "panelled_extended_radiation_source",
                        &tss::
                            extendedRadiationSourceModelSettingsWithOccultationMap,
                        py::arg("panel_radiosity_settings"),
                        py::arg("number_of_panels_per_ring"),
                        py::arg("original_source_occulting_bodies") =
                            std::map<std::string, std::vector<std::string>>(),
"");


                    ///////////////////////////////////////////////////////////
                    ///////////   RADIATION TARGET MODELS
                    ///////////////////////////////////////////////////////////

                    py::class_<tss::RadiationPressureTargetModelSettings,
                               std::shared_ptr<
                                   tss::RadiationPressureTargetModelSettings>>(
                        m, "RadiationPressureTargetModelSettings",
"");


                    m.def(
                        "cannonball_radiation_target",
                        &tss::
                            cannonballRadiationPressureTargetModelSettingsWithOccultationMap,
                        py::arg("reference_area"),
                        py::arg("radiation_pressure_coefficient"),
                        py::arg("per_source_occulting_bodies") =
                            std::map<std::string, std::vector<std::string>>(),
"");

                    m.def(
                        "panelled_radiation_target",
                        &tss::
                            paneledRadiationPressureTargetModelSettingsWithOccultationMap,
                        py::arg("source_to_target_occulting_bodies") =
                            std::map<std::string, std::vector<std::string>>(),
"");
                }

            }  // namespace radiation_pressure
        }  // namespace environment_setup
    }  // namespace numerical_simulation
}  // namespace tudatpy
