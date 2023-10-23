/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_radiation_pressure_setup.h"

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

namespace tudatpy {
namespace numerical_simulation {
namespace environment_setup {
namespace radiation_pressure {

    void expose_radiation_pressure_setup(py::module &m) {

        /////////////////////////////////////////////////////////////////////////////
        // createRadiationPressureInterface.h
        /////////////////////////////////////////////////////////////////////////////
        py::enum_<tss::RadiationPressureType>(m, "RadiationPressureType",
                                              get_docstring("RadiationPressureType").c_str())
                .value("cannonball_radiation_pressure_interface",
                        tss::RadiationPressureType::cannon_ball_radiation_pressure_interface,
                        get_docstring("RadiationPressureType.cannonball_radiation_pressure_interface").c_str())
//                .value("panelled_radiation_pressure_interface",
//                       tss::RadiationPressureType::panelled_radiation_pressure_interface,
//                       get_docstring("RadiationPressureType.panelled_radiation_pressure_interface").c_str())
//                .value("solar_sailing_radiation_pressure_interface",
//                       tss::RadiationPressureType::solar_sailing_radiation_pressure_interface,
//                       get_docstring("RadiationPressureType.solar_sailing_radiation_pressure_interface").c_str())
                .export_values();


        py::class_<tss::RadiationPressureInterfaceSettings,
                std::shared_ptr<tss::RadiationPressureInterfaceSettings>>(
                m, "RadiationPressureInterfaceSettings",
                get_docstring("RadiationPressureInterfaceSettings").c_str());
//            .def(py::init<const tss::RadiationPressureType, const std::string &,
//                 const std::vector<std::string>>(),
//                 py::arg("radiation_pressure_type"), py::arg("source_body"),
//                 py::arg("occulting_bodies") = std::vector<std::string>());

        py::class_<tss::CannonBallRadiationPressureInterfaceSettings,
                std::shared_ptr<tss::CannonBallRadiationPressureInterfaceSettings>,
                tss::RadiationPressureInterfaceSettings>(
                m, "CannonBallRadiationPressureInterfaceSettings",
                get_docstring("CannonBallRadiationPressureInterfaceSettings").c_str());

//            .def(py::init<const std::string &, const double, const double,
//                 const std::vector<std::string> &>(),
//                 py::arg("source_body"), py::arg("area"),
//                 py::arg("radiation_pressure_coefficient"),
//                 py::arg("occulting_bodies") = std::vector<std::string>())

        m.def("cannonball",
              py::overload_cast<const std::string &, const double, const double, const std::vector<std::string> &>(
                      &tss::cannonBallRadiationPressureSettings),
              py::arg("source_body"),
              py::arg("reference_area"),
              py::arg("radiation_pressure_coefficient"),
              py::arg("occulting_bodies") = std::vector<std::string>(),
              get_docstring("cannonball").c_str()
              );

        m.def("panelled",
              &tss::panelledRadiationPressureInterfaceSettings,
              py::arg("source_body"),
              py::arg("emissivities"),
              py::arg("areas"),
              py::arg("diffusion_coefficients"),
              py::arg("surface_normals_in_body_fixed_frame"),
              py::arg("occulting_bodies") = std::vector<std::string>(),
              get_docstring("panelled").c_str()
              );


///////////////////////////////////////////////////////////
///////////   ENUMS
///////////////////////////////////////////////////////////

        enum class SecondDegreeZonalPeriodicSurfacePropertyDistributionModel
        {
            custom,
            albedo_knocke, /**< Knocke Earth albedo model: Knocke, Philip et al. "Earth radiation pressure effects on satellites." Astrodynamics Conference. American Institute of Aeronautics and Astronautics, 1988. */
            emissivity_knocke /**< Knocke Earth emissivity model: Knocke, Philip et al. "Earth radiation pressure effects on satellites." Astrodynamics Conference. American Institute of Aeronautics and Astronautics, 1988. */
        };

        py::enum_<tss::SecondDegreeZonalPeriodicSurfacePropertyDistributionModel>(m, "SecondDegreeZonalPeriodicSurfacePropertyDistributionModel")
            .value("custom", tss::SecondDegreeZonalPeriodicSurfacePropertyDistributionModel::custom)
            .value("albedo_knocke", tss::SecondDegreeZonalPeriodicSurfacePropertyDistributionModel::albedo_knocke)
            .value("emissivity_knocke", tss::SecondDegreeZonalPeriodicSurfacePropertyDistributionModel::emissivity_knocke)
            .export_values();

///////////////////////////////////////////////////////////
///////////   LUMINOSITY MODELS
///////////////////////////////////////////////////////////

        py::class_<tss::LuminosityModelSettings,
            std::shared_ptr<tss::LuminosityModelSettings>>(
            m, "LuminosityModelSettings",
            get_docstring("LuminosityModelSettings").c_str());

        m.def("constant_luminosity",
              &tss::constantLuminosityModelSettings,
              py::arg("luminosity"),
              get_docstring("constant_luminosity").c_str()
        );

        m.def("irradiance_based_luminosity",
              &tss::irradianceBasedLuminosityModelSettings,
              py::arg("irradiance_function_at_reference_distance"),
              py::arg("irradiance_based_luminosity"),
              get_docstring("constant_luminosity").c_str()
        );


///////////////////////////////////////////////////////////
///////////   SURFACE PROPERTY MODELS
///////////////////////////////////////////////////////////


        py::class_<tss::SurfacePropertyDistributionSettings,
            std::shared_ptr<tss::SurfacePropertyDistributionSettings>>(
            m, "SurfacePropertyDistributionSettings",
            get_docstring("SurfacePropertyDistributionSettings").c_str());


        m.def("constant_surface_property_distribution",
              &tss::constantSurfacePropertyDistributionSettings,
              py::arg("constant_value"),
              get_docstring("constant_surface_property_distribution").c_str()
        );

        m.def("spherical_harmonic_surface_property_distribution",
              py::overload_cast< const Eigen::MatrixXd&, const Eigen::MatrixXd& >(
                  &tss::sphericalHarmonicsSurfacePropertyDistributionSettings ),
              py::arg("cosine_coefficients"),
              py::arg("sine_coefficients"),
              get_docstring("constant_surface_property_distribution").c_str()
        );

        m.def("predefined_spherical_harmonic_surface_property_distribution",
              py::overload_cast< tss::SphericalHarmonicsSurfacePropertyDistributionModel >(
                  &tss::sphericalHarmonicsSurfacePropertyDistributionSettings ),
              py::arg("predefined_model"),
              get_docstring("predefined_spherical_harmonic_surface_property_distribution").c_str()
        );


        m.def("knocke_type_surface_property_distribution",
                  &tss::manualSecondDegreeZonalPeriodicSurfacePropertyDistributionSettings,
              py::arg("constant_contribution"),
              py::arg("constant_degree_one_contribution"),
              py::arg("cosine_periodic_degree_one_contribution"),
              py::arg("sine_periodic_degree_one_contribution"),
              py::arg("constant_degree_two_contribution"),
              py::arg("reference_epoch"),
              py::arg("period"),
              get_docstring("knocke_type_property_distribution").c_str()
        );


        m.def("predefined_knocke_type_surface_property_distribution",
              py::overload_cast< tss::SecondDegreeZonalPeriodicSurfacePropertyDistributionModel >(
                  &tss::secondDegreeZonalPeriodicSurfacePropertyDistributionSettings ),
              py::arg("predefined_model"),
              get_docstring("predefined_knocke_type_property_distribution").c_str()
        );


///////////////////////////////////////////////////////////
///////////   PANEL RADIOSITY MODELS
///////////////////////////////////////////////////////////


        py::class_<tss::PanelRadiosityModelSettings,
            std::shared_ptr<tss::PanelRadiosityModelSettings>>(
            m, "PanelRadiosityModelSettings",
            get_docstring("PanelRadiosityModelSettings").c_str());

        m.def("constant_panel_radiosity",
              &tss::constantPanelRadiosityModelSettings,
              py::arg("constant_radiosoty"),
              get_docstring("constant_panel_radiosity").c_str()
        );

        m.def("albedo_constant_surface_radiosity",
              py::overload_cast< double, const std::string& >( &tss::albedoPanelRadiosityModelSettings ),
              py::arg("constant_albedo"),
              py::arg("original_source_name"),
              get_docstring("constant_albedo_panel_radiosity").c_str()
        );

        m.def("albedo_variable_surface_radiosity",
              &tss::albedoPanelRadiosityModelSettingsGeneric,
              py::arg("albedo_distribution_model"),
              py::arg("original_source_name"),
              get_docstring("constant_albedo_panel_radiosity").c_str()
        );

        m.def("thermal_emission_blackbody_constant_emissivity",
              py::overload_cast< double, const std::string& >( &tss::delayedThermalPanelRadiosityModelSettings ),
              py::arg("constant_emissivity"),
              py::arg("original_source_name"),
              get_docstring("thermal_emission_delayed_radiosity_constant_emissivity").c_str()
        );

        m.def("thermal_emission_blackbody_variable_emissivity",
              &tss::delayedThermalPanelRadiosityModelSettingsGeneric,
              py::arg("emissivity_distribution_model"),
              py::arg("original_source_name"),
              get_docstring("thermal_emission_delayed_radiosity_constant_emissivity").c_str()
        );

        m.def("thermal_emission_angle_based_radiosity",
              &tss::angleBasedThermalPanelRadiosityModelSettings,
              py::arg("minimum_temperature"),
              py::arg("maximum_temperature"),
              py::arg("constant_emissivity"),
              py::arg("original_source_name"),
              get_docstring("thermal_emission_delayed_radiosity_constant_emissivity").c_str()
        );


///////////////////////////////////////////////////////////
///////////   PANEL REFLECTION MODELS
///////////////////////////////////////////////////////////

        py::class_<tss::BodyPanelReflectionLawSettings,
            std::shared_ptr<tss::BodyPanelReflectionLawSettings>>(
            m, "BodyPanelReflectionLawSettings",
            get_docstring("BodyPanelReflectionLawSettings").c_str());


        m.def("specular_diffuse_body_panel_reflection",
              &tss::specularDiffuseBodyPanelReflectionLawSettings,
              py::arg("specular_reflectivity"),
              py::arg("diffuse_reflectivity"),
              py::arg("with_instantaneous_reradiation"),
            get_docstring("specular_diffuse_body_panel_reflection").c_str()
        );


///////////////////////////////////////////////////////////
///////////   RADIATION SOURCE MODELS
///////////////////////////////////////////////////////////


        py::class_<tss::RadiationSourceModelSettings,
            std::shared_ptr<tss::RadiationSourceModelSettings>>(
            m, "RadiationSourceModelSettings",
            get_docstring("RadiationSourceModelSettings").c_str());



        m.def("isotropic_radiation_source",
              &tss::isotropicPointRadiationSourceModelSettings,
              py::arg("luminosity_model"),
              get_docstring("isotropic_radiation_source").c_str()
              );

        m.def("panelled_extended_radiation_source",
              &tss::extendedRadiationSourceModelSettingsWithOccultationMap,
              py::arg("panel_radiosity_settings"),
              py::arg("number_of_panels_per_ring"),
              py::arg("original_source_occulting_bodies" ) = std::map<std::string, std::vector<std::string>>(),
              get_docstring("panelled_extended_radiation_source").c_str()
              );


///////////////////////////////////////////////////////////
///////////   RADIATION TARGET MODELS
///////////////////////////////////////////////////////////

        py::class_<tss::RadiationPressureTargetModelSettings,
            std::shared_ptr<tss::RadiationPressureTargetModelSettings>>(
            m, "RadiationPressureTargetModelSettings",
            get_docstring("RadiationPressureTargetModelSettings").c_str());


        m.def("cannonball_radiation_target",
                  &tss::cannonballRadiationPressureTargetModelSettingsWithOccultationMap,
              py::arg("reference_area"),
              py::arg("radiation_pressure_coefficient"),
              py::arg("per_source_occulting_bodies") = std::map<std::string, std::vector<std::string>>( ),
              get_docstring("cannonball_radiation_target").c_str()
        );

        m.def("panelled_radiation_target",
              &tss::paneledRadiationPressureTargetModelSettingsWithOccultationMap,
              get_docstring("panelled_radiation_target").c_str()
        );

    }

}// namespace radiation_pressure
}// namespace environment_setup
}// namespace numerical_simulation
}// namespace tudatpy
