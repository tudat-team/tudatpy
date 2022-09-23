/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_gravity_field_setup.h"

#include "tudatpy/docstrings.h"
#include <tudat/simulation/environment_setup.h>
#include <tudat/astro/reference_frames/referenceFrameTransformations.h>

//#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;

namespace tudat
{

namespace simulation_setup
{

inline std::shared_ptr< GravityFieldSettings > fromFileSphericalHarmonicsGravityFieldSettings(
        const std::string& filePath,
        const std::string& associatedReferenceFrame,
        const int maximumDegree,
        const int maximumOrder,
        const int gravitationalParameterIndex,
        const int referenceRadiusIndex )
{

    return std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >(
                filePath, associatedReferenceFrame, maximumDegree,  maximumOrder,
                gravitationalParameterIndex, referenceRadiusIndex );
}


inline std::shared_ptr< GravityFieldSettings > predefinedSphericalHarmonic(
        const SphericalHarmonicsModel sphericalHarmonicsModel,
        const int maximumDegree = -1 )
{
    return std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( sphericalHarmonicsModel, maximumDegree );
}


}

}

namespace tudatpy {
namespace numerical_simulation {
namespace environment_setup {
namespace gravity_field {

void expose_gravity_field_setup(py::module &m) {
    /////////////////////////////////////////////////////////////////////////////
    // createGravityField.h
    /////////////////////////////////////////////////////////////////////////////
    py::enum_<tss::GravityFieldType>(m, "GravityFieldType",
                                     get_docstring("GravityFieldType").c_str())
            .value("central_gravity", tss::GravityFieldType::central, get_docstring("GravityFieldType.central_gravity").c_str())
            .value("central_spice_gravity", tss::GravityFieldType::central_spice, get_docstring("GravityFieldType.central_spice_gravity").c_str())
            .value("spherical_harmonic_gravity", tss::GravityFieldType::spherical_harmonic, get_docstring("GravityFieldType.spherical_harmonic_gravity").c_str())
            .export_values();

    py::enum_<tss::SphericalHarmonicsModel>(m, "PredefinedSphericalHarmonicsModel",
                                            get_docstring("PredefinedSphericalHarmonicsModel").c_str())
//            .value("custom_model", tss::SphericalHarmonicsModel::customModel, get_docstring("SphericalHarmonicsModel.custom_model").c_str())
            .value("egm96", tss::SphericalHarmonicsModel::egm96, get_docstring("PredefinedSphericalHarmonicsModel.egm96").c_str())
            .value("ggm02c", tss::SphericalHarmonicsModel::ggm02c, get_docstring("PredefinedSphericalHarmonicsModel.ggm02c").c_str())
            .value("ggm02s", tss::SphericalHarmonicsModel::ggm02s, get_docstring("PredefinedSphericalHarmonicsModel.ggm02s").c_str())
            .value("glgm3150", tss::SphericalHarmonicsModel::glgm3150, get_docstring("PredefinedSphericalHarmonicsModel.glgm3150").c_str())
            .value("lpe200", tss::SphericalHarmonicsModel::lpe200, get_docstring("PredefinedSphericalHarmonicsModel.lpe200").c_str())
            .value("jgmro120d", tss::SphericalHarmonicsModel::jgmro120d, get_docstring("PredefinedSphericalHarmonicsModel.jgmro120d").c_str())
            .value("jgmess160a", tss::SphericalHarmonicsModel::jgmess160a, get_docstring("PredefinedSphericalHarmonicsModel.jgmess160a").c_str())
            .value("shgj180u", tss::SphericalHarmonicsModel::shgj180u, get_docstring("PredefinedSphericalHarmonicsModel.shgj180u").c_str())
            .export_values();

    py::class_<tss::GravityFieldSettings, std::shared_ptr<tss::GravityFieldSettings>>(
                m, "GravityFieldSettings",
                get_docstring("GravityFieldSettings").c_str())
            //            .def(py::init<const tss::GravityFieldType>(),
            //                 py::arg("gravity_field_type"))
            .def_property_readonly("gravity_field_type", &tss::GravityFieldSettings::getGravityFieldType,
                                   get_docstring("GravityFieldSettings.gravity_field_type").c_str());


    py::class_<tss::CentralGravityFieldSettings, std::shared_ptr<tss::CentralGravityFieldSettings>,
            tss::GravityFieldSettings>(m, "CentralGravityFieldSettings",
                                       get_docstring("CentralGravityFieldSettings").c_str())
            //            .def(py::init<double>(), py::arg("gravitational_parameter") )
            .def_property("gravitational_parameter", &tss::CentralGravityFieldSettings::getGravitationalParameter,
                          &tss::CentralGravityFieldSettings::resetGravitationalParameter,
                          get_docstring("CentralGravityFieldSettings.gravitational_parameter").c_str());


    py::class_<tss::SphericalHarmonicsGravityFieldSettings, std::shared_ptr<tss::SphericalHarmonicsGravityFieldSettings>,
            tss::GravityFieldSettings>(m, "SphericalHarmonicsGravityFieldSettings",
                                       get_docstring("SphericalHarmonicsGravityFieldSettings").c_str())
            //            .def(py::init<const double, const double, const Eigen::MatrixXd, const Eigen::MatrixXd, const std::string&>(),
            //                 py::arg("gravitational_parameter"), py::arg("reference_radius"), py::arg("cosine_coefficients"),
            //                 py::arg("sine_coefficients"), py::arg("associated_reference_frame"))
            .def_property("gravitational_parameter",
                          &tss::SphericalHarmonicsGravityFieldSettings::getGravitationalParameter,
                          &tss::SphericalHarmonicsGravityFieldSettings::resetGravitationalParameter,
                          get_docstring("SphericalHarmonicsGravityFieldSettings.gravitational_parameter").c_str())
            .def_property("normalized_cosine_coefficients",
                          &tss::SphericalHarmonicsGravityFieldSettings::getCosineCoefficients,
                          &tss::SphericalHarmonicsGravityFieldSettings::resetCosineCoefficients,
                          get_docstring("SphericalHarmonicsGravityFieldSettings.normalized_cosine_coefficients").c_str())
            .def_property("normalized_sine_coefficients",
                          &tss::SphericalHarmonicsGravityFieldSettings::getSineCoefficients,
                          &tss::SphericalHarmonicsGravityFieldSettings::resetSineCoefficients,
                          get_docstring("SphericalHarmonicsGravityFieldSettings.normalized_sine_coefficients").c_str())
            .def_property("associated_reference_frame",
                          &tss::SphericalHarmonicsGravityFieldSettings::getAssociatedReferenceFrame,
                          &tss::SphericalHarmonicsGravityFieldSettings::resetAssociatedReferenceFrame,
                          get_docstring("SphericalHarmonicsGravityFieldSettings.associated_reference_frame").c_str())
            .def_property("create_time_dependent_field",
                          &tss::SphericalHarmonicsGravityFieldSettings::getCreateTimeDependentField,
                          &tss::SphericalHarmonicsGravityFieldSettings::setCreateTimeDependentField,
                          get_docstring("SphericalHarmonicsGravityFieldSettings.create_time_dependent_field").c_str())
            .def_property_readonly("reference_radius",
                                   &tss::SphericalHarmonicsGravityFieldSettings::getReferenceRadius,
                                   get_docstring("SphericalHarmonicsGravityFieldSettings.reference_radius").c_str());


    py::class_<tss::FromFileSphericalHarmonicsGravityFieldSettings, std::shared_ptr<tss::FromFileSphericalHarmonicsGravityFieldSettings>,
            tss::SphericalHarmonicsGravityFieldSettings>(m, "FromFileSphericalHarmonicsGravityFieldSettings",
                                                         get_docstring("FromFileSphericalHarmonicsGravityFieldSettings").c_str());


    m.def("central",
          &tss::centralGravitySettings,
          py::arg("gravitational_parameter"),
          get_docstring("central").c_str()
          );

    m.def("central_spice",
          &tss::centralGravityFromSpiceSettings,
          get_docstring("central_spice").c_str()
          );

    m.def("spherical_harmonic",
          py::overload_cast< const double,
          const double,
          const Eigen::MatrixXd,
          const Eigen::MatrixXd,
          const std::string& >( &tss::sphericalHarmonicsGravitySettings ),
          py::arg("gravitational_parameter"),
          py::arg("reference_radius"),
          py::arg("normalized_cosine_coefficients"),
          py::arg("normalized_sine_coefficients"),
          py::arg("associated_reference_frame"),
          get_docstring("spherical_harmonic").c_str()
          );

    m.def("spherical_harmonic_triaxial_body",
          &tss::createHomogeneousTriAxialEllipsoidGravitySettings,
          py::arg("axis_a"),
          py::arg("axis_b"),
          py::arg("axis_c"),
          py::arg("density"),
          py::arg("maximum_degree"),
          py::arg("maximum_order"),
          py::arg("associated_reference_frame"),
          get_docstring("spherical_harmonic_triaxial_body").c_str()
          );

    m.def("predefined_spherical_harmonic",
          tss::fromFileSphericalHarmonicsGravityFieldSettings,
          py::arg("file"),
          py::arg("associated_reference_frame"),
          py::arg("maximum_degree"),
          py::arg("maximum_order"),
          py::arg("gravitational_parameter_index"),
          py::arg("reference_radius_index"),
          get_docstring("predefined_spherical_harmonic").c_str()
          );

    m.def("from_file_spherical_harmonic",
          tss::predefinedSphericalHarmonic,
          py::arg("predefined_model"),
          py::arg("maximum_degree") = -1,
          get_docstring("from_file_spherical_harmonic").c_str()
          );
}

}// namespace gravity_field
}// namespace environment_setup
}// namespace numerical_simulation
}// namespace tudatpy
