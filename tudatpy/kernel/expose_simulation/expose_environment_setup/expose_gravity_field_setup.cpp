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
//#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;

namespace tudatpy {
namespace simulation {
namespace environment_setup {


    void expose_gravity_field_setup(py::module &m) {
        /////////////////////////////////////////////////////////////////////////////
        // createGravityField.h
        /////////////////////////////////////////////////////////////////////////////
        py::enum_<tss::GravityFieldType>(m, "GravityFieldType", "<no doc>")
                .value("central_gravity", tss::GravityFieldType::central)
                .value("central_spice_gravity", tss::GravityFieldType::central_spice)
                .value("spherical_harmonic_gravity", tss::GravityFieldType::spherical_harmonic)
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
//            .def(py::init<const tss::GravityFieldType>(),
//                 py::arg("gravity_field_type"))
                .def_property_readonly("gravity_field_type", &tss::GravityFieldSettings::getGravityFieldType);

        py::class_<tss::CentralGravityFieldSettings, std::shared_ptr<tss::CentralGravityFieldSettings>,
                tss::GravityFieldSettings>(m, "CentralGravityFieldSettings", "<no doc>")
//            .def(py::init<double>(), py::arg("gravitational_parameter") )
                .def_property("gravitational_parameter", &tss::CentralGravityFieldSettings::getGravitationalParameter,
                              &tss::CentralGravityFieldSettings::resetGravitationalParameter);


        py::class_<tss::SphericalHarmonicsGravityFieldSettings, std::shared_ptr<tss::SphericalHarmonicsGravityFieldSettings>,
                tss::GravityFieldSettings>(m, "SphericalHarmonicsGravityFieldSettings", "<no doc>")
//            .def(py::init<const double, const double, const Eigen::MatrixXd, const Eigen::MatrixXd, const std::string&>(),
//                 py::arg("gravitational_parameter"), py::arg("reference_radius"), py::arg("cosine_coefficients"),
//                 py::arg("sine_coefficients"), py::arg("associated_reference_frame"))
                .def_property("gravitational_parameter",
                              &tss::SphericalHarmonicsGravityFieldSettings::getGravitationalParameter,
                              &tss::SphericalHarmonicsGravityFieldSettings::resetGravitationalParameter)
                .def_property("normalized_cosine_coefficients",
                              &tss::SphericalHarmonicsGravityFieldSettings::getCosineCoefficients,
                              &tss::SphericalHarmonicsGravityFieldSettings::resetCosineCoefficients)
                .def_property("normalized_sine_coefficients",
                              &tss::SphericalHarmonicsGravityFieldSettings::getSineCoefficients,
                              &tss::SphericalHarmonicsGravityFieldSettings::resetSineCoefficients)
                .def_property("associated_reference_frame",
                              &tss::SphericalHarmonicsGravityFieldSettings::getAssociatedReferenceFrame,
                              &tss::SphericalHarmonicsGravityFieldSettings::resetAssociatedReferenceFrame)
                .def_property("create_time_dependent_field",
                              &tss::SphericalHarmonicsGravityFieldSettings::getCreateTimeDependentField,
                              &tss::SphericalHarmonicsGravityFieldSettings::setCreateTimeDependentField)
                .def_property_readonly("reference_radius",
                                       &tss::SphericalHarmonicsGravityFieldSettings::getReferenceRadius);

        py::class_<tss::FromFileSphericalHarmonicsGravityFieldSettings, std::shared_ptr<tss::FromFileSphericalHarmonicsGravityFieldSettings>,
                tss::SphericalHarmonicsGravityFieldSettings>(m, "FromFileSphericalHarmonicsGravityFieldSettings",
                                                             "<no doc>");


        m.def("central",
              &tss::centralGravitySettings,
              py::arg("gravitational_parameter"),
              get_docstring("central").c_str());

        m.def("central_spice",
              &tss::centralGravityFromSpiceSettings),
                get_docstring("central_spice").c_str();

        m.def("spherical_harmonic",
              &tss::sphericalHarmonicsGravitySettings,
              py::arg("gravitational_parameter"),
              py::arg("reference_radius"),
              py::arg("normalized_cosine_coefficients"),
              py::arg("normalized_sine_coefficients"),
              py::arg("associated_reference_frame"),
              get_docstring("spherical_harmonic").c_str());

        m.def("spherical_harmonic_triaxial_body",
              &tss::createHomogeneousTriAxialEllipsoidGravitySettings,
              py::arg("axis_a"),
              py::arg("axis_b"),
              py::arg("axis_c"),
              py::arg("density"),
              py::arg("maximum_degree"),
              py::arg("maximum_order"),
              py::arg("associated_reference_frame"),
              get_docstring("spherical_harmonic_triaxial_body").c_str());
    }


}// namespace environment_setup
}// namespace simulation
}// namespace tudatpy