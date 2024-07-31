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

#include "docstrings.h"
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
namespace tpc = tudat::physical_constants;

namespace tudat
{

namespace simulation_setup
{

inline std::shared_ptr< GravityFieldSettings > fromFileSphericalHarmonicsGravityFieldSettings(
        const std::string& filePath,
        const int maximumDegree,
        const int maximumOrder,
        const std::string& associatedReferenceFrame = "",
        const int gravitationalParameterIndex = 0,
        const int referenceRadiusIndex = 1 )
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

inline std::shared_ptr< GravityFieldSettings > createHomogeneousTriAxialEllipsoidGravitySettingsDeprecated(
        const double axisA, const double axisB, const double axisC, const double ellipsoidDensity,
        const int maximumDegree, const int maximumOrder,
        const std::string& associatedReferenceFrame,
        const double gravitationalConstant = tpc::GRAVITATIONAL_CONSTANT )
{
    static bool isWarningPrinted = false;
    if( isWarningPrinted == false )
    {
        tudat::utilities::printDeprecationWarning( "tudatpy.numerical_simulation.environment_setup.gravity_field.spherical_harmonic_triaxial_body",
                             "tudatpy.numerical_simulation.environment_setup.gravity_field.sh_triaxial_ellipsoid_from_density");
        isWarningPrinted = true;
    }

    return createHomogeneousTriAxialEllipsoidGravitySettings(
                axisA, axisB, axisC, ellipsoidDensity, maximumDegree, maximumOrder, associatedReferenceFrame, gravitationalConstant );

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
            .value("polyhedron_gravity", tss::GravityFieldType::polyhedron, get_docstring("GravityFieldType.polyhedron_gravity").c_str())
            .value("ring_gravity", tss::GravityFieldType::one_dimensional_ring, get_docstring("GravityFieldType.ring_gravity").c_str())
        .export_values();


    py::enum_<tss::SphericalHarmonicsModel>(m, "PredefinedSphericalHarmonicsModel",
                                            get_docstring("PredefinedSphericalHarmonicsModel").c_str())
            .value("egm96", tss::SphericalHarmonicsModel::egm96, get_docstring("PredefinedSphericalHarmonicsModel.egm96").c_str())
            .value("ggm02c", tss::SphericalHarmonicsModel::ggm02c, get_docstring("PredefinedSphericalHarmonicsModel.ggm02c").c_str())
            .value("ggm02s", tss::SphericalHarmonicsModel::ggm02s, get_docstring("PredefinedSphericalHarmonicsModel.ggm02s").c_str())
            .value("goco05c", tss::SphericalHarmonicsModel::ggm02s, get_docstring("PredefinedSphericalHarmonicsModel.goco05c").c_str())
            .value("glgm3150", tss::SphericalHarmonicsModel::glgm3150, get_docstring("PredefinedSphericalHarmonicsModel.glgm3150").c_str())
            .value("lpe200", tss::SphericalHarmonicsModel::lpe200, get_docstring("PredefinedSphericalHarmonicsModel.lpe200").c_str())
            .value("gggrx1200", tss::SphericalHarmonicsModel::lpe200, get_docstring("PredefinedSphericalHarmonicsModel.ggglpe200rx1200").c_str())
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
            .def_property("scaled_mean_moment_of_inertia",
                         &tss::SphericalHarmonicsGravityFieldSettings::getScaledMeanMomentOfInertia,
                         &tss::SphericalHarmonicsGravityFieldSettings::setScaledMeanMomentOfInertia,
                         get_docstring("SphericalHarmonicsGravityFieldSettings.scaled_mean_moment_of_inertia").c_str())
            .def_property_readonly("reference_radius",
                                   &tss::SphericalHarmonicsGravityFieldSettings::getReferenceRadius,
                                   get_docstring("SphericalHarmonicsGravityFieldSettings.reference_radius").c_str());


    py::class_<tss::FromFileSphericalHarmonicsGravityFieldSettings, std::shared_ptr<tss::FromFileSphericalHarmonicsGravityFieldSettings>,
            tss::SphericalHarmonicsGravityFieldSettings>(m, "FromFileSphericalHarmonicsGravityFieldSettings",
                                                         get_docstring("FromFileSphericalHarmonicsGravityFieldSettings").c_str());


    py::class_<tss::PolyhedronGravityFieldSettings, std::shared_ptr<tss::PolyhedronGravityFieldSettings>,
            tss::GravityFieldSettings>(m, "PolyhedronGravityFieldSettings",
                                       get_docstring("PolyhedronGravityFieldSettings").c_str())
            .def_property ("gravitational_parameter",
                           &tss::PolyhedronGravityFieldSettings::getGravitationalParameter,
                           &tss::PolyhedronGravityFieldSettings::resetGravitationalParameter,
                           get_docstring("PolyhedronGravityFieldSettings.gravitational_parameter").c_str())
            .def_property ("density",
                           &tss::PolyhedronGravityFieldSettings::getDensity,
                           &tss::PolyhedronGravityFieldSettings::resetDensity,
                           get_docstring("PolyhedronGravityFieldSettings.density").c_str())
            .def_property("associated_reference_frame",
                          &tss::PolyhedronGravityFieldSettings::getAssociatedReferenceFrame,
                          &tss::PolyhedronGravityFieldSettings::resetAssociatedReferenceFrame,
                          get_docstring("PolyhedronGravityFieldSettings.associated_reference_frame").c_str())
            .def_property_readonly ("vertices_coordinates",
                           &tss::PolyhedronGravityFieldSettings::getVerticesCoordinates,
                           get_docstring("PolyhedronGravityFieldSettings.vertices_coordinates").c_str())
            .def_property_readonly ("vertices_defining_each_facet",
                           &tss::PolyhedronGravityFieldSettings::getVerticesDefiningEachFacet,
                           get_docstring("PolyhedronGravityFieldSettings.vertices_defining_each_facet").c_str());


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

    m.def("from_file_spherical_harmonic",
          tss::fromFileSphericalHarmonicsGravityFieldSettings,
          py::arg("file"),
          py::arg("maximum_degree"),
          py::arg("maximum_order"),
          py::arg("associated_reference_frame") = "",
          py::arg("gravitational_parameter_index") = 0,
          py::arg("reference_radius_index") = 1,
          get_docstring("from_file_spherical_harmonic").c_str()
          );

    m.def("predefined_spherical_harmonic",
          tss::predefinedSphericalHarmonic,
          py::arg("predefined_model"),
          py::arg("maximum_degree") = -1,
          get_docstring("predefined_spherical_harmonic").c_str()
          );


    m.def("polyhedron_from_mu",
          py::overload_cast< const double,
          const Eigen::MatrixXd,
          const Eigen::MatrixXi,
          const std::string&,
          const double >( &tss::polyhedronGravitySettingsFromMu ),
          py::arg("gravitational_parameter"),
          py::arg("vertices_coordinates"),
          py::arg("vertices_defining_each_facet"),
          py::arg("associated_reference_frame"),
          py::arg("gravitational_constant") = tpc::GRAVITATIONAL_CONSTANT,
          get_docstring("polyhedron_from_mu").c_str()
          );

    m.def("polyhedron_from_density",
          py::overload_cast<
          const double,
          const Eigen::MatrixXd,
          const Eigen::MatrixXi,
          const std::string&,
          const double >( &tss::polyhedronGravitySettings ),
          py::arg("density"),
          py::arg("vertices_coordinates"),
          py::arg("vertices_defining_each_facet"),
          py::arg("associated_reference_frame"),
          py::arg("gravitational_constant") = tpc::GRAVITATIONAL_CONSTANT,
          get_docstring("polyhedron_from_density").c_str()
          );

    // Triaxial ellipsoid: overload 1
    m.def("sh_triaxial_ellipsoid_from_density",
          py::overload_cast< const double, const double, const double, const double, const int, const int,
          const std::string&, const double >(&tss::createHomogeneousTriAxialEllipsoidGravitySettings),
          py::arg("axis_a"),
          py::arg("axis_b"),
          py::arg("axis_c"),
          py::arg("density"),
          py::arg("maximum_degree"),
          py::arg("maximum_order"),
          py::arg("associated_reference_frame"),
          py::arg("gravitational_constant") = tudat::physical_constants::GRAVITATIONAL_CONSTANT,
          get_docstring("sh_triaxial_ellipsoid_from_density").c_str()
          );

    // Triaxial ellipsoid: overload 2
    m.def("sh_triaxial_ellipsoid_from_gravitational_parameter",
          py::overload_cast< const double, const double, const double, const int, const int,
          const std::string&, const double >(&tss::createHomogeneousTriAxialEllipsoidGravitySettings),
          py::arg("axis_a"),
          py::arg("axis_b"),
          py::arg("axis_c"),
          py::arg("maximum_degree"),
          py::arg("maximum_order"),
          py::arg("associated_reference_frame"),
          py::arg("gravitational_parameter"),
          get_docstring("sh_triaxial_ellipsoid_from_gravitational_parameter").c_str()
          );

    m.def("spherical_harmonic_triaxial_body",
          py::overload_cast< const double, const double, const double, const double, const int, const int,
          const std::string&, const double >(&tss::createHomogeneousTriAxialEllipsoidGravitySettingsDeprecated),
          py::arg("axis_a"),
              py::arg("axis_b"),
              py::arg("axis_c"),
              py::arg("density"),
              py::arg("maximum_degree"),
              py::arg("maximum_order"),
              py::arg("associated_reference_frame"),
              py::arg("gravitational_constant") = tudat::physical_constants::GRAVITATIONAL_CONSTANT
              );

    m.def("ring_model",
          &tss::ringGravitySettings,
          py::arg("gravitational_parameter"),
          py::arg("ring_radius"),
          py::arg("associated_reference_frame"),
          py::arg("elliptic_integral_s_from_d_and_b"),
          get_docstring("ring_model").c_str()
    );
}

}// namespace gravity_field
}// namespace environment_setup
}// namespace numerical_simulation
}// namespace tudatpy
