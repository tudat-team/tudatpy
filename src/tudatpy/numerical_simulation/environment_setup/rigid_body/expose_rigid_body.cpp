/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


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

namespace tudatpy {
namespace numerical_simulation {
namespace environment_setup {
namespace rigid_body {

PYBIND11_MODULE(expose_rigid_body, m) {


//    py::enum_<tss::RigidBodyPropertiesType>(m, "RigidBodyPropertiesType",
//        .export_values();
//
    py::class_<tss::RigidBodyPropertiesSettings, std::shared_ptr<tss::RigidBodyPropertiesSettings>>(
        m, "RigidBodyPropertiesSettings",
"")
        .def_property_readonly("body_mass_property_type", &tss::RigidBodyPropertiesSettings::getRigidBodyPropertiesType,
"");


    m.def("constant_rigid_body_properties",
          tss::constantRigidBodyPropertiesSettings,
          py::arg("mass"),
          py::arg("center_of_mass") = Eigen::Vector3d::Constant( TUDAT_NAN ),
          py::arg("inertia_tensor") = Eigen::Matrix3d::Constant( TUDAT_NAN ),
""
    );

    m.def("custom_time_dependent_rigid_body_properties",
          tss::fromFunctionRigidBodyPropertiesSettings,
          py::arg("mass_function"),
          py::arg("center_of_mass_function") = nullptr,
          py::arg("inertia_tensor_function") = nullptr,
""
    );

    m.def("custom_mass_dependent_rigid_body_properties",
          tss::massDependentMassDistributionSettings,
          py::arg("mass"),
          py::arg("center_of_mass_function") = nullptr,
          py::arg("inertia_tensor_function") = nullptr,
""
    );


}

}// namespace rigid_bpdy
}// namespace environment_setup
}// namespace numerical_simulation
}// namespace tudatpy
