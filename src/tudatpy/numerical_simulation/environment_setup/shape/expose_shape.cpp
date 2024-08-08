/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


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
namespace shape {

    PYBIND11_MODULE(expose_shape, m) {

        py::class_<tss::BodyShapeSettings,
                std::shared_ptr<tss::BodyShapeSettings>>(m, "BodyShapeSettings",
                                                         get_docstring("BodyShapeSettings").c_str());

        py::class_<tss::SphericalBodyShapeSettings,
                std::shared_ptr<tss::SphericalBodyShapeSettings>,
                tss::BodyShapeSettings>(m, "SphericalBodyShapeSettings",
                                        get_docstring("SphericalBodyShapeSettings").c_str())
                .def_property("radius", &tss::SphericalBodyShapeSettings::getRadius,
                              &tss::SphericalBodyShapeSettings::resetRadius,
                              get_docstring("SphericalBodyShapeSettings.radius").c_str());


        py::class_<tss::OblateSphericalBodyShapeSettings,
                std::shared_ptr<tss::OblateSphericalBodyShapeSettings>,
                tss::BodyShapeSettings>(m, "OblateSphericalBodyShapeSettings",
                                        get_docstring("OblateSphericalBodyShapeSettings").c_str())
                .def_property("equatorial_radius", &tss::OblateSphericalBodyShapeSettings::getEquatorialRadius,
                              &tss::OblateSphericalBodyShapeSettings::resetEquatorialRadius,
                              get_docstring("OblateSphericalBodyShapeSettings.equatorial_radius").c_str())
                .def_property("flattening", &tss::OblateSphericalBodyShapeSettings::getFlattening,
                              &tss::OblateSphericalBodyShapeSettings::resetFlattening,
                              get_docstring("OblateSphericalBodyShapeSettings.flattening").c_str());

        py::class_<tss::PolyhedronBodyShapeSettings,
                std::shared_ptr<tss::PolyhedronBodyShapeSettings>,
                tss::BodyShapeSettings>(m, "PolyhedronBodyShapeSettings",
                                        get_docstring("PolyhedronBodyShapeSettings").c_str())
                .def_property("vertices_coordinates", &tss::PolyhedronBodyShapeSettings::getVerticesCoordinates,
                              &tss::PolyhedronBodyShapeSettings::resetVerticesCoordinates,
                              get_docstring("PolyhedronBodyShapeSettings.vertices_coordinates").c_str())
                .def_property("vertices_defining_each_facet", &tss::PolyhedronBodyShapeSettings::getVerticesDefiningEachFacet,
                              &tss::PolyhedronBodyShapeSettings::resetVerticesDefiningEachFacet,
                              get_docstring("PolyhedronBodyShapeSettings.vertices_defining_each_facet").c_str())
                .def_property("compute_altitude_with_sign", &tss::PolyhedronBodyShapeSettings::getComputeAltitudeWithSign,
                              &tss::PolyhedronBodyShapeSettings::resetComputeAltitudeWithSign,
                              get_docstring("PolyhedronBodyShapeSettings.compute_altitude_with_sign").c_str())
                .def_property("just_compute_distance_to_vertices", &tss::PolyhedronBodyShapeSettings::getJustComputeDistanceToVertices,
                              &tss::PolyhedronBodyShapeSettings::resetJustComputeDistanceToVertices,
                              get_docstring("PolyhedronBodyShapeSettings.just_compute_distance_to_vertices").c_str());

        py::class_<tss::HybridBodyShapeSettings,
                std::shared_ptr<tss::HybridBodyShapeSettings>,
                tss::BodyShapeSettings>(m, "HybridBodyShapeSettings",
                                        get_docstring("HybridBodyShapeSettings").c_str())
                .def_property("low_resolution_body_shape_settings", &tss::HybridBodyShapeSettings::getLowResolutionBodyShapeSettings,
                              &tss::HybridBodyShapeSettings::resetLowResolutionBodyShapeSettings,
                              get_docstring("HybridBodyShapeSettings.vertices_coordinates").c_str())
                .def_property("high_resolution_body_shape_settings", &tss::HybridBodyShapeSettings::getHighResolutionBodyShapeSettings,
                              &tss::HybridBodyShapeSettings::resetHighResolutionBodyShapeSettings,
                              get_docstring("HybridBodyShapeSettings.vertices_defining_each_facet").c_str())
                .def_property("switchover_altitude", &tss::HybridBodyShapeSettings::getSwitchoverAltitude,
                              &tss::HybridBodyShapeSettings::resetSwitchoverAltitude,
                              get_docstring("HybridBodyShapeSettings.compute_altitude_with_sign").c_str());


        m.def("spherical",
              &tss::sphericalBodyShapeSettings,
              py::arg("radius"),
              get_docstring("spherical").c_str());

        m.def("spherical_spice",
              &tss::fromSpiceSphericalBodyShapeSettings,
              get_docstring("spherical_spice").c_str());

        m.def("oblate_spherical",
              &tss::oblateSphericalBodyShapeSettings,
              py::arg("equatorial_radius"),
              py::arg("flattening"),
              get_docstring("oblate_spherical").c_str());

        m.def("polyhedron",
              &tss::polyhedronBodyShapeSettings,
              py::arg("vertices_coordinates"),
              py::arg("vertices_defining_each_facet"),
              py::arg("compute_altitude_with_sign") = true,
              py::arg("just_compute_distance_to_vertices") = false,
              get_docstring("polyhedron").c_str());

        m.def("hybrid",
              &tss::hybridBodyShapeSettings,
              py::arg("low_resolution_body_shape_settings"),
              py::arg("high_resolution_body_shape_settings"),
              py::arg("switchover_altitude"),
              get_docstring("hybrid").c_str());

    }

}// namespace shape
}// namespace environment_setup
}// namespace numerical_simulation
}// namespace tudatpy
