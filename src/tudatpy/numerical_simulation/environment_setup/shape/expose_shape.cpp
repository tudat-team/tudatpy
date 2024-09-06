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
R"doc(Base class for providing settings for body shape model.

	Functional (base) class for settings of body shape models that require no information in addition to their type.
	Body shape model settings requiring additional information must be defined using an object derived from this class.

)doc");

        py::class_<tss::SphericalBodyShapeSettings,
                std::shared_ptr<tss::SphericalBodyShapeSettings>,
                tss::BodyShapeSettings>(m, "SphericalBodyShapeSettings",
R"doc(Class for defining model settings of a strictly spherical body shape.

	`BodyShapeSettings` derived class for strictly spherical body shape model settings.
)doc")
                .def_property("radius", &tss::SphericalBodyShapeSettings::getRadius,
                              &tss::SphericalBodyShapeSettings::resetRadius,
R"doc(Radius specifying spherical body shape.
	)doc");


        py::class_<tss::OblateSphericalBodyShapeSettings,
                std::shared_ptr<tss::OblateSphericalBodyShapeSettings>,
                tss::BodyShapeSettings>(m, "OblateSphericalBodyShapeSettings",
R"doc(Class for defining model settings of a oblate spherical body shape.

	`BodyShapeSettings` derived class for oblate spherical body shape model settings.
)doc")
                .def_property("equatorial_radius", &tss::OblateSphericalBodyShapeSettings::getEquatorialRadius,
                              &tss::OblateSphericalBodyShapeSettings::resetEquatorialRadius,
R"doc(Equatorial radius of the oblate spherical body shape.
	)doc")
                .def_property("flattening", &tss::OblateSphericalBodyShapeSettings::getFlattening,
                              &tss::OblateSphericalBodyShapeSettings::resetFlattening,
R"doc(Flattening of spheroid shape model.
	)doc");

        py::class_<tss::PolyhedronBodyShapeSettings,
                std::shared_ptr<tss::PolyhedronBodyShapeSettings>,
                tss::BodyShapeSettings>(m, "PolyhedronBodyShapeSettings",
R"doc(Class for defining model settings of a polyhedron body shape.

	`BodyShapeSettings` derived class for polyhedron body shape model settings.
)doc")
                .def_property("vertices_coordinates", &tss::PolyhedronBodyShapeSettings::getVerticesCoordinates,
                              &tss::PolyhedronBodyShapeSettings::resetVerticesCoordinates,
R"doc(Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
row per vertex, 3 columns).

	)doc")
                .def_property("vertices_defining_each_facet", &tss::PolyhedronBodyShapeSettings::getVerticesDefiningEachFacet,
                              &tss::PolyhedronBodyShapeSettings::resetVerticesDefiningEachFacet,
R"doc(Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
when seen from the outside of the polyhedron.

	)doc")
                .def_property("compute_altitude_with_sign", &tss::PolyhedronBodyShapeSettings::getComputeAltitudeWithSign,
                              &tss::PolyhedronBodyShapeSettings::resetComputeAltitudeWithSign,
R"doc(Flag indicating whether the altitude should be computed with sign (i.e. >0 if above surface, <0 otherwise) or
having always a positive value. If the the sign of the altitude is not relevant, then setting it to *false* is
recommended, as it reduces the CPU time for computing the altitude.

	)doc")
                .def_property("just_compute_distance_to_vertices", &tss::PolyhedronBodyShapeSettings::getJustComputeDistanceToVertices,
                              &tss::PolyhedronBodyShapeSettings::resetJustComputeDistanceToVertices,
R"doc(Flag indicating whether the altitude should be computed just with respect to the polyhedron vertices (if flag
is set to *true*) or to all polyhedron features (vertices, facets and edges; happens if flag is set to
*false*). Depending on the application, it might be useful to set the flag to *true* for medium to high
altitudes, as it allows significantly reducing the CPU time (the resulting altitude errors depend on the
resolution of the used polyhedron and altitude itself).


description: |
  Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one

  description: |
    Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one

description: |
  Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
  the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise

  description: |
    Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
    the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise

description: |
  Flag indicating whether the altitude should be computed with sign (i.e. >0 if above surface, <0 otherwise) or

  description: |
    Flag indicating whether the altitude should be computed with sign (i.e. >0 if above surface, <0 otherwise) or
    having always a positive value. If the the sign of the altitude is not relevant, then setting it to *false* is

description: |
  Flag indicating whether the altitude should be computed just with respect to the polyhedron vertices (if flag
  is set to *true*) or to all polyhedron features (vertices, facets and edges; happens if flag is set to

  description: |
    Flag indicating whether the altitude should be computed just with respect to the polyhedron vertices (if flag
    is set to *true*) or to all polyhedron features (vertices, facets and edges; happens if flag is set to
    *false*). Depending on the application, it might be useful to set the flag to *true* for medium to high
    altitudes, as it allows significantly reducing the CPU time (the resulting altitude errors depend on the

	)doc");

        py::class_<tss::HybridBodyShapeSettings,
                std::shared_ptr<tss::HybridBodyShapeSettings>,
                tss::BodyShapeSettings>(m, "HybridBodyShapeSettings",
R"doc(Class for defining model settings of a hybrid body shape.

	`BodyShapeSettings` derived class for hybrid body shape model settings.
)doc")
                .def_property("low_resolution_body_shape_settings", &tss::HybridBodyShapeSettings::getLowResolutionBodyShapeSettings,
                              &tss::HybridBodyShapeSettings::resetLowResolutionBodyShapeSettings,
"")
                .def_property("high_resolution_body_shape_settings", &tss::HybridBodyShapeSettings::getHighResolutionBodyShapeSettings,
                              &tss::HybridBodyShapeSettings::resetHighResolutionBodyShapeSettings,
"")
                .def_property("switchover_altitude", &tss::HybridBodyShapeSettings::getSwitchoverAltitude,
                              &tss::HybridBodyShapeSettings::resetSwitchoverAltitude,
"");


        m.def("spherical",
              &tss::sphericalBodyShapeSettings,
              py::arg("radius"),
R"doc(Factory function for creating spherical body shape model settings.

	Factory function for settings object, defining strictly spherical body shape model entirely from single radius parameter.


	:param radius:
		Radius specifying spherical body shape.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.shape.BodyShapeModel` derived :class:`~tudatpy.numerical_simulation.environment_setup.shape.SphericalBodyShapeSettings` class
)doc");

        m.def("spherical_spice",
              &tss::fromSpiceSphericalBodyShapeSettings,
R"doc(Factory function for creating spherical body shape model settings entirely from spice.

	Factory function for settings object, defining spherical body shape model entirely from spice parameters.

	:return:
		Instance of :class:`~tudatpy.numerical_simulation.environment_setup.shape.BodyShapeModel` class
)doc");

        m.def("oblate_spherical",
              &tss::oblateSphericalBodyShapeSettings,
              py::arg("equatorial_radius"),
              py::arg("flattening"),
R"doc(Factory function for creating oblate spherical body shape model settings.

	Factory function for settings object, defining oblate spherical body shape model from equatorial radius and flattening parameter.


	:param equatorial_radius:
		Equatorial radius specifying oblate spherical body shape.
	:param flattening:
		Flattening parameter specifying oblate spherical body shape.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.shape.BodyShapeModel` derived :class:`~tudatpy.numerical_simulation.environment_setup.shape.OblateSphericalBodyShapeSettings` class
)doc");

        m.def("polyhedron",
              &tss::polyhedronBodyShapeSettings,
              py::arg("vertices_coordinates"),
              py::arg("vertices_defining_each_facet"),
              py::arg("compute_altitude_with_sign") = true,
              py::arg("just_compute_distance_to_vertices") = false,
R"doc(Factory function for creating a polyhedron body shape model settings.

	Factory function for settings object, defining a polyhedron shape model.

	Note 1: The evaluation of the altitude with a polyhedron model tends to be computationally expensive. To reduce the
	computational time, it might be useful to instead define a hybrid shape model (see
	:func:`~tudatpy.numerical_simulation.environment_setup.shape.hybrid`), which allows using a high-resolution
	polyhedron (with a large number of facets) at low altitudes and a low-resolution one (with smaller number of facets)
	at high-altitudes.

	Note 2: If the goal of using the shape model is only to detect collisions with the surface and not to explicitly
	obtain the altitude, it is instead recommended to use the Laplacian of the gravitational potential (see
	:func:`~tudatpy.numerical_simulation.environment_setup.dependent_variable.yaml.gravity_field_laplacian_of_potential`).
	This allows reducing the computational time, but is only valid if the same polyhedron model that is used to define
	the gravitational acceleration should also be used to detect the impacts.


	:param vertices_coordinates:
		Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
		row per vertex, 3 columns).

	:param vertices_defining_each_facet:
		Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
		the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
		when seen from the outside of the polyhedron.

	:param compute_altitude_with_sign:
		Flag indicating whether the altitude should be computed with sign (i.e. >0 if above surface, <0 otherwise) or
		having always a positive value. If the the sign of the altitude is not relevant, then setting it to *false* is
		recommended, as it reduces the CPU time.

	:param just_compute_distance_to_vertices:
		Flag indicating whether the altitude should be computed just with respect to the polyhedron vertices (if flag
		is set to *true*) or to all polyhedron features (vertices, facets and edges; happens if flag is set to
		*false*). Depending on the application, it might be useful to set the flag to *true* for medium to high
		altitudes, as it allows significantly reducing the CPU time (the resulting altitude errors depend on the
		resolution of the used polyhedron and altitude itself).

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.shape.BodyShapeModel` derived
		:class:`~tudatpy.numerical_simulation.environment_setup.shape.PolyhedronBodyShapeSettings` class

)doc");

        m.def("hybrid",
              &tss::hybridBodyShapeSettings,
              py::arg("low_resolution_body_shape_settings"),
              py::arg("high_resolution_body_shape_settings"),
              py::arg("switchover_altitude"),
R"doc(Factory function for creating hybrid body shape model settings.

	Factory function for settings object, defining a hybrid shape model.

	The hybrid shape model is constituded by two shape models: a low-resolution model which is used at high altitudes
	(above the switchover altitude) and a high-resolution model used at low altitudes (below the switchover altitude).
	In each computation of the altitude, the altitude is first computed with the low-resolution model. The
	low-resolution altitude is then compared to the switchover altitude to decide whether to compute the high-resolution
	altitude.

	The hybrid shape model is useful when the evaluation of the high-resolution model is computationally expensive
	(e.g. polyhedron model).


	:param low_resolution_body_shape_settings:
		Settings of the shape model that is to be used to compute the altitude at high altitudes (above the switchover
		altitude).

	:param high_resolution_body_shape_settings:
		Settings of the shape model that is to be used to compute the altitude at low altitudes (below the switchover
		altitude).

	:param switchover_altitude:
		Altitude at which the model used to compute the altitude is changed. The high-resolution model is used for
		altitudes below the switchover altitude, the low-resolution model for altitudes above it.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.shape.BodyShapeModel` derived
		:class:`~tudatpy.numerical_simulation.environment_setup.shape.HybridBodyShapeSettings` class

)doc");

    }

}// namespace shape
}// namespace environment_setup
}// namespace numerical_simulation
}// namespace tudatpy
