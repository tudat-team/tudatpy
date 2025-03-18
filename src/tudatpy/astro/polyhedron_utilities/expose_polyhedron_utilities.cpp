/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <tudat/astro/basic_astro/polyhedronFuntions.h>

namespace tba = tudat::basic_astrodynamics;

namespace py = pybind11;

namespace tudatpy
{

namespace astro
{
namespace polyhedron_utilities
{

PYBIND11_MODULE( expose_polyhedron_utilities, m )
{
    m.def( "surface_area",
           &tba::computePolyhedronSurfaceArea,
           py::arg( "vertices_coordinates" ),
           py::arg( "vertices_defining_each_facet" ),
           R"doc(

 Computes the surface area of a polyhedron [1]_.


 Parameters
 ----------
 vertices_coordinates : numpy.ndarray
     Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
     row per vertex, 3 columns).

 vertices_defining_each_facet : numpy.ndarray
     Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
     the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
     when seen from the outside of the polyhedron.

 Returns
 -------
 float
     Surface area.






     )doc" );

    m.def( "volume",
           &tba::computePolyhedronVolume,
           py::arg( "vertices_coordinates" ),
           py::arg( "vertices_defining_each_facet" ),
           R"doc(

 Computes the volume of a polyhedron [1]_.


 Parameters
 ----------
 vertices_coordinates : numpy.ndarray
     Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
     row per vertex, 3 columns).

 vertices_defining_each_facet : numpy.ndarray
     Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
     the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
     when seen from the outside of the polyhedron.

 Returns
 -------
 float
     Volume.






     )doc" );

    m.def( "centroid",
           &tba::computePolyhedronCentroidPosition,
           py::arg( "vertices_coordinates" ),
           py::arg( "vertices_defining_each_facet" ),
           R"doc(

 Computes the position of the centroid of a polyhedron [1]_.


 Parameters
 ----------
 vertices_coordinates : numpy.ndarray
     Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
     row per vertex, 3 columns).

 vertices_defining_each_facet : numpy.ndarray
     Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
     the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
     when seen from the outside of the polyhedron.

 Returns
 -------
 numpy.ndarray
     Position of the centroid.






     )doc" );

    m.def( "modify_centroid",
           &tba::modifyPolyhedronCentroidPosition,
           py::arg( "vertices_coordinates" ),
           py::arg( "vertices_defining_each_facet" ),
           py::arg( "desired_centroid" ),
           R"doc(

 Modifies vertex coordinates of the polyhedron based on the desired position of the centroid.

 Modifies the coordinates of the polyhedron vertices, such that the centroid of the modified polyhedron coincides
 with the specified position. The centroid is computed according to Dobrovolskis [1]_.


 Parameters
 ----------
 vertices_coordinates : numpy.ndarray
     Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
     row per vertex, 3 columns).

 vertices_defining_each_facet : numpy.ndarray
     Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
     the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
     when seen from the outside of the polyhedron.

 desired_centroid : numpy.ndarray
     Desired position of the centroid.

 Returns
 -------
 numpy.ndarray
     Vertices coordinates of the modified polyhedron, which has the specified centroid position.






     )doc" );

    m.def( "inertia_tensor_from_density",
           py::overload_cast< const Eigen::MatrixXd &, const Eigen::MatrixXi &, const double >(
                   &tba::computePolyhedronInertiaTensor ),
           py::arg( "vertices_coordinates" ),
           py::arg( "vertices_defining_each_facet" ),
           py::arg( "density" ),
           R"doc(

 Compute the inertia tensor of a polyhedron, from the density.

 Computes the inertia tensor of a polyhedron, according to Dobrovolskis [1]_.

 The mass distribution is defined using the density of the polyhedron. To instead use the gravitational
 parameter see :func:`~tudatpy.astro.polyhedron_utilities.inertia_tensor_from_gravitational_parameter`.


 Parameters
 ----------
 vertices_coordinates : numpy.ndarray
     Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
     row per vertex, 3 columns).

 vertices_defining_each_facet : numpy.ndarray
     Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
     the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
     when seen from the outside of the polyhedron.

 density : float
     Density of the polyhedron

 Returns
 -------
 numpy.ndarray
     Inertia tensor.






     )doc" );

    m.def( "inertia_tensor_from_gravitational_parameter",
           py::overload_cast< const Eigen::MatrixXd &,
                              const Eigen::MatrixXi &,
                              const double,
                              const double >( &tba::computePolyhedronInertiaTensor ),
           py::arg( "vertices_coordinates" ),
           py::arg( "vertices_defining_each_facet" ),
           py::arg( "gravitational_parameter" ),
           py::arg( "gravitational_constant" ),
           R"doc(

 Compute the inertia tensor of a polyhedron, from the gravitational parameter.

 Computes the inertia tensor of a polyhedron, according to Dobrovolskis [1]_.

 The mass distribution is defined using the gravitational parameter of the polyhedron. To instead use the density
 see :func:`~tudatpy.astro.polyhedron_utilities.inertia_tensor_from_density`.


 Parameters
 ----------
 vertices_coordinates : numpy.ndarray
     Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
     row per vertex, 3 columns).

 vertices_defining_each_facet : numpy.ndarray
     Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
     the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
     when seen from the outside of the polyhedron.

 gravitational_parameter : float
     Gravitational parameter :math:`\mu` of gravity field.
 Returns
 -------
 numpy.ndarray
     Inertia tensor.






     )doc" );
}

}  // namespace polyhedron_utilities
}  // namespace astro
}  // namespace tudatpy
