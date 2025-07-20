import numpy
import pybind11_stubgen.typing_ext
import typing
__all__ = ['centroid', 'inertia_tensor_from_density', 'inertia_tensor_from_gravitational_parameter', 'modify_centroid', 'surface_area', 'volume']

def centroid(vertices_coordinates: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], vertices_defining_each_facet: typing.Annotated[numpy.ndarray, numpy.int32, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]:
    """Computes the position of the centroid of a polyhedron :cite:p:`Dobrovolskis1996`.
    
    
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
        Position of the centroid."""

def inertia_tensor_from_density(vertices_coordinates: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], vertices_defining_each_facet: typing.Annotated[numpy.ndarray, numpy.int32, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], density: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
    """Compute the inertia tensor of a polyhedron, from the density.
    
    Computes the inertia tensor of a polyhedron, according to :cite:t:`Dobrovolskis1996`.
    
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
        Inertia tensor."""

def inertia_tensor_from_gravitational_parameter(vertices_coordinates: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], vertices_defining_each_facet: typing.Annotated[numpy.ndarray, numpy.int32, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], gravitational_parameter: float, gravitational_constant: float) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 3)]:
    """Compute the inertia tensor of a polyhedron, from the gravitational parameter.
    
    Computes the inertia tensor of a polyhedron, according to :cite:t:`Dobrovolskis1996`.
    
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
        Gravitational parameter :math:`\\mu` of gravity field.
    Returns
    -------
    numpy.ndarray
        Inertia tensor."""

def modify_centroid(vertices_coordinates: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], vertices_defining_each_facet: typing.Annotated[numpy.ndarray, numpy.int32, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], desired_centroid: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.FixedSize(3, 1)]) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
    """Modifies vertex coordinates of the polyhedron based on the desired position of the centroid.
    
    Modifies the coordinates of the polyhedron vertices, such that the centroid of the modified polyhedron coincides
    with the specified position. The centroid is computed according to :cite:t:`Dobrovolskis1996`.
    
    
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
        Vertices coordinates of the modified polyhedron, which has the specified centroid position."""

def surface_area(vertices_coordinates: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], vertices_defining_each_facet: typing.Annotated[numpy.ndarray, numpy.int32, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]) -> float:
    """Computes the surface area of a polyhedron :cite:p:`Dobrovolskis1996`.
    
    
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
        Surface area."""

def volume(vertices_coordinates: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], vertices_defining_each_facet: typing.Annotated[numpy.ndarray, numpy.int32, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]) -> float:
    """Computes the volume of a polyhedron :cite:p:`Dobrovolskis1996`.
    
    
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
        Volume."""