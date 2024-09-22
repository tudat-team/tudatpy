import typing
import numpy
__all__ = ['centroid', 'inertia_tensor_from_density', 'inertia_tensor_from_gravitational_parameter', 'modify_centroid', 'surface_area', 'volume']

def centroid(vertices_coordinates: numpy.ndarray, vertices_defining_each_facet: numpy.ndarray) -> numpy.ndarray:
    """Computes the position of the centroid of a polyhedron.
	
	Computes the position of the centroid of a polyhedron, according to Dobrovolskis [1]_.
	
	
	:param vertices_coordinates:
			Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
			row per vertex, 3 columns).
	
	:param vertices_defining_each_facet:
			Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
			the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
			when seen from the outside of the polyhedron.
	
	:return:
			Position of the centroid.
	"""

def inertia_tensor_from_density(vertices_coordinates: numpy.ndarray, vertices_defining_each_facet: numpy.ndarray, density: float) -> numpy.ndarray:
    """Compute the inertia tensor of a polyhedron, from the density.
	
	Computes the inertia tensor of a polyhedron, according to Dobrovolskis [1]_.
	
	The mass distribution is defined using the density of the polyhedron. To instead use the gravitational
	parameter see :func:`~tudatpy.astro.polyhedron_utilities.inertia_tensor_from_gravitational_parameter`.
	
	
	:param vertices_coordinates:
			Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
			row per vertex, 3 columns).
	
	:param vertices_defining_each_facet:
			Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
			the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
			when seen from the outside of the polyhedron.
	
	:param density:
			Density of the polyhedron
	
	:return:
			Inertia tensor.
	"""

def inertia_tensor_from_gravitational_parameter(vertices_coordinates: numpy.ndarray, vertices_defining_each_facet: numpy.ndarray, gravitational_parameter: float, gravitational_constant: float) -> numpy.ndarray:
    """Compute the inertia tensor of a polyhedron, from the gravitational parameter.
	
	Computes the inertia tensor of a polyhedron, according to Dobrovolskis [1]_.
	
	The mass distribution is defined using the gravitational parameter of the polyhedron. To instead use the density
	see :func:`~tudatpy.astro.polyhedron_utilities.inertia_tensor_from_density`.
	
	
	:param vertices_coordinates:
			Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
			row per vertex, 3 columns).
	
	:param vertices_defining_each_facet:
			Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
			the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
			when seen from the outside of the polyhedron.
	
	:param gravitational_parameter:
			Gravitational parameter :math:`\\mu` of gravity field.
	:return:
			Inertia tensor.
	"""

def modify_centroid(vertices_coordinates: numpy.ndarray, vertices_defining_each_facet: numpy.ndarray, desired_centroid: numpy.ndarray) -> numpy.ndarray:
    """Modifies the position of the centroid of the polyhedron.
	
	Modifies the coordinates of the polyhedron vertices, such that the centroid of the modified polyhedron coincides
	with the specified position. The centroid is computed according to Dobrovolskis [1]_.
	
	
	:param vertices_coordinates:
			Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
			row per vertex, 3 columns).
	
	:param vertices_defining_each_facet:
			Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
			the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
			when seen from the outside of the polyhedron.
	
	:param desired_centroid:
			Desired position of the centroid.
	
	:return:
			Vertices coordinates of the modified polyhedron, which has the specified centroid position.
	"""

def surface_area(vertices_coordinates: numpy.ndarray, vertices_defining_each_facet: numpy.ndarray) -> float:
    """Computes the surface area of a polyhedron.
	
	Computes the surface area of a polyhedron, according to Dobrovolskis [1]_.
	
	
	:param vertices_coordinates:
			Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
			row per vertex, 3 columns).
	
	:param vertices_defining_each_facet:
			Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
			the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
			when seen from the outside of the polyhedron.
	
	:return:
			Surface area.
	"""

def volume(vertices_coordinates: numpy.ndarray, vertices_defining_each_facet: numpy.ndarray) -> float:
    """Computes the volume of a polyhedron.
	
	Computes the volume of a polyhedron, according to Dobrovolskis [1]_.
	
	
	:param vertices_coordinates:
			Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
			row per vertex, 3 columns).
	
	:param vertices_defining_each_facet:
			Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
			the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
			when seen from the outside of the polyhedron.
	
	:return:
			Volume.
	"""