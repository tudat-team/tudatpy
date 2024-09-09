import typing
import numpy
__all__ = ['BodyShapeSettings', 'HybridBodyShapeSettings', 'OblateSphericalBodyShapeSettings', 'PolyhedronBodyShapeSettings', 'SphericalBodyShapeSettings', 'hybrid', 'oblate_spherical', 'polyhedron', 'spherical', 'spherical_spice']

class BodyShapeSettings:
    """Base class for providing settings for body shape model.
	
	Functional (base) class for settings of body shape models that require no information in addition to their type.
	Body shape model settings requiring additional information must be defined using an object derived from this class.
	"""

class HybridBodyShapeSettings(BodyShapeSettings):
    """Class for defining model settings of a hybrid body shape.
	
	`BodyShapeSettings` derived class for hybrid body shape model settings.
	"""
    high_resolution_body_shape_settings: BodyShapeSettings
    low_resolution_body_shape_settings: BodyShapeSettings
    switchover_altitude: float

class OblateSphericalBodyShapeSettings(BodyShapeSettings):
    """Class for defining model settings of a oblate spherical body shape.
	
	`BodyShapeSettings` derived class for oblate spherical body shape model settings.
	"""

    @property
    def equatorial_radius(self) -> float:
        """
        Equatorial radius of the oblate spherical body shape.
        	
        """

    @equatorial_radius.setter
    def equatorial_radius(self, arg1: float) -> None:
        ...

    @property
    def flattening(self) -> float:
        """
        Flattening of spheroid shape model.
        	
        """

    @flattening.setter
    def flattening(self, arg1: float) -> None:
        ...

class PolyhedronBodyShapeSettings(BodyShapeSettings):
    """Class for defining model settings of a polyhedron body shape.
	
	`BodyShapeSettings` derived class for polyhedron body shape model settings.
	"""

    @property
    def compute_altitude_with_sign(self) -> bool:
        """
        Flag indicating whether the altitude should be computed with sign (i.e. >0 if above surface, <0 otherwise) or
        having always a positive value. If the the sign of the altitude is not relevant, then setting it to *false* is
        recommended, as it reduces the CPU time for computing the altitude.
        """

    @compute_altitude_with_sign.setter
    def compute_altitude_with_sign(self, arg1: bool) -> None:
        ...

    @property
    def just_compute_distance_to_vertices(self) -> bool:
        """
        Flag indicating whether the altitude should be computed just with respect to the polyhedron vertices (if flag
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
        """

    @just_compute_distance_to_vertices.setter
    def just_compute_distance_to_vertices(self, arg1: bool) -> None:
        ...

    @property
    def vertices_coordinates(self) -> numpy.ndarray:
        """
        Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
        row per vertex, 3 columns).
        """

    @vertices_coordinates.setter
    def vertices_coordinates(self, arg1: numpy.ndarray) -> None:
        ...

    @property
    def vertices_defining_each_facet(self) -> numpy.ndarray:
        """
        Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
        the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
        when seen from the outside of the polyhedron.
        """

    @vertices_defining_each_facet.setter
    def vertices_defining_each_facet(self, arg1: numpy.ndarray) -> None:
        ...

class SphericalBodyShapeSettings(BodyShapeSettings):
    """Class for defining model settings of a strictly spherical body shape.
	
	`BodyShapeSettings` derived class for strictly spherical body shape model settings.
	"""

    @property
    def radius(self) -> float:
        """
        Radius specifying spherical body shape.
        	
        """

    @radius.setter
    def radius(self, arg1: float) -> None:
        ...

def hybrid(low_resolution_body_shape_settings: BodyShapeSettings, high_resolution_body_shape_settings: BodyShapeSettings, switchover_altitude: float) -> BodyShapeSettings:
    """Factory function for creating hybrid body shape model settings.
	
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
	"""

def oblate_spherical(equatorial_radius: float, flattening: float) -> BodyShapeSettings:
    """Factory function for creating oblate spherical body shape model settings.
	
	Factory function for settings object, defining oblate spherical body shape model from equatorial radius and flattening parameter.
	
	
	:param equatorial_radius:
			Equatorial radius specifying oblate spherical body shape.
	:param flattening:
			Flattening parameter specifying oblate spherical body shape.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.shape.BodyShapeModel` derived :class:`~tudatpy.numerical_simulation.environment_setup.shape.OblateSphericalBodyShapeSettings` class
	"""

def polyhedron(vertices_coordinates: numpy.ndarray, vertices_defining_each_facet: numpy.ndarray, compute_altitude_with_sign: bool=True, just_compute_distance_to_vertices: bool=False) -> BodyShapeSettings:
    """Factory function for creating a polyhedron body shape model settings.
	
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
	"""

def spherical(radius: float) -> BodyShapeSettings:
    """Factory function for creating spherical body shape model settings.
	
	Factory function for settings object, defining strictly spherical body shape model entirely from single radius parameter.
	
	
	:param radius:
			Radius specifying spherical body shape.
	:return:
			Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.shape.BodyShapeModel` derived :class:`~tudatpy.numerical_simulation.environment_setup.shape.SphericalBodyShapeSettings` class
	"""

def spherical_spice() -> BodyShapeSettings:
    """Factory function for creating spherical body shape model settings entirely from spice.
	
	Factory function for settings object, defining spherical body shape model entirely from spice parameters.
	
	:return:
			Instance of :class:`~tudatpy.numerical_simulation.environment_setup.shape.BodyShapeModel` class
	"""