import numpy
import pybind11_stubgen.typing_ext
import typing
__all__ = ['BodyShapeSettings', 'HybridBodyShapeSettings', 'OblateSphericalBodyShapeSettings', 'PolyhedronBodyShapeSettings', 'SphericalBodyShapeSettings', 'hybrid', 'oblate_spherical', 'oblate_spherical_spice', 'polyhedron', 'spherical', 'spherical_spice']

class BodyShapeSettings:
    """Base class for providing settings for body shape model.
    
    Functional (base) class for settings of body shape models that require no information in addition to their type.
    Body shape model settings requiring additional information must be defined using an object derived from this class."""

class HybridBodyShapeSettings(BodyShapeSettings):
    """Class for defining model settings of a hybrid body shape.
    
    `BodyShapeSettings` derived class for hybrid body shape model settings."""

    @property
    def high_resolution_body_shape_settings(self) -> BodyShapeSettings:
        """
        No documentation found.
        """

    @high_resolution_body_shape_settings.setter
    def high_resolution_body_shape_settings(self, arg1: BodyShapeSettings) -> None:
        ...

    @property
    def low_resolution_body_shape_settings(self) -> BodyShapeSettings:
        """
        No documentation found.
        """

    @low_resolution_body_shape_settings.setter
    def low_resolution_body_shape_settings(self, arg1: BodyShapeSettings) -> None:
        ...

    @property
    def switchover_altitude(self) -> float:
        """
        No documentation found.
        """

    @switchover_altitude.setter
    def switchover_altitude(self, arg1: float) -> None:
        ...

class OblateSphericalBodyShapeSettings(BodyShapeSettings):
    """Class for defining model settings of a oblate spherical body shape.
    
    `BodyShapeSettings` derived class for oblate spherical body shape model settings."""

    @property
    def equatorial_radius(self) -> float:
        """
                 **read-only**
        
                 Equatorial radius of the oblate spherical body shape.
        
                 :type: float
        """

    @equatorial_radius.setter
    def equatorial_radius(self, arg1: float) -> None:
        ...

    @property
    def flattening(self) -> float:
        """
                 **read-only**
        
                 Flattening of spheroid shape model.
        
                 :type: float
        """

    @flattening.setter
    def flattening(self, arg1: float) -> None:
        ...

class PolyhedronBodyShapeSettings(BodyShapeSettings):
    """Class for defining model settings of a polyhedron body shape.
    
    `BodyShapeSettings` derived class for polyhedron body shape model settings."""

    @property
    def compute_altitude_with_sign(self) -> bool:
        """
                 Flag indicating whether the altitude should be computed with sign (i.e. >0 if above surface, <0 otherwise) or
                 having always a positive value. If the the sign of the altitude is not relevant, then setting it to *false* is
                 recommended, as it reduces the CPU time for computing the altitude.
        
        
                 :type: bool, default=True
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
        
        
                 :type: bool, default=False
        """

    @just_compute_distance_to_vertices.setter
    def just_compute_distance_to_vertices(self, arg1: bool) -> None:
        ...

    @property
    def vertices_coordinates(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
                 row per vertex, 3 columns).
        
        
                 :type: numpy.ndarray
        """

    @vertices_coordinates.setter
    def vertices_coordinates(self, arg1: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]) -> None:
        ...

    @property
    def vertices_defining_each_facet(self) -> typing.Annotated[numpy.ndarray, numpy.int32, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
                 the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
                 when seen from the outside of the polyhedron.
        
        
                 :type: numpy.ndarray
        """

    @vertices_defining_each_facet.setter
    def vertices_defining_each_facet(self, arg1: typing.Annotated[numpy.ndarray, numpy.int32, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]) -> None:
        ...

class SphericalBodyShapeSettings(BodyShapeSettings):
    """Class for defining model settings of a strictly spherical body shape.
    
    `BodyShapeSettings` derived class for strictly spherical body shape model settings."""

    @property
    def radius(self) -> float:
        """
                 **read-only**
        
                 Radius specifying spherical body shape.
        
                 :type: float
        """

    @radius.setter
    def radius(self, arg1: float) -> None:
        ...

def hybrid(low_resolution_body_shape_settings: BodyShapeSettings, high_resolution_body_shape_settings: BodyShapeSettings, switchover_altitude: float) -> BodyShapeSettings:
    """Function for creating hybrid body shape model settings.
    
    Function for settings object, defining a hybrid shape model.
    
    The hybrid shape model is constituted by two shape models: a low-resolution model which is used at high altitudes
    (above the switchover altitude) and a high-resolution model used at low altitudes (below the switchover altitude).
    In each computation of the altitude, the altitude is first computed with the low-resolution model. The
    low-resolution altitude is then compared to the switchover altitude to decide whether to compute the high-resolution
    altitude.
    
    The hybrid shape model is useful when the evaluation of the high-resolution model is computationally expensive
    (e.g. polyhedron model).
    
    
    Parameters
    ----------
    low_resolution_body_shape_settings : BodyShapeSettings
        Settings of the shape model that is to be used to compute the altitude at high altitudes (above the switchover
        altitude).
    
    high_resolution_body_shape_settings : BodyShapeSettings
        Settings of the shape model that is to be used to compute the altitude at low altitudes (below the switchover
        altitude).
    
    switchover_altitude : float
        Altitude at which the model used to compute the altitude is changed. The high-resolution model is used for
        altitudes below the switchover altitude, the low-resolution model for altitudes above it.
    
    Returns
    -------
    HybridBodyShapeSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.shape.BodyShapeSettings` derived
        :class:`~tudatpy.dynamics.environment_setup.shape.HybridBodyShapeSettings` class"""

def oblate_spherical(equatorial_radius: float, flattening: float) -> BodyShapeSettings:
    """Function for creating oblate spherical body shape model settings.
    
    Function for settings object, defining oblate spherical body shape model from equatorial radius and flattening parameter.
    
    
    Parameters
    ----------
    equatorial_radius : float
        Equatorial radius specifying oblate spherical body shape.
    flattening : float
        Flattening parameter specifying oblate spherical body shape.
    Returns
    -------
    OblateSphericalBodyShapeSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.shape.BodyShapeSettings` derived :class:`~tudatpy.dynamics.environment_setup.shape.OblateSphericalBodyShapeSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create a :class:`~tudatpy.dynamics.environment_setup.shape.BodyShapeSettings` using a perfectly oblate spherical shape model:
    
    .. code-block:: python
    
       # define parameters describing oblate spherical model
       body_radius = 6378.0E3
       body_flattening = 1.0 / 300.0
       # create shape model settings
       body_settings.get( "Earth" ).shape_settings = environment_setup.shape.oblate_spherical( body_radius, body_flattening )"""

def oblate_spherical_spice() -> BodyShapeSettings:
    """No documentation found."""

def polyhedron(vertices_coordinates: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], vertices_defining_each_facet: typing.Annotated[numpy.ndarray, numpy.int32, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], compute_altitude_with_sign: bool=True, just_compute_distance_to_vertices: bool=False) -> BodyShapeSettings:
    """Function for creating a polyhedron body shape model settings.
    
    Function for settings object, defining a polyhedron shape model.
    
    Note 1: The evaluation of the altitude with a polyhedron model tends to be computationally expensive. To reduce the
    computational time, it might be useful to instead define a hybrid shape model (see
    :func:`~tudatpy.dynamics.environment_setup.shape.hybrid`), which allows using a high-resolution
    polyhedron (with a large number of facets) at low altitudes and a low-resolution one (with smaller number of facets)
    at high-altitudes.
    
    Note 2: If the goal of using the shape model is only to detect collisions with the surface and not to explicitly
    obtain the altitude, it is instead recommended to use the Laplacian of the gravitational potential (see
    :func:`~tudatpy.dynamics.propagation_setup.dependent_variable.gravity_field_laplacian_of_potential`).
    This allows reducing the computational time, but is only valid if the same polyhedron model that is used to define
    the gravitational acceleration should also be used to detect the impacts.
    
    
    Parameters
    ----------
    vertices_coordinates : numpy.ndarray
        Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
        row per vertex, 3 columns).
    
    vertices_defining_each_facet : numpy.ndarray
        Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
        the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
        when seen from the outside of the polyhedron.
    
    compute_altitude_with_sign : bool, default=True
        Flag indicating whether the altitude should be computed with sign (i.e. >0 if above surface, <0 otherwise) or
        having always a positive value. If the the sign of the altitude is not relevant, then setting it to *false* is
        recommended, as it reduces the CPU time.
    
    just_compute_distance_to_vertices : bool, default=False
        Flag indicating whether the altitude should be computed just with respect to the polyhedron vertices (if flag
        is set to *true*) or to all polyhedron features (vertices, facets and edges; happens if flag is set to
        *false*). Depending on the application, it might be useful to set the flag to *true* for medium to high
        altitudes, as it allows significantly reducing the CPU time (the resulting altitude errors depend on the
        resolution of the used polyhedron and altitude itself).
    
    Returns
    -------
    PolyhedronBodyShapeSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.shape.BodyShapeSettings` derived
        :class:`~tudatpy.dynamics.environment_setup.shape.PolyhedronBodyShapeSettings` class"""

def spherical(radius: float) -> BodyShapeSettings:
    """Function for creating spherical body shape model settings.
    
    Function for settings object, defining strictly spherical body shape model entirely from single radius parameter.
    
    
    Parameters
    ----------
    radius : float
        Radius specifying spherical body shape.
    Returns
    -------
    SphericalBodyShapeSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.shape.BodyShapeSettings` derived :class:`~tudatpy.dynamics.environment_setup.shape.SphericalBodyShapeSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create a :class:`~tudatpy.dynamics.environment_setup.shape.BodyShapeSettings` using a perfectly spherical shape model:
    
    .. code-block:: python
    
       # define parameters describing perfectly spherical model
       body_radius = 6378.0E3
       # create shape model settings
       body_settings.get( "Earth" ).shape_settings = environment_setup.shape.spherical( body_radius )"""

def spherical_spice() -> BodyShapeSettings:
    """Function for creating spherical body shape model settings entirely from spice.
    
    Function for settings object, defining spherical body shape model entirely from spice parameters.
    
    Returns
    -------
    BodyShapeSettings
        Instance of :class:`~tudatpy.dynamics.environment_setup.shape.BodyShapeSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create a :class:`~tudatpy.dynamics.environment_setup.shape.BodyShapeSettings` using a perfectly spherical shape model and data from Spice:
    
    .. code-block:: python
    
       # create shape model settings
       body_settings.get( "Earth" ).shape_settings = environment_setup.shape.spherical_spice( )"""