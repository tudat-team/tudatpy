import numpy
import pybind11_stubgen.typing_ext
import typing
__all__ = ['CentralGravityFieldSettings', 'FromFileSphericalHarmonicsGravityFieldSettings', 'GravityFieldSettings', 'GravityFieldType', 'PolyhedronGravityFieldSettings', 'PredefinedSphericalHarmonicsModel', 'SphericalHarmonicsGravityFieldSettings', 'central', 'central_gravity', 'central_spice', 'central_spice_gravity', 'egm96', 'from_file_spherical_harmonic', 'gggrx1200', 'ggm02c', 'ggm02s', 'glgm3150', 'goco05c', 'jgmess160a', 'jgmro120d', 'lpe200', 'polyhedron_from_density', 'polyhedron_from_mu', 'polyhedron_gravity', 'predefined_spherical_harmonic', 'ring_gravity', 'ring_model', 'sh_triaxial_ellipsoid_from_density', 'sh_triaxial_ellipsoid_from_gravitational_parameter', 'shgj180u', 'spherical_harmonic', 'spherical_harmonic_gravity', 'spherical_harmonic_triaxial_body']

class CentralGravityFieldSettings(GravityFieldSettings):
    """`GravityFieldSettings` derived class defining settings of point mass gravity field.
    
    Derived class of `GravityFieldSettings` for central gravity fields, which are defined by a single gravitational parameter."""

    @property
    def gravitational_parameter(self) -> float:
        """
                 Gravitational parameter of central gravity field.
        
                 :type: float
        """

    @gravitational_parameter.setter
    def gravitational_parameter(self, arg1: float) -> None:
        ...

class FromFileSphericalHarmonicsGravityFieldSettings(SphericalHarmonicsGravityFieldSettings):
    """No documentation found."""

class GravityFieldSettings:
    """Base class for providing settings for automatic gravity field model creation.
    
    This class is a functional base class for settings of gravity field models that require no information in addition to their type.
    Gravity field model classes requiring additional information must be created using an object derived from this class."""

    @property
    def gravity_field_type(self) -> GravityFieldType:
        """
                 **read-only**
        
                 Type of gravity field model that is to be created.
        
                 :type: GravityFieldType
        """

class GravityFieldType:
    """Enumeration of gravity field types.
    
             Enumeration of gravity field types supported by tudat.
    
    
    
    
    
          
    
    Members:
    
      central_gravity : 
          
    
      central_spice_gravity : 
          
    
      spherical_harmonic_gravity : 
          
    
      polyhedron_gravity : 
          
    
      ring_gravity : No documentation found."""
    __members__: typing.ClassVar[dict[str, GravityFieldType]]
    central_gravity: typing.ClassVar[GravityFieldType]
    central_spice_gravity: typing.ClassVar[GravityFieldType]
    polyhedron_gravity: typing.ClassVar[GravityFieldType]
    ring_gravity: typing.ClassVar[GravityFieldType]
    spherical_harmonic_gravity: typing.ClassVar[GravityFieldType]

    def __eq__(self, other: typing.Any) -> bool:
        ...

    def __getstate__(self) -> int:
        ...

    def __hash__(self) -> int:
        ...

    def __index__(self) -> int:
        ...

    def __init__(self, value: int) -> None:
        ...

    def __int__(self) -> int:
        ...

    def __ne__(self, other: typing.Any) -> bool:
        ...

    def __repr__(self) -> str:
        ...

    def __setstate__(self, state: int) -> None:
        ...

    def __str__(self) -> str:
        ...

    @property
    def name(self) -> str:
        ...

    @property
    def value(self) -> int:
        ...

class PolyhedronGravityFieldSettings(GravityFieldSettings):
    """`GravityFieldSettings` derived class defining settings of a polyhedron gravity field representation.
    
    Derived class of `GravityFieldSettings` for gravity fields, which are defined by a polyhedron gravity field representation."""

    @property
    def associated_reference_frame(self) -> str:
        """
                 Identifier for body-fixed reference frame with which the vertices coordinates are associated.
        
                 :type: str
        """

    @associated_reference_frame.setter
    def associated_reference_frame(self, arg1: str) -> None:
        ...

    @property
    def density(self) -> float:
        """
                 Density of the polyhedron.
        
                 :type: float
        """

    @density.setter
    def density(self, arg1: float) -> None:
        ...

    @property
    def gravitational_parameter(self) -> float:
        """
                 Gravitational parameter of gravity field.
        
                 :type: float
        """

    @gravitational_parameter.setter
    def gravitational_parameter(self, arg1: float) -> None:
        ...

    @property
    def vertices_coordinates(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
                 row per vertex, 3 columns).
        
        
                 :type: numpy.ndarray
        """

    @property
    def vertices_defining_each_facet(self) -> typing.Annotated[numpy.ndarray, numpy.int32, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
                 the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
                 when seen from the outside of the polyhedron.
        
        
                 :type: numpy.ndarray
        """

class PredefinedSphericalHarmonicsModel:
    """Enumeration of predefined spherical harmonics models.
    
    Enumeration of predefined spherical harmonics models supported by tudat, for which thee coefficient files are automatically available (downloaded from
    `here <https://github.com/tudat-team/tudat-resources/tree/master/resource/gravity_models>`_). The directory where these files are stored can be
    extracted using the :func:`~tudatpy.data.get_gravity_models_path` function.
                     
            
    
    Members:
    
      egm96 : 
    
    Coefficients for EGM96 Earth gravity field up to degree and order 200, (see `link <https://cddis.gsfc.nasa.gov/926/egm96/egm96.html>`__ )
    
    
    
      ggm02c : 
    
    Coefficients for the combined GGM02 Earth gravity field up to degree and order 200, (see `link <https://www2.csr.utexas.edu/grace/gravity/ggm02/>`__ )
    
                        
    
      ggm02s : 
    
    Coefficients for the GRACE-only GGM02 Earth gravity field up to degree and order 160, (see `link <https://www2.csr.utexas.edu/grace/gravity/ggm02/>`__ )
    
    
    
      goco05c : 
    
    Coefficients for the GOCO05c combined Earth gravity field up to degree and order 719, (see `link <https://www2.csr.utexas.edu/grace/gravity/ggm02/>`__ )
    
        
    
      glgm3150 : 
    
    Coefficients for the GLGM3150 Moon gravity field up to degree and order 150, (see `link <https://pds.nasa.gov/ds-view/pds/viewProfile.jsp?dsid=LP-L-RSS-5-GLGM3/GRAVITY-V1.0>`__ )
    
    
    
      lpe200 : 
    
    Coefficients for the LPE200 Moon gravity field up to degree and order 200, (see `link <https://pds.nasa.gov/ds-view/pds/viewProfile.jsp?dsid=LP-L-RSS-5-GLGM3/GRAVITY-V1.0>`__ )
    
        
    
      gggrx1200 : 
    
    Coefficients for the GRGM1200A Moon gravity field up to degree and order 1199, (see `link <https://pgda.gsfc.nasa.gov/products/50>`__ )
    
    
    
      jgmro120d : 
    
    Coefficients for the MRO120D Moon gravity field up to degree and order 120, (see `link <https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/>`__ )
    
    
    
      jgmess160a : 
    
    Coefficients for the MESS160A Moon gravity field up to degree and order 160, (see `link <https://pds-geosciences.wustl.edu/messenger/mess-h-rss_mla-5-sdp-v1/messrs_1001/data/shadr/jgmess_160a_sha.lbl>`__ )
    
    
    
      shgj180u : 
    
    Coefficients for the SHGJ180U Moon gravity field up to degree and order 180, (see `link <https://pds-geosciences.wustl.edu/mgn/mgn-v-rss-5-gravity-l2-v1/mg_5201/gravity/shgj120u.lbl>`__ )"""
    __members__: typing.ClassVar[dict[str, PredefinedSphericalHarmonicsModel]]
    egm96: typing.ClassVar[PredefinedSphericalHarmonicsModel]
    gggrx1200: typing.ClassVar[PredefinedSphericalHarmonicsModel]
    ggm02c: typing.ClassVar[PredefinedSphericalHarmonicsModel]
    ggm02s: typing.ClassVar[PredefinedSphericalHarmonicsModel]
    glgm3150: typing.ClassVar[PredefinedSphericalHarmonicsModel]
    goco05c: typing.ClassVar[PredefinedSphericalHarmonicsModel]
    jgmess160a: typing.ClassVar[PredefinedSphericalHarmonicsModel]
    jgmro120d: typing.ClassVar[PredefinedSphericalHarmonicsModel]
    lpe200: typing.ClassVar[PredefinedSphericalHarmonicsModel]
    shgj180u: typing.ClassVar[PredefinedSphericalHarmonicsModel]

    def __eq__(self, other: typing.Any) -> bool:
        ...

    def __getstate__(self) -> int:
        ...

    def __hash__(self) -> int:
        ...

    def __index__(self) -> int:
        ...

    def __init__(self, value: int) -> None:
        ...

    def __int__(self) -> int:
        ...

    def __ne__(self, other: typing.Any) -> bool:
        ...

    def __repr__(self) -> str:
        ...

    def __setstate__(self, state: int) -> None:
        ...

    def __str__(self) -> str:
        ...

    @property
    def name(self) -> str:
        ...

    @property
    def value(self) -> int:
        ...

class SphericalHarmonicsGravityFieldSettings(GravityFieldSettings):
    """`GravityFieldSettings` derived class defining settings of spherical harmonic gravity field representation.
    
    Derived class of `GravityFieldSettings` for gravity fields, which are defined by a spherical harmonic gravity field representation."""

    @property
    def associated_reference_frame(self) -> str:
        """
                 Identifier for body-fixed reference frame with which the coefficients are associated.
        
                 :type: str
        """

    @associated_reference_frame.setter
    def associated_reference_frame(self, arg1: str) -> None:
        ...

    @property
    def create_time_dependent_field(self) -> bool:
        """
                 Boolean that denotes whether the field should be created as time-dependent (even if no variations are imposed initially).
        
                 :type: bool
        """

    @create_time_dependent_field.setter
    def create_time_dependent_field(self, arg1: bool) -> None:
        ...

    @property
    def gravitational_parameter(self) -> float:
        """
                 Gravitational parameter of gravity field.
        
                 :type: float
        """

    @gravitational_parameter.setter
    def gravitational_parameter(self, arg1: float) -> None:
        ...

    @property
    def normalized_cosine_coefficients(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 Cosine spherical harmonic coefficients (geodesy normalized). Entry (i,j) denotes coefficient at degree i and order j.
        
                 :type: numpy.ndarray
        """

    @normalized_cosine_coefficients.setter
    def normalized_cosine_coefficients(self, arg1: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]) -> None:
        ...

    @property
    def normalized_sine_coefficients(self) -> typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]:
        """
                 Sine spherical harmonic coefficients (geodesy normalized). Entry (i,j) denotes coefficient at degree i and order j.
        
                 :type: numpy.ndarray
        """

    @normalized_sine_coefficients.setter
    def normalized_sine_coefficients(self, arg1: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')]) -> None:
        ...

    @property
    def reference_radius(self) -> float:
        """
                 **read-only**
        
                 Reference radius of spherical harmonic field expansion.
        
                 :type: float
        """

    @property
    def scaled_mean_moment_of_inertia(self) -> float:
        """
                 Value of the scaled mean moment of inertia :math:`I_{xx}+I_{yy}+I_{zz}/(MR^{2})`. This value does not influence the gravity field itself,
                 but together with the degree 2 gravity field coefficients defines the body's inertia tensor.
        
        
                 :type: float
        """

    @scaled_mean_moment_of_inertia.setter
    def scaled_mean_moment_of_inertia(self, arg1: float) -> None:
        ...

def central(gravitational_parameter: float) -> GravityFieldSettings:
    """Function for central gravity field settings object.
    
    Function for settings object, defining a point-mass gravity field model with user-defined gravitational parameter :math:`\\mu`. The gravitational potential is the defined as:
    
    .. math::
       U(\\mathbf{r})=\\frac{\\mu}{||\\mathbf{r}||}
    
    with :math:`\\mathbf{r}` the position vector measured from the body's center of mass.
    
    
    Parameters
    ----------
    gravitational_parameter : float
        Gravitational parameter defining the point-mass gravity field.
    Returns
    -------
    CentralGravityFieldSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.dynamics.environment_setup.gravity_field.CentralGravityFieldSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.gravity_field.GravityFieldSettings` for Earth using a simple central gravity field model:
    
    .. code-block:: python
    
       # define parameters describing central gravity model
       gravitational_parameter = 3.986e14
       # create gravity field settings
       body_settings.get( "Earth" ).gravity_field_settings = environment_setup.gravity_field.central( gravitational_parameter )"""

def central_spice(body_name_to_use: str='') -> GravityFieldSettings:
    """Function to create central gravity field settings from Spice settings.
    
    Function for settings object, defining a point-mass gravity field model. This function provides the same model as :func:`~tudatpy.dynamics.environment_setup.gravity_field.central`), but with gravitational parameter :math:`\\mu` from Spice.
    
    Parameters
    ----------
    body_name_to_use : str, default = ""
        Body from which Spice gravitational paramerer is queried (if empty, it uses the name of the body to which the settings are assigned, see :func:`~tudatpy.dynamics.environment_setup.ephemeris.direct_spice` for example of analogous functionality for spice ephemeris).
    
    Returns
    -------
    GravityFieldSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.gravity_field.GravityFieldSettings` class of gravity field type ``central_spice``
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.gravity_field.GravityFieldSettings` for Earth using a simple central gravity field model and data from Spice:
    
    .. code-block:: python
    
       # create gravity field settings
       body_settings.get( "Earth" ).gravity_field_settings = environment_setup.gravity_field.central_spice( )"""

def from_file_spherical_harmonic(file: str, maximum_degree: int, maximum_order: int, associated_reference_frame: str='', gravitational_parameter_index: int=0, reference_radius_index: int=1) -> GravityFieldSettings:
    """Function to load a custom spherical harmonics gravity field settings from a file.
    
    Function to load a custom spherical harmonics gravity field settings from a file. The file should contain **fully normalized** spherical harmonic coefficients.
    The associated gravitational parameter and reference radius should be given in m^3/s^2 and m, respectively. The file format should be the same as that used for the files
    in the directories `here <https://github.com/tudat-team/tudat-resources/tree/master/resource/gravity_models>`__. Specifically, the file should contain
    
    - The first line should be a series of text blocks (typically numerical data). Two of these blocks (by default the first and second one) should be the gravitational parameter and reference radius, respectively. The text block should be separated by spaces, tabs and/or commas
    - Each subsequent line should contain a set of spherical harmonic coefficients (first ordered in ascending order by degree, then in ascending order by order), where the first, second, third and fourth value of the line should be: degree :math:`l`, order :math:`m`, normalized cosine coefficient :math:`\\bar{C}_{lm}`, normalized sine coefficient :math:`\\bar{S}_{lm}`. Additional entries (for instance with coefficient uncertainties) are ignored.
    
    .. warning::
       
       The function expects exponents to be indicated by either "e" or "E". If the exponent is indicated by "d" or "D", the exponent will not be parsed and the function will read wrong coefficients! 
    
    The following example shows the first lines of a file with correct format:
    
    .. code-block:: text
    
      0.3986004415E+15 0.6378136300E+07
       0    0  1.000000000000E+00  0.000000000000E+00  0.00000E+00  0.00000E+00
       1    0  0.000000000000E+00  0.000000000000E+00  0.00000E+00  0.00000E+00
       1    1  0.000000000000E+00  0.000000000000E+00  0.00000E+00  0.00000E+00
       2    0 -4.841693259705E-04  0.000000000000E+00  4.68460E-11  0.00000E+00
       2    1 -2.189810040712E-10  1.467451636117E-09  7.75160E-12  7.81670E-12
       2    2  2.439349093502E-06 -1.400284857733E-06  7.80670E-12  7.80760E-12
       ...
    
    
    Parameters
    ----------
    file : str
        Full file path and name where th gravity field file is located
    maximum_degree : int
        Maximum degree of the coefficients that are to be loaded
    maximum_order : int
        Maximum order of the coefficients that are to be loaded
    associated_reference_frame : str, default = ""
        Name of the body-fixed reference frame to which the gravity field is to be fixed. If left empty, this reference frame will automatically be set to the body-fixed frame defined by this body's rotation (see :ref:`rotation_model` for specifying rotation models).
    gravitational_parameter_index : int, default = 0
        Index of the values in the file header (first line of file) that contains the gravitational parameter
    reference_radius_index : int, default = 1
        Index of the values in the file header (first line of file) that contains the reference radius
    Returns
    -------
    SphericalHarmonicsGravityFieldSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.dynamics.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.gravity_field.GravityFieldSettings` for the Moon using GGGRX spherical harmonics gravity model, up to degree and order 300:
    
    .. code-block:: python
    
      # Create the gravity field settings for Moon with Spherical Harmonics loaded in from a Spherical Harmonics file
      body_settings.get("Moon").gravity_field_settings = environment_setup.gravity_field.from_file_spherical_harmonic(
          r"...\\.tudat\\resource\\gravity_models\\Moon\\gggrx_1200l_sha.tab", 300, 300)"""

def polyhedron_from_density(density: float, vertices_coordinates: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], vertices_defining_each_facet: typing.Annotated[numpy.ndarray, numpy.int32, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], associated_reference_frame: str, gravitational_constant: float=6.67259e-11) -> GravityFieldSettings:
    """Function for creating a polyhedron gravity field settings object, using the density.
    
    Function for settings object, defining a gravity field model through a polyhedron.
    The ``associated_reference_frame`` must be the same frame ID as the target frame of the body’s rotation model.
    It represents the frame in which the polyhedron field is defined.
    
    The gravitational potential, acceleration, Laplacian of potential and Hessian of potential are computed according
    to :cite:t:`werner1997`.
    
    This function uses the density to define the gravity field. To instead use the
    gravitational parameter see :func:`~tudatpy.astro.gravitation.polyhedron_from_mu`.
    
    
    Parameters
    ----------
    density : float, default=TUDAT_NAN
        Density of the polyhedron.
    
    vertices_coordinates : numpy.ndarray
        Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
        row per vertex, 3 columns).
    
    vertices_defining_each_facet : numpy.ndarray
        Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
        the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
        when seen from the outside of the polyhedron.
    
    associated_reference_frame : str
        Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
    gravitational_constant : float, default=GRAVITATIONAL_CONSTANT
        Newton's gravitational constant G, used to computed the gravitational parameter
    
    Returns
    -------
    PolyhedronGravityFieldSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.dynamics.environment_setup.gravity_field.PolyhedronGravityFieldSettings` class"""

def polyhedron_from_mu(gravitational_parameter: float, vertices_coordinates: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], vertices_defining_each_facet: typing.Annotated[numpy.ndarray, numpy.int32, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], associated_reference_frame: str, gravitational_constant: float=6.67259e-11) -> GravityFieldSettings:
    """Function for creating a polyhedron gravity field settings object, using the gravitational parameter.
    
    Function for settings object, defining a gravity field model through a polyhedron.
    The ``associated_reference_frame`` must be the same frame ID as the target frame of the body’s rotation model.
    It represents the frame in which the polyhedron field is defined.
    
    The gravitational potential, acceleration, Laplacian of potential and Hessian of potential are computed according
    to :cite:t:`werner1997`.
    
    This function uses the gravitational parameter to define the gravity field. To instead use the density
    constant see :func:`~tudatpy.astro.gravitation.polyhedron_from_density`. Since both models tend to be computationally intensive,
    it is recommended to use polyhedra with the lowest number of facets that allows meeting the desired accuracy. The number of facets of a polyhedron
    model can be reduced using any mesh processing software, for example `PyMeshLab <https://pymeshlab.readthedocs.io/en/latest/>`__.
    Additionally, different functions to process a polyhedron are available in `Polyhedron utilities <https://py.api.tudat.space/en/latest/polyhedron_utilities.html>`__.
    
    
    Parameters
    ----------
    gravitational_parameter : float
        Gravitational parameter :math:`\\mu` of gravity field.
    vertices_coordinates : numpy.ndarray
        Cartesian coordinates of each polyhedron vertex. Entry (i,j) denotes vertex i, coordinate j (one
        row per vertex, 3 columns).
    
    vertices_defining_each_facet : numpy.ndarray
        Index (0 based) of the vertices constituting each facet. Entry (i,j) denotes facet i, and the jth vertex of
        the facet (one row per facet, 3 columns). In each row, the vertices' indices should be ordered counterclockwise
        when seen from the outside of the polyhedron.
    
    associated_reference_frame : str
        Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
    gravitational_constant : float, default=GRAVITATIONAL_CONSTANT
        Newton's gravitational constant G, used to computed the density
    
    Returns
    -------
    PolyhedronGravityFieldSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.dynamics.environment_setup.gravity_field.PolyhedronGravityFieldSettings` class"""

def predefined_spherical_harmonic(predefined_model: PredefinedSphericalHarmonicsModel, maximum_degree: int=-1) -> GravityFieldSettings:
    """Function for spherical harmonics gravity field settings of a predefined model.
    
    Function for spherical harmonics gravity field settings of a predefined model
    
    
    Parameters
    ----------
    predefined_model : PredefinedSphericalHarmonicsModel
        Identified for gravity field model that is to be loaded
    maximum_degree : int, default = -1
        Maximum degree and order to which the coefficients are to be loaded. If value is negative, all coefficients for the specified gravity field are loaded
    Returns
    -------
    SphericalHarmonicsGravityFieldSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.dynamics.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.gravity_field.GravityFieldSettings` for Earth using EGM96 spherical harmonics gravity model, up to degree and order 32:
    
    .. code-block:: python
    
      # Create the gravity field settings for Earth with Spherical Harmonics from a triaxial ellipsoid
      body_settings.get( "Earth" ).gravity_field_settings = environment_setup.gravity_field.predefined_spherical_harmonic(
          environment_setup.gravity_field.egm96, 32 )"""

def ring_model(gravitational_parameter: float, ring_radius: float, associated_reference_frame: str, elliptic_integral_s_from_d_and_b: bool) -> GravityFieldSettings:
    """No documentation found."""

def sh_triaxial_ellipsoid_from_density(axis_a: float, axis_b: float, axis_c: float, density: float, maximum_degree: int, maximum_order: int, associated_reference_frame: str, gravitational_constant: float=6.67259e-11) -> SphericalHarmonicsGravityFieldSettings:
    """Function for spherical harmonics gravity field settings object from triaxial ellipsoid parameters, using the density to define the mass distribution.
    
    Function for settings object, defining a gravity field model through spherical harmonic expansion of a homogeneous triaxial ellipsoid, same as :class:`~tudatpy.dynamics.environment_setup.gravity_field.spherical_harmonic`
    The constant mass distribution in the specified ellipsoid shape is expanded to obtain a spherical harmonic coefficient representation.
    Gravity fields from this setting object are expressed in normalized spherical harmonic coefficients.
    The constant mass distribution is defined by the density and gravitational constant (optional).
    The body-fixed x-, y- and z- axes are assumed to be along the A-, B- and C- axes.
    This function implements the models of (see :cite:t:`balmino1994`).
    
    
    Parameters
    ----------
    axis_a : float
        Dimension of largest axis of triaxial ellipsoid.
    axis_b : float
        Dimension of intermediate axis of triaxial ellipsoid.
    axis_c : float
        Dimension of smallest axis of triaxial ellipsoid.
    density : float
        Density of ellipsoid.
    maximum_degree : int
        Maximum degree of spherical harmonics expansion.
    maximum_order : int
        Maximum order of spherical harmonics expansion.
    associated_reference_frame : str
        Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
    gravitational_constant : float, default=physical_constants::GRAVITATIONAL_CONSTANT
        Gravitational constant G of the gravity field.
    Returns
    -------
    SphericalHarmonicsGravityFieldSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.dynamics.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.gravity_field.GravityFieldSettings` for Earth using the expansion of a homogeneous triaxial ellipsoid into a spherical harmonics gravity model:
    
    .. code-block:: python
    
      # Create the gravity field settings for Earth with Spherical Harmonics from a triaxial ellipsoid
      body_settings.get( "Earth" ).gravity_field_settings = environment_setup.gravity_field.spherical_harmonic_triaxial_ellipsoid_from_density(
          axis_a=6378171.88,
          axis_b=6378102.03,
          axis_c=6356752.24,
          density=5520,
          maximum_degree=50,
          maximum_order=50,
          associated_reference_frame="IAU_Earth" )"""

def sh_triaxial_ellipsoid_from_gravitational_parameter(axis_a: float, axis_b: float, axis_c: float, maximum_degree: int, maximum_order: int, associated_reference_frame: str, gravitational_parameter: float) -> SphericalHarmonicsGravityFieldSettings:
    """Function for spherical harmonics gravity field settings object from triaxial ellipsoid parameters, using the gravitational parameter to define the mass distribution..
    
    Function for settings object, defining a gravity field model through spherical harmonic expansion of a homogeneous triaxial ellipsoid, same as :class:`~tudatpy.dynamics.environment_setup.gravity_field.spherical_harmonic`
    The constant mass distribution in the specified ellipsoid shape is expanded to obtain a spherical harmonic coefficient representation.
    Gravity fields from this setting object are expressed in normalized spherical harmonic coefficients.
    The constant mass distribution is defined by the gravitational parameter.
    The body-fixed x-, y- and z- axes are assumed to be along the A-, B- and C- axes.
    This function implements the models of (see :cite:t:`balmino1994`).
    
    
    Parameters
    ----------
    axis_a : float
        Dimension of largest axis of triaxial ellipsoid.
    axis_b : float
        Dimension of intermediate axis of triaxial ellipsoid.
    axis_c : float
        Dimension of smallest axis of triaxial ellipsoid.
    maximum_degree : int
        Maximum degree of spherical harmonics expansion.
    maximum_order : int
        Maximum order of spherical harmonics expansion.
    associated_reference_frame : str
        Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
    gravitational_parameter : float
        Gravitational parameter :math:`\\mu` of gravity field.
    Returns
    -------
    SphericalHarmonicsGravityFieldSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.dynamics.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class"""

def spherical_harmonic(gravitational_parameter: float, reference_radius: float, normalized_cosine_coefficients: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], normalized_sine_coefficients: typing.Annotated[numpy.ndarray, numpy.float64, pybind11_stubgen.typing_ext.DynamicSize('m', 'n')], associated_reference_frame: str) -> GravityFieldSettings:
    """Function for creating a spherical harmonics gravity field settings object.
    
    Function for settings object, defining a gravity field model through spherical harmonic expansion.
    The ``associated_reference_frame`` must be the same frame ID as the target frame of the body’s rotation model.
    It represents the frame in which the spherical harmonic field is defined.
    
    The gravitational potential is the defined as:
    
    .. math::
       U(\\mathbf{r})=\\sum_{l=0}^{l_{max}}\\sum_{m=0}^{l}\\mu\\left(\\frac{{R}^{l}}{r^{l+1}}\\right)\\bar{P}_{lm}(\\sin\\phi)\\left(\\bar{C}_{lm}\\cos m\\theta+\\bar{S}_{lm}\\sin m\\theta\\right)
    
    with :math:`\\mathbf{r}` the position vector of the evaluation point, measured from the body's center of mass. The angles :math:`\\phi` and :math:`\\theta` are the body-fixed latitude and longitude of the evaluation point, and :math:`\\bar{P}_{lm}` is the associated Legendre polynomial (at degree/order :math:`l/m`).
    
    For the spherical harmonic gravity field (including other spherical harmonic functions), the normalized mean moment of inertia must be set by the user, to allow an inertia tensor to be computed. This is done using the :attr:`~tudatpy.dynamics.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings.scaled_mean_moment_of_inertia` attribute of the :class:`~tudatpy.dynamics.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class, as in the example below
    
    .. code-block:: python
    
      # Add gravity field model settings to body of spherical harmonic type
      body_settings.get( "Mars" ).gravity_field = ...
    
      # Add setting for moment of inertia
      body_settings.get( "Mars" ).gravity_field.scaled_mean_moment_of_inertia = 0.365
    
    This code snippet will automatically create a rigid body properties for Mars, with the inertia tensor computed from this value of 0.365 and the degree 2 gravity field coefficients. Note that, if gravity field variations are used for the body, time-variability of the degree 1- and 2- coefficients will be reflected in time-variability of the body's center of mass and inertia tensor.
    
    Note: Spherical harmonic coefficients used for this environment model must *always* be fully normalized.
    To normalize un-normalized spherical harmonic coefficients, see :func:`~tudatpy.astro.gravitation.normalize_spherical_harmonic_coefficients`.
    
    
    Parameters
    ----------
    gravitational_parameter : float
        Gravitational parameter :math:`\\mu` of gravity field.
    reference_radius : float
        Reference radius :math:`R` of spherical harmonic field expansion.
    normalized_cosine_coefficients : numpy.ndarray
        Cosine spherical harmonic coefficients (geodesy normalized). Entry (i,j) denotes coefficient :math:`\\bar{C}_{ij}` at degree i and order j.
        As such, note that entry (0,0) of cosine coefficients should be equal to 1.
    
    normalized_sine_coefficients : numpy.ndarray
        Sine spherical harmonic coefficients (geodesy normalized). Entry (i,j) denotes coefficient :math:`\\bar{S}_{ij}`
        at degree i and order j.
    
    associated_reference_frame : str
        Identifier for body-fixed reference frame with which the spherical harmonics coefficients are associated.
    Returns
    -------
    SphericalHarmonicsGravityFieldSettings
        Instance of the :class:`~tudatpy.dynamics.environment_setup.gravity_field.GravityFieldSettings` derived :class:`~tudatpy.dynamics.environment_setup.gravity_field.SphericalHarmonicsGravityFieldSettings` class
    
    
    
    
    
    Examples
    --------
    In this example, we create :class:`~tudatpy.dynamics.environment_setup.gravity_field.GravityFieldSettings` for Earth using a spherical harmonics gravity model:
    
    .. code-block:: python
    
      # Define the spherical harmonics gravity model
      gravitational_parameter = 3986004.415E+8
      reference_radius = 6378136.3
      # Normalized coefficients taken from https://cddis.nasa.gov/archive/egm96/general_info/egm96_to360.ascii
      # The above file is described in https://cddis.nasa.gov/archive/egm96/general_info/readme.egm96
      normalized_cosine_coefficients = [
          [1,                   0,                   0,                   0],
          [0,                   0,                   0,                   0],
          [-0.484165371736E-03, -0.186987635955E-09, 0.243914352398E-05,  0],
          [0.957254173792E-06,  0.202998882184E-05,  0.904627768605E-06,  0.721072657057E-06]
      ]
      normalized_sine_coefficients = [
          [0,                   0,                   0,                   0],
          [0,                   0,                   0,                   0],
          [0,                   0.119528012031E-08,  -0.140016683654E-05, 0],
          [0,                   0.248513158716E-06,  -0.619025944205E-06, 0.141435626958E-05]
      ]
      associated_reference_frame = "IAU_Earth"
      # Create the gravity field settings and add them to the body "Earth"
      body_settings.get( "Earth" ).gravity_field_settings = environment_setup.gravity_field.spherical_harmonic(
          gravitational_parameter,
          reference_radius,
          normalized_cosine_coefficients,
          normalized_sine_coefficients,
          associated_reference_frame )"""

def spherical_harmonic_triaxial_body(axis_a: float, axis_b: float, axis_c: float, density: float, maximum_degree: int, maximum_order: int, associated_reference_frame: str, gravitational_constant: float=6.67259e-11) -> GravityFieldSettings:
    ...
central_gravity: GravityFieldType
central_spice_gravity: GravityFieldType
egm96: PredefinedSphericalHarmonicsModel
gggrx1200: PredefinedSphericalHarmonicsModel
ggm02c: PredefinedSphericalHarmonicsModel
ggm02s: PredefinedSphericalHarmonicsModel
glgm3150: PredefinedSphericalHarmonicsModel
goco05c: PredefinedSphericalHarmonicsModel
jgmess160a: PredefinedSphericalHarmonicsModel
jgmro120d: PredefinedSphericalHarmonicsModel
lpe200: PredefinedSphericalHarmonicsModel
polyhedron_gravity: GravityFieldType
ring_gravity: GravityFieldType
shgj180u: PredefinedSphericalHarmonicsModel
spherical_harmonic_gravity: GravityFieldType